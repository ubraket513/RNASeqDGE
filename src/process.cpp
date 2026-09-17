#include "rnaseq/process.hpp"
#include <cerrno>
#include <csignal>
#include <cstring>
#include <cstdlib>
#include <string_view>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <set>
#include <sys/prctl.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>
extern char** environ;
namespace rnaseq {
namespace {
volatile sig_atomic_t interrupted = 0;
void interrupt_handler(int signal) { interrupted = signal; }
struct Signals {
    struct sigaction old_int{}, old_term{};
    Signals() {
        interrupted = 0;
        struct sigaction action{}; action.sa_handler = interrupt_handler; sigemptyset(&action.sa_mask);
        sigaction(SIGINT, &action, &old_int); sigaction(SIGTERM, &action, &old_term);
    }
    ~Signals() { sigaction(SIGINT, &old_int, nullptr); sigaction(SIGTERM, &old_term, nullptr); }
};
}
namespace {
// Freeze parents before walking their children, closing the spawn/kill race.
// Tools may establish their own process groups, so group signalling alone is
// insufficient. Linux subreaping lets this driver reap orphaned descendants.
void stop_tree(pid_t pid) {
    kill(pid, SIGSTOP);
    std::error_code error;
    const auto tasks = std::filesystem::path("/proc") / std::to_string(pid) / "task";
    for(const auto& task : std::filesystem::directory_iterator(tasks,error)) {
        std::ifstream children(task.path()/"children"); pid_t child;
        while(children>>child)stop_tree(child);
    }
    kill(pid,SIGKILL);
}
void cancel_children() {
    for(;;) {
        std::ifstream children("/proc/self/task/"+std::to_string(getpid())+"/children");
        pid_t pid; while(children>>pid)stop_tree(pid);
        int status;
        const auto result=waitpid(-1,&status,WNOHANG);
        if(result<0 && errno==ECHILD)return;
        const timespec delay{0,10000000}; nanosleep(&delay,nullptr);
    }
}
struct Subreaper {
    int previous=0;
    Subreaper() {
        if(prctl(PR_GET_CHILD_SUBREAPER,&previous)<0 || prctl(PR_SET_CHILD_SUBREAPER,1)<0)
            throw std::runtime_error("cannot enable worker descendant reaping");
    }
    ~Subreaper(){prctl(PR_SET_CHILD_SUBREAPER,previous);}
};
}
void run_jobs(std::size_t count, std::size_t concurrency, const std::function<void(std::size_t)>& action) {
    if(!concurrency)throw std::runtime_error("job concurrency must be positive");
    Signals signals; Subreaper subreaper;
    std::set<pid_t> active;
    std::size_t next=0;
    std::cout.flush(); std::cerr.flush();
    try {
        while(next<count || !active.empty()) {
            if(interrupted)throw ProcessError("local jobs interrupted",128+interrupted);
            while(next<count && active.size()<concurrency && !interrupted) {
                const auto task=next++;
                const pid_t child=fork();
                if(child<0)throw std::runtime_error("worker fork failed");
                if(child==0) {
                    setpgid(0,0); signal(SIGINT,SIG_DFL);signal(SIGTERM,SIG_DFL);
                    int code=0;
                    try {action(task);}
                    catch(const ProcessError& e){std::cerr<<e.what()<<'\n';code=e.status;}
                    catch(const std::exception& e){std::cerr<<e.what()<<'\n';code=1;}
                    std::cout.flush();std::cerr.flush();_exit(code);
                }
                setpgid(child,child);active.insert(child);
            }
            for(auto it=active.begin();it!=active.end();) {
                int status=0;const auto result=waitpid(*it,&status,WNOHANG);
                if(result==0 || (result<0 && errno==EINTR)){++it;continue;}
                if(result<0)throw std::runtime_error("worker waitpid failed");
                it=active.erase(it);
                const int code=WIFEXITED(status)?WEXITSTATUS(status):128+WTERMSIG(status);
                if(code)throw ProcessError("local alignment worker failed with exit "+std::to_string(code),code);
            }
            if(!active.empty()){const timespec delay{0,10000000};nanosleep(&delay,nullptr);}
        }
        if(interrupted)throw ProcessError("local jobs interrupted",128+interrupted);
    }catch(...){cancel_children();throw;}
}
void run_process(const std::vector<std::string>& argv, const std::filesystem::path& log,
                 const std::filesystem::path& stdout_log) {
    if (argv.empty() || argv[0].empty() || argv[0][0] != '/') throw std::runtime_error("executable must be absolute");
    std::vector<char*> arguments;
    for (const auto& argument : argv) {
        if (argument.find('\0') != std::string::npos) throw std::runtime_error("NUL in argument");
        arguments.push_back(const_cast<char*>(argument.c_str()));
    }
    arguments.push_back(nullptr);
    // Build a private environment before fork. Upstream /usr/bin/env wrappers
    // must find interpreters staged beside the explicitly selected executable.
    const auto directory = std::filesystem::path(argv[0]).parent_path().string();
    if (directory.find(':') != std::string::npos)
        throw std::runtime_error("executable directory cannot contain a PATH separator");
    const char* inherited_path = std::getenv("PATH");
    std::vector<std::string> environment;
    for (char** entry = environ; *entry; ++entry)
        if (!std::string_view(*entry).starts_with("PATH=")) environment.emplace_back(*entry);
    environment.push_back("PATH=" + directory +
                          (inherited_path && *inherited_path ? ":" + std::string(inherited_path) : ""));
    std::vector<char*> envp;
    for (auto& entry : environment) envp.push_back(entry.data());
    envp.push_back(nullptr);
    int fd = open(log.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC, 0600);
    if (fd < 0) throw std::runtime_error("cannot open process log: " + log.string());
    if (fd < 3) { int copy = fcntl(fd, F_DUPFD_CLOEXEC, 3); close(fd); fd = copy; }
    if (fd < 0) throw std::runtime_error("cannot duplicate process log");
    int output_fd = fd;
    if (!stdout_log.empty()) {
        output_fd = open(stdout_log.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC, 0600);
        if (output_fd >= 0 && output_fd < 3) { int copy = fcntl(output_fd, F_DUPFD_CLOEXEC, 3); close(output_fd); output_fd = copy; }
        if (output_fd < 0) { close(fd); throw std::runtime_error("cannot open stdout log: " + stdout_log.string()); }
    }
    Signals signals;
    pid_t child = fork();
    if (child < 0) { if (output_fd != fd) close(output_fd); close(fd); throw std::runtime_error("fork failed"); }
    if (child == 0) {
        setpgid(0, 0);
        signal(SIGINT, SIG_DFL); signal(SIGTERM, SIG_DFL);
        if (dup2(output_fd, STDOUT_FILENO) < 0 || dup2(fd, STDERR_FILENO) < 0) _exit(126);
        if (output_fd != fd) close(output_fd);
        close(fd);
        execve(arguments[0], arguments.data(), envp.data());
        _exit(127);
    }
    if (output_fd != fd) close(output_fd);
    close(fd); setpgid(child, child);
    int status = 0;
    while (true) {
        if (interrupted) kill(-child, SIGKILL);
        const auto waited = waitpid(child, &status, WNOHANG);
        if (waited == 0) {
            // Bounded polling prevents a signal arriving just before waitpid
            // from leaving the parent blocked on an uncooperative child.
            const timespec delay{0, 10000000}; nanosleep(&delay, nullptr); continue;
        }
        if (waited == child) break;
        if (waited < 0 && errno == EINTR) continue;
        kill(-child, SIGKILL);
        while (waitpid(child, &status, 0) < 0 && errno == EINTR) {}
        throw std::runtime_error("waitpid failed");
    }
    int code = interrupted ? 128 + interrupted : WIFEXITED(status) ? WEXITSTATUS(status) : 128 + WTERMSIG(status);
    if (code) throw ProcessError(argv[0] + " failed with exit " + std::to_string(code) + "; log: " + log.string(), code);
}
}
