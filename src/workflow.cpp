#include "rnaseq/workflow.hpp"
#include "rnaseq/process.hpp"
#include "workflow/internal.hpp"
#include <charconv>
#include <cstdlib>

namespace rnaseq {
using namespace workflow_detail;
void workflow_command(int argc, char **argv) {
    const std::string command = argv[1];
    if (command != "workflow-plan" && command != "workflow-local" && command != "workflow-submit" &&
        command != "workflow-task")
        fail("unknown workflow command");
    if (command == "workflow-task" ? argc != 4 : argc != 3)
        fail("usage: rnaseq workflow-plan|workflow-local|workflow-submit CONFIG; workflow-task "
             "CONFIG index|array|finish");
    Workflow w(argv[2], command);
    if (command == "workflow-plan") {
        w.plan();
        return;
    }
    for (const auto &name : {"OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"})
        setenv(name, "1", 1);
    fs::create_directory(w.run);
    Lock lock(w.run / "workflow.lock", command == "workflow-task", command == "workflow-task");
    w.snapshot(command != "workflow-task");
    fs::create_directory(w.run / "align");
    if (command == "workflow-submit") {
        w.submit();
        return;
    }
    if (command == "workflow-local") {
        if (fs::exists(w.run / "submission"))
            fail("run has a submission record; use scheduler tasks, not a local driver");
        w.index();
        run_jobs(w.runs.rows.size(), w.local_jobs, [&](std::size_t i) { w.alignment(i); });
        w.finish();
        return;
    }
    if (command != "workflow-task")
        fail("unknown workflow command");
    const std::string task = argv[3];
    if (task == "index")
        w.index();
    else if (task == "finish")
        w.finish();
    else if (task == "array") {
        const char *value = getenv("SLURM_ARRAY_TASK_ID");
        if (!value)
            fail("SLURM_ARRAY_TASK_ID is required");
        std::string s = value;
        unsigned long i = 0;
        auto [end, ec] = std::from_chars(s.data(), s.data() + s.size(), i);
        if (ec != std::errc() || end != s.data() + s.size())
            fail("invalid SLURM_ARRAY_TASK_ID");
        w.alignment(i);
    } else
        fail("unknown workflow task");
}
} // namespace rnaseq
