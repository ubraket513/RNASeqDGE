#include "rnaseq/workflow.hpp"
#include "rnaseq/alignment.hpp"
#include "rnaseq/counts.hpp"
#include "rnaseq/file_hash.hpp"
#include "rnaseq/manifests.hpp"
#include "rnaseq/process.hpp"
#include "rnaseq/table.hpp"
#include <algorithm>
#include <charconv>
#include <chrono>
#include <memory>
#include <limits>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <fcntl.h>
#include <sys/file.h>
#include <sys/resource.h>
#include <sched.h>
#include <unistd.h>
namespace rnaseq { namespace {
namespace fs = std::filesystem;
using Settings = std::map<std::string,std::string>;
void fail(const std::string& s) { throw std::runtime_error("workflow: " + s); }
std::string read(const fs::path& p) {
    std::ifstream f(p); if (!f) fail("cannot read " + p.string());
    return {std::istreambuf_iterator<char>(f), {}};
}
void write(const fs::path& p,const std::string& s) {
    std::ofstream f(p); f << s; f.close(); if (!f) fail("cannot write " + p.string());
}
void atomic_write(const fs::path& p,const std::string& s) {
    const auto tmp=p.string()+".tmp-"+std::to_string(getpid()); write(tmp,s); fs::rename(tmp,p);
}
void clean(const std::string& s) { if (s.find_first_of("\t\r\n")!=std::string::npos) fail("control character in path/value"); }
fs::path resolve(const fs::path& base,const std::string& s) { clean(s); return fs::weakly_canonical(base / fs::path(s)); }
int number(const std::string& s,const std::string& name) {
    int n=0; auto [end,ec]=std::from_chars(s.data(),s.data()+s.size(),n);
    if(ec!=std::errc() || end!=s.data()+s.size() || n<1) fail(name+" must be a positive integer");
    return n;
}
std::size_t col(const Table& t,const std::string& name) {
    auto i=std::find(t.header.begin(),t.header.end(),name); if(i==t.header.end()) fail("missing column "+name); return i-t.header.begin();
}
std::string table_text(const Table& t) {
    std::ostringstream o; auto row=[&](const auto& r){bool first=true;for(const auto& v:r){if(!first)o<<'\t'; first=false; clean(v); o<<v;}o<<'\n';};
    row(t.header); for(const auto& r:t.rows)row(r); return o.str();
}
struct Lock {
    int fd=-1;
    Lock(const fs::path& p,bool shared=false,bool wait=false) {
        fd=open(p.c_str(),O_RDWR|O_CREAT|O_CLOEXEC,0600);
        if(fd<0)fail("cannot open lock "+p.string());
        if(flock(fd,(shared?LOCK_SH:LOCK_EX)|(wait?0:LOCK_NB))!=0){close(fd);fd=-1;fail("busy lock "+p.string());}
    }
    ~Lock(){if(fd>=0)close(fd);} Lock(const Lock&)=delete;
};
void quarantine(const fs::path& p) {
    if(!fs::exists(p))return;
    std::string suffix=".invalid-"+std::to_string(getpid()); int n=0;
    auto target=p.string()+suffix; while(fs::exists(target))target=p.string()+suffix+"-"+std::to_string(++n);
    fs::rename(p,target); std::cout<<"quarantine "<<target<<'\n';
}
// Output manifests include every regular file; additions, removals and changes invalidate.
std::string inventory(const fs::path& dir) {
    std::vector<fs::path> files;
    for(const auto& e:fs::recursive_directory_iterator(dir)) {
        if(e.is_symlink())fail("symlink in stage output "+e.path().string());
        if(e.is_regular_file() && e.path()!=dir/"STAGE.tsv" && !(e.path().parent_path()==dir && e.path().filename().string().find("STAGE.tsv.tmp-")==0))files.push_back(e.path());
    }
    std::sort(files.begin(),files.end()); std::string s;
    for(const auto& p:files)s+=sha256_file(p)+"\t"+p.lexically_relative(dir).string()+"\n";
    if(files.empty())fail("empty stage "+dir.string());
    return s;
}
bool valid(const fs::path& dir,const std::string& identity) {
    try{return fs::is_regular_file(dir/"STAGE.tsv") && read(dir/"STAGE.tsv")==identity+inventory(dir);}catch(const std::exception&){return false;}
}
std::string shell_quote(const std::string& value) { std::string s="'";for(char c:value)s+=c=='\''?"'\\''":std::string(1,c);return s+"'"; }
struct Workflow {
    Settings c; fs::path config,run,exe,fasta,gtf; Table runs,refs; std::string source_identity;
    std::map<std::string,std::string> snapshots;
    std::vector<fs::path> sources;
    int threads=1,workers=1,local_jobs=1,local_cpus=1,local_mem_mb=0,local_job_mem_mb=0;
    bool memory_fallback=true,slurm_resources=false;
    explicit Workflow(const fs::path& path,const std::string& command) {
        config=fs::canonical(path); const auto base=config.parent_path();
        auto t=read_table(config.string()); if(t.header!=std::vector<std::string>{"key","value"})fail("config header must be key/value");
        const std::set<std::string> allowed{"samples","runs","references","analysis","contrasts","genes","annotation","bin_dir","rscript","r_script","tool_lock","r_lock","run_dir","index_cache","backend","threads","workers","star_sa_bases","star_chr_bits","slurm_bin_dir","slurm_cpus","slurm_mem_mb","slurm_time","slurm_partition","slurm_concurrency","local_jobs","local_cpus","local_mem_mb","local_job_mem_mb"};
        for(const auto& r:t.rows){if(!allowed.count(r[0])||r[1].empty()||!c.emplace(r[0],r[1]).second)fail("unknown, empty or duplicate config key: "+r[0]);clean(r[1]);}
        for(const auto& k:{"samples","runs","references","analysis","contrasts","genes","bin_dir","rscript","r_script","tool_lock","r_lock","run_dir"})if(!c.count(k))fail(std::string("required config key: ")+k);
        for(const auto& [k,v]:Settings{{"backend","hisat2"},{"threads","1"},{"workers","1"},{"star_sa_bases","14"},{"star_chr_bits","18"},{"slurm_concurrency","1"}})if(!c.count(k))c[k]=v;
        if(c.at("backend")!="hisat2"&&c.at("backend")!="star")fail("backend must be hisat2 or star");
        threads=number(c.at("threads"),"threads");workers=number(c.at("workers"),"workers");number(c.at("slurm_concurrency"),"slurm_concurrency");
        if(number(c.at("star_sa_bases"),"star_sa_bases")>14||number(c.at("star_chr_bits"),"star_chr_bits")>18)fail("STAR sizing exceeds allowed bounds");
        slurm_resources=command=="workflow-submit" || (command=="workflow-plan" && c.count("slurm_cpus"));
        // Local settings still have a strict schema even when submitting, but
        // the login host's affinity/allocation does not constrain remote jobs.
        for(const auto& key:{"local_jobs","local_cpus","local_mem_mb","local_job_mem_mb"})
            if(c.count(key))number(c.at(key),key);
        if(slurm_resources) {
            if(!c.count("slurm_cpus"))fail("submission requires slurm_cpus");
            local_cpus=number(c.at("slurm_cpus"),"slurm_cpus");
            local_mem_mb=c.count("slurm_mem_mb")?number(c.at("slurm_mem_mb"),"slurm_mem_mb"):0;
        }else {
            if(const char* budget=getenv("SLURM_CPUS_PER_TASK"))if(std::max(threads,workers)>number(budget,"SLURM_CPUS_PER_TASK"))fail("threads/workers exceeds SLURM_CPUS_PER_TASK");
            cpu_set_t affinity; CPU_ZERO(&affinity);
            local_cpus=sched_getaffinity(0,sizeof(affinity),&affinity)==0?CPU_COUNT(&affinity):1;
            if(c.count("local_cpus"))local_cpus=std::min(local_cpus,number(c.at("local_cpus"),"local_cpus"));
            if(const char* budget=getenv("SLURM_CPUS_PER_TASK"))local_cpus=std::min(local_cpus,number(budget,"SLURM_CPUS_PER_TASK"));
            if(std::max(threads,workers)>local_cpus)fail("threads/workers exceeds resolved local CPU budget");
            local_jobs=c.count("local_jobs")?number(c.at("local_jobs"),"local_jobs"):local_cpus/threads;
            local_jobs=std::min(local_jobs,local_cpus/threads);
            if(c.count("local_mem_mb"))local_mem_mb=number(c.at("local_mem_mb"),"local_mem_mb");
            if(c.count("local_job_mem_mb"))local_job_mem_mb=number(c.at("local_job_mem_mb"),"local_job_mem_mb");
            if(const char* memory=getenv("SLURM_MEM_PER_NODE")) {
                const int allocated=number(memory,"SLURM_MEM_PER_NODE");
                local_mem_mb=local_mem_mb?std::min(local_mem_mb,allocated):allocated;
            }
            if(const char* memory=getenv("SLURM_MEM_PER_CPU")) {
                const auto allocated=std::min<long long>(std::numeric_limits<int>::max(),
                    static_cast<long long>(number(memory,"SLURM_MEM_PER_CPU"))*local_cpus);
                local_mem_mb=local_mem_mb?std::min(local_mem_mb,static_cast<int>(allocated)):static_cast<int>(allocated);
            }
            memory_fallback=!local_mem_mb || !local_job_mem_mb;
            if(memory_fallback)local_jobs=1;
            else {
                if(local_job_mem_mb>local_mem_mb)fail("local_job_mem_mb exceeds resolved memory budget");
                local_jobs=std::min(local_jobs,local_mem_mb/local_job_mem_mb);
            }
        }
        if(c.count("slurm_cpus")&&number(c.at("slurm_cpus"),"slurm_cpus")<std::max(threads,workers))fail("slurm_cpus smaller than threads/workers");
        if(c.count("slurm_mem_mb"))number(c.at("slurm_mem_mb"),"slurm_mem_mb");
        if(c.count("slurm_time") && c.at("slurm_time").find_first_not_of("0123456789:-")!=std::string::npos)fail("slurm_time must use numeric Slurm time syntax");
        for(const auto& k:{"samples","runs","references","analysis","contrasts","genes","annotation","bin_dir","rscript","r_script","tool_lock","r_lock","run_dir","index_cache","slurm_bin_dir"})if(c.count(k))c[k]=resolve(base,c.at(k)).string();
        run=c.at("run_dir"); exe=fs::canonical("/proc/self/exe");
        if(!fs::is_directory(run.parent_path()))fail("run_dir parent must exist");
        if(c.count("index_cache")&&!fs::is_directory(fs::path(c.at("index_cache")).parent_path()))fail("index_cache parent must exist");
        if(c.at("bin_dir").find(':')!=std::string::npos)fail("bin_dir contains PATH separator");
        validate_bundle(c.at("samples"),c.at("runs"),c.at("references"),c.at("analysis"),c.at("contrasts"));
        runs=read_table(c.at("runs")); refs=read_table(c.at("references"));
        for(auto& r:runs.rows)for(const auto& key:{"fastq_1","fastq_2"}) {
            auto& v=r[col(runs,key)];if(v.empty())continue;v=resolve(fs::path(c.at("runs")).parent_path(),v).string();
            if(fs::path(v).extension()==".gz")fail("P4 requires uncompressed FASTQ");
            sources.emplace_back(v);
        }
        for(auto& r:refs.rows) {
            auto& v=r[col(refs,"path")];v=resolve(fs::path(c.at("references")).parent_path(),v).string();sources.emplace_back(v);
            if(r[col(refs,"role")]=="genome")fasta=v;else gtf=v;
        }
        auto genes=read_table(c.at("genes")); if(genes.header!=std::vector<std::string>{"gene_id"}||genes.rows.empty())fail("genes must be nonempty gene_id table");
        std::set<std::string> ids;for(const auto& r:genes.rows)if(r[0].empty()||!ids.insert(r[0]).second)fail("empty/duplicate gene_id");
        if(c.count("annotation"))validate_annotation(c.at("genes"),c.at("annotation"));
        sources.push_back(config);sources.push_back(exe);
        for(const auto& k:{"samples","runs","references","analysis","contrasts","genes","annotation","rscript","r_script","tool_lock","r_lock"})if(c.count(k))sources.emplace_back(c.at(k));
        for(const auto& k:{"samples","analysis","contrasts","genes","annotation","tool_lock","r_lock"})if(c.count(k))snapshots[std::string(k)+".tsv"]=read(c.at(k));
        snapshots["r-preflight.R"]="stopifnot(getRversion() == '4.5.3')\npins <- c(DESeq2='1.50.2', apeglm='1.32.0', BiocParallel='1.44.0', pheatmap='1.0.13')\nfor (p in names(pins)) stopifnot(as.character(packageVersion(p)) == pins[[p]])\ncat('Pinned R package versions verified\\n')\n";
        snapshots["runs.tsv"]=table_text(runs); snapshots["references.tsv"]=table_text(refs);
        std::string resolved="key\tvalue\n";for(const auto& [k,v]:c)resolved+=k+"\t"+v+"\n";snapshots["config.tsv"]=resolved;
        const auto bin=fs::path(c.at("bin_dir"));
        for(const auto& name:std::vector<std::string>{c.at("backend")=="star"?"STAR":"hisat2-align-s",c.at("backend")=="star"?"STAR":"hisat2-build-s","samtools","featureCounts"})if(access((bin/name).c_str(),X_OK)!=0)fail("missing executable "+(bin/name).string());
        if(access(c.at("rscript").c_str(),X_OK)!=0)fail("rscript not executable");
        // Include wrapper companions/interpreters in the selected tool directory.
        for(const auto& e:fs::directory_iterator(bin))if(e.is_regular_file())sources.push_back(e.path());
        std::sort(sources.begin(),sources.end());sources.erase(std::unique(sources.begin(),sources.end()),sources.end());
        source_identity="sha256\tpath\n";
        for(const auto& p:sources){if(!fs::is_regular_file(p)||access(p.c_str(),R_OK)!=0)fail("unreadable input "+p.string());source_identity+=sha256_file(p)+"\t"+p.string()+"\n";}
    }
    void verify_sources() const {
        std::string current="sha256\tpath\n";for(const auto& p:sources)current+=sha256_file(p)+"\t"+p.string()+"\n";
        if(current!=source_identity)fail("inputs changed during execution; use a new run directory");
    }
    void snapshot(bool create) {
        const auto dir=run/"snapshot";
        if(fs::exists(dir)) {
            if(!valid(dir,"snapshot-v1\n")||read(dir/"sources.tsv")!=source_identity)fail("snapshot/source/config/tool hash changed; use a new run directory");
            for(const auto& [name,data]:snapshots)if(read(dir/name)!=data)fail("snapshot differs; use a new run directory");
            return;
        }
        if(!create)fail("snapshot is missing; prepare with workflow-local or workflow-submit");
        auto tmp=run/("snapshot.partial-"+std::to_string(getpid()));quarantine(tmp);fs::create_directory(tmp);
        for(const auto& [name,data]:snapshots)write(tmp/name,data);
        write(tmp/"sources.tsv",source_identity);
        verify_sources();atomic_write(tmp/"STAGE.tsv","snapshot-v1\n"+inventory(tmp));fs::rename(tmp,dir);
    }
    std::string base_id()const{return "workflow-v1\n"+sha256_file(run/"snapshot"/"STAGE.tsv")+"\n";}
    void stage(const fs::path& dir,const std::string& id,const std::function<void()>& action) {
        using Clock=std::chrono::steady_clock;
        const auto started=Clock::now();
        auto elapsed=[](auto start){return std::chrono::duration<double>(Clock::now()-start).count();};
        rusage self_before{},child_before{};
        getrusage(RUSAGE_SELF,&self_before);getrusage(RUSAGE_CHILDREN,&child_before);
        double validation=0,execution=0,sources_time=0,publication=0;
        const auto profile=[&](const std::string& outcome) {
            rusage self{},child{};getrusage(RUSAGE_SELF,&self);getrusage(RUSAGE_CHILDREN,&child);
            auto seconds=[](timeval t){return t.tv_sec+t.tv_usec/1e6;};
            std::ostringstream o;
            o<<"stage\tstatus\tvalidation_seconds\texecution_seconds\tsource_validation_seconds\tpublication_hash_seconds\twall_seconds\tuser_cpu_seconds\tsystem_cpu_seconds\tprocess_peak_rss_kb\tchildren_peak_rss_kb\tlocal_jobs\tlocal_cpus\tlocal_mem_mb\tlocal_job_mem_mb\n";
            o<<dir.string()<<'\t'<<outcome<<'\t'<<validation<<'\t'<<execution<<'\t'<<sources_time<<'\t'<<publication<<'\t'<<elapsed(started)<<'\t'
             <<seconds(self.ru_utime)+seconds(child.ru_utime)-seconds(self_before.ru_utime)-seconds(child_before.ru_utime)<<'\t'
             <<seconds(self.ru_stime)+seconds(child.ru_stime)-seconds(self_before.ru_stime)-seconds(child_before.ru_stime)<<'\t'
             <<self.ru_maxrss<<'\t'<<child.ru_maxrss<<'\t'<<local_jobs<<'\t'<<local_cpus<<'\t'<<local_mem_mb<<'\t'<<local_job_mem_mb<<'\n';
            fs::create_directories(run/"profiles");
            atomic_write(run/"profiles"/(dir.filename().string()+"-"+std::to_string(getpid())+"-"+std::to_string(started.time_since_epoch().count())+".tsv"),o.str());
        };
        try {
            Lock lock(dir.string()+".lock");
            auto start=Clock::now();const bool reusable=valid(dir,id);validation=elapsed(start);
            if(reusable){std::cout<<"reuse "<<dir<<'\n';profile("reused");return;}
            quarantine(dir);std::cout<<"run "<<dir<<'\n';
            start=Clock::now();
            try {action();}catch(...){execution=elapsed(start);throw;}
            execution=elapsed(start);
            start=Clock::now();verify_sources();sources_time=elapsed(start);
            write(dir/"generation.txt",std::to_string(std::chrono::system_clock::now().time_since_epoch().count())+"\n");
            start=Clock::now();atomic_write(dir/"STAGE.tsv",id+inventory(dir));publication=elapsed(start);
            profile("executed");
        }catch(...){try{profile("failed");}catch(...){}throw;}
    }
    std::string cache_key()const {
        // SHA256 of complete index identity is produced via the ordinary file hasher.
        std::string s="genome-only-v1\n"+c.at("backend")+"\n"+sha256_file(fasta)+"\n"+sha256_file(gtf)+"\n";
        for(const auto& k:{"threads","star_sa_bases","star_chr_bits","tool_lock"})s+=c.at(k)+"\n";
        s+=sha256_file(c.at("tool_lock"))+"\n"+sha256_file(exe)+"\n";
        std::vector<fs::path> tools;
        for(const auto& e:fs::directory_iterator(c.at("bin_dir")))if(e.is_regular_file())tools.push_back(e.path());
        std::sort(tools.begin(),tools.end());
        for(const auto& tool:tools)s+=tool.filename().string()+"\t"+sha256_file(tool)+"\n";
        return s;
    }
    fs::path index_path() {
        if(!c.count("index_cache"))return run/"index";
        const auto key=cache_key(); const auto temp=run/("cache-key-"+std::to_string(getpid()));write(temp,key);auto hash=sha256_file(temp);fs::remove(temp);
        return fs::path(c.at("index_cache"))/hash;
    }
    void align(const std::string& cmd,const fs::path& out,const std::vector<std::string>& extra) {
        std::vector<std::string> args{"rnaseq",cmd,"--backend",c.at("backend"),"--bin-dir",c.at("bin_dir"),"--fasta",fasta.string(),"--gtf",gtf.string(),"--threads",c.at("threads"),"--output",out.string()};
        args.insert(args.end(),extra.begin(),extra.end());std::vector<char*> ptrs;for(auto& a:args)ptrs.push_back(a.data());alignment_command(ptrs.size(),ptrs.data());
    }
    void index() {
        if(c.count("index_cache"))fs::create_directories(c.at("index_cache"));
        const auto dir=index_path();
        const auto id=c.count("index_cache")?"cache-v1\n"+cache_key():base_id()+"index\n";
        stage(dir,id,[&]{std::vector<std::string> extra;if(c.at("backend")=="star")extra={"--star-sa-bases",c.at("star_sa_bases"),"--star-chr-bits",c.at("star_chr_bits")};align("index",dir,extra);});
    }
    std::string index_id() { return c.count("index_cache")?"cache-v1\n"+cache_key():base_id()+"index\n"; }
    void alignment(std::size_t i) {
        if(i>=runs.rows.size())fail("array task outside immutable run mapping");
        const auto idx=index_path();
        // Shared cache lock holds the checked index stable throughout alignment.
        Lock index_lock(idx.string()+".lock",true);
        if(!valid(idx,index_id()))fail("index incomplete or corrupt; rerun index stage");
        const auto& row=runs.rows[i]; const auto name=row[col(runs,"run_id")];const auto dir=run/"align"/name;
        const auto id=base_id()+"align\t"+name+"\n"+sha256_file(idx/"STAGE.tsv")+"\n";
        stage(dir,id,[&]{std::vector<std::string> extra{"--index",idx.string(),"--reads1",row[col(runs,"fastq_1")],"--layout",row[col(runs,"layout")],"--strandedness",row[col(runs,"strandedness")]};
            if(!row[col(runs,"fastq_2")].empty())extra.insert(extra.end(),{"--reads2",row[col(runs,"fastq_2")]});
            align("align-count",dir,extra);});
    }
    void finish() {
        const auto idx=index_path();Lock index_lock(idx.string()+".lock",true);if(!valid(idx,index_id()))fail("invalid index before merge");
        std::vector<std::unique_ptr<Lock>> alignment_locks;
        std::string merge_id=base_id()+"merge\n",inputs="run_id\tformat\tpath\tcolumn\n";
        for(const auto& row:runs.rows) {
            const auto name=row[col(runs,"run_id")];const auto dir=run/"align"/name;
            alignment_locks.push_back(std::make_unique<Lock>(dir.string()+".lock",true));
            const auto id=base_id()+"align\t"+name+"\n"+sha256_file(idx/"STAGE.tsv")+"\n";
            if(!valid(dir,id))fail("incomplete or corrupt alignment "+name);
            merge_id+=sha256_file(dir/"STAGE.tsv")+"\n";
            std::ifstream f(dir/"counts.txt");TsvReader reader(f,(dir/"counts.txt").string(),true);std::vector<std::string> header;
            if(!reader.next(header)||header.size()!=7||header[0]!="Geneid")fail("featureCounts header must contain exactly one BAM column");
            inputs+=name+"\tfeaturecounts\t"+(dir/"counts.txt").string()+"\t"+header[6]+"\n";
        }
        const auto snap=run/"snapshot",merged=run/"merge";
        stage(merged,merge_id,[&]{fs::create_directory(merged);write(merged/"inputs.tsv",inputs);merge_counts(snap/"samples.tsv",snap/"runs.tsv",snap/"genes.tsv",merged/"inputs.tsv",merged/"counts.tsv",c.count("annotation")?snap/"annotation.tsv":fs::path{});});
        const auto out=run/"analysis";
        stage(out,base_id()+"R\n"+sha256_file(merged/"STAGE.tsv")+"\n",[&]{
            // R owns its own atomic publication; wrapper logs live outside its destination.
            std::vector<std::string> args{c.at("rscript"),"--vanilla",c.at("r_script"),"--counts",(merged/"counts.tsv").string(),"--samples",(snap/"samples.tsv").string(),"--analysis",(snap/"analysis.tsv").string(),"--contrasts",(snap/"contrasts.tsv").string(),"--out",out.string(),"--workers",c.at("workers")};
            if(c.count("annotation"))args.insert(args.end(),{"--annotation",(snap/"annotation.tsv").string()});
            for(const auto& name:{"R_HOME","R_LIBS","LD_LIBRARY_PATH","LD_PRELOAD"})unsetenv(name);
            for(const auto& name:{"R_LIBS_USER","R_LIBS_SITE"})setenv(name,"/dev/null",1);
            const auto preflight_log=run/("R-preflight-"+std::to_string(getpid())+".log");
            quarantine(preflight_log);
            run_process({c.at("rscript"),"--vanilla",(snap/"r-preflight.R").string()},preflight_log);
            const auto log=run/("R-"+std::to_string(getpid())+".log");quarantine(log);run_process(args,log);
            if(!fs::is_directory(out))fail("R returned success without output directory");
            fs::copy_file(log,out/"workflow-R.log");
            fs::copy_file(preflight_log,out/"workflow-R-preflight.log");
        });
    }
    void plan()const {
        std::cout<<"DAG: index -> alignment/counting["<<runs.rows.size()<<"] -> sample merge -> offline R\n";
        for(const auto& [k,v]:c)std::cout<<k<<'\t'<<v<<'\n';
        std::cout<<"reference_fasta\t"<<fasta<<"\nreference_gtf\t"<<gtf<<"\nBLAS/OpenMP: 1\n";
        if(slurm_resources) {
            std::cout<<"resource_scope\tslurm_requested\nresolved_task_cpus\t"<<local_cpus
                     <<"\nresolved_task_mem_mb\t"<<local_mem_mb<<'\n';
        }else {
            std::cout<<"resource_scope\tlocal_allocation\n";
            std::cout<<"local_resolved_cpus\t"<<local_cpus<<"\nlocal_resolved_mem_mb\t"<<local_mem_mb
                     <<"\nlocal_job_mem_mb\t"<<local_job_mem_mb<<"\nlocal_resolved_jobs\t"<<local_jobs
                     <<"\nlocal_memory_policy\t"<<(memory_fallback?"serial fallback: total/per-job memory estimate missing":"reservation estimates; not a hard RSS limit")<<'\n';
        }
        std::cout<<table_text(runs);
    }
    void submit() {
        if(getenv("SLURM_JOB_ID"))fail("nested submission is forbidden");
        for(const auto& k:{"slurm_bin_dir","slurm_cpus","slurm_mem_mb","slurm_time"})if(!c.count(k))fail(std::string("submission requires ")+k);
        const auto bin=fs::path(c.at("slurm_bin_dir"));for(const auto& name:{"sbatch","scancel"})if(access((bin/name).c_str(),X_OK)!=0)fail(std::string("missing scheduler executable ")+name);
        if(fs::exists(run/"submission"))fail("submission already attempted; inspect recorded jobs before using a new run directory");
        fs::create_directory(run/"submission");const auto dir=run/"submission";
        for(const auto& task:{"index","array","finish"})write(dir/(std::string(task)+".sh"),"#!/bin/sh\nset -eu\nexec "+shell_quote(exe.string())+" workflow-task "+shell_quote(config.string())+" "+task+"\n");
        std::vector<std::string> accepted;
        try {
            for(const auto& task:{"index","array","finish"}) {
                std::vector<std::string> args{(bin/"sbatch").string(),"--parsable","--cpus-per-task="+c.at("slurm_cpus"),"--mem="+c.at("slurm_mem_mb"),"--time="+c.at("slurm_time"),"--output="+(dir/(std::string(task)+"-%A_%a.log")).string()};
                if(c.count("slurm_partition"))args.push_back("--partition="+c.at("slurm_partition"));
                if(!accepted.empty())args.push_back("--dependency=afterok:"+accepted.back());
                if(std::string(task)=="array")args.push_back("--array=0-"+std::to_string(runs.rows.size()-1)+"%"+c.at("slurm_concurrency"));
                args.push_back((dir/(std::string(task)+".sh")).string());
                std::string record;for(const auto& a:args)record+=a+"\n";write(dir/(std::string(task)+".argv"),record);
                const auto log=dir/(std::string(task)+".submit.stdout");
                run_process(args,dir/(std::string(task)+".submit.stderr"),log);
                std::string id=read(log);while(!id.empty()&&(id.back()=='\n'||id.back()=='\r'))id.pop_back();const auto semicolon=id.find(';');if(semicolon!=std::string::npos)id=id.substr(0,semicolon);
                if(id.empty()||id.find_first_not_of("0123456789")!=std::string::npos)fail("sbatch did not return a parsable numeric job ID; inspect "+log.string());
                accepted.push_back(id);std::string ids;for(const auto& job:accepted)ids+=job+"\n";atomic_write(dir/"accepted.txt",ids);std::cout<<"accepted "<<task<<" job "<<id<<'\n';
            }
            atomic_write(dir/"SUBMITTED","jobs accepted; live execution is not yet verified\n");
        }catch(...) {
            for(const auto& id:accepted)try{run_process({(bin/"scancel").string(),id},dir/("cancel-"+id+".log"));}catch(const std::exception& e){std::cerr<<"rollback cancellation failed: "<<e.what()<<'\n';}
            write(dir/"FAILED","submission failed; inspect accepted IDs and cancellation logs\n");throw;
        }
    }
};
} // namespace
void workflow_command(int argc,char** argv) {
    const std::string command=argv[1];
    if(command!="workflow-plan" && command!="workflow-local" && command!="workflow-submit" && command!="workflow-task")fail("unknown workflow command");
    if(command=="workflow-task" ? argc!=4 : argc!=3)fail("usage: rnaseq workflow-plan|workflow-local|workflow-submit CONFIG; workflow-task CONFIG index|array|finish");
    Workflow w(argv[2],command);if(command=="workflow-plan"){w.plan();return;}
    for(const auto& name:{"OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS","NUMEXPR_NUM_THREADS"})setenv(name,"1",1);
    fs::create_directory(w.run);Lock lock(w.run/"workflow.lock",command=="workflow-task",command=="workflow-task");
    w.snapshot(command!="workflow-task");fs::create_directory(w.run/"align");
    if(command=="workflow-submit"){w.submit();return;}
    if(command=="workflow-local") {if(fs::exists(w.run/"submission"))fail("run has a submission record; use scheduler tasks, not a local driver");w.index();run_jobs(w.runs.rows.size(),w.local_jobs,[&](std::size_t i){w.alignment(i);});w.finish();return;}
    if(command!="workflow-task")fail("unknown workflow command");
    const std::string task=argv[3];if(task=="index")w.index();else if(task=="finish")w.finish();else if(task=="array") {
        const char* value=getenv("SLURM_ARRAY_TASK_ID");if(!value)fail("SLURM_ARRAY_TASK_ID is required");std::string s=value;unsigned long i=0;auto [end,ec]=std::from_chars(s.data(),s.data()+s.size(),i);if(ec!=std::errc()||end!=s.data()+s.size())fail("invalid SLURM_ARRAY_TASK_ID");w.alignment(i);
    }else fail("unknown workflow task");
}
} // namespace rnaseq
