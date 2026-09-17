#!/usr/bin/env Rscript
source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])),if(basename(dirname(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])))=='integration')'../../tools/lib_helpers.R' else 'lib_helpers.R'))
packages <- function(kind) read_json(file.path(ROOT,paste0('config/',kind,'-linux-64.provenance.json')))$packages
check_records <- function(prefix, expected) {
 records<-function(rows) {
  result<-list(); fields<-c('name','version','build','sha256')
  for(row in rows) { assert(all(vapply(fields,function(k)!is.null(row[[k]])&&nzchar(row[[k]]),TRUE))&&!row$name %in% names(result),'incomplete or duplicate package record'); result[[row$name]]<-row[fields] }
  result[sort(names(result))]
 }
 actual<-lapply(list.files(file.path(prefix,'conda-meta'),'[.]json$',full.names=TRUE),read_json)
 # Explicit MD5 installs may omit SHA256 from conda-meta and record a legacy
 # tar.bz2 cache suffix even for .conda archives. Verify the actual retained
 # archive against the pinned SHA256 instead of weakening the identity check.
 for(i in seq_along(actual)) if(is.null(actual[[i]]$sha256) || !nzchar(actual[[i]]$sha256)) {
  row<-actual[[i]]; match<-Filter(function(w)identical(w$name,row$name),expected)
  assert(length(match)==1L,'unexpected installed package');wanted<-match[[1]]
  assert(identical(row$url,wanted$url)&&identical(row$md5,wanted$md5),'installed package URL/MD5 mismatch')
  cached<-row$package_tarball_full_path
  assert(!is.null(cached)&&nzchar(cached),'missing retained package archive location')
  candidates<-unique(c(cached,file.path(dirname(cached),basename(wanted$url))))
  candidates<-candidates[file.exists(candidates)&!dir.exists(candidates)]
  assert(length(candidates)>0L,paste('SHA256 absent; retain pinned package archive for',row$name))
  actual[[i]]$sha256<-sha256(candidates[1])
 }
 assert(identical(records(expected),records(actual)),paste(prefix,'package lock mismatch'))
}
check_lock <- function(path, expected) {
 lines<-trimws(readLines(path)); lines<-lines[nzchar(lines)&!startsWith(lines,'#')]
 wanted<-vapply(expected,function(row)paste0(row$url,'#',row$md5),'')
 assert(length(lines)>0L && lines[1]=='@EXPLICIT' && identical(sort(lines[-1]),sort(wanted)),paste(path,'explicit URLs/checksums do not match provenance'))
}
stage_toolchain <- function(kind,prefix,mamba='mamba') {
 assert(Sys.info()[['sysname']]=='Linux' && Sys.info()[['machine']] %in% c('x86_64','amd64'),'locks support Linux x86_64 only')
 prefix<-new_destination(prefix); lock<-file.path(ROOT,paste0('config/',kind,'-linux-64.explicit.txt')); check_lock(lock,packages(kind))
 manager<-Sys.which(mamba); assert(nzchar(manager),'mamba executable not found')
 cat(jsonlite::toJSON(list(online_staging_argv=c(manager,'create','--yes','--prefix',prefix,'--file',lock)),auto_unbox=TRUE),'\n')
 with_environment(clean_environment(),run(c(manager,'create','--yes','--prefix',prefix,'--file',lock)))
 check_records(prefix,packages(kind)); cat('Staged and package-verified:',prefix,'\n')
}
preflight_toolchain <- function(alignment,rprefix,destination,r_kind='runtime-r') {
 alignment<-absolute(alignment); rprefix<-absolute(rprefix); destination<-new_destination(destination)
 check_lock(file.path(ROOT,paste0('config/',r_kind,'-linux-64.explicit.txt')),packages(r_kind)); check_records(rprefix,packages(r_kind))
 if(dir.exists(file.path(alignment,'conda-meta'))) {
  check_lock(file.path(ROOT,'config/alignment-linux-64.explicit.txt'),packages('alignment')); check_records(alignment,packages('alignment'))
 } else {
  # Runtime bundle schema is validated independently of historical conda records.
  assert(file.exists(file.path(alignment,'COMPLETE')),'incomplete native runtime bundle')
  bundle<-table_read(file.path(alignment,'provenance/files.tsv')); assert(nrow(bundle)>0L,'empty native runtime manifest')
  for(item in rows_list(bundle)) { assert(!grepl('(^/|(^|/)[.][.](/|$))',item$path),'unsafe runtime manifest path'); assert(sha256(file.path(alignment,item$path))==item$sha256,'runtime file hash mismatch') }
  local({ old<-getwd(); on.exit(setwd(old)); setwd(alignment); run(c('sha256sum','-c','SHA256SUMS')) })
  assert(!any(grepl('(^|/)(python[^/]*|snakemake)$',list.files(alignment,recursive=TRUE,all.files=TRUE))),'Python/Snakemake found in runtime bundle')
 }
 manifest<-file.path(ROOT,'vendor/toolchain/manifest.json'); evidence<-read_json(manifest)
 for(item in evidence$retained_files) assert(sha256(file.path(ROOT,item$path))==item$sha256,paste('retained recipe/license hash mismatch:',item$path))
 scratch<-new_directory('.preflight-',dirname(destination)); okay<-FALSE
 on.exit(if(!okay)cat('Failed preflight diagnostics retained at',scratch,'\n',file=stderr()))
 probes<-list(STAR=c('--version','^2[.]7[.]10b$'),'hisat2-align-s'=c('--version','(^|/)hisat2-align-s version 2[.]2[.]1$'),'hisat2-build-s'=c('--version','(^|/)hisat2-build-s version 2[.]2[.]1$'),samtools=c('--version','^samtools 1[.]17$'),featureCounts=c('-v','^featureCounts v2[.]0[.]6$'))
 observed<-with_environment(c(clean_environment(),list(PATH=paste(file.path(alignment,'bin'),Sys.getenv('PATH'),sep=':'))),{
  values<-list()
  for(name in names(probes)) { binary<-file.path(alignment,'bin',name); output<-run(c(binary,probes[[name]][1]),timeout=30)$output
   writeLines(output,file.path(scratch,paste0(name,'.version.txt'))); assert(any(grepl(probes[[name]][2],output)),paste(name,'pinned version probe failed'))
   values[[length(values)+1L]]<-list(name=name,executable=binary,sha256=sha256(binary),version_output=paste(output,collapse='\n')) }
  run(c(file.path(rprefix,'bin/Rscript'),'--vanilla',file.path(ROOT,'tools/preflight_r.R'),scratch),log=file.path(scratch,'R.log'),timeout=180); values
 })
 write_json(list(platform=as.list(Sys.info()),alignment_prefix=alignment,r_prefix=rprefix,tools=observed,tool_manifest_sha256=sha256(manifest),r_lock_sha256=sha256(file.path(ROOT,paste0('config/',r_kind,'-linux-64.explicit.txt'))),scope='Installed records, native runtime hashes, exact versions, R computation and PNG; no fetch or install.'),file.path(scratch,'verified.json'))
 assert(!exists_any(destination),'output appeared during preflight'); assert(file.rename(scratch,destination),'cannot publish preflight'); okay<-TRUE
 cat('PASS offline toolchain and R preflight:',destination,'\n')
}
main <- function() {
 argv<-commandArgs(TRUE); assert(length(argv)>0L,'Usage: toolchain.R stage KIND --prefix PATH | verify --kind KIND --prefix PATH | preflight --out PATH')
 command<-argv[1]; argv<-argv[-1]
 if(command=='stage') { assert(length(argv)>0L,'missing kind'); kind<-argv[1]; args<-parse_cli(list(prefix=NULL,mamba='mamba'),required='prefix',args=argv[-1]); stage_toolchain(kind,args$prefix,args$mamba) }
 else if(command=='verify') { args<-parse_cli(list(kind=NULL,prefix=NULL),required=c('kind','prefix'),args=argv); check_lock(file.path(ROOT,paste0('config/',args$kind,'-linux-64.explicit.txt')),packages(args$kind)); check_records(args$prefix,packages(args$kind)); cat('PASS package lock and installed records\n') }
 else if(command=='preflight') { args<-parse_cli(list(alignment=file.path(ROOT,'.deps/runtime-tools'),'r-prefix'=file.path(ROOT,'.deps/runtime-r'),'r-kind'='runtime-r',out=NULL),required='out',args=argv); preflight_toolchain(args$alignment,args[['r-prefix']],args$out,args[['r-kind']]) }
 else stop('unknown command')
}
if(sys.nframe()==0L)main()
