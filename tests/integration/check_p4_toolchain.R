#!/usr/bin/env Rscript
source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])),'../../tools/toolchain.R'))
root<-new_directory('toolchain-tests-')
fails<-function(expr)stopifnot(inherits(tryCatch({force(expr);NULL},error=identity),'error'))
dir.create(file.path(root,'conda-meta')); record<-list(name='pkg',version='1',build='a',sha256=strrep('a',64)); path<-file.path(root,'conda-meta/pkg.json')
write_json(record,path); check_records(root,list(record))
for(key in c('version','build','sha256')) {changed<-record; changed[[key]]<-'changed';write_json(changed,path);fails(check_records(root,list(record)))}
write_json(record,path); extra<-record;extra$name<-'extra';write_json(extra,file.path(root,'conda-meta/extra.json'));fails(check_records(root,list(record)))
lock<-file.path(root,'lock'); expected<-list(list(url='https://example.org/pkg.conda',md5=strrep('a',32)))
writeLines(c('@EXPLICIT',paste0(expected[[1]]$url,'#',expected[[1]]$md5)),lock);check_lock(lock,expected)
writeLines(c('@EXPLICIT',paste0(expected[[1]]$url,'#',strrep('b',32))),lock);fails(check_lock(lock,expected))
fails(new_destination(root));file.symlink(file.path(root,'absent'),file.path(root,'link'));fails(new_destination(file.path(root,'link')));stopifnot(new_destination(file.path(root,'new'))==file.path(root,'new'))
stub<-file.path(root,'mamba');writeLines(c('#!/bin/sh','printf "%s\\n" "$@" > "$P7_ARGS"','exit 7'),stub);Sys.chmod(stub,'0755')
with_environment(list(P7_ARGS=file.path(root,'args')),fails(stage_toolchain('alignment',file.path(root,'space ; literal'),stub)))
argv<-readLines(file.path(root,'args'));stopifnot(argv[match('--prefix',argv)+1L]==file.path(root,'space ; literal'),!file.exists(file.path(root,'space ; literal')))
with_environment(list(R_LIBS_USER='/untrusted',LD_LIBRARY_PATH='/wrong',OMP_NUM_THREADS='999'),with_environment(clean_environment(),stopifnot(Sys.getenv('R_LIBS_USER')=='/dev/null',Sys.getenv('LD_LIBRARY_PATH')=='',Sys.getenv('OMP_NUM_THREADS')=='1')))
# MD5-explicit mamba records can omit SHA256 and retain a legacy cache suffix.
# A verified cached archive supplies SHA256; absent or modified archives fail.
unlink(file.path(root,'conda-meta/extra.json'))
archive<-file.path(root,'pkg.conda');writeLines('pinned archive content',archive)
pinned<-c(record,list(url='https://example.org/pkg.conda',md5=unname(tools::md5sum(archive))))
pinned$sha256<-sha256(archive);installed<-pinned;installed$sha256<-NULL
installed$package_tarball_full_path<-file.path(root,'pkg.tar.bz2');write_json(installed,path)
check_records(root,list(pinned))
writeLines('modified archive',archive);fails(check_records(root,list(pinned)))
unlink(archive);fails(check_records(root,list(pinned)))
unlink(root,recursive=TRUE);cat('PASS exact package/lock records, absent output, literal failed staging, clean environment\n')
