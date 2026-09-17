#!/usr/bin/env Rscript
source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])),'lib_helpers.R'))
# Parse only the source subset used by these five pinned rendered recipes.
# Unknown syntax fails; literal full recipes are always retained unchanged.
recipe_source<-function(path) {
 lines<-readLines(path);start<-which(lines=='source:');assert(length(start)==1L,'missing recipe source')
 end<-which(seq_along(lines)>start & grepl('^[^ #]',lines))[1]; if(is.na(end))end<-length(lines)+1L
 source<-list();for(line in lines[seq.int(start+1L,end-1L)]) {
  if(!nzchar(trimws(line)))next
  if(grepl('^  (url|sha256): ',line)){key<-sub('^  ([^:]+):.*','\\1',line);source[[key]]<-sub('^  [^:]+: ','',line)}
  else if(line=='  patches:')source$patches<-list()
  else if(startsWith(line,'    - '))source$patches<-c(source$patches,list(sub('^    - ', '',line)))
  else stop('unsupported source YAML: ',line)
 }
 assert(!is.null(source$url)&&grepl('^[0-9a-f]{64}$',source$sha256),'invalid source declaration');source
}
main<-function(){
 args<-parse_cli(list(prefix=NULL,out=file.path(ROOT,'vendor/toolchain')),required='prefix',positional='prefix');prefix<-absolute(args$prefix);target<-new_destination(args$out);dir.create(target)
 records<-lapply(sort(list.files(file.path(prefix,'conda-meta'),'[.]json$',full.names=TRUE)),read_json)
 manifest<-list(format=1,platform='linux-64',distribution='Pinned conda binaries; recipes retained as build provenance, not rebuilt locally.',transitive_dependencies='config/alignment-linux-64.provenance.json',source_hash_scope='Upstream archive hashes declared by retained recipes; upstream archives not downloaded by this recorder.',packages=list(),retained_files=list())
 for(name in c('star','hisat2','subread','samtools','htslib')) {
  matches<-Filter(function(r)r$name==name,records);assert(length(matches)==1L,'missing/duplicate package');record<-matches[[1]]
  assert(sha256(record$package_tarball_full_path)==record$sha256,'package archive hash mismatch')
  info<-file.path(record$extracted_package_dir,'info');dest<-file.path(target,name);dir.create(dest)
  assert(file.copy(file.path(info,'recipe'),dest,recursive=TRUE),'cannot retain recipe')
  if(dir.exists(file.path(info,'licenses')))assert(file.copy(file.path(info,'licenses'),dest,recursive=TRUE),'cannot retain licenses')
  for(filename in c('about.json','hash_input.json','index.json','git'))if(file.exists(file.path(info,filename)))copy_files(file.path(info,filename),dest)
  item<-record[c('name','version','build','url','sha256','license','depends')];item$source<-recipe_source(file.path(info,'recipe/meta.yaml'))
  item$recipe<-file.path('vendor/toolchain',name,'recipe/meta.yaml');item$build_script<-file.path('vendor/toolchain',name,'recipe/build.sh');item$build_flags_scope<-'Literal upstream recipe scripts and conda_build_config.yaml retained; no local rebuild claim.'
  manifest$packages[[length(manifest$packages)+1L]]<-item
 }
 for(path in sort(list.files(target,recursive=TRUE,all.files=TRUE,full.names=TRUE)))if(!dir.exists(path))manifest$retained_files[[length(manifest$retained_files)+1L]]<-list(path=file.path('vendor/toolchain',substring(path,nchar(target)+2L)),sha256=sha256(path))
 write_json(manifest,file.path(target,'manifest.json'));cat('Recorded',length(manifest$packages),'verified package archives and',length(manifest$retained_files),'provenance files\n')
}
if(sys.nframe()==0L)main()
