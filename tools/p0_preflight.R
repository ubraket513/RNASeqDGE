#!/usr/bin/env Rscript
# Read-only inventory of the supported native/R toolchain, never installation.
source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])),'lib_helpers.R'))
flags<-list('g++'='--version',make='--version',Rscript='--version',STAR='--version','hisat2-align-s'='--version','hisat2-build-s'='--version',featureCounts='-v',samtools='--version','fasterq-dump'='--version',sbatch='--version')
report<-list(platform=as.list(Sys.info()),R=R.version.string,tools=list())
for(name in names(flags)){path<-Sys.which(name);report$tools[[name]]<-if(!nzchar(path))list(path=NULL)else c(list(path=unname(path)),tryCatch(run(c(path,flags[[name]]),timeout=30,check=FALSE),error=function(e)list(error=conditionMessage(e))))}
report$r_packages<-setNames(lapply(c('DESeq2','apeglm','BiocParallel','pheatmap','jsonlite'),function(name)if(requireNamespace(name,quietly=TRUE))as.character(packageVersion(name))else 'MISSING'),c('DESeq2','apeglm','BiocParallel','pheatmap','jsonlite'))
cat(jsonlite::toJSON(report,auto_unbox=TRUE,pretty=TRUE,null='null'),'\n')
