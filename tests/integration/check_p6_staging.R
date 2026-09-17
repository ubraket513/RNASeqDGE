#!/usr/bin/env Rscript
source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])),'../../tools/stage_gse80336_local.R'))
check_copy<-function(text,n=1L){incoming<-textConnection(text);result<-character();outgoing<-textConnection('result','w',local=TRUE);on.exit({close(incoming);close(outgoing)});copy_fastq_prefix(incoming,outgoing,n);result}
stopifnot(identical(check_copy('@a\nAC\n+\n!!\n@b\nGG\n+\n!!'),c('@a','AC','+','!!')))
for(data in c('@a\nAC\n+\n!','@a\nAC\n+','bad\nAC\n+\n!!'))stopifnot(inherits(tryCatch({check_copy(data);NULL},error=identity),'error'))
cat('PASS FASTQ prefix stops at complete record and rejects delimiters, truncation, quality lengths\n')
