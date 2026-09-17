#!/usr/bin/env Rscript
source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])),'../../tools/stage_gse80336_local.R'))
check_copy<-function(text,n=1L){incoming<-textConnection(text);result<-character();outgoing<-textConnection('result','w',local=TRUE);on.exit({close(incoming);close(outgoing)});copy_fastq_prefix(incoming,outgoing,n);result}
stopifnot(identical(check_copy('@a\nAC\n+\n!!\n@b\nGG\n+\n!!'),c('@a','AC','+','!!')))
for(data in c('@a\nAC\n+\n!','@a\nAC\n+','bad\nAC\n+\n!!'))stopifnot(inherits(tryCatch({check_copy(data);NULL},error=identity),'error'))
cat('PASS FASTQ prefix stops at complete record and rejects delimiters, truncation, quality lengths\n')
# Every-fifth selections must never be relabeled as contiguous archive prefixes.
original<-list(sampling='first N complete FASTQ records in archive order; nonrandom prefix sampling may be biased',reads=list(list(prefix_records=50000)))
validate_read_stage(original)
thinned<-original;thinned$sampling<-'every fifth record of first50k archive prefix; nonrandom'
stopifnot(inherits(tryCatch({validate_read_stage(thinned);NULL},error=identity),'error'))
mislabelled<-original;mislabelled$reads[[1]]$selection<-'zero-based records 0,5,...,49995 of 50k archive prefix'
stopifnot(inherits(tryCatch({validate_read_stage(mislabelled);NULL},error=identity),'error'))
