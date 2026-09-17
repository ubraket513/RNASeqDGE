#!/usr/bin/env Rscript
source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])),'../../tools/lib_helpers.R'))
f <- file.path(ROOT,'tests/fixtures/p0'); rows<-function(name)table_read(file.path(f,name))
references<-rows('references.tsv'); stopifnot(identical(references$role,c('genome','annotation')))
for(i in seq_len(nrow(references))) stopifnot(sha256(file.path(f,references$path[i]))==references$sha256[i])
runs<-rows('runs.tsv'); inputs<-rows('count_inputs.tsv')
stopifnot(identical(runs[c('run_id','sample_id')],inputs[c('run_id','sample_id')]))
for(s in unique(runs$sample_id)) stopifnot(nrow(unique(runs[runs$sample_id==s,c('layout','strandedness')]))==1L)
stopifnot(identical(nzchar(runs$fastq_2),runs$layout=='paired'))
reads<-list()
for(path in unique(c(runs$fastq_1,runs$fastq_2[nzchar(runs$fastq_2)]))) {
 lines<-readLines(file.path(f,path)); stopifnot(length(lines)%%4L==0L)
 for(i in seq(1L,length(lines),4L)) { stopifnot(startsWith(lines[i],'@'),lines[i+2L]=='+',nchar(lines[i+1L])==75L,nchar(lines[i+3L])==75L); reads[[substring(lines[i],2L)]]<-lines[i+1L] }
}
sequence<-paste(readLines(file.path(f,'reference.fa'))[-1],collapse=''); origins<-rows('read_origins.tsv')
stopifnot(setequal(names(reads),origins$read_id),substr(sequence,301,302)=='GT',substr(sequence,499,500)=='AG')
for(i in seq_len(nrow(origins))) { intervals<-strsplit(origins$intervals_1based_inclusive[i],',')[[1]];
 expected<-paste(vapply(intervals,function(x){p<-as.integer(strsplit(x,'-')[[1]]); substr(sequence,p[1],p[2])},''),collapse='')
 if(origins$orientation[i]=='-') expected<-reverse_complement(expected)
 stopifnot(reads[[origins$read_id[i]]]==expected)
}
samples<-rows('samples.tsv'); annotation<-rows('annotation.tsv'); genes<-annotation$gene_id
stopifnot(!anyDuplicated(samples$sample_id),!anyDuplicated(genes),!anyDuplicated(inputs$run_id),setequal(samples$sample_id,inputs$sample_id))
totals<-matrix(0,nrow=length(genes),ncol=nrow(samples),dimnames=list(genes,samples$sample_id))
for(i in seq_len(nrow(inputs))) { x<-rows(inputs$path[i]); stopifnot(nrow(x)==length(genes),setequal(x$gene_id,genes),all(grepl('^[0-9]+$',x$count))); totals[x$gene_id,inputs$sample_id[i]]<-totals[x$gene_id,inputs$sample_id[i]]+as.numeric(x$count) }
expected<-rows('expected/counts.tsv'); stopifnot(identical(names(expected),c('gene_id',samples$sample_id)),identical(expected$gene_id,sort(genes)))
stopifnot(all(as.matrix(sapply(expected[-1],as.numeric))==totals[expected$gene_id,,drop=FALSE]))
stopifnot(sum(samples$condition=='control')==2L,sum(samples$condition=='treated')==2L,
 identical(unname(unlist(rows('contrasts.tsv'))),c('treated_vs_control','condition','treated','control')),
 annotation$gene_symbol[1]==annotation$gene_symbol[2],annotation$gene_symbol[3]=='',
 identical(kv_read(file.path(f,'analysis.tsv')),list(version='1',design_terms='condition',alpha='0.05',filter='zero_total',shrinkage='apeglm')))
cat('PASS independent P0 sequence, arithmetic, annotation and contrast oracles\n')
