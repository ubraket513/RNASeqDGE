#!/usr/bin/env Rscript
source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])),'../../tools/lib_helpers.R'))
args<-parse_cli(list(exe=file.path(ROOT,'build/rnaseq'),'bin-dir'=file.path(ROOT,'.deps/alignment/bin'),output=NULL))
f<-file.path(ROOT,'tests/fixtures/p0'); out<-if(is.null(args$output))new_directory('p4-real-r-',file.path(ROOT,'tests/output')) else args$output
dir.create(out,recursive=TRUE,showWarnings=FALSE); exe<-absolute(args$exe); tools<-absolute(args[['bin-dir']])
common<-c('--bin-dir',tools,'--fasta',file.path(f,'reference.fa'),'--gtf',file.path(f,'reference.gtf')); evidence<-list()
expected<-list(single=c('gene_A.1'=2,'gene_B.2'=1,gene_C=0,gene_zero=0),paired=c('gene_A.1'=1,'gene_B.2'=0,gene_C=0,gene_zero=0))
for(backend in c('hisat2','star')) for(threads in 1:2) {
 index<-file.path(out,paste0(backend,'-',threads,'-index'))
 run(c(exe,'index','--backend',backend,common,'--threads',threads,'--output',index,if(backend=='star')c('--star-sa-bases',3,'--star-chr-bits',10)),timeout=180)
 for(layout in c('single','paired')) {
  dest<-file.path(out,paste0(backend,'-',threads,'-',layout)); reads<-c('--reads1',file.path(f,'reads',if(layout=='single')'single.fastq' else 'paired_1.fastq'))
  if(layout=='paired')reads<-c(reads,'--reads2',file.path(f,'reads/paired_2.fastq'))
  run(c(exe,'align-count','--backend',backend,common,'--threads',threads,'--index',index,'--layout',layout,'--strandedness','unstranded',reads,'--output',dest),timeout=180)
  counts<-read.delim(file.path(dest,'counts.txt'),comment.char='#',check.names=FALSE); observed<-setNames(counts[[ncol(counts)]],counts[[1]])
  stopifnot(setequal(names(observed),names(expected[[layout]])),all(observed[names(expected[[layout]])]==expected[[layout]]))
  summary<-table_read(file.path(dest,'counts.txt.summary')); values<-setNames(as.integer(summary[[2]]),summary[[1]])
  stopifnot(values[['Assigned']]==if(layout=='single')3 else 1,values[['Unassigned_Ambiguity']]==if(layout=='single')1 else 0,sum(values)==if(layout=='single')4 else 1)
  records<-run(c(file.path(tools,'samtools'),'view',file.path(dest,'aligned.bam')))$output
  if(layout=='single') {junction<-records[startsWith(records,'junction_plus\t')]; stopifnot(length(junction)==1L,strsplit(junction,'\t')[[1]][6]=='37M200N38M')}
  stopifnot(file.exists(file.path(dest,'COMPLETE')))
  evidence[[length(evidence)+1L]]<-list(backend=backend,threads=threads,layout=layout,counts=as.list(observed),summary=as.list(values))
  cat('PASS',backend,'threads=',threads,layout,'counts, summary, junction, BAM\n')
 }
}
write_json(evidence,file.path(out,'verified.json')); cat('Evidence:',out,'\n')
