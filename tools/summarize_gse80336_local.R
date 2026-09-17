#!/usr/bin/env Rscript
# Summarize completed reduced workflows without claiming whole-genome parity.
source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])),'lib_helpers.R'))
correlation<-function(x,y)if(length(x)<2L||length(unique(x))<2L||length(unique(y))<2L)NULL else unname(cor(x,y))
main<-function(){
 args<-parse_cli(list(directory=NULL,'expected-repeats'='3',out=NULL),required='directory',positional='directory');root<-absolute(args$directory);records<-read_json(file.path(root,'runs.json'))
 wanted<-if(is.null(records$requested_replicates))as.integer(args[['expected-repeats']])else records$requested_replicates
 assert(wanted>=1L,'invalid replicate count');expected<-as.vector(outer(c('hisat2','star'),seq_len(wanted),paste,sep='-'))
 observed<-vapply(records$runs,function(r)paste(r$backend,r[['repeat']],sep='-'),'')
 assert(setequal(observed,expected)&&length(observed)==length(expected),'incomplete or duplicate backend/repetition pairs')
 records$requested_replicates<-wanted;records$complete<-TRUE
 for(i in seq_along(records$runs)){
  r<-records$runs[[i]];assert(r$status==0L,'refusing comparison with failed workflows');name<-paste(r$backend,r[['repeat']],sep='-');path<-file.path(root,name)
  files<-list.files(path,recursive=TRUE,full.names=TRUE,all.files=TRUE);r$run_directory_bytes_excluding_shared_index<-sum(file.info(files)$size,na.rm=TRUE);r$resource_usage_by_attempt<-list()
  timing<-if(is.null(r$timing_files))paste0(name,'.time.txt')else unlist(r$timing_files)
  for(filename in timing){lines<-trimws(readLines(file.path(root,filename)));lines<-lines[grepl(': ',lines,fixed=TRUE)];keys<-sub(': .*','',lines);values<-substring(lines,nchar(keys)+3L);r$resource_usage_by_attempt[[filename]]<-as.list(setNames(values,keys))}
  r$assignment<-list();for(summary in sort(Sys.glob(file.path(path,'align/*/counts.txt.summary')))){x<-table_read(summary);r$assignment[[basename(dirname(summary))]]<-as.list(setNames(as.numeric(x[[2]]),x$Status))}
  r$analysis<-kv_read(file.path(path,'analysis/contrasts/bipolar_vs_control/summary.tsv'));records$runs[[i]]<-r
 }
 data<-setNames(lapply(c('hisat2','star'),function(b)table_read(file.path(root,paste0(b,'-1/merge/counts.tsv')))),c('hisat2','star'));h<-data$hisat2;s<-data$star
 assert(identical(h$gene_id,s$gene_id)&&identical(names(h),names(s)),'incomparable gene/sample identity or order');ids<-names(h)[-1]
 hm<-sapply(h[-1],as.numeric);sm<-sapply(s[-1],as.numeric);hv<-as.vector(t(hm));sv<-as.vector(t(sm))
 records$counts_comparison<-list(genes=nrow(h),samples=length(ids),differing_cells=sum(hv!=sv),total_absolute_difference=sum(abs(hv-sv)),hisat2_total=sum(hv),star_total=sum(sv),pearson=correlation(hv,sv))
 for(backend in c('hisat2','star')){first<-Filter(function(r)r$backend==backend&&r[['repeat']]==1L,records$runs)[[1]];assert(sum(vapply(first$assignment,function(x)x$Assigned,0))==if(backend=='hisat2')sum(hv)else sum(sv),'merged counts disagree with Assigned totals')}
 differences<-rowSums(abs(hm-sm));order<-order(-differences,h$gene_id,decreasing=c(FALSE,TRUE),method='radix')[seq_len(min(10L,nrow(h)))]
 records$counts_comparison$largest_gene_differences<-lapply(order,function(i)list(gene_id=h$gene_id[i],summed_absolute_difference=differences[i]))
 config<-kv_read(file.path(root,'hisat2-1.tsv'));source<-read_json(file.path(dirname(config$runs),'provenance.json'));records$mapping_primary_records<-list()
 for(backend in c('hisat2','star')){
  records$mapping_primary_records[[backend]]<-list()
  for(read in source$reads){run_id<-read$run_id;total<-read$prefix_records;metrics<-list(input_reads=total,mapped=0L,unique=0L,multimapped=0L,spliced=0L)
   # Capture to a temporary regular file so subprocess exit status is verified
   # before streaming SAM records in bounded chunks.
   sam<-tempfile('mapping-',tmpdir=root);log<-tempfile('mapping-log-',tmpdir=root)
   run(c(file.path(config$bin_dir,'samtools'),'view','-o',sam,file.path(root,paste0(backend,'-1/align'),run_id,'aligned.bam')),log=log)
   con<-file(sam,'r');tryCatch(repeat{lines<-readLines(con,n=10000L);if(!length(lines))break;for(line in lines){fields<-strsplit(line,'\t')[[1]];flag<-as.integer(fields[2]);if(bitwAnd(flag,4L+256L+2048L)!=0L)next
    metrics$mapped<-metrics$mapped+1L;nh<-fields[grepl('^NH:i:',fields)];assert(length(nh)==1L,'missing/duplicate NH tag');nh<-as.integer(sub('^NH:i:','',nh));assert(!is.na(nh),'invalid NH tag');key<-if(nh==1L)'unique'else'multimapped';metrics[[key]]<-metrics[[key]]+1L;metrics$spliced<-metrics$spliced+as.integer(grepl('N',fields[6],fixed=TRUE))
   }},finally=close(con));unlink(c(sam,log));assert(metrics$mapped<=total,'invalid BAM mapping total');metrics$mapping_rate<-metrics$mapped/total;records$mapping_primary_records[[backend]][[run_id]]<-metrics
  }
 }
 result<-setNames(lapply(c('hisat2','star'),function(b){x<-table_read(file.path(root,paste0(b,'-1/analysis/contrasts/bipolar_vs_control/results.tsv')));rownames(x)<-x$gene_id;x}),c('hisat2','star'))
 genes<-sort(intersect(result$hisat2$gene_id,result$star$gene_id));finite<-function(x)is.finite(suppressWarnings(as.numeric(x)))
 genes<-genes[finite(result$hisat2[genes,'log2FoldChange'])&finite(result$star[genes,'log2FoldChange'])]
 lfc<-lapply(result,function(x)as.numeric(x[genes,'log2FoldChange']));deg<-lapply(result,function(x){value<-suppressWarnings(as.numeric(x$padj));sort(x$gene_id[is.finite(value)&value<.05])});union<-union(deg$hisat2,deg$star)
 records$statistics_comparison<-list(common_finite_lfc_genes=length(genes),lfc_pearson=correlation(lfc$hisat2,lfc$star),median_absolute_lfc_difference=if(length(genes))median(abs(lfc$hisat2-lfc$star))else NULL,lfc_sign_disagreements=sum((lfc$hisat2>0)!=(lfc$star>0)&lfc$hisat2!=0&lfc$star!=0),deg=lapply(deg,as.list),deg_jaccard=if(length(union))length(intersect(deg$hisat2,deg$star))/length(union)else NULL)
 records$repeat_equality<-list()
 for(backend in c('hisat2','star')){repeats<-vapply(Filter(function(r)r$backend==backend,records$runs),function(r)r[['repeat']],0);records$repeat_equality[[backend]]<-list()
  for(artifact in c('merge/counts.tsv','analysis/contrasts/bipolar_vs_control/results.tsv','analysis/contrasts/bipolar_vs_control/shrunken.tsv')){hashes<-vapply(repeats,function(n)sha256(file.path(root,paste0(backend,'-',n),artifact)),'');records$repeat_equality[[backend]][[artifact]]<-list(repeats=length(repeats),identical=if(length(repeats)>1L)length(unique(hashes))==1L else NULL,assessment=if(length(repeats)>1L)'checked'else'not assessed: single completed run')}
 }
 records$limitations<-as.list(c('Restricted chromosome reference can misassign reads from omitted chromosomes.','Prefix sampling is not random and preserves sequencing-order bias.','Low-depth DEG/LFC differences are diagnostic, not biological findings.','OS cache uncontrolled; first run includes indexing, later runs reuse index.','Timing includes provenance hashes and R; these historical runs do not have independent per-stage resource profiling.'))
 target<-if(is.null(args$out))file.path(root,'comparison.json')else args$out;write_json(records,target);cat(jsonlite::toJSON(records[c('counts_comparison','statistics_comparison','repeat_equality')],auto_unbox=TRUE,pretty=TRUE,null='null'),'\n')
}
if(sys.nframe()==0L)main()
