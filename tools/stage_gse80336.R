#!/usr/bin/env Rscript
source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])),'lib_helpers.R'))
metadata<-function(lines){
 samples<-list();current<-NULL; labels<-c('age (years)'='age',Sex='sex','postmortem interval (hours)'='pmi',rin='rin')
 for(line in lines){
  if(startsWith(line,'^SAMPLE = ')){current<-sub('^\\^SAMPLE = ','',line);assert(!current %in% names(samples),'duplicate GEO sample');samples[[current]]<-list(sample_id=current)}
  else if(!is.null(current)&&grepl(' = ',line,fixed=TRUE)){
   key<-sub(' = .*','',line);value<-substring(line,nchar(key)+4L)
   if(key=='!Sample_title'){assert(grepl('^(control|bipolar)-[0-9]+$',value),'unexpected study title');parts<-strsplit(value,'-')[[1]];samples[[current]]<-c(samples[[current]],list(condition=parts[1],title=value,count_column=paste0(if(parts[1]=='control')'C_'else'BD_',parts[2])))}
   else if(key=='!Sample_relation'&&grepl('SRX[0-9]+',value))samples[[current]]$experiment<-regmatches(value,regexpr('SRX[0-9]+',value))
   else if(key=='!Sample_characteristics_ch1'){label<-sub(': .*','',value);assert(label %in% names(labels),'unknown GEO covariate');samples[[current]][[labels[[label]]]]<-substring(value,nchar(label)+3L)}
  }
 }
 assert(length(samples)==36L&&any(grepl('PRJNA318642',lines))&&any(grepl('SRP073382',lines)),'wrong study or incomplete samples')
 for(sample in samples)assert(setequal(names(sample),c('sample_id','condition','title','count_column','experiment','age','sex','pmi','rin')),'incomplete GEO metadata');samples
}
main<-function(){
 args<-parse_cli(list(source=NULL,out=NULL),required=c('source','out'));source<-absolute(args$source);destination<-new_destination(args$out)
 info<-read_json(file.path(ROOT,'config/gse80336-sources.json'))
 for(name in names(info$files))assert(sha256(file.path(source,name))==info$files[[name]]$sha256,paste('source checksum mismatch:',name))
 con<-gzfile(file.path(source,'GSE80336_family.soft.gz'),'rt');samples<-metadata(readLines(con));close(con)
 ena<-table_read(file.path(source,'ena-runs.tsv'));experiments<-vapply(samples,function(s)s$experiment,'')
 assert(nrow(ena)==36L&&!anyDuplicated(ena$experiment_accession)&&setequal(ena$experiment_accession,experiments),'GEO/ENA mapping must be one-to-one')
 con<-gzfile(file.path(source,'GSE80336_Counts.txt.gz'),'rt');rows<-table_read(con);close(con)
 expected<-c('Ensembl_ID','GeneSymbol','Biotype','Chromosome');columns<-vapply(samples,function(s)s$count_column,'')
 assert(identical(names(rows)[1:4],expected)&&setequal(names(rows)[-(1:4)],columns),'count columns do not match cohort')
 rows<-rows[order(rows$Ensembl_ID,method='radix'),,drop=FALSE]
 assert(nrow(rows)==47886L&&!anyDuplicated(rows$Ensembl_ID),'unexpected gene identities')
 for(column in columns)assert(all(grepl('^[0-9]+$',rows[[column]]))&&all(as.numeric(rows[[column]])<=2147483647),'counts must be raw compatible integers')
 tmp<-new_directory('.gse80336-',dirname(destination));mapping<-list()
 for(sid in sort(names(samples))){s<-samples[[sid]];r<-ena[ena$experiment_accession==s$experiment,,drop=FALSE]
  assert(r$study_accession=='PRJNA318642'&&r$library_layout=='SINGLE'&&r$library_strategy=='RNA-Seq','unexpected ENA study/layout/strategy')
  mapping[[sid]]<-data.frame(sample_id=sid,title=s$title,count_column=s$count_column,experiment_id=s$experiment,run_id=r$run_accession,layout='single',strandedness='reverse',fastq_ftp=r$fastq_ftp,fastq_md5=r$fastq_md5,fastq_bytes=r$fastq_bytes,read_count=r$read_count)
 }
 table_write(do.call(rbind,mapping),file.path(tmp,'sample_run_mapping.tsv'))
 cohorts<-list(full36=sort(names(samples)),without_C28=sort(setdiff(names(samples),'GSM2124750')),subset6=c('GSM2124739','GSM2124742','GSM2124744','GSM2124760','GSM2124763','GSM2124766'))
 for(name in names(cohorts)){ids<-cohorts[[name]];out<-file.path(tmp,name);dir.create(out)
  header<-c('sample_id','condition','age','sex','pmi','rin');sample_rows<-do.call(rbind,lapply(ids,function(s)as.data.frame(samples[[s]][header],stringsAsFactors=FALSE)))
  assert(setequal(sample_rows$condition,c('control','bipolar')),'cohort needs both conditions');table_write(sample_rows,file.path(out,'samples.tsv'))
  count_table<-rows[c('Ensembl_ID',vapply(ids,function(s)samples[[s]]$count_column,''))];names(count_table)<-c('gene_id',ids);table_write(count_table,file.path(out,'counts.tsv'))
  annotation<-rows[expected];names(annotation)<-c('gene_id','gene_symbol','biotype','chromosome');table_write(annotation,file.path(out,'annotation.tsv'))
  kv_write(list(version='1',design_terms='condition',alpha='0.05',filter='zero_total',shrinkage='apeglm',type.age='numeric',type.sex='categorical',type.pmi='numeric',type.rin='numeric'),file.path(out,'analysis.tsv'))
  table_write(data.frame(contrast_id='bipolar_vs_control',factor='condition',numerator='bipolar',denominator='control'),file.path(out,'contrasts.tsv'))
 }
 info$genes<-nrow(rows);info$cohort_sizes<-lapply(cohorts,length);info$compressed_fastq_bytes<-sum(as.numeric(ena$fastq_bytes));info$design<-'~ condition; covariates retained as metadata, not included in report-model replication';info$source_count_scope<-'Published TopHat/HTSeq processed counts, including upstream symbol collapse; not new reference-gene featureCounts output';info$subset_selection<-'Three per condition, explicitly selected with similar ages and matched sex distribution, not first N runs';info$direction<-'bipolar/control; reverse of displayed report results contrast'
 write_json(info,file.path(tmp,'provenance.json'));assert(!exists_any(destination)&&file.rename(tmp,destination),'cannot publish stage');cat('Staged',nrow(rows),'genes and',length(samples),'mapped samples:',destination,'\n')
}
if(sys.nframe()==0L)main()
