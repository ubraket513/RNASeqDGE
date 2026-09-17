#!/usr/bin/env Rscript
source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])),'lib_helpers.R'))
generate_p0 <- function(out) {
 dir.create(out,recursive=TRUE,showWarnings=FALSE); dir.create(file.path(out,'reads'),showWarnings=FALSE)
 tmp<-tempfile(); on.exit(unlink(tmp))
 bytes<-unlist(lapply(0:99,function(i) {writeBin(charToRaw(paste0('RNASeqDGE-P0-',i)),tmp); h<-sha256(tmp); strtoi(substring(h,seq(1,63,2),seq(2,64,2)),16L)}))
 sequence<-paste(c('A','C','G','T')[bytes%%4L+1L][1:3000],collapse='')
 substr(sequence,301,302)<-'GT'; substr(sequence,499,500)<-'AG'
 writeLines(c('>chrSynthetic',substring(sequence,seq(1,3000,60),seq(60,3000,60))),file.path(out,'reference.fa'))
 genes<-c('gene_A.1','gene_A.1','gene_B.2','gene_C','gene_zero'); strands<-c('+','+','-','+','+')
 starts<-c(101,501,1001,1101,2001); ends<-c(300,700,1300,1400,2300)
 writeLines(sprintf('chrSynthetic\tP0\texon\t%d\t%d\t.\t%s\t.\tgene_id "%s"; transcript_id "%s.tx";',starts,ends,strands,genes,genes),file.path(out,'reference.gtf'))
 fastq<-function(name,ids,seqs) writeLines(unlist(Map(function(id,s)c(paste0('@',id),s,'+',strrep('I',nchar(s))),ids,seqs)),file.path(out,'reads',name))
 fastq('single.fastq',c('exon_plus','junction_plus','exon_minus','ambiguous_overlap'),c(substr(sequence,151,225),paste0(substr(sequence,264,300),substr(sequence,501,538)),reverse_complement(substr(sequence,1001,1075)),substr(sequence,1151,1225)))
 fastq('paired_1.fastq','fragment_A/1',substr(sequence,151,225)); fastq('paired_2.fastq','fragment_A/2',reverse_complement(substr(sequence,576,650)))
 table_write(data.frame(read_id=c('exon_plus','junction_plus','exon_minus','ambiguous_overlap','fragment_A/1','fragment_A/2'),
 intervals_1based_inclusive=c('151-225','264-300,501-538','1001-1075','1151-1225','151-225','576-650'),orientation=c('+','+','-','+','+','-'),
 case=c('exonic_gene_A.1','spliced_gene_A.1','exonic_gene_B.2','gene_B.2_and_gene_C_unstranded','paired_gene_A.1','paired_gene_A.1')),file.path(out,'read_origins.tsv'))
 paired<-c(FALSE,TRUE,TRUE,FALSE,FALSE)
 table_write(data.frame(run_id=c('run_t2','run_c1a','run_c1b','run_t1','run_c2'),sample_id=c('treated_2','control_1','control_1','treated_1','control_2'),fastq_1=ifelse(paired,'reads/paired_1.fastq','reads/single.fastq'),fastq_2=ifelse(paired,'reads/paired_2.fastq',''),layout=ifelse(paired,'paired','single'),strandedness='unstranded'),file.path(out,'runs.tsv'))
 names<-c('reference.fa','reference.gtf'); table_write(data.frame(role=c('genome','annotation'),path=names,sha256=vapply(file.path(out,names),sha256,''),source='synthetic:tools/generate_p0_sequences.R',release='p0-v1'),file.path(out,'references.tsv'))
}
if(sys.nframe()==0L) {args<-parse_cli(list(output=file.path(ROOT,'tests/fixtures/p0'))); generate_p0(args$output)}
