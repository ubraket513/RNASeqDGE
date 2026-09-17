#!/usr/bin/env Rscript
source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])),'../../tools/lib_helpers.R'))
args<-parse_cli(list('bin-dir'=file.path(ROOT,'.deps/alignment/bin'),rscript=file.path(ROOT,'.deps/p0/bin/Rscript')))
work<-new_directory('p5-real-r-',file.path(ROOT,'tests/output')); settings<-setNames(as.list(file.path(ROOT,'tests/fixtures/p0',paste0(c('samples','runs','references','analysis','contrasts'),'.tsv'))),c('samples','runs','references','analysis','contrasts'))
settings<-c(settings,list(genes=file.path(ROOT,'tests/fixtures/p2/genes.tsv'),bin_dir=absolute(args[['bin-dir']]),rscript=absolute(args$rscript),r_script=file.path(ROOT,'run_deg_analysis_offline.R'),tool_lock=file.path(ROOT,'config/alignment-linux-64.explicit.txt'),r_lock=file.path(ROOT,'config/p0-linux-64.explicit.txt'),run_dir=file.path(work,'run')))
config<-file.path(work,'workflow.tsv'); kv_write(settings,config)
result<-with_environment(c(clean_environment(),list(PATH='/usr/bin:/bin')),run(c(file.path(ROOT,'build/rnaseq'),'workflow-local',config),log=file.path(work,'workflow.log'),check=FALSE))
stopifnot(result$status!=0L); dest<-settings$run_dir
stopifnot(file.exists(file.path(dest,'merge/STAGE.tsv')),!file.exists(file.path(dest,'analysis/STAGE.tsv')))
x<-table_read(file.path(dest,'merge/counts.tsv')); rownames(x)<-x$gene_id
for(sample in c('treated_2','control_1','treated_1','control_2')) stopifnot(as.integer(x['gene_A.1',sample])==2L,as.integer(x['gene_B.2',sample])==if(sample=='control_1')0L else 1L,as.integer(x['gene_C',sample])==0L,as.integer(x['gene_zero',sample])==0L)
logs<-list.files(dest,'^R-.*[.]log$',full.names=TRUE); logs<-logs[!startsWith(basename(logs),'R-preflight-')]; stopifnot(length(logs)==1L)
stopifnot(any(grepl('dispersion|too few|fewer',readLines(logs))))
write_json(list(alignment='HISAT2 2.2.1 native, samtools 1.17, featureCounts 2.0.6',run_count=5,sample_count=4,counts='independent P0 SE and PE expectations passed',R='real tiny-data failure retained; no analysis success manifest',exit_status=result$status),file.path(work,'verified.json'))
cat('PASS real P5 alignment/count merge and honest R failure:',work,'\n')
