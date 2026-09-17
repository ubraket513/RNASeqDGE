#!/usr/bin/env Rscript
# P0's independent origin oracle is retained in the broader native real-tool gate.
# This entry point now exercises both supported backends and both worker budgets.
source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])),'check_p4_alignment.R'))
