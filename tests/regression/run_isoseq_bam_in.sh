#!/bin/bash
# * Run pbtranscript on the rat_bax1 data using the current 
#   build taking bam as input.

#NAME=run_isoseq_bam_in
INDIR=/pbi/dept/secondary/siv/testdata/pbtranscript-unittest/data/regression
DSDIR=/pbi/dept/secondary/siv/testdata/SA3-DS/rat/2450619/0005/Analysis_Results
FLNC=$INDIR/isoseq_flnc.fasta
NFL=$INDIR/isoseq_nfl.fasta
BAM=$DSDIR/m131018_081703_42161_c100585152550000001823088404281404_s1_p0.1.subreadset.xml
CCS=$DSDIR/m131018_081703_42161_c100585152550000001823088404281404_s1_p0.1.consensusreadset.xml
#BAM_OUTDIR=`pwd`/$NAME
BAM_OUTDIR=$1
OUTFA=$BAM_OUTDIR/cluster_out.fasta

echo Running pbtranscript cluster with --quiver
rm -rf $BAM_OUTDIR $OUTFA
pbtranscript cluster $FLNC $OUTFA -d $BAM_OUTDIR --bas_fofn $BAM --ccs_fofn $CCS --quiver --nfl_fa $NFL --use_finer_qv --blasr_nproc 8 --quiver_nproc 8 || (echo pbtranscript FAILED WITH CODE $? && exit 1)
