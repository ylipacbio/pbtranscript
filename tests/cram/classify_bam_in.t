#Test pbtranscript classify, bam | dataset input.

  $ . $TESTDIR/setup.sh

  $ D=$SIVDATDIR/rat_downsampled
  $ MOVIE=m131018_081703_42161_c100585152550000001823088404281404_s1_p0
  $ CCSBAM=$D/$MOVIE.1.ccs.bam
  $ ONAME=output
  $ BAM_OD=$OUTDIR/test_classify_bam
  $ BAM_OFA=$BAM_OD/$ONAME.fasta
  $ rm -rf $BAM_OFA $BAM_OD 

  $ pbtranscript classify $CCSBAM $BAM_OFA -d $BAM_OD && echo $?
  0

  $ CCSDS=$D/$MOVIE.1.consensusreadset.xml
  $ DS_OD=$OUTDIR/test_classify_dataset
  $ DS_OFA=$DS_OD/$ONAME.fasta
  $ rm -rf $DS_OFA $DS_OD 

  $ pbtranscript classify $CCSDS $DS_OFA -d $DS_OD && echo $?
  0

  $ diff $BAM_OFA $DS_OFA

  $ diff $BAM_OD/$ONAME.classify_summary.txt $DS_OD/$ONAME.classify_summary.txt

  $ diff $BAM_OD/$ONAME.primer_info.csv $DS_OD/$ONAME.primer_info.csv

