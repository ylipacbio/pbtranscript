#Test pbtranscript cluster, bam input.

  $ . $TESTDIR/setup.sh

  $ D=$SIVDATDIR/sa3
  $ BAS=$D/bam.subreadset.xml
  $ CCS=$D/tinyccsbam.consensusreadset.xml # use a smaller dataset other than $D/ccsbam.consensusset.xml
  $ NFL=$D/nfl.contigset.xml
  $ FLNC=$D/flnc.contigset.xml

  $ CLASSIFYOUTDIR=$OUTDIR/test_classify
  $ OUTFA=$OUTDIR/classify_out_1.fasta
  $ OUTXML=$OUTDIR/classify_out_1.contigset.xml
  $ rm -rf $OUTFA $OUTXML $OUTCSV $OUTSUMMARY $CLASSIFYOUTDIR
  $ pbtranscript classify $CCS $OUTXML -d $CLASSIFYOUTDIR --cpus 8
