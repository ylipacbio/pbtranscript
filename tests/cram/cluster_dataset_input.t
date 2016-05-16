#Test pbtranscript cluster, bam input.

  $ . $TESTDIR/setup.sh

  $ D=$SIVDATDIR/sa3
  $ BAS=$D/bam.subreadset.xml
  $ CCS=$D/ccsbam.consensusreadset.xml
  $ NFL=$D/nfl.contigset.xml
  $ FLNC=$D/flnc.contigset.xml

  $ OFA=$OUTDIR/cluster_bam_out.contigset.xml
  $ OD=$OUTDIR/test_cluster_dataset

  $ rm -rf $OFA $OD && mkdir -p $OD
  $ pbtranscript cluster $FLNC $OFA -d $OD --bas_fofn $BAS --ccs_fofn $CCS 
  $ wc $OFA > /dev/null

# Test pbtranscript cluster, dataset input, using finer qvs, quiver.
  $ rm -rf $OFA $OD && mkdir -p $OD
  $ pbtranscript cluster $FLNC $OFA -d $OD --bas_fofn $BAS --ccs_fofn $CCS  --use_finer_qv --quiver --nfl $NFL
  $ wc $OFA > /dev/null
