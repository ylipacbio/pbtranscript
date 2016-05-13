Test if --tmp_dir works well in pbtranscript cluster, ice_quiver and ice_polish.

Set up
  $ . $TESTDIR/setup.sh

  $ BASFOFN=$SIVDATDIR/bam.fofn
  $ CCSFOFN=$SIVDATDIR/ccsbam.fofn
  $ NFLFA=$SIVDATDIR/nfl.fasta
  $ FLNCFA=$SIVDATDIR/flnc.fasta
  $ PBTRANSCRIPT=pbtranscript
  $ ICE_QUIVER_PY=ice_quiver.py
  $ ICE_PARTIAL_PY=ice_partial.py
  $ ICE_POLISH_PY=ice_polish.py


##################################################
# Run a cluster job with --tmp_dir without quiver
  $ mkdir -p $OUTDIR/test_tmp_dir
  $ XDIR=$OUTDIR/test_tmp_dir/cluster_no_quiver
  $ XFA=$OUTDIR/output.fasta
  $ TDIR=$OUTDIR/test_tmp_dir/cluster_no_quiver_tmp
  $ rm -rf $XDIR $XFA $TDIR && mkdir -p $XDIR

  $ $PBTRANSCRIPT --quiet cluster $FLNCFA $XFA -d $XDIR --nfl_fa $NFLFA --bas_fofn $BASFOFN --ccs_fofn $CCSFOFN --tmp_dir=$TDIR && echo $?
  0

#Whether tmp_dir exists and has the right cluster dir
  $ ls $TDIR/0/c0/g_consensus.fasta > /dev/null && echo $?
  0

  $ $ICE_PARTIAL_PY split $XDIR $NFLFA 1
  $ $ICE_PARTIAL_PY i $XDIR 0
  $ $ICE_PARTIAL_PY merge $XDIR 1

# ice_quiver, no --tmp_dir, 1 chunk
  $ $ICE_QUIVER_PY i $XDIR 1 0 --bas_fofn $BASFOFN --tmp_dir=$TDIR && echo $?
  0

# ice_quiver, --tmp_dir, 1 chunk
  $ $ICE_QUIVER_PY i $XDIR 1 0 --bas_fofn $BASFOFN && echo $?
  0

# ice_quiver, --tmp_dir, 2 chunk
  $ $ICE_QUIVER_PY i $XDIR 2 0 --bas_fofn $BASFOFN && echo $?
  0
  $ $ICE_QUIVER_PY i $XDIR 2 1 --bas_fofn $BASFOFN && echo $?
  0
# merge and postprocess
  $ $ICE_QUIVER_PY merge $XDIR 1
  $ $ICE_QUIVER_PY postprocess $XDIR
  $ grep ">" $XDIR/all_quivered_hq.100_30_0.99.fasta | wc -l 
  7

##################################################
# Run a cluster job with --tmp_dir quiver
  $ XDIR=$OUTDIR/test_tmp_dir/cluster_quiver
  $ XFA=$OUTDIR/output.fasta
  $ TDIR=$OUTDIR/test_tmp_dir/cluster_quiver_tmp
  $ rm -rf $XDIR $XFA $TDIR && mkdir -p $XDIR

  $ $PBTRANSCRIPT --quiet cluster $FLNCFA $XFA -d $XDIR --nfl_fa $NFLFA --bas_fofn $BASFOFN --ccs_fofn $CCSFOFN --tmp_dir=$TDIR --quiver && echo $?
  0

  $ grep ">" $XDIR/all_quivered_hq.100_30_0.99.fasta | wc -l 
  7

#Whether tmp_dir exists and has the right cluster dir
  $ ls $TDIR/0/c0/g_consensus.fasta > /dev/null && echo $?
  0


##################################################
#Run a cluster job no --tmp_dir, no --quiver
  $ XDIR=$OUTDIR/test_tmp_dir/polish
  $ XFA=$OUTDIR/output.fasta
  $ TDIR=$OUTDIR/test_tmp_dir/polish_tmp
  $ rm -rf $XDIR $XFA $TDIR && mkdir -p $XDIR

  $ $PBTRANSCRIPT --quiet cluster $FLNCFA $XFA -d $XDIR --nfl_fa $NFLFA --bas_fofn $BASFOFN --ccs_fofn $CCSFOFN && echo $?
  0

#Remove root_dir/tmp
  $ rm -rf $XDIR/tmp
  $ $ICE_POLISH_PY $XDIR $NFLFA --tmp_dir $TDIR --bas_fofn $BASFOFN && echo $?
  0
  $ grep ">" $XDIR/all_quivered_hq.100_30_0.99.fasta | wc -l 
  7
