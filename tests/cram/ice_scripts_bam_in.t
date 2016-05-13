Set up
  $ . $TESTDIR/setup.sh

# Run Ice iterative clustering without Quiver polishing
  $ BASFOFN=$SIVDATDIR/bam.fofn
  $ CCSFOFN=$SIVDATDIR/ccsbam.fofn
  $ NFLFA=$SIVDATDIR/nfl.fasta
  $ FLNCFA=$SIVDATDIR/flnc.fasta


# Run a cluster job without quiver
  $ XDIR=$OUTDIR/test_ice_scripts_bam_in
  $ XFA=$ODIR.fasta
  $ rm -rf $XDIR $XFA && mkdir -p $XDIR

  $ pbtranscript --quiet cluster $FLNCFA $XFA -d $XDIR --nfl_fa $NFLFA --bas_fofn $BASFOFN --ccs_fofn $CCSFOFN

# Test ice_partial.py all without sge
  $ ODIR=$OUTDIR/test_ice_partial_all_bam_in
  $ REFFA=$ODIR/output/final.consensus.fasta
  $ OUTPICKLE=$ODIR/nfl_all.pickle
  $ rm -rf $ODIR $OFA && cp -r $XDIR $ODIR

  $ ice_partial.py all $NFLFA $REFFA $OUTPICKLE --ccs_fofn $CCSFOFN --root_dir $ODIR

# Test ice_partial.py -> split -> i -> merge;
  $ ROOTDIR=$OUTDIR/test_ice_partial_split_merge_bam_in
  $ rm -rf $ROOTDIR && cp -r $XDIR $ROOTDIR

# Test ice_partial.py split
  $ N=3
  $ ice_partial.py split $ROOTDIR $NFLFA $N

# Test ice_partial.py i
  $ ice_partial.py i $ROOTDIR 0 --ccs_fofn $CCSFOFN
  $ ice_partial.py i $ROOTDIR 1 --ccs_fofn $CCSFOFN
  $ ice_partial.py i $ROOTDIR 2 --ccs_fofn $CCSFOFN

# Test ice_partial.py merge
  $ ice_partial.py merge $ROOTDIR $N

# Test ice_quiver.py all  on $ODIR
  $ ice_quiver.py all $ODIR --bas_fofn $BASFOFN 

# Test ice_quiver.py -> i -> merge -> postprocess on $ROOTDIR
# Test ice_quiver.py i
  $ ice_quiver.py i $ROOTDIR 2 0 --bas_fofn $BASFOFN
  $ ice_quiver.py i $ROOTDIR 2 1 --bas_fofn $BASFOFN

# Test ice_quiver.py merge
  $ ice_quiver.py merge $ROOTDIR 2

# Test ice_quiver.py postprocess
  $ ice_quiver.py postprocess $ROOTDIR

# Test whether 'ice_partial.py all' -> 'ice_quiver.py all' can produce identical
# results with 'ice_partial.py split,i,merge' -> 'ice_quiver.py i,merge,postprocess'.
  $ diff $ODIR/all_quivered_hq.100_30_0.99.fasta $ROOTDIR/all_quivered_hq.100_30_0.99.fasta
  $ grep -c '>' $ROOTDIR/all_quivered_hq.100_30_0.99.fasta | grep -c '0'
  0
  [1]
