Test ice_daligner

Set up
  $ . $TESTDIR/setup.sh

  $ ice_daligner_out_dir=$OUTDIR/test_ice_daligner
  $ query=$DATDIR/test_daligner_against_ref/test_daligner_query.fasta 
  $ target=$DATDIR/test_daligner_against_ref/test_daligner_target.fasta 

  $ oquery=$ice_daligner_out_dir/test_daligner_query.fasta 
  $ otarget=$ice_daligner_out_dir/test_daligner_target.fasta 

  $ rm -rf $ice_daligner_out_dir && mkdir -p $ice_daligner_out_dir
  $ cp $query $oquery && cp $target $otarget

  $ ice_daligner.py $oquery $otarget $ice_daligner_out_dir && echo $?
  0

  $ ls $ice_daligner_out_dir/*.out | wc -l
  8
