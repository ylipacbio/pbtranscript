Test ice_daligner

Set up
  $ . $TESTDIR/setup.sh

  $ out_dir=$OUTDIR/test_ice_daligner
  $ query=$DATDIR/test_daligner_against_ref/test_daligner_query.fasta 
  $ target=$DATDIR/test_daligner_against_ref/test_daligner_target.fasta 

  $ oquery=$out_dir/test_daligner_query.fasta 
  $ otarget=$out_dir/test_daligner_target.fasta 

  $ rm -rf $out_dir && mkdir -p $out_dir
  $ cp $query $oquery && cp $target $otarget

  $ ice_daligner.py $oquery $otarget $out_dir && echo $?
  0

  $ ls $out_dir/*.out | wc -l
  8
