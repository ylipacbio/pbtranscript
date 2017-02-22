Test chain_samples.py

Set up
  $ . $TESTDIR/setup.sh

  $ cfg_fn=$SIVDATDIR/chain_samples/sample.config
  $ chain_samples.py $cfg_fn norm_nfl --max_fuzzy_junction=5 1>/dev/null 2>/dev/null && echo $?
  0

  $ std_gff_fn=$SIVSTDDIR/chain_samples/all_samples.chained.gff
  $ std_count_fn=$SIVSTDDIR/chain_samples/all_samples.chained_count.txt
  $ std_ids_fn=$SIVSTDDIR/chain_samples/all_samples.chained_ids.txt
  $ diff all_samples.chained.gff $std_gff_fn
  $ diff all_samples.chained_count.txt $std_count_fn
  $ diff all_samples.chained_ids.txt $std_ids_fn
