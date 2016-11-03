Test ice_partial

Set up
  $ . $TESTDIR/setup.sh

  $ subreads=/pbi/dept/secondary/siv/testdata/SA3-Sequel/rc0/315/3150353/r54086_20160831_010819/4_D01_tiny/tiny_flea_isoseq.subreadset.xml
  $ src_tasks_dir=$SIVDATDIR/test_ice_entries/tasks
  $ nfl=$src_tasks_dir/pbcoretools.tasks.gather_contigset-3/file.contigset.xml

  $ out_dir=$OUTDIR/test_ice_partial_1

  $ rm -rf $out_dir && mkdir -p $out_dir
  $ cp -r $src_tasks_dir/pbtranscript.tasks.separate_flnc-0/3to4kb_part0/cluster_out/* $out_dir/
  $ consensus_fa=$out_dir/output/final.consensus.fasta
  $ out_pickle=$out_dir/out_ice_partial.pickle

test ice_partial all
  $ ice_partial.py --verbose all $nfl $consensus_fa $out_pickle --root_dir $out_dir 1>/dev/null 2>/dev/null && echo $?
  0
  $ ls $out_pickle 1>/dev/null 2>/dev/null && echo $?
  0

test ice_partial split|i|merge
  $ rm -rf $out_dir && mkdir -p $out_dir
  $ cp -r $src_tasks_dir/pbtranscript.tasks.separate_flnc-0/3to4kb_part0/cluster_out/* $out_dir/

  $ out_pickle=$out_dir/output/map_noFL/nfl.all.partial_uc.pickle
  $ ice_partial.py --verbose split $out_dir $nfl 3  1>/dev/null 2>/dev/null && echo $?
  0
  $ ice_partial.py --verbose i $out_dir   0   1>/dev/null 2>/dev/null && echo $?
  0
  $ ice_partial.py --verbose i $out_dir   1   1>/dev/null 2>/dev/null && echo $?
  0
  $ ice_partial.py --verbose i $out_dir   2   1>/dev/null 2>/dev/null && echo $?
  0
  $ ice_partial.py --verbose merge $out_dir 3 1>/dev/null 2>/dev/null && echo $?
  0
  $ ls $out_pickle 1>/dev/null 2>/dev/null && echo $?
  0
