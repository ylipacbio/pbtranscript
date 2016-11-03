Test ice_quiver

Set up
  $ . $TESTDIR/setup.sh

  $ subreads=/pbi/dept/secondary/siv/testdata/SA3-Sequel/rc0/315/3150353/r54086_20160831_010819/4_D01_tiny/tiny_flea_isoseq.subreadset.xml
  $ src_tasks_dir=$SIVDATDIR/test_ice_entries/tasks
  $ nfl=$src_tasks_dir/pbcoretools.tasks.gather_contigset-3/file.contigset.xml

  $ out_dir=$OUTDIR/test_ice_quiver_1

  $ rm -rf $out_dir && mkdir -p $out_dir
  $ cp -r $src_tasks_dir/pbtranscript.tasks.separate_flnc-0/3to4kb_part0/cluster_out/* $out_dir/
  $ ice_quiver.py --verbose all $out_dir --bas_fofn=$subreads 1>/dev/null 2>/dev/null & echo $?
  0
  $ ls $out_dir/all_quivered_hq.100_30_0.99.fasta $out_dir/all_quivered_hq.100_30_0.99.fastq $out_dir/all_quivered_lq.fasta $out_dir/all_quivered_lq.fastq 2>&1 > /dev/null && echo $?
  0

  $ rm -rf $out_dir && mkdir -p $out_dir
  $ cp -r $src_tasks_dir/pbtranscript.tasks.separate_flnc-0/3to4kb_part0/cluster_out/* $out_dir/
  $ ice_quiver.py --verbose i   $out_dir 3 0 --bas_fofn=$subreads 1>/dev/null 2>/dev/null & echo $?
  0
  $ ice_quiver.py --verbose i   $out_dir 3 1 --bas_fofn=$subreads 1>/dev/null 2>/dev/null & echo $?
  0
  $ ice_quiver.py --verbose i   $out_dir 3 2 --bas_fofn=$subreads 1>/dev/null 2>/dev/null & echo $?
  0
  $ ice_quiver.py --verbose merge   $out_dir 3 1>/dev/null 2>/dev/null & echo $?
  0
  $ ice_quiver.py --verbose postprocess   $out_dir 1>/dev/null 2>/dev/null & echo $?
  0
  $ ls $out_dir/all_quivered_hq.100_30_0.99.fasta $out_dir/all_quivered_hq.100_30_0.99.fastq $out_dir/all_quivered_lq.fasta $out_dir/all_quivered_lq.fastq 2>&1 > /dev/null && echo $?
  0
