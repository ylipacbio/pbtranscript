Test ice_combine_cluster_bins

Set up
  $ . $TESTDIR/setup.sh

  $ src_tasks_dir=$SIVDATDIR/test_ice_entries/tasks
  $ out_dir=$OUTDIR/test_ice_combine_cluster_bins

  $ rm -rf $out_dir && mkdir -p $out_dir
  $ cp -r $src_tasks_dir/pbtranscript.tasks.separate_flnc-0/* $out_dir/

  $ tmp_dir=$out_dir/tmp
  $ ice_combine_cluster_bins.py --verbose $out_dir/separate_flnc.pickle $tmp_dir --cluster_bin_dirs=$out_dir/3to4kb_part0/cluster_out,$out_dir/4to5kb_part0/cluster_out 1>/dev/null 2>/dev/null && echo $?
  0

  $ ls $tmp_dir/all.cluster_report.csv $tmp_dir/all.cluster_summary.json $tmp_dir/all.consensus_isoforms.fasta $tmp_dir/all.polished_hq.fasta $tmp_dir/all.polished_lq.fasta $tmp_dir/all.polished_hq.fastq $tmp_dir/all.polished_lq.fastq 1>/dev/null 2>/dev/null && echo $?
  0
