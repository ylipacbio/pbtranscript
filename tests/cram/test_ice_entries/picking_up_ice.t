Test picking_up_ice.py

How to create testing dataset which is a halted `pbtranscript cluster` job?
Call the cmd below and then Ctrl + c to halt the job.
$ pbtranscript --verbose cluster isoseq_flnc.fasta out.fasta --bas_fofn=$subreads -d=output --flnc_reads_per_split=3

How to test picking_up_ice.py?
$ picking_up_ice.py output/output/input.split_001.fasta.pickle --root_dir=output


Set up
  $ . $TESTDIR/setup.sh

  $ subreads=/pbi/dept/secondary/siv/testdata/SA3-Sequel/rc0/315/3150353/r54086_20160831_010819/4_D01_tiny/tiny_flea_isoseq.subreadset.xml
  $ flnc=$SIVDATDIR/test_ice_entries/test_picking_up_ice/isoseq_flnc.fasta
  $ out_dir=$OUTDIR/test_picking_up_ice

  $ rm -rf $out_dir && mkdir -p $out_dir
  $ pbtranscript --verbose cluster $flnc $out_dir/out.fasta --bas_fofn=$subreads -d=$out_dir --flnc_reads_per_split=4 1>/dev/null 2>/dev/null && echo $?
  0

Delete intermediate ICE output pickle file for split_002
  $ rm $out_dir/output/input.split_002.fasta.consensus.fasta  $out_dir/output/input.split_002.fasta.pickle && echo $?
  0
 
Restart from ICE pickle file for split_001
  $ picking_up_ice.py $out_dir/output/input.split_001.fasta.pickle --root_dir=$out_dir --flnc=$flnc 1>/dev/null 2>/dev/null && echo $?
  0

  $ ls $out_dir/output/input.split_002.fasta.consensus.fasta  $out_dir/output/input.split_002.fasta.pickle 1>/dev/null && echo $?
  0
