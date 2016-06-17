#Test map_isoforms.py and collapse_isoforms.py.

  $ . $TESTDIR/setup.sh

  $ D=$SIVDATDIR/test_collapsing
  $ GMAP_OUTPUT=$D/gmap-output.sam
  $ SORTED_GMAP_OUTPUT=$D/sorted-gmap-output.sam
  $ ISOFORM_FA=$D/gmap-input.fasta
  $ ISOFORM_FQ=$D/gmap-input.fastq
  $ ISOFORM_FQ_DS=$D/gmap-input.fastq.contigset.xml
  $ ISOFORM_FA_DS=$D/gmap-input.fasta.contigset.xml

  $ OD=$OUTDIR/cram_test_collapsing

  $ rm -rf $OD && mkdir -p $OD

# Test map_isoforms
  $ map_isoforms.py --quiet $ISOFORM_FA $OD/gmap-output.sam --gmap_name=SIRV --gmap_db=/pbi/dept/secondary/siv/testdata/pbtranscript-unittest/data/gmap_db --gmap_nproc=10 && echo $?
  0

  $ map_isoforms.py --quiet $ISOFORM_FA $OD/gmap-output.sam --gmap_ds=$SIVDATDIR/test_map_isoforms/sirv.gmapreferenceset.xml --gmap_nproc=10 && echo $?
  0

# Test collapse_isoforms
  $ collapse_isoforms.py --quiet $ISOFORM_FA $SORTED_GMAP_OUTPUT $OD/fa_in.fasta && echo $?
  0

  $ collapse_isoforms.py --quiet $ISOFORM_FA $SORTED_GMAP_OUTPUT $OD/fa_in.contigset.xml && echo $?
  0

  $ collapse_isoforms.py --quiet $ISOFORM_FQ $SORTED_GMAP_OUTPUT $OD/fq_in.fasta && echo $?
  0

  $ collapse_isoforms.py --quiet $ISOFORM_FQ $SORTED_GMAP_OUTPUT $OD/fq_in.contigset.xml && echo $?
  0

  $ collapse_isoforms.py --quiet $ISOFORM_FA_DS $SORTED_GMAP_OUTPUT $OD/fa_ds_in.fasta && echo $?
  0

  $ collapse_isoforms.py --quiet $ISOFORM_FA_DS $SORTED_GMAP_OUTPUT $OD/fa_ds_in.contigset.xml && echo $?
  0

  $ collapse_isoforms.py --quiet $ISOFORM_FQ_DS $SORTED_GMAP_OUTPUT $OD/fq_ds_in.fastq && echo $?
  0

  $ collapse_isoforms.py --quiet $ISOFORM_FQ_DS $SORTED_GMAP_OUTPUT $OD/fq_ds_in.contigset.xml & echo $?
  0
