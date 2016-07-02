#Test map_isoforms_to_genome.py and collapse_mapped_isoforms.py.

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
  $ map_isoforms_to_genome.py --quiet $ISOFORM_FA $OD/gmap-output.sam --gmap_name=gmap_db --gmap_db=$SIVDATDIR/gmap-referenceset-root-dir/SIRV/ --gmap_nproc=10 && echo $?
  0

  $ map_isoforms_to_genome.py --quiet $ISOFORM_FA $OD/gmap-output.sam --gmap_ds=$SIVDATDIR/gmap-referenceset-root-dir/SIRV/gmapreferenceset.xml --gmap_nproc=10 && echo $?
  0

# Test collapse_mapped_isoforms
  $ collapse_mapped_isoforms.py --quiet $ISOFORM_FA $SORTED_GMAP_OUTPUT $OD/fa_in && echo $?
  0

  $ collapse_mapped_isoforms.py --quiet $ISOFORM_FQ $SORTED_GMAP_OUTPUT $OD/fq_in --collapsed_isoforms=$OD/output_fq_in.fastq && echo $?
  0
  $ cat $OD/output_fq_in.fastq | wc -l
  240

  $ collapse_mapped_isoforms.py --quiet $ISOFORM_FA_DS $SORTED_GMAP_OUTPUT $OD/fa_ds_in && echo $?
  0

  $ collapse_mapped_isoforms.py --quiet $ISOFORM_FQ_DS $SORTED_GMAP_OUTPUT $OD/fq_ds_in --collapsed_isoforms=$OD/output_fq_ds_in.fastq & echo $?
  0
# Strangely, $OD/output_fq_ds_in.fastq (a soft link output) can not be accessed within this cram test, even with sleep 10.
