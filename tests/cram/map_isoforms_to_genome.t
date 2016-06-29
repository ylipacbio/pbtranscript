#Test make_abundance.py

  $ . $TESTDIR/setup.sh

  $ D=$SIVDATDIR/test_make_abundance
  $ in_hq=$D/combined/all.polished_hq.fastq
  $ in_ha=$D/combined/all.polished_hq.fasta
  $ gmap_db=$SIVDATDIR/gmap-referenceset-root-dir/SIRV
  $ gmap_name=gmap_db
  $ gmap_ds=$SIVDATDIR/gmap-referenceset-root-dir/SIRV/gmapreferenceset.xml

  $ OD=$OUTDIR/cram_test_map_isoforms_to_genome

  $ rm -rf $OD && mkdir -p $OD

# Test map_isoforms_to_genome
  $ out_sam=$OD/output1.sam
  $ map_isoforms_to_genome.py --quiet $in_hq $out_sam --gmap_ds $gmap_ds && echo $?
  0

  $ out_sam=$OD/output2.sam
  $ map_isoforms_to_genome.py --quiet $in_hq $out_sam --gmap_db $gmap_db --gmap_name gmap_db && echo $?
  0
