#Test post_mapping_to_genome.py

  $ . $TESTDIR/setup.sh

  $ D=$SIVDATDIR/test_make_abundance
  $ in_hq=$D/combined/all.polished_hq.fastq
  $ in_ha=$D/combined/all.polished_hq.fasta
  $ in_pickle=$D/combined/all.hq_lq_pre_dict.pickle
  $ in_sam=$D/sorted_gmap_alignments.sam

  $ OD=$OUTDIR/cram_test_post_mapping_to_genome

  $ rm -rf $OD && mkdir -p $OD

# Test post_mapping_to_genome.py
# Test fasta input
  $ O_fa=$OD/output_mapped.fasta
  $ O_fq=$OD/output_mapped.fastq
  $ O_gff=$OD/output_mapped.gff
  $ O_abundance=$OD/output_mapped.abundance.txt
  $ O_group=$OD/output_mapped.group.txt
  $ O_read_stat=$OD/output_mapped.read_stat.txt

  $ post_mapping_to_genome.py --quiet $in_ha $in_sam $in_pickle $O_fa $O_gff --abundance_fn $O_abundance --group_fn $O_group --read_stat_fn $O_read_stat && echo $?
  0

  $ cat $O_abundance | grep -v '^#' |head -2
  pbid\tcount_fl\tcount_nfl\tcount_nfl_amb\tnorm_fl\tnorm_nfl\tnorm_nfl_amb (esc)
  PB.1.1\t30\t30\t30.83\t6.3667e-03\t3.0367e-03\t3.0637e-03 (esc)

  $ head -1 $O_group
  PB.1.1\ti0_HQ_sample18ba5d|c11/f3p0/458,i0_HQ_sample18ba5d|c1353/f2p0/462,i0_HQ_sample18ba5d|c152/f3p0/462,i0_HQ_sample18ba5d|c1543/f8p1/465,i0_HQ_sample18ba5d|c38/f2p0/460,i0_HQ_sample18ba5d|c563/f2p0/473,i0_HQ_sample18ba5d|c633/f2p1/462,i0_HQ_sample18ba5d|c81/f3p0/462,i0_HQ_sample18ba5d|c140/f3p0/462,i0_HQ_sample18ba5d|c204/f2p0/422 (esc)

  $ head -2 $O_read_stat
  id\tlength\tis_fl\tstat\tpbid (esc)
  m54006_160328_233933/17957574/30_1837_CCS\t1807\tY\tunique\tPB.5.4 (esc)


# Test fasta input  
  $ post_mapping_to_genome.py --quiet $in_hq $in_sam $in_pickle $O_fq $O_gff && echo $?
  0

  $ cat $O_fq | wc -l
  140

  $ cat $O_gff | grep -v '^##' |wc -l
  244

  $ cat $O_gff | grep -v '^##' |head -1
  SIRV1\tPacBio\ttranscript\t10713\t11643\t.\t+\t.\tgene_id "PB.1"; transcript_id "PB.1.1"; (esc)

  $ rm -f $OD/output_mapped.*
