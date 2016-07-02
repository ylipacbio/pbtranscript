#Test filter_collapsed_isoforms.py

  $ . $TESTDIR/setup.sh

  $ D=$SIVDATDIR/test_filtering
  $ in_rep_fq=$D/in.rep.fastq
  $ in_gff=$D/in.gff
  $ in_abundance=$D/in.abundance.txt
  $ in_group=$D/in.group.txt

  $ OD=$OUTDIR/cram_test_filtering

  $ rm -rf $OD && mkdir -p $OD

# Test filter_collapsed_isoforms.py --min_count 20
  $ out_rep_fq=$OD/filter_by_count.min_count_20.rep.fastq
  $ out_gff=$OD/filter_by_count.min_count_20.gff
  $ out_abundance=$OD/filter_by_count.min_count_20.abundance.txt
  $ out_group=$OD/filter_by_count.min_count_20.group.txt

  $ filter_collapsed_isoforms.py --quiet $in_rep_fq $out_rep_fq --min_count 20 && echo $?
  0

  $ grep -v '^#' $out_abundance | cut -f 1
  pbid
  PB.2.5
  PB.5.1
  PB.7.1
  PB.10.2
  PB.10.42
  PB.12.1

  $ cut -f 9 $out_gff |cut -f 2 -d ';'|uniq
   transcript_id "PB.2.5"
   transcript_id "PB.5.1"
   transcript_id "PB.7.1"
   transcript_id "PB.10.2"
   transcript_id "PB.10.42"
   transcript_id "PB.12.1"

  $ cut -f 1 -d '|' $out_rep_fq |grep '^@'
  @PB.2.5
  @PB.5.1
  @PB.7.1
  @PB.10.2
  @PB.10.42
  @PB.12.1


# Test filter_collapsed_isoforms.py --min_count 0 --no_filter_subsets
  $ out_rep_fq=$OD/filter_by_count.min_count_0.no_subsets.rep.fastq
  $ out_gff=$OD/filter_by_count.min_count_0.no_subsets.gff
  $ out_abundance=$OD/filter_by_count.min_count_0.no_subsets.abundance.txt
  $ out_group=$OD/filter_by_count.min_count_0.no_subsets.group.txt

  $ filter_collapsed_isoforms.py --quiet $in_rep_fq $out_rep_fq --min_count 20 --no_filter_subsets && echo $?
  0

  $ grep -v '^#' $out_abundance | cut -f 1
  pbid
  PB.2.5
  PB.5.1
  PB.7.1
  PB.10.2
  PB.10.42
  PB.12.1

  $ cut -f 9 $out_gff |cut -f 2 -d ';'|uniq
   transcript_id "PB.2.5"
   transcript_id "PB.5.1"
   transcript_id "PB.7.1"
   transcript_id "PB.10.2"
   transcript_id "PB.10.42"
   transcript_id "PB.12.1"

  $ cut -f 1 -d '|' $out_rep_fq |grep '^@'
  @PB.2.5
  @PB.5.1
  @PB.7.1
  @PB.10.2
  @PB.10.42
  @PB.12.1
