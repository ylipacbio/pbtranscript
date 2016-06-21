#Test make_abundance.py

  $ . $TESTDIR/setup.sh

  $ D=$SIVDATDIR/test_make_abundance
  $ GROUP=$D/group.txt
  $ PICKLE=$D/combined/all.hq_lq_pre_dict.pickle

  $ OD=$OUTDIR/cram_test_counting

  $ rm -rf $OD && mkdir -p $OD

# Test make_abundance.py
  $ O_READ_STAT=$OD/output.read_stat.txt
  $ O_ABUNDANCE=$OD/output.abundance.txt
  $ make_abundance.py --quiet $GROUP $PICKLE $O_READ_STAT $O_ABUNDANCE && echo $?
  0

  $ cat $O_READ_STAT |wc -l
  10416

  $ head -2 $O_READ_STAT 
  id	length	is_fl	stat	pbid
  m54006_160328_233933/17957574/30_1837_CCS	1807	Y	unique	PB.5.4

  $ cat $O_ABUNDANCE | grep -v ^# | wc -l
  39
