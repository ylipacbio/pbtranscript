Test pbtranscript modes other than 'cluster'

Set up
  $ . $TESTDIR/setup.sh

#Test pbtranscript classify
  $ OUTFA=$OUTDIR/classify_out_1.fasta
  $ OUTCSV=$OUTDIR/classify_out_1.primer_info.csv 
  $ OUTSUMMARY=$OUTDIR/classify_out_1.classify_summary.txt 

  $ CLASSIFYOUTDIR=$OUTDIR/test_classify
  $ rm -rf $OUTFA ${OUTFA}.fai $OUTCSV $OUTSUMMARY $CLASSIFYOUTDIR
  $ pbtranscript classify $DATDIR/test.fasta $OUTFA -d $CLASSIFYOUTDIR

  $ ls $OUTFA
  *classify_out_1.fasta (glob)

  $ ls $OUTCSV
  *classify_out_1.primer_info.csv (glob)


  $ OUTFA=$OUTDIR/classify_out_2.fasta 
  $ OUTCSV=$OUTDIR/classify_out_2.csv 
  $ OUTSUMMARY=$OUTDIR/classify_out_2.txt 
  $ CLASSIFYOUTDIR=$OUTDIR/test_classify
  $ rm -rf $OUTFA $OUTCSV $OUTSUMMARY $CLASSIFYOUTDIR
  $ pbtranscript classify $DATDIR/reads_of_insert.fasta $OUTFA -d $CLASSIFYOUTDIR --report $OUTCSV --summary $OUTSUMMARY

  $ ls $OUTFA
  *classify_out_2.fasta (glob)

  $ ls $OUTCSV
  *classify_out_2.csv (glob)

  $ cat $OUTSUMMARY
  Number of consensus reads=22
  Number of five prime reads=11
  Number of three prime reads=13
  Number of poly-A reads=13
  Number of filtered short reads=0
  Number of non-full-length reads=13
  Number of full-length reads=9
  Number of full-length non-chimeric reads=9
  Number of full-length non-chimeric bases=35564
  Mean full-length non-chimeric read length=3951

  $ rm -rf $OUTTXT $OUTFA
  $ INFA=$DATDIR/test_subset.fasta 
  $ OUTFA=$OUTDIR/test_subset.fasta
  $ OUTTXT=$OUTDIR/test_subset.txt 
  $ pbtranscript subset $INFA $OUTFA --FL --nonChimeric

  $ rm -rf $OUTTXT
  $ pbtranscript subset $INFA $OUTTXT --FL --nonChimeric --printReadLengthOnly
  $ cat $OUTTXT
  3889
  4194
  3950
  4409
  4106
  3784
 
  $ rm -rf $OUTTXT
  $ pbtranscript subset $INFA $OUTTXT --chimeric --printReadLengthOnly
  $ wc -l $OUTTXT | cut -f 1 -d ' ' 
  5

# Test classify with --detect_chimera_nfl
  $ rm -rf $OUTFA $OUTCSV $OUTSUMMARY $CLASSIFYOUTDIR
  $ pbtranscript classify $DATDIR/test.fasta $OUTFA -d $CLASSIFYOUTDIR --detect_chimera_nfl

  $ ls $CLASSIFYOUTDIR/nflnc.fasta
  *nflnc.fasta (glob)
  $ ls $CLASSIFYOUTDIR/nflc.fasta
  *nflc.fasta (glob)


# Test customized primers and ignore_polyA
  $ CLASSIFYOUTDIR=$OUTDIR/test_custom_primer_out
  $ OUTFA=$OUTDIR/test_custom_primer.fasta
  $ rm -rf $CLASSIFYOUTDIR $OUTFA
  $ pbtranscript classify $DATDIR/test_custom_primer.fasta $OUTFA -d $CLASSIFYOUTDIR --primer $DATDIR/customized_primer.fasta --ignore_polyA
  $ grep '>' $OUTFA |wc -l 
  1
