gmap_db_dir=/pbi/dept/secondary/siv/testdata/pbtranscript-unittest/data/gmap_db/
gmap_db_name=SIRV
fq=gmap-input.fastq
gmap -D $gmap_db_dir -d $gmap_db_name -n 0 -t 15 -z sense_force --cross-species -f samse $fq > gmap-output.sam
