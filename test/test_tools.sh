# /bin/bash

tmpdir=_tmp_test_tools
LOG=$tmpdir/_test_tools.log
SEP='--------------------------------------------------------------------------------'

mkdir -p $tmpdir
mkdir -p $tmpdir/FASTQs
mkdir -p $tmpdir/db

nreads=1000000

echo '\n' | tee -a $LOG
echo `date` | tee -a $LOG
echo '\n' | tee -a $LOG

echo ' - downloading Hi-C experiments 1/3' | tee -a $LOG
             #sample ID #no quality header #sequence header #split reads #quality filter #quality filter #stop download after 1M
fastq-dump -A SRR4433970 -DQ '+' --defline-seq '@$ac.$si' --split-files --qual-filter --qual-filter-1 -X $nreads
mv  SRR4433970_1.fastq $tmpdir/FASTQs/SRR4433970a_1.fastq
mv  SRR4433970_2.fastq $tmpdir/FASTQs/SRR4433970a_2.fastq

echo ' - downloading Hi-C experiments 2/3' | tee -a $LOG
             #sample ID #no quality header #sequence header #split reads #quality filter #quality filter #start download after 1M #stop downnload after 2M
fastq-dump -A SRR4433970 -DQ '+' --defline-seq '@$ac.$si' --split-files --qual-filter --qual-filter-1 -N $nreads -X $((nreads/2))
mv  SRR4433970_1.fastq $tmpdir/FASTQs/SRR4433970b_1.fastq
mv  SRR4433970_2.fastq $tmpdir/FASTQs/SRR4433970b_2.fastq

echo ' - downloading Hi-C experiments 3/3' | tee -a $LOG
fastq-dump -A SRR4433971 -DQ '+' --defline-seq '@$ac.$si' --split-files -E  --qual-filter-1 -X $((nreads*2))
mv  SRR4433971_1.fastq $tmpdir/FASTQs/SRR4433971_1.fastq
mv  SRR4433971_2.fastq $tmpdir/FASTQs/SRR4433971_2.fastq

echo '\n' | tee -a $LOG
echo `date` | tee -a $LOG
echo '\n' | tee -a $LOG

echo ' - downloading Yeast genome' | tee -a $LOG
wget http://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz

echo ''
echo ' - uncompressing Yeast genome' | tee -a $LOG
tar xzvf chromFa.tar.gz --to-stdout > $tmpdir/db/yeast.fa
rm -f chromFa.tar.gz

echo ''
echo ' - indexing' | tee -a $LOG
gem-indexer -t 8 -i $tmpdir/db/yeast.fa -o $tmpdir/db/yeast

echo '\n' | tee -a $LOG
echo `date` | tee -a $LOG
echo '\n' | tee -a $LOG

START=$(date +%s);

echo '' | tee -a $LOG
echo 'Mapping\n' | tee -a $LOG
echo $SEP"\n   $ " tadbit map $tmpdir/rep1 --fastq $tmpdir/FASTQs/SRR4433970a_1.fastq --read 1 --index $tmpdir/db/yeast.gem --renz Sau3AI | tee -a $LOG
tadbit map $tmpdir/rep1 --fastq $tmpdir/FASTQs/SRR4433970a_1.fastq --read 1 --index $tmpdir/db/yeast.gem --renz Sau3AI 2>> $LOG

echo $SEP"\n   $ " tadbit map $tmpdir/rep1 --fastq $tmpdir/FASTQs/SRR4433970a_2.fastq --read 2 --index $tmpdir/db/yeast.gem --renz Sau3AI | tee -a $LOG
tadbit map $tmpdir/rep1 --fastq $tmpdir/FASTQs/SRR4433970a_2.fastq --read 2 --index $tmpdir/db/yeast.gem --renz Sau3AI 2>> $LOG

echo $SEP"\n   $ " tadbit map $tmpdir/rep1 --fastq $tmpdir/FASTQs/SRR4433970b_1.fastq --read 1 --index $tmpdir/db/yeast.gem --renz Sau3AI | tee -a $LOG
tadbit map $tmpdir/rep1 --fastq $tmpdir/FASTQs/SRR4433970b_1.fastq --read 1 --index $tmpdir/db/yeast.gem --renz Sau3AI 2>> $LOG

echo $SEP"\n   $ " tadbit map $tmpdir/rep1 --fastq $tmpdir/FASTQs/SRR4433970b_2.fastq --read 2 --index $tmpdir/db/yeast.gem --renz Sau3AI | tee -a $LOG
tadbit map $tmpdir/rep1 --fastq $tmpdir/FASTQs/SRR4433970b_2.fastq --read 2 --index $tmpdir/db/yeast.gem --renz Sau3AI 2>> $LOG

echo $SEP"\n   $ " tadbit map $tmpdir/rep2 --fastq $tmpdir/FASTQs/SRR4433971_1.fastq  --read 1 --index $tmpdir/db/yeast.gem --renz HindIII | tee -a $LOG
tadbit map $tmpdir/rep2 --fastq $tmpdir/FASTQs/SRR4433971_1.fastq  --read 1 --index $tmpdir/db/yeast.gem --renz HindIII 2>> $LOG

echo $SEP"\n   $ " tadbit map $tmpdir/rep2 --fastq $tmpdir/FASTQs/SRR4433971_2.fastq  --read 2 --index $tmpdir/db/yeast.gem --renz HindIII | tee -a $LOG
tadbit map $tmpdir/rep2 --fastq $tmpdir/FASTQs/SRR4433971_2.fastq  --read 2 --index $tmpdir/db/yeast.gem --renz HindIII 2>> $LOG

echo '\n' | tee -a $LOG
echo `date` | tee -a $LOG
echo '\n' | tee -a $LOG

echo '' | tee -a $LOG
echo 'Parsing 1/2\n' | tee -a $LOG
echo $SEP"\n   $ " tadbit parse $tmpdir/rep1 --genome $tmpdir/db/yeast.fa | tee -a $LOG
tadbit parse $tmpdir/rep1 --genome $tmpdir/db/yeast.fa 2>> $LOG
echo 'Parsing 2/2\n'
echo $SEP"\n   $ " tadbit parse $tmpdir/rep2 --genome $tmpdir/db/yeast.fa | tee -a $LOG
tadbit parse $tmpdir/rep2 --genome $tmpdir/db/yeast.fa 2>> $LOG

echo '' | tee -a $LOG
echo 'Filtering 1/2\n' | tee -a $LOG
echo $SEP"\n   $ " tadbit filter $tmpdir/rep1 | tee -a $LOG
tadbit filter $tmpdir/rep1  2>> $LOG
echo 'Filtering 2/2\n' | tee -a $LOG
echo $SEP"\n   $ " tadbit filter $tmpdir/rep2 | tee -a $LOG
tadbit filter $tmpdir/rep2 2>> $LOG

echo '' | tee -a $LOG
echo 'Normalizing 1/2\n' | tee -a $LOG
echo $SEP"\n   $ " tadbit normalize $tmpdir/rep1 -r 100000 | tee -a $LOG
tadbit normalize $tmpdir/rep1 -r 100000 2>> $LOG
echo 'Normalizing 2/2\n'
echo $SEP"\n   $ " tadbit normalize $tmpdir/rep2 -r 100000 --min_count 100 | tee -a $LOG
tadbit normalize $tmpdir/rep2 -r 100000 --min_count 100 2>> $LOG

echo '' | tee -a $LOG
echo 'Merging\n' | tee -a $LOG
echo $SEP"\n   $ " tadbit merge $tmpdir/both -w1 $tmpdir/rep1 -w2 $tmpdir/rep2 -r 100000 --norm | tee -a $LOG
tadbit merge $tmpdir/both -w1 $tmpdir/rep1 -w2 $tmpdir/rep2 -r 100000 --norm 2>> $LOG

echo '' | tee -a $LOG
echo 'Normalizing\n' | tee -a $LOG
echo $SEP"\n   $ " tadbit normalize $tmpdir/both -r 100000 | tee -a $LOG
tadbit normalize $tmpdir/both -r 100000 2>> $LOG
echo $SEP"\n   $ " tadbit normalize $tmpdir/both -r 10000 | tee -a $LOG
tadbit normalize $tmpdir/both -r 10000 2>> $LOG

echo '' | tee -a $LOG
echo 'Binning\n' | tee -a $LOG
echo $SEP"\n   $ " tadbit bin $tmpdir/both -r 100000 -c chrII --norm raw norm decay | tee -a $LOG
tadbit bin $tmpdir/both -r 100000 -c chrII --norm raw norm decay 2>> $LOG
echo $SEP"\n   $ " tadbit bin $tmpdir/both -r 100000 -c chrIII -c2 chrVII:100000-1000000 | tee -a $LOG
tadbit bin $tmpdir/both -r 100000 -c chrIII -c2 chrVII:100000-1000000 2>> $LOG
echo $SEP"\n   $ " tadbit bin $tmpdir/both -r 100000 --plot --matrix | tee -a $LOG
tadbit bin $tmpdir/both -r 100000 --plot --matrix 2>> $LOG

echo '' | tee -a $LOG
echo 'TADs/compartments\n' | tee -a $LOG
echo $SEP"\n   $ " tadbit segment $tmpdir/both -r 10000 --fasta $tmpdir/db/yeast.fa | tee -a $LOG
tadbit segment $tmpdir/both -r 10000 --fasta $tmpdir/db/yeast.fa 2>> $LOG

echo '\n' | tee -a $LOG
echo `date` | tee -a $LOG
echo '\n' | tee -a $LOG


END=$(date +%s)

echo $((END-START)) | awk '{print "TADbit done in: " int($1/60)"m "int($1%60)"s"}' | tee -a $LOG

errors=`grep -ic error $LOG`

echo "\n\n -> Found" $errors "errors\n" | tee -a $LOG


if [ $errors -eq 0]
then
    echo 'Cleanning temporary directory'
    rm -rf $tmpdir
else
    echo '  ==>> Check LOG in: ' $LOG 'for details'
fi

echo "Done." | tee -a $LOG
