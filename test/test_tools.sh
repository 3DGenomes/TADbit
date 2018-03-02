# /bin/sh

tmpdir=_tmp_test_tools

mkdir -p $tmpdir
mkdir -p $tmpdir/FASTQs
mkdir -p $tmpdir/db

echo '\n'
date
echo '\n'

echo ' - downloading Hi-C experiments 1/3'
fastq-dump -A SRR4433970 -DQ '+' --defline-seq '@$ac.$si' --split-files -E  --qual-filter-1 -X 500000
mv  SRR4433970_1.fastq $tmpdir/FASTQs/SRR4433970a_1.fastq
mv  SRR4433970_2.fastq $tmpdir/FASTQs/SRR4433970a_2.fastq

echo ' - downloading Hi-C experiments 2/3'
fastq-dump -A SRR4433970 -DQ '+' --defline-seq '@$ac.$si' --split-files -E  --qual-filter-1 -N 500000 -X 1000000
mv  SRR4433970_1.fastq $tmpdir/FASTQs/SRR4433970b_1.fastq
mv  SRR4433970_2.fastq $tmpdir/FASTQs/SRR4433970b_2.fastq

echo ' - downloading Hi-C experiments 3/3'
fastq-dump -A SRR4433971 -DQ '+' --defline-seq '@$ac.$si' --split-files -E  --qual-filter-1 -X 1000000
mv  SRR4433971_1.fastq $tmpdir/FASTQs/SRR4433971_1.fastq
mv  SRR4433971_2.fastq $tmpdir/FASTQs/SRR4433971_2.fastq

echo '\n'
date
echo '\n'

echo ' - downloading Yeast genome'
wget http://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz

echo ''
echo ' - uncompressing Yeast genome'
tar xzvf chromFa.tar.gz --to-stdout > $tmpdir/db/yeast.fa

echo ''
echo ' - indexing'
gem-indexer -t 8 -i $tmpdir/db/yeast.fa -o $tmpdir/db/yeast

echo '\n'
date
echo '\n'

echo ''
echo 'Mapping'
echo tadbit map $tmpdir/rep1 --fastq $tmpdir/FASTQs/SRR4433970a_1.fastq --read 1 --index $tmpdir/db/yeast.gem --renz Sau3AI
tadbit map $tmpdir/rep1 --fastq $tmpdir/FASTQs/SRR4433970a_1.fastq --read 1 --index $tmpdir/db/yeast.gem --renz Sau3AI

echo tadbit map $tmpdir/rep1 --fastq $tmpdir/FASTQs/SRR4433970a_2.fastq --read 2 --index $tmpdir/db/yeast.gem --renz Sau3AI
tadbit map $tmpdir/rep1 --fastq $tmpdir/FASTQs/SRR4433970a_2.fastq --read 2 --index $tmpdir/db/yeast.gem --renz Sau3AI

echo tadbit map $tmpdir/rep1 --fastq $tmpdir/FASTQs/SRR4433970b_1.fastq --read 1 --index $tmpdir/db/yeast.gem --renz Sau3AI
tadbit map $tmpdir/rep1 --fastq $tmpdir/FASTQs/SRR4433970b_1.fastq --read 1 --index $tmpdir/db/yeast.gem --renz Sau3AI

echo tadbit map $tmpdir/rep1 --fastq $tmpdir/FASTQs/SRR4433970b_2.fastq --read 2 --index $tmpdir/db/yeast.gem --renz Sau3AI
tadbit map $tmpdir/rep1 --fastq $tmpdir/FASTQs/SRR4433970b_2.fastq --read 2 --index $tmpdir/db/yeast.gem --renz Sau3AI

echo tadbit map $tmpdir/rep2 --fastq $tmpdir/FASTQs/SRR4433971_1.fastq  --read 1 --index $tmpdir/db/yeast.gem --renz HindIII
tadbit map $tmpdir/rep2 --fastq $tmpdir/FASTQs/SRR4433971_1.fastq  --read 1 --index $tmpdir/db/yeast.gem --renz HindIII

echo tadbit map $tmpdir/rep2 --fastq $tmpdir/FASTQs/SRR4433971_2.fastq  --read 2 --index $tmpdir/db/yeast.gem --renz HindIII
tadbit map $tmpdir/rep2 --fastq $tmpdir/FASTQs/SRR4433971_2.fastq  --read 2 --index $tmpdir/db/yeast.gem --renz HindIII

echo '\n'
date
echo '\n'

echo ''
echo 'Parsing 1/2'
tadbit parse $tmpdir/rep1 --genome $tmpdir/db/yeast.fa &
echo 'Parsing 2/2'
tadbit parse $tmpdir/rep2 --genome $tmpdir/db/yeast.fa

echo ''
echo 'Filtering 1/2'
tadbit filter $tmpdir/rep1 &
echo 'Filtering 2/2'
tadbit filter $tmpdir/rep2 

echo ''
echo 'Normalizing 1/2'
tadbit normalize $tmpdir/rep1 -r 100000 &
echo 'Normalizing 2/2'
tadbit normalize $tmpdir/rep2 -r 100000 --min_count 100

echo 'Merging'
tadbit merge $tmpdir/both -w1 $tmpdir/rep1 -w2 $tmpdir/rep2 -r 100000 --norm

echo 'Normalizing'
tadbit normalize $tmpdir/both -r 100000
tadbit normalize $tmpdir/both -r 10000

echo 'TADs/compartments'
tadbit segment $tmpdir/both -r 10000 --fasta $tmpdir/db/yeast.fa
