# Short reads

Variables:

```shell
DIR_BASE=/lizardfs/guarracino/keane_mouse_pangenome
RUN_ODGI=/home/guarracino/tools/odgi/bin/odgi-34f006f31c3f6b35a1eb8d58a4edb1c458583de3
```

## Download

```shell
mkdir -p $DIR_BASE/short_reads

# Obtain FASTQ urls from ENA accessions
rm accessions.urls.tsv; cut -f 2 $DIR_BASE/data/accessions.tsv | sed '1d' | while read ACC; do echo $ACC; wget -c https://www.ebi.ac.uk/ena/portal/api/filereport\?accession\="${ACC}"\&result\=read_run\&fastq_ftp,fastq_md5 -O - | sed '1d' | cut -f 2,4 | awk '{split($1,a,";"); split($2,b,";"); print(a[1],b[1],"\n"a[2],b[2],"\n"a[3],b[3])}' | grep '_' | awk -v acc=$ACC -v OFS='\t' '{print(acc,$1,$2)}' >> accessions.urls.tsv; done

cd /scratch
cut -f 2 /lizardfs/guarracino/keane_mouse_pangenome/short_reads/accessions.urls.tsv | parallel -j 8 'x=$(echo {} | rev | cut -f 1 -d / | rev); wget -c {}; mv $x /lizardfs/guarracino/keane_mouse_pangenome/short_reads; echo got {}'
```

Statistics (slow):

```shell
# Create output file
echo -e "filename\tnumber.of.reads\tread.length" > output.tsv

# Loop through .fastq.gz files
for file in *_1.fastq.gz
do
  echo $file

  # Count number of reads (4 lines per read in a FASTQ file)
  count=$(echo $(($(zcat $file | wc -l) / 4)))
  
  # Determine read length (assuming all reads are the same length)
  # get second line of the first read and count its length
  length=$(zcat $file | head -n 2 | tail -n 1 | wc -c)

  # Write to output file
  echo -e "$file\t$count\t$length" >> output.tsv
done
```

Subsample to get the same coverage in all samples:

```shell
mm39_ref_size=2723431143
target_coverage=10

# Loop through .fastq.gz files
ls *_1.fastq.gz -lht | rev | cut -f 1 -d ' ' | rev | cut -f 1 -d '_'  | while read ERR; do
  if [ ! -f ${ERR}_1.cov${target_coverage}X.fastq.gz ]; then
    # Determine read length (assuming all reads are the same length)
    length=$(zcat ${ERR}_1.fastq.gz | head -n 2 | tail -n 1 | wc -c)

    # Count number of reads needed to get `target_coverage`
    count=$(echo "$mm39_ref_size * $target_coverage / $length" | bc)

    # Write to output file
    echo "$ERR\t$length\t$count"

    seqtk sample -s 17 ${ERR}_1.fastq.gz $count | pigz -9 -c > ${ERR}_1.cov${target_coverage}X.fastq.gz
    seqtk sample -s 17 ${ERR}_2.fastq.gz $count | pigz -9 -c > ${ERR}_2.cov${target_coverage}X.fastq.gz
  fi
done
```

## Seq-to-pangenome alignment

Indexing:

```shell
cd $DIR_BASE/short_reads

GFA=$DIR_BASE/chr1-19+XY.ref+pan.p95.gfa

# Work on /scratch (of octopus11)
sbatch -p workers -w octopus11 -c 48 --job-name mice --wrap "hostname; cd /scratch; \time -v $RUN_ODGI paths -i $GFA -f > $DIR_BASE/chr1-19+XY.ref+pan.p95.fa; \time -v bwa index $DIR_BASE/chr1-19+XY.ref+pan.p95.fa"
```

Alignment, injection and packing:

```shell
cd $DIR_BASE/short_reads
mkdir -p $DIR_BASE/short_reads/alignments

GFA=$DIR_BASE/graphs/chr1-19+XY.ref+pan.p95.gfa
FASTA=$DIR_BASE/graphs/chr1-19+XY.ref+pan.p95.fa
target_coverage=10

sed '1d' $DIR_BASE/short_reads/accessions.tsv | while read -r STRAIN ACC; do
  STRAIN=$(echo $STRAIN | tr '/' '_')
  PACK=$DIR_BASE/short_reads/alignments/$STRAIN.cov${target_coverage}X.pack
  if [ ! -f $PACK ]; then
    echo $STRAIN $ACC;
    r1=$DIR_BASE/short_reads/${ACC}_1.cov${target_coverage}X.fastq.gz;
    r2=$DIR_BASE/short_reads/${ACC}_2.cov${target_coverage}X.fastq.gz;

    sbatch -p workers -x octopus07 -c 48 --job-name mice --wrap "hostname; cd /scratch; \time -v bwa mem -t 48 $FASTA $r1 $r2 -m 5 > $STRAIN.cov${target_coverage}X.sam; \time -v samtools view -b $STRAIN.cov${target_coverage}X.sam -@ 48 > $STRAIN.cov${target_coverage}X.bam; pigz -9 $STRAIN.cov${target_coverage}X.sam; \time -v gfainject --gfa $GFA --bam $STRAIN.cov${target_coverage}X.bam > $STRAIN.cov${target_coverage}X.gfainject.gaf; \time -v gafpack --graph $GFA --alignments $STRAIN.cov${target_coverage}X.gfainject.gaf -c > $STRAIN.cov${target_coverage}X.pack; pigz -9 $STRAIN.cov${target_coverage}X.gfainject.gaf; pigz ; mv $STRAIN.cov${target_coverage}X.* $DIR_BASE/short_reads/alignments/"
  fi
done
```

Prepare the coverage matrix:

```shell

# Get samples
grep '^##' *.pack -m 1 | cut -f 3 -d ':' | tr -d '[:blank:]' | rev | cut -f 1 -d '/' | rev | sed "s/.cov10X.gfainject.gaf//g" > samples.txt

# Remove rows with the same value in all columns
paste *.pack | sed '1,2d' | awk '{for(i=2; i<=NF; i++) {if ($i != $1) {print; break}}}' | pigz -9 > matrix.tsv.gz

# if we want to remove something
#grep '^##' `ls *.pack | grep -v 'CAST_EiJ\|JF1_MsJ\|LEWES_EiJ'` -m 1 | cut -f 3 -d ':' | tr -d '[:blank:]' | rev | cut -f 1 -d '/' | rev | sed "s/.cov10X.gfainject.gaf//g" > samples.txt
#paste `ls *.pack | grep -v 'CAST_EiJ\|JF1_MsJ\|LEWES_EiJ'` | sed '1,2d' | awk '{for(i=2; i<=NF; i++) {if ($i != $1) {print; break}}}' | pigz -9 > matrix.tsv.gz

# Subsample rows
zcat matrix.tsv.gz | shuf -n 50000000 | pigz -9 > matrix.50M.tsv.gz
```

Compute PCA:

```shell
Rscript pca.R
```
