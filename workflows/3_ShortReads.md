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

## Seq-to-pangenome alignment

Indexing:

```shell
cd $DIR_BASE/short_reads

GFA=$DIR_BASE/chr1-19+XY.ref+pan.p95.gfa

# Work on /scratch (of octopus11)
sbatch -p workers -w octopus11 -c 48 --job-name mice --wrap "hostname; cd /scratch; \time -v $RUN_ODGI paths -i $GFA -f > $DIR_BASE/chr1-19+XY.ref+pan.p95.fa; \time -v bwa index $DIR_BASE/chr1-19+XY.ref+pan.p95.fa"
```

Alignment:

```shell
cd $DIR_BASE/short_reads

GFA=$DIR_BASE/chr1-19+XY.ref+pan.p95.gfa
FASTA=$DIR_BASE/chr1-19+XY.ref+pan.p95.fa

sed '1d' $DIR_BASE/data/accessions.tsv | head -n 1 | while read -r STRAIN ACC; do
    STRAIN=$(echo $STRAIN | tr '/' '_')
    echo $STRAIN $ACC;
    r1=$DIR_BASE/short_reads/${ACC}_1.fastq.gz;
    r2=$DIR_BASE/short_reads/${ACC}_2.fastq.gz;

    # Work on /scratch (of octopus11)
    sbatch -p workers -c 48 --job-name mice --wrap "hostname; cd /scratch; \time -v bwa mem -t 48 $FASTA $r1 $r2 > $STRAIN.sam; samtools view -b $STRAIN.sam -@ 48 > $STRAIN.bam; gfainject --gfa $GFA --bam $STRAIN.bam > $STRAIN.gfainject.gaf; gafpack --graph $GFA --alignments $STRAIN.gfainject.gaf > $STRAIN.pack; mv $STRAIN.* $DIR_BASE/short_reads/"
done
```
