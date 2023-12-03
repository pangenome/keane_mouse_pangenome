# Pangenome

Variables:

```shell
DIR_BASE=/lizardfs/guarracino/keane_mouse_pangenome
RUN_PGGB=/home/guarracino/tools/pggb/pggb-de5303e24d3e5594a5a2c9bdeb49aba420b24b0c
RUN_VCFBUB=/home/guarracino/tools/vcfbub/target/release/vcfbub-d92ae50da4926612f9bcf5d3ea8506a78a348734
RUN_ODGI=/home/guarracino/tools/odgi/bin/odgi-34f006f31c3f6b35a1eb8d58a4edb1c458583de3
RUN_VCFWAVE=/gnu/store/hkbkw85gjvsyqvx9vv9bw0ynmad989ag-vcflib-1.0.3+a36dbe9-11/bin/vcfwave
```

## Pangenome building

Run `pggb`:

```shell
mkdir -p $DIR_BASE/graphs/
cd $DIR_BASE/graphs/

( seq 1 19; echo X; echo Y) | while read i; do
    sbatch -p workers -x octopus02 -c 48 --job-name mice --wrap "hostname; $RUN_PGGB -i $DIR_BASE/partitioning/chr$i.ref+pan.fa.gz -p 95 -s 50000 -n 18 -o $DIR_BASE/graphs/chr$i.ref+pan.p95 -V mm39:#:100000 -t 48 -D /scratch;"
done

# C57BL_6NJ assembly has 2 chrM's contigs (C57BL_6NJ#1#21 and C57BL_6NJ#1#MT)
i=M
$RUN_PGGB -i $DIR_BASE/partitioning/chr$i.ref+pan.fa.gz -p 95 -s 1000 -n 19 -o $DIR_BASE/graphs/chr$i.ref+pan.p95.s1000.n19 -V mm39:#:100000 -t 48 -D /scratch;
```

Compress VCF files:

```shell
cd $DIR_BASE/graphs/

ls */*.vcf | while read v; do
    echo $v
    bgzip $v -@ 48
done
```

Merge all VCF files, except chrY's one because it misses a few samples (AKR_J and NOD_ShiLtJ):

```shell
bcftools concat -f <(ls $DIR_BASE/graphs/chr*.ref+pan.p95*/*.decomposed.vcf | grep -v chrY) | bcftools sort -T /scratch/bcftools.XXXXX | bgzip -@ 48 -c > chr1-19+X.ref+pan.p95.decomposed.vcf.gz
```

Squeeze graphs:

```shell
ls $DIR_BASE/graphs/*p95*/*.og | sort -k 1,1 -V > $DIR_BASE/graphs/graphs_to_squeeze.p95.txt

sbatch -p highmem -c 48 --job-name mice --wrap "hostname; cd /scratch; $RUN_ODGI squeeze -f $DIR_BASE/graphs/graphs_to_squeeze.p95.txt -o - -t 1 -P | $RUN_ODGI view -i - -g > chr1-19+XY.ref+pan.p95.gfa; mv chr1-19+XY.ref+pan.p95.gfa $DIR_BASE/graphs/"
```

Compress graphs:

```shell
sbatch -p headnode -c 1 --job-name mice --wrap "hostname; zstd -12 $DIR_BASE/graphs/chr1-19+XY.ref+pan.p95.gfa"

ls $DIR_BASE/graphs/*p95*/*.gfa | while read GFA; do
    echo $GFA
    zstd -12 $GFA
done
```

Statistics:

```shell
cd $DIR_BASE/graphs/

# Graphs
echo -e "chromosome\tlength\tnodes\tedges\tpaths\tsteps" > graph.statistics.tsv
(seq 1 19; echo X; echo Y; echo M) | while read i; do
  $RUN_ODGI stats -i chr$i.ref+pan.p95/chr$i.ref+pan.fa.gz.7976304.417fcdf.4cea6f5.smooth.final.og -S | grep '^#' -v | awk -v OFS='\t' -v chr=chr$i '{print(chr,$0)}'
done >> graph.statistics.tsv

# Variants
REFERENCES_FASTA_GZ=$DIR_BASE/assemblies/mm39.fasta.gz

echo -e "sample\tchromosome\tnum.samples\tnum.records\tno.alts\tsnps\tmnps\tindels\tothers\tmultiallelic.sites\tmultiallelic.snp.sites" > variant.stats.tsv
bcftools query -l chr1.ref+pan.p95/chr1.ref+pan.fa.gz.7976304.417fcdf.4cea6f5.smooth.final.mm39.decomposed.vcf.gz | while read SAMPLE; do
    echo $SAMPLE

    (seq 1 19; echo X; echo Y; echo M) | while read i; do
        VCF=chr$i.ref+pan.p95*/chr*.ref+pan.fa.gz.*.smooth.final.mm39.decomposed.vcf.gz
        echo $VCF
        bcftools view -a -s ${SAMPLE} -Ou ${VCF} | \
            bcftools norm -f ${REFERENCES_FASTA_GZ} -c s -m - -Ou | \
            bcftools view -e 'GT="ref" | GT~"\."' -f 'PASS,.' -Ou | \
            bcftools sort -T /scratch/bcftools-sort.XXXXXX -Ou | \
            bcftools norm -d exact -Oz | \
            bcftools stats | grep '^SN' | cut -f 4 | paste - - - - - - - - - | awk -v OFS='\t' -v sample=$SAMPLE -v chr=chr$i '{print(sample,chr,$0)}' >> variant.stats.tsv
    done
done
```
