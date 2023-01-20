# Pangenome

Variables:

```shell
DIR_BASE=/lizardfs/guarracino/keane_mouse_pangenome
RUN_PGGB=/home/guarracino/tools/pggb/pggb-288a395abf4a9f4755375633093f8ac3af59a081
RUN_VCFBUB=/home/guarracino/tools/vcfbub/target/release/vcfbub-26a1f0cb216a423f8547c4ad0e0ce38cb9d324b9
RUN_ODGI=/home/guarracino/tools/odgi/bin/odgi-f483f9ed5a514a531fbd64833d49cd931ea59943
RUN_VCFWAVE=/gnu/store/hkbkw85gjvsyqvx9vv9bw0ynmad989ag-vcflib-1.0.3+a36dbe9-11/bin/vcfwave
```

## Pangenome building

Run `pggb`:

```shell
mkdir -p $DIR_BASE/graphs/
cd $DIR_BASE/graphs/

( seq 1 19; echo X; echo Y; echo M ) | while read i; do
    sbatch -p workers -c 48 --job-name mice --wrap "hostname; $RUN_PGGB -i $DIR_BASE/partitioning/chr$i.ref+pan.fa.gz -p 95 -s 50000 -n 18 -o $DIR_BASE/graphs/chr$i.ref+pan.p95 -V mm39:# -D /scratch;"
done
```

Compress VCF files:

```shell
cd $DIR_BASE/graphs/

ls */*.vcf | while read v; do
    echo $v
    bgzip $v -@ 48
done
```

Decompose VCF files:

```shell
ls */*.vcf.gz | while read v; do
    echo $v;
    prefix=${v%.vcf.gz};

    sbatch -p workers -w octopus03,octopus11 -c 1 --job-name mice --wrap "hostname; $RUN_VCFBUB -l 0 -a 100000 --input $v | $RUN_VCFWAVE -I 1000 -t 1 > $prefix.decomposed.tmp.vcf"
done
```

Finish decomposition:

```shell
ls */*.vcf.gz | while read v; do
    echo $v;
    prefix=${v%.vcf.gz};
    bcftools annotate -x INFO/TYPE $prefix.decomposed.tmp.vcf  | awk '$5 != "."' | bgzip -@ 48 -c > $prefix.decomposed.vcf.gz
    #rm $prefix.decomposed.tmp.vcf
done
```

Merge all VCF files:

```shell
bcftools concat $DIR_BASE/graphs/chr*.ref+pan.p95/*.decomposed.vcf.gz | bcftools sort -T /scratch/bcftools.XXXXX | bgzip -@ 24 -c > chrALL.ref+pan.p95.vcf.gz
```

Squeeze graphs:

```shell
ls $DIR_BASE/graphs/*p95/*.og | sort -k 1,1 -V > $DIR_BASE/graphs/graphs_to_squeeze.p95.txt

sbatch -p headnode -c 1 --job-name mice --wrap "$RUN_ODGI squeeze -f $DIR_BASE/graphs/graphs_to_squeeze.p95.txt -o - -t 1 -P | $RUN_ODGI view -i - -g > $DIR_BASE/graphs/chrALL.ref+pan.p95.gfa"
zstd chrALL.ref+pan.p95.gfa
```