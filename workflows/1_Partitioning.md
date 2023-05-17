# Partitioning


Variables:

```shell
DIR_BASE=/lizardfs/guarracino/keane_mouse_pangenome
RUN_WFMASH=/home/guarracino/tools/wfmash/build/bin/wfmash-8ba3c53f327731ca515abd1ef32179f15acb9732
REFERENCES_FASTA_GZ=$DIR_BASE/assemblies/mm39.fasta.gz
```

## Partitioning by chromosome

```shell
mkdir -p $DIR_BASE/partitioning
cd $DIR_BASE/partitioning



# Unplaced files have all contigs
ls $DIR_BASE/assemblies/*fasta.gz | grep mm39 -v | grep unplaced | while read FASTA_GZ; do
  HAPLOTYPE=$(basename $FASTA_GZ .fasta.gz | cut -f 1,2 -d '.');
  echo $HAPLOTYPE

  PAF=$DIR_BASE/partitioning/$HAPLOTYPE.vs.ref.paf
  sbatch -p workers -c 6 --wrap "$RUN_WFMASH -t 6 -m -N -s 50k -l 150k -p 90 -H 0.001 $REFERENCES_FASTA_GZ $FASTA_GZ > $PAF"
done
```

Collect unmapped contigs and remap them in split mode:

```shell
ls $DIR_BASE/assemblies/*fasta.gz | grep mm39 -v | grep unplaced | while read FASTA_GZ; do
  HAPLOTYPE=$(basename $FASTA_GZ .fasta.gz | cut -f 1,2 -d '.');
  echo $HAPLOTYPE

  UNALIGNED=$DIR_BASE/partitioning/$HAPLOTYPE.unaligned
  
  PAF=$DIR_BASE/partitioning/$HAPLOTYPE.vs.ref.paf
  comm -23 <(cut -f 1 $FASTA_GZ.fai | sort) <(cut -f 1 $PAF | sort) > $UNALIGNED.txt
  if [[ $(wc -l $UNALIGNED.txt | cut -f 1 -d\ ) != 0 ]];
  then
    samtools faidx $FASTA_GZ $(tr '\n' ' ' < $UNALIGNED.txt) > $UNALIGNED.fa
    samtools faidx $UNALIGNED.fa
    sbatch -p workers -c 6 --wrap "$RUN_WFMASH -t 6 -m -s 50k -l 150k -p 90 -H 0.001 $REFERENCES_FASTA_GZ $UNALIGNED.fa > $UNALIGNED.split.vs.ref.paf"
  fi
done
```

Collect our best mapping for each of our attempted split rescues:

```shell
ls *.unaligned.split.vs.ref.paf | while read PAF; do
  cat $PAF | awk -v OFS='\t' '{ print $1,$11,$0 }' | sort -n -r -k 1,2 | \
    awk -v OFS='\t' '$1 != last { print($3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15); last = $1; }'
done > rescues.paf
```

Collect partitioned contigs:

```shell
( seq 1 19; echo X; echo Y; echo M ) | while read i; do cat *paf | grep -P "\tmm39#1#chr$i\t" | cut -f 1 | sort | uniq | awk '{print("^"$0"$")}' > chr$i.contigs.txt; done
```

Check if contigs are in multiple chromosomes. If any, align them (not just map) and take the best alignment target for each one:

```shell
sort chr*.contigs.txt | uniq -c | sort -k 1,1n | awk '$1 > 1' | tr -s ' ' | cut -f 3 -d ' ' > chrALL.contigs_in_multuple_chromosomes.txt

# Prepare the FASTA with all multi-chromosome contigs
( seq 1 19; echo X; echo Y; echo M ) | while read i; do
    echo chr$i

    rm chrALL.contigs_in_multuple_chromosomes.fa*
    ls ../assemblies/*.fasta.gz | grep unplaced | while read PATH_FASTA; do
      echo $PATH_FASTA

      samtools faidx $PATH_FASTA $(cut -f 1 $PATH_FASTA.fai | grep -f chrALL.contigs_in_multuple_chromosomes.txt) >> chrALL.contigs_in_multuple_chromosomes.fa
    done
    bgzip -@ 48 chrALL.contigs_in_multuple_chromosomes.fa
    samtools faidx chrALL.contigs_in_multuple_chromosomes.fa.gz
done

# Align all multi-chromosome contigs
$RUN_WFMASH -t 48 -N -s 50k -l 150k -p 90 -H 0.001 $REFERENCES_FASTA_GZ chrALL.contigs_in_multuple_chromosomes.fa.gz > chrALL.contigs_in_multuple_chromosomes.vs.ref.paf

# Compute alignment_length * gap-compressed-identity, sort by that and take the target of the best alignment for each query (each contig)
sed 's/gi:f://g' chrALL.contigs_in_multuple_chromosomes.vs.ref.paf | awk -v OFS= '{print($0,$11*$13)}' | sort -k 1,1 -k 21,21nr | awk -v OFS='\t' '$1 != last { print($0,$11*$13); last = $1; }' > chrALL.contigs_in_multuple_chromosomes.vs.ref.best_alignment.paf
```

Collect partitioned contigs again:


```shell
( seq 1 19; echo X; echo Y; echo M ) | while read i; do
  cat *paf | grep -P "\tmm39#1#chr$i\t" | cut -f 1 | sort | uniq | grep -f chrALL.contigs_in_multuple_chromosomes.txt -v | awk '{print("^"$0"$")}' > chr$i.contigs.fixed.txt;

  grep -P "\tmm39#1#chr$i\t" chrALL.contigs_in_multuple_chromosomes.vs.ref.best_alignment.paf | cut -f 1 | awk '{print("^"$0"$")}' >> chr$i.contigs.fixed.txt
done

# It has to be 0!
sort chr*.contigs.fixed.txt | uniq -c | sort -k 1,1n | awk '$1 > 1' | wc -l

( seq 1 19; echo X; echo Y; echo M ) | while read i; do
    echo chr$i

    rm chr$i.ref+pan.fa*
    samtools faidx $REFERENCES_FASTA_GZ "mm39#1#chr$i" > chr$i.ref+pan.fa
    ls ../assemblies/*.fasta.gz | grep unplaced | while read PATH_FASTA; do
      echo $PATH_FASTA

      samtools faidx $PATH_FASTA $(cut -f 1 $PATH_FASTA.fai | grep -f chr$i.contigs.fixed.txt) >> chr$i.ref+pan.fa
    done
    bgzip -@ 48 chr$i.ref+pan.fa
    samtools faidx chr$i.ref+pan.fa.gz
done
```

Note that a few contigs are partitioned into different chromosomes than specified in their names.

```shell
( seq 1 19; echo X; echo Y; echo M ) | while read i; do
    echo chr$i
    
    grep -f <(cat chr$i.contigs.txt | grep chr | grep chr$i -v | sed 's/.$//g') chr$i.ref+pan.fa.gz.fai
done | column -t
```

Output:

```bash
chr1
JF1_MsJ#1#chr2_unloc_1      511612   2358240973  60  61
JF1_MsJ#1#chr2_unloc_2      506306   2358761136  60  61
NZO_HlLtJ#1#chr2_unloc_2    1576622  2964375130  60  61
chr2
NOD_ShiLtJ#1#chrX_unloc_4   67075    2523559164  60  61
chr3
chr4
DBA_2J#1#chrX_unloc_11      61397    1377461601  60  61
NOD_ShiLtJ#1#chrX_unloc_16  16818    2156782899  60  61
NOD_ShiLtJ#1#chr5_unloc_5   39975    2156800025  60  61
chr5
chr6
WSB_EiJ#1#chr2_unloc_4      24119    2535927868  60  61
chr7
NZO_HlLtJ#1#chrX_unloc_7    37713    1938935089  60  61
chr8
JF1_MsJ#1#chrX_unloc_12     754125   1406147958  60  61
NOD_ShiLtJ#1#chrX_unloc_19  14419    1679865717  60  61
chr9
chr10
chr11
LP_J#1#chr2_unloc_11        586699   1462228748  60  61
LP_J#1#chr2_unloc_3         1413569  1462825247  60  61
LP_J#1#chr5_unloc_4         541241   1464262397  60  61
LP_J#1#chr9_unloc_4         662196   1464812680  60  61
LP_J#1#chr12_unloc_4        902544   1465485935  60  61
LP_J#1#chr12_unloc_1        3804029  1466403544  60  61
NOD_ShiLtJ#1#chrX_unloc_14  20678    1597906455  60  61
chr12
NOD_ShiLtJ#1#chr7_unloc_3   32211    1528602130  60  61
chr13
chr14
NOD_ShiLtJ#1#chrX_unloc_8   38559    1564140118  60  61
chr15
chr16
DBA_2J#1#chr19_unloc_1      43338    967573023   60  61
NOD_ShiLtJ#1#chrX_unloc_11  28376    1260426462  60  61
NOD_ShiLtJ#1#chrX_unloc_17  15663    1260455339  60  61
chr17
DBA_2J#1#chrX_unloc_12      57835    839023626   60  61
NOD_ShiLtJ#1#chrX_unloc_20  13272    1222962526  60  61
chr18
chr19
chrX
JF1_MsJ#1#chr9_unloc_1      448480   1911635942  60  61
NOD_ShiLtJ#1#chr4_unloc_3   11884    2140564223  60  61
WSB_EiJ#1#chr2_unloc_2      67262    2772619412  60  61
chrY
LP_J#1#chr9_unloc_3         710282   96527849    60  61
LP_J#1#chr18_unloc_5        317819   97249992    60  61
WSB_EiJ#1#chr4_unloc_1      169623   110516964   60  61
chrM
```

Statistics:

```shell
cd $DIR_BASE/partitioning

echo -e "chromosome\tstrain\tnum.contigs" > assembly.statistics.tsv
(seq 1 19; echo X; echo Y; echo M) | while read i; do
  cut -f 1 chr$i.ref+pan.fa.gz.fai -d '#' | uniq -c | awk -v OFS='\t' -v chr=chr$i '{print(chr,$2,$1)}';
done >> assembly.statistics.tsv
```
