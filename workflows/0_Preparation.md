# Preparation

## PGGB

```shell
mkdir -p ~/tools/
cd ~/tools/

git clone --recursive https://github.com/waveygang/wfmash
cd wfmash
git checkout master
git pull
git checkout cb0ce952a9bec3f2c8c78b98679375e5275e05db
git submodule update --init --recursive
cmake -H. -DCMAKE_BUILD_TYPE=Release -Bbuild && cmake --build build -- -j $(nproc)
mv build/bin/wfmash build/bin/wfmash-cb0ce952a9bec3f2c8c78b98679375e5275e05db
cd ..

clone --recursive https://github.com/ekg/seqwish
cd seqwish
git checkout master
git pull
git checkout f362f6f5ea89dbb6a0072a8b8ba215e663301d33
git submodule update --init --recursive
cmake -H. -DCMAKE_BUILD_TYPE=Release -DEXTRA_FLAGS='-march=native' -Bbuild && cmake --build build -- -j $(nproc)
mv bin/seqwish bin/seqwish-f362f6f5ea89dbb6a0072a8b8ba215e663301d33
cd ..

git clone --recursive https://github.com/pangenome/smoothxg
cd smoothxg
git checkout master
git pull
git checkout c12f2d2685e566fe04868fd4749e544eb5a6bc37
git submodule update --init --recursive
cmake -H. -DCMAKE_BUILD_TYPE=Release -Bbuild && cmake --build build -- -j $(nproc)
mv bin/smoothxg bin/smoothxg-c12f2d2685e566fe04868fd4749e544eb5a6bc37
cd ..

git clone --recursive https://github.com/pangenome/odgi.git
cd odgi
git checkout master
git pull
git checkout f483f9ed5a514a531fbd64833d49cd931ea59943
git submodule update --init --recursive
cmake -H. -Bbuild && cmake --build build -- -j 48
mv bin/odgi bin/odgi-f483f9ed5a514a531fbd64833d49cd931ea59943
cd ..

git clone --recursive https://github.com/pangenome/pggb.git
cd pggb
git checkout master
git pull
git checkout 288a395abf4a9f4755375633093f8ac3af59a081
sed 's,"$fmt" wfmash,"$fmt" ~/tools/wfmash/build/bin/wfmash-cb0ce952a9bec3f2c8c78b98679375e5275e05db,g' pggb -i
sed 's,"$fmt" seqwish,"$fmt" ~/tools/seqwish/bin/seqwish-f362f6f5ea89dbb6a0072a8b8ba215e663301d33,g' pggb -i
sed 's,"$fmt" smoothxg,"$fmt" ~/tools/smoothxg/bin/smoothxg-c12f2d2685e566fe04868fd4749e544eb5a6bc37,g' pggb -i
sed 's,"$fmt" odgi,"$fmt" ~/tools/odgi/bin/odgi-f483f9ed5a514a531fbd64833d49cd931ea59943,g' pggb -i
mv pggb pggb-288a395abf4a9f4755375633093f8ac3af59a081
cd ..
```

## Assemblies

```shell
mkdir -p /lizardfs/guarracino/keane_mouse_pangenome/assemblies

# TODO: to upload the REL-2205-Assembly.urls.txt file

sbatch -p workers -c 48 --wrap "cd /lizardfs/guarracino/keane_mouse_pangenome/assemblies; (echo https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz; cat ../data/REL-2205-Assembly.urls.txt) | parallel -j 4 'wget -q {} && echo got {}'"

# Apply PanSN-spec to the mm39 reference. The assemblies already follow PanSN-spec.
zcat mm39.fa.gz | sed 's/>chr/>mm39#1#chr/g' | bgzip -@ 48 -c > mm39.fasta.gz
samtools faidx mm39.fasta.gz
rm mm39.fa.gz
```

### Partitioning by chromosome

```shell
mkdir -p /lizardfs/guarracino/keane_mouse_pangenome/partitioning
cd /lizardfs/guarracino/keane_mouse_pangenome/partitioning

RUN_WFMASH=/home/guarracino/tools/wfmash/build/bin/wfmash-cb0ce952a9bec3f2c8c78b98679375e5275e05db
REFERENCES_FASTA_GZ=/lizardfs/guarracino/keane_mouse_pangenome/assemblies/mm39.fasta.gz

ls /lizardfs/guarracino/keane_mouse_pangenome/assemblies/*fasta.gz | grep mm39 -v | grep unplaced | while read FASTA_GZ; do
  HAPLOTYPE=$(basename $FASTA_GZ .fasta.gz | cut -f 1,2 -d '.');
  echo $HAPLOTYPE

  PAF=/lizardfs/guarracino/keane_mouse_pangenome/partitioning/$HAPLOTYPE.vs.ref.paf
  sbatch -p workers -c 6 --wrap "$RUN_WFMASH -t 6 -m -N -s 50k -l 150k -p 90 -H 0.001 $REFERENCES_FASTA_GZ $FASTA_GZ > $PAF"
done
```

Collect unmapped contigs and remap them in split mode:

```shell
ls /lizardfs/guarracino/keane_mouse_pangenome/assemblies/*fasta.gz | grep mm39 -v | grep unplaced | while read FASTA_GZ; do
  HAPLOTYPE=$(basename $FASTA_GZ .fasta.gz | cut -f 1,2 -d '.');
  echo $HAPLOTYPE

  UNALIGNED=/lizardfs/guarracino/keane_mouse_pangenome/partitioning/$HAPLOTYPE.unaligned
  
  PAF=/lizardfs/guarracino/keane_mouse_pangenome/partitioning/$HAPLOTYPE.vs.ref.paf
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
( seq 1 19; echo X; echo Y; echo M ) | while read i; do cat *paf | grep -P "mm39#1#chr$i\t" | cut -f 1 | sort | uniq | awk '{print($0"$")}' > chr$i.contigs.txt; done

( seq 1 19; echo X; echo Y; echo M ) | while read i; do
    echo chr$i

    rm chr$i.ref+pan.fa*
    samtools faidx $REFERENCES_FASTA_GZ "mm39#1#chr$i" > chr$i.ref+pan.fa
    ls ../assemblies/*.fasta.gz | grep unplaced | while read PATH_FASTA; do
      echo $PATH_FASTA

      samtools faidx $PATH_FASTA $(cut -f 1 $PATH_FASTA.fai | grep -f chr$i.contigs.txt) >> chr$i.ref+pan.fa
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
JF1_MsJ#1#chr2_unloc_1      511612   2357066079  60  61
JF1_MsJ#1#chr2_unloc_2      506306   2357586242  60  61
NZO_HlLtJ#1#chr2_unloc_2    1576622  2962088283  60  61
chr2
NOD_ShiLtJ#1#chrX_unloc_4   67075    2524067938  60  61
chr3
chr4
DBA_2J#1#chrX_unloc_11      61397    1377888279  60  61
NOD_ShiLtJ#1#chrX_unloc_16  16818    2156953427  60  61
chr5
chr6
WSB_EiJ#1#chr2_unloc_4      24119    2536728415  60  61
chr7
NZO_HlLtJ#1#chrX_unloc_7    37713    1939635765  60  61
chr8
NOD_ShiLtJ#1#chrX_unloc_19  14419    1681775376  60  61
chr9
chr10
chr11
chr12
NOD_ShiLtJ#1#chr7_unloc_3   32211    1529595179  60  61
chr13
chr14
NOD_ShiLtJ#1#chrX_unloc_8   38559    1564399476  60  61
chr15
chr16
NOD_ShiLtJ#1#chrX_unloc_11  28376    1260054872  60  61
NOD_ShiLtJ#1#chrX_unloc_17  15663    1260083749  60  61
chr17
NOD_ShiLtJ#1#chrX_unloc_20  13272    1220929652  60  61
chr18
chr19
chrX
JF1_MsJ#1#chr9_unloc_1      448480   1913034795  60  61
NOD_ShiLtJ#1#chr4_unloc_3   11884    2142058411  60  61
WSB_EiJ#1#chr2_unloc_2      67262    2774246432  60  61
chrY
LP_J#1#chr9_unloc_3         710282   96910817    60  61
LP_J#1#chr18_unloc_5        317819   97632960    60  61
WSB_EiJ#1#chr4_unloc_1      169623   111149347   60  61
chrM
```
