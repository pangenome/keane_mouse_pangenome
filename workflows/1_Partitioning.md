# Partitioning

Variables:

```shell
DIR_BASE=/lizardfs/guarracino/keane_mouse_pangenome
REF_FASTA=$DIR_BASE/assemblies/mm39.fasta.gz

PGGB=/home/guarracino/tools/pggb/pggb-13482bd06359a7ad8e3d3e0dd6eb6d9399f26046
ODGI=/home/guarracino/tools/odgi/bin/odgi-861b1c04f5622c5bb916a161c1abe812c213f1a5
WFMASH=/home/guarracino/tools/wfmash/build/bin/wfmash-0b191bb84ffdfd257354c1aa82a7f1e13dc536d0 
```

conda create --prefix /lizardfs/guarracino/condatools/liftoff/1.6.3/ -c conda-forge -c bioconda liftoff=1.6.3 -y
conda activate /lizardfs/guarracino/condatools/liftoff/1.6.3
liftoff assemblies/129S1_SvImJ.fa assemblies/mm39.fasta -g data/loci_of_interest.gff
conda deactivate

```shell
mkdir -p $DIR_BASE/pangenome/loci_of_interest
cd $DIR_BASE/pangenome/loci_of_interest

ls $DIR_BASE/assemblies/*fa.gz | while read FASTA; do
    SAMPLE=$(basename $FASTA .fa.gz)
    echo $SAMPLE

    sbatch -c 24 -p workers --job-name wfmash-$SAMPLE --wrap "$WFMASH $REF_FASTA $FASTA -t 24 -p 90 -m > $DIR_BASE/pangenome/loci_of_interest/$SAMPLE.vs.ref.mappings.paf"
done

cd $DIR_BASE/pangenome/loci_of_interest
cat $DIR_BASE/data/loci_of_interest.bed | while read -r chrom start end; do
    name=${chrom}_${start}_${end}
    echo $chrom $start $end $name

    bedtools intersect \
        -a <(cat $DIR_BASE/pangenome/loci_of_interest/*s.paf | awk -v OFS='\t' '{print($6,$8,$9,$1,".", $5)}') \
        -b <(bedtools slop -i <(echo "$chrom\t$start\t$end") \
        -g <(cut -f 1,2 $DIR_BASE/assemblies/mm39.fasta.gz.fai ) -b 100) | cut -f 4 | sort | uniq > $DIR_BASE/pangenome/loci_of_interest/$name.contigs.txt


    samtools faidx $REF_FASTA $chrom >> $DIR_BASE/pangenome/loci_of_interest/$name.contigs.fa
    ls $DIR_BASE/assemblies/*fa.gz | while read FASTA; do
        SAMPLE=$(basename $FASTA .fa.gz)
        echo $SAMPLE

        samtools faidx $FASTA $(grep -Ff $DIR_BASE/pangenome/loci_of_interest/$name.contigs.txt $FASTA.fai | cut -f 1) >> $DIR_BASE/pangenome/loci_of_interest/$name.contigs.fa
    done

    bgzip -f -@ 48 -l 9 $DIR_BASE/pangenome/loci_of_interest/$name.contigs.fa && samtools faidx $DIR_BASE/pangenome/loci_of_interest/$name.contigs.fa.gz
    rm $DIR_BASE/pangenome/loci_of_interest/$name.contigs.txt

    sbatch -c 48 -p workers --job-name pggb-$name --wrap "$PGGB -i $DIR_BASE/pangenome/loci_of_interest/$name.contigs.fa.gz -o $DIR_BASE/pangenome/loci_of_interest/pggb.$name -s 10k -p 98 -t 48 -D /scratch/$name"
done

cd $DIR_BASE/pangenome/loci_of_interest
cat $DIR_BASE/data/loci_of_interest.bed | while read -r chrom start end; do
    name=${chrom}_${start}_${end}
    echo $chrom $start $end $name
    
    mkdir $DIR_BASE/pangenome/loci_of_interest/pggb.$name/extractions
    cd $DIR_BASE/pangenome/loci_of_interest/pggb.$name/extractions
    
    NAME=$name
    
    if [[ ! -f $NAME.og ]]; then
        odgi extract \
            -i ../$name.*.smooth.final.og \
            -o - \
            -r $chrom:$start-$end \
            -t 48 -P | odgi sort -i - -o - -O -Y -t 48 -x 500 -P --temp-dir /scratch | \
            odgi flip -i - -o $NAME.og
    fi

    odgi view -i $NAME.og -g | sed 's/_inv//g' > $NAME.gfa
    #odgi view -i $NAME.og -g > $NAME.gfa

    odgi viz -i $NAME.gfa -o $NAME.1.N-bases.png -N --max-num-of-characters 64 -x 2000
    odgi viz -i $NAME.gfa -o $NAME.2.orientation.png -z --max-num-of-characters 64 -x 2000
    odgi viz -i $NAME.gfa -o $NAME.3.depth.png -m --max-num-of-characters 64 -x 2000
    odgi viz -i $NAME.gfa -o $NAME.4.by-sample.png -s '#' --max-num-of-characters 64 -x 2000

    odgi layout -i $NAME.gfa -o $NAME.lay -t 48 -x 500 -G 20 -I 10000 -l 1000 --temp-dir /scratch -P
    odgi draw -i $NAME.gfa -c $NAME.lay -p $NAME.5.2D.png

    #guix install bandage
    Bandage image $NAME.gfa $NAME.6.2D.bandage.png --height 1000 --width 1000

    cd ../..
done

# Run PGGB locally
cd $DIR_BASE/pangenome/loci_of_interest
cat $DIR_BASE/data/loci_of_interest.bed | grep chr1 | while read -r chrom start end; do
    name=${chrom}_${start}_${end}
    echo $chrom $start $end $name
    
    cd $DIR_BASE/pangenome/loci_of_interest/pggb.$name/extractions

    NAME=$name

    if [[ ! -f $NAME.fa.gz ]]; then
        # In the GFA the name are without the `_inv` suffix
        odgi paths -i $NAME.gfa -f | bgzip -@ 48 -l 9 > $NAME.fa.gz
        samtools faidx $NAME.fa.gz
    fi

    #$PGGB -i $NAME.fa.gz -o pggb.$NAME.1 -D /scratch/$name
    $PGGB -i $NAME.fa.gz -o pggb.$NAME.2 -p 95 -D /scratch/$name
    $PGGB -i $NAME.fa.gz -o pggb.$NAME.3 -p 95 -k 0 -D /scratch/$name
    $PGGB -i $NAME.fa.gz -o pggb.$NAME.4 -p 95 -s 1k -k 47 -D /scratch/$name
    $PGGB -i $NAME.fa.gz -o pggb.$NAME.5 -p 95 -s 1k -k 47 -G 1400,1800,2200 -D /scratch/$name
    $PGGB -i $NAME.fa.gz -o pggb.$NAME.6 -p 80 -D /scratch/$name
    $PGGB -i $NAME.fa.gz -o pggb.$NAME.6 -p 80 -s 1k -k 47 -D /scratch/$name
    $PGGB -i $NAME.fa.gz -o pggb.$NAME.5 -p 95 -s 1k -k 47 -G 1400,1800,2200 -D /scratch/$name
    $PGGB -i $NAME.fa.gz -o pggb.$NAME.7 -p 80 -s 10k -D /scratch/$name

    ls pggb.$NAME.*/*.og | while read GRAPH; do
        echo $GRAPH
        $ODGI sort -i $GRAPH -o - -p gYs -x 300 -t 48 --temp-dir /scratch/ -H <($ODGI paths -i $GRAPH -L | grep mm39) | $ODGI viz -i - -o $(echo $GRAPH | sed 's/.og//g').sort-by-ref.png -N -x 2000 -a 10
    done

    cd ../..
done



odgi similarity -i chr1_146713676_146736261_Fah.gfa -D '#' -p 1 > chr1_146713676_146736261_Fah.similarity.tsv



NAME=chr1_146713676_146736261_Fah
cd /lizardfs/guarracino/rat/proteome/specific_loci/pggb.Fah/extractions
odgi viz -i $NAME.og -o $NAME.4.by-sample.wide.png -s '#' -p $NAME.path_names.txt -x 50000


vg convert -g $NAME.gfa -x > $NAME.xg
vg viz -x $NAME.xg --out $NAME.svg -X 500 -C
```
