# Pangenome building

Run `pggb`:

```shell
DIR_BASE=/lizardfs/guarracino/keane_mouse_pangenome
RUN_PGGB=/home/guarracino/tools/pggb/pggb-288a395abf4a9f4755375633093f8ac3af59a081

mkdir -p $DIR_BASE/graphs/

( seq 1 19; echo X; echo Y; echo M ) | while read i; do
    sbatch -p workers -c 48 --job-name mice --wrap "hostname; $RUN_PGGB -i $DIR_BASE/partitioning/chr$i.ref+pan.fa.gz -p 95 -s 50000 -n 18 -o $DIR_BASE/graphs/chr$i.ref+pan.fa -V mm39:# -D /scratch;"
done
```
