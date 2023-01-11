# Pangenome building

Run `pggb`:

```shell
RUN_PGGB=/home/guarracino/tools/pggb/pggb-288a395abf4a9f4755375633093f8ac3af59a081

mkdir -p /lizardfs/guarracino/keane_mouse_pangenome/graphs/

( seq 1 19; echo X; echo Y; echo M ) | while read i; do
    sbatch -p workers -c 48 --job-name mice --wrap "hostname; $RUN_PGGB -i /lizardfs/guarracino/keane_mouse_pangenome/partitioning/chr$i.ref+pan.fa.gz -p 95 -s 50000 -n 18 -o /lizardfs/guarracino/keane_mouse_pangenome/graphs/chr$i.ref+pan.fa -V mm39:# -D /scratch;"
done
```
