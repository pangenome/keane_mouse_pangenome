# Preparation

Variables:

```shell
DIR_BASE=/lizardfs/guarracino/keane_mouse_pangenome
DATASETS=/lizardfs/guarracino/tools/ncbi_datasets/datasets
PGGB=/home/guarracino/tools/pggb/pggb-b26786adce800bbf2664e41e4c95d57cb3ee048e
ODGI=/home/guarracino/tools/odgi/bin/odgi-fbdb4d23a512cdcbede21c7a453f62415a2ca7e6
```

## Tools

```shell
mkdir -p ~/tools/
cd ~/tools/

git clone --recursive https://github.com/waveygang/wfmash
cd wfmash
git checkout master && git pull && git submodule update --init --recursive
git checkout ecdfb95e660eda46c3c51312a7db54394270271a
cmake -H. -Bbuild && cmake --build build -- -j 48
cp build/bin/wfmash build/bin/wfmash-ecdfb95e660eda46c3c51312a7db54394270271a
cd ..

git clone --recursive https://github.com/ekg/seqwish.git
cd seqwish
git checkout master && git pull && git submodule update --init --recursive
git checkout f44b402f0c2e02988d431d9b2e5eba9727cf93a9
rm -rf build/
cmake -H. -Bbuild && cmake --build build -- -j 48
cp bin/seqwish bin/seqwish-f44b402f0c2e02988d431d9b2e5eba9727cf93a9
cd ..

git clone --recursive https://github.com/pangenome/smoothxg.git
cd smoothxg
git checkout master && git pull && git submodule update --init --recursive
git checkout aaa0b283e13ca57c4e6e4e67a03451925f5342f1
rm -rf build/
cmake -H. -Bbuild && cmake --build build -- -j 48
cp bin/smoothxg bin/smoothxg-aaa0b283e13ca57c4e6e4e67a03451925f5342f1
cd ..

git clone --recursive https://github.com/pangenome/odgi.git
cd odgi
git checkout master && git pull && git submodule update --init --recursive
git checkout fbdb4d23a512cdcbede21c7a453f62415a2ca7e6
rm -rf build/
cmake -H. -Bbuild && cmake --build build -- -j 48
cp bin/odgi bin/odgi-fbdb4d23a512cdcbede21c7a453f62415a2ca7e6
cd ..

clone https://github.com/marschall-lab/GFAffix.git
cd GFAffix
git checkout main && git pull && git submodule update --init --recursive
git checkout d630eb7d9827340f5f292e57cb3cb5e31e6f86f0
env -i bash -c 'PATH=:/usr/local/bin:/usr/bin:/bin ~/.cargo/bin/cargo build --release'
cp target/release/gfaffix target/release/gfaffix-d630eb7d9827340f5f292e57cb3cb5e31e6f86f0
cd ..

git clone --recursive https://github.com/pangenome/pggb.git
cd pggb
git checkout master && git pull && git submodule update --init --recursive
git checkout b26786adce800bbf2664e41e4c95d57cb3ee048e
cp pggb pggb-x
sed 's,"$fmt" wfmash,"$fmt" ~/tools/wfmash/build/bin/wfmash-ecdfb95e660eda46c3c51312a7db54394270271a,g' pggb-x -i
sed 's,$(wfmash,$(~/tools/wfmash/build/bin/wfmash-ecdfb95e660eda46c3c51312a7db54394270271a,g' pggb-x -i
sed 's,"$fmt" seqwish,"$fmt" ~/tools/seqwish/bin/seqwish-f44b402f0c2e02988d431d9b2e5eba9727cf93a9,g' pggb-x -i
sed 's,"$fmt" smoothxg,"$fmt" ~/tools/smoothxg/bin/smoothxg-aaa0b283e13ca57c4e6e4e67a03451925f5342f1,g' pggb-x -i
sed 's,"$fmt" odgi,"$fmt" ~/tools/odgi/bin/odgi-fbdb4d23a512cdcbede21c7a453f62415a2ca7e6,g' pggb-x -i
sed 's,"$fmt" gfaffix,"$fmt" ~/tools/GFAffix/target/release/gfaffix-d630eb7d9827340f5f292e57cb3cb5e31e6f86f0,g' pggb-x -i
mv pggb-x pggb-b26786adce800bbf2664e41e4c95d57cb3ee048e
cp partition-before-pggb partition-before-pggb-x
sed 's,"$fmt" wfmash,"$fmt" ~/tools/wfmash/build/bin/wfmash-ecdfb95e660eda46c3c51312a7db54394270271a,g' partition-before-pggb-x -i
sed 's,$(wfmash,$(~/tools/wfmash/build/bin/wfmash-ecdfb95e660eda46c3c51312a7db54394270271a,g' partition-before-pggb-x -i
sed 's,"$fmt" seqwish,"$fmt" ~/tools/seqwish/bin/seqwish-f44b402f0c2e02988d431d9b2e5eba9727cf93a9,g' partition-before-pggb-x -i
sed 's,"$fmt" smoothxg,"$fmt" ~/tools/smoothxg/bin/smoothxg-aaa0b283e13ca57c4e6e4e67a03451925f5342f1,g' partition-before-pggb-x -i
sed 's,"$fmt" odgi,"$fmt" ~/tools/odgi/bin/odgi-fbdb4d23a512cdcbede21c7a453f62415a2ca7e6,g' partition-before-pggb-x -i
sed 's,"$fmt" gfaffix,"$fmt" ~/tools/GFAffix/target/release/gfaffix-d630eb7d9827340f5f292e57cb3cb5e31e6f86f0,g' partition-before-pggb-x -i
mv partition-before-pggb-x partition-before-pggb-b26786adce800bbf2664e41e4c95d57cb3ee048e
cd ..

wget -c https://github.com/vgteam/vg/releases/download/v1.52.0/vg
chmod +x vg

git clone --recursice https://github.com/pangenome/vcfbub.git
cd vcfbub
git checkout d92ae50da4926612f9bcf5d3ea8506a78a348734 # version 0.1.0
cargo build --release

/gnu/store/79r54wk4p3705dk89jg9hidyvf4754jp-vcflib-1.0.3+fdcdaad-10/bin/vcfwave
```

## Assemblies

Download:

```shell
mkdir -p $DIR_BASE/assemblies

sed '1d' $DIR_BASE/data/assemblies.tsv | while read STRAIN ACC; do
  STRAIN=$(echo $STRAIN | tr '/' '_')
  echo $STRAIN $ACC
  $DATASETS download genome accession --inputfile <(echo $ACC) --filename $STRAIN.zip
done
ls *zip | while read f; do unzip -o $f; done
mv ncbi_dataset/*/*/*fna .
rm ncbi_dataset -rf README.md *zip

wget -c https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
```

PanSN-spec:

```shell
ls *fna | while read f; do
  ACC=$(basename $f .fna | cut -f 1,2 -d '_')
  STRAIN=$(echo $ACC | sed -f <(cat $DIR_BASE/data/assemblies.tsv | tr '/' '_' | awk 'BEGIN{FS=OFS="\t"} {print "s/" $2 "/" $1 "/g"}'))
  echo $STRAIN $ACC
  sed "s/>/>$STRAIN#1#/g" $f | cut -f 1 -d ' ' | bgzip -@ 48 -l 9 > $STRAIN.fa.gz && samtools faidx $STRAIN.fa.gz
done

zcat mm39.fa.gz | sed 's/>chr/>mm39#1#chr/g' | bgzip -@ 48 -c > mm39.fasta.gz
samtools faidx mm39.fasta.gz
rm mm39.fa.gz
```
