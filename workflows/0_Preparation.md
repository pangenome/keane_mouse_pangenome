# Preparation

Variables:

```shell
DIR_BASE=/lizardfs/guarracino/keane_mouse_pangenome
```

## Tools

```shell
mkdir -p ~/tools/
cd ~/tools/

git clone --recursive https://github.com/waveygang/wfmash
cd wfmash
git checkout master && git pull && git submodule update --init --recursive
git checkout 8ba3c53f327731ca515abd1ef32179f15acb9732
cmake -H. -Bbuild && cmake --build build -- -j 48
cp build/bin/wfmash build/bin/wfmash-8ba3c53f327731ca515abd1ef32179f15acb9732
cd ..

git clone --recursive https://github.com/ekg/seqwish.git
cd seqwish
git checkout master && git pull && git submodule update --init --recursive
git checkout d9e7ab59e73258f57875f2a060437735a460475e
rm -rf build/
cmake -H. -Bbuild && cmake --build build -- -j 48
cp bin/seqwish bin/seqwish-d9e7ab59e73258f57875f2a060437735a460475e
cd ..

git clone --recursive https://github.com/pangenome/smoothxg.git
cd smoothxg
git checkout master && git pull && git submodule update --init --recursive
git checkout 956eb75644522bb2b96b4cca44b7bafa9cf02f4a
rm -rf build/
cmake -H. -Bbuild && cmake --build build -- -j 48
cp bin/smoothxg bin/smoothxg-956eb75644522bb2b96b4cca44b7bafa9cf02f4a
cd ..

git clone --recursive https://github.com/pangenome/odgi.git
cd odgi
git checkout master && git pull && git submodule update --init --recursive
git checkout 34f006f31c3f6b35a1eb8d58a4edb1c458583de3
rm -rf build/
cmake -H. -Bbuild && cmake --build build -- -j 48
cp bin/odgi bin/odgi-34f006f31c3f6b35a1eb8d58a4edb1c458583de3
cd ..

git clone --recursive https://github.com/pangenome/pggb.git
cd pggb
git checkout master && git pull && git submodule update --init --recursive
git checkout de5303e24d3e5594a5a2c9bdeb49aba420b24b0c
cp pggb pggb-x
sed 's,"$fmt" wfmash,"$fmt" ~/tools/wfmash/build/bin/wfmash-8ba3c53f327731ca515abd1ef32179f15acb9732,g' pggb-x -i
sed 's,$(wfmash,$(~/tools/wfmash/build/bin/wfmash-8ba3c53f327731ca515abd1ef32179f15acb9732,g' pggb-x -i
sed 's,"$fmt" seqwish,"$fmt" ~/tools/seqwish/bin/seqwish-d9e7ab59e73258f57875f2a060437735a460475e,g' pggb-x -i
sed 's,"$fmt" smoothxg,"$fmt" ~/tools/smoothxg/bin/smoothxg-956eb75644522bb2b96b4cca44b7bafa9cf02f4a,g' pggb-x -i
sed 's,"$fmt" odgi,"$fmt" ~/tools/odgi/bin/odgi-34f006f31c3f6b35a1eb8d58a4edb1c458583de3,g' pggb -i
mv pggb-x pggb-de5303e24d3e5594a5a2c9bdeb49aba420b24b0c
cd ..

git clone https://github.com/pangenome/vcfbub
cd vcfbub
git checkout master
git pull
git checkout 26a1f0cb216a423f8547c4ad0e0ce38cb9d324b9
git submodule update --init --recursive
cargo build --release
mv target/release/vcfbub target/release/vcfbub-26a1f0cb216a423f8547c4ad0e0ce38cb9d324b9
cd ..
```

## Assemblies

```shell
mkdir -p $DIR_BASE/assemblies

sbatch -p workers -c 48 --wrap "cd $DIR_BASE/assemblies; (echo https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz; cat ../data/REL-2205-Assembly.urls.txt) | parallel -j 4 'wget -q {} && echo got {}'"

# Apply PanSN-spec to the mm39 reference. The assemblies already follow PanSN-spec.
zcat mm39.fa.gz | sed 's/>chr/>mm39#1#chr/g' | bgzip -@ 48 -c > mm39.fasta.gz
samtools faidx mm39.fasta.gz
rm mm39.fa.gz
```
