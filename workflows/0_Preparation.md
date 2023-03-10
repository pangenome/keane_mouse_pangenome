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
