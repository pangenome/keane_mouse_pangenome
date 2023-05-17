# Short reads

## Download

```shell
# Obtain FASTQ urls from ENA accessions
rm accessions.urls.tsv; cut -f 2 accessions.tsv | sed '1d' | while read ACC; do echo $ACC; wget -c https://www.ebi.ac.uk/ena/portal/api/filereport\?accession\="${ACC}"\&result\=read_run\&fastq_ftp,fastq_md5 -O - | sed '1d' | cut -f 2,4 | awk '{split($1,a,";"); split($2,b,";"); print(a[1],b[1],"\n"a[2],b[2],"\n"a[3],b[3])}' | grep '_' | awk -v acc=$ACC -v OFS='\t' '{print(acc,$1,$2)}' >> accessions.urls.tsv; done

cut -f 2 accessions.urls.tsv | parallel -j 2 'wget -q {} && echo got {}'
```
