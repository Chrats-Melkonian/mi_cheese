#!/usr/bin/env bash

# assuming you have models organized into subfolders as shown below, the following script loops through the media compositions in media_db.tsv and simulates each community/subdfolder in each media present in the file
#models/
#├── milk_aer
#│   ├── CHCC10675.xml
#│   ├── CHCC4895.xml
#│   ├── CHCC5614.xml
#│   ├── CHCC6086.xml
#├── milk_ana
#│   ├── CHCC10675.xml
#│   ├── CHCC4895.xml
#│   ├── CHCC5614.xml
#│   ├── CHCC6086.xml
#├── mm_aer
#│   ├── CHCC10675.xml
#│   ├── CHCC4895.xml
#│   ├── CHCC5614.xml
#│   ├── CHCC6086.xml
#├── mm_ana
#│   ├── CHCC10675.xml
#│   ├── CHCC4895.xml
#│   ├── CHCC5614.xml
#│   ├── CHCC6086.xml
#├── no_gf
#│   ├── CHCC10675.xml
#│   ├── CHCC4895.xml
#│   ├── CHCC5614.xml
#│   ├── CHCC6086.xml

while read modelSet;do
    while read simMedia;do
        echo smetana --verbose --detailed -o ${simMedia}_${modelSet} --mediadb media_db.tsv -m $simMedia models/$modelSet/*.xml;
    done< <(less media_db.tsv |cut -f1|tail -n +2|sort|uniq)
done< <(ls models/)
