ONT-TU-analyse
==============

This is the Python script used to analyse microsynth ont reads from Nashlab MoClo plasmids.

[![build passing](https://img.shields.io/badge/build-passing-brightgreen.svg?style=flat)](https://github.com/lorenzwidmer/Ont-TU-analyse)  [![version 0.1](https://img.shields.io/badge/version-0.1-blue.svg?style=flat)](https://github.com/lorenzwidmer/Ont-TU-analyse) [![contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/lorenzwidmer/Ont-TU-analyse) [![awesomeness 110%](https://img.shields.io/badge/awesomeness-110%25-red.svg?style=flat)](https://www.buymeacoffee.com/ntIOpc7) 

---
<a href="https://www.buymeacoffee.com/ntIOpc7" target="_blank"><img src="https://www.buymeacoffee.com/assets/img/custom_images/orange_img.png" alt="Buy Me A Coffee" style="height: 41px !important;width: 174px !important;" ></a>

Installation
------------

```shell
# clone repo
git clone git@github.com:lorenzwidmer/Ont-TU-analyse.git 250725_example_run
cd 250725_example_run

# create python virtual environment
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt 
```

Usage
-----

```shell
# convert parts exported from Benchling to fasta format. Teh script strips adapters (only the sequence between two BsaI sites is kept) and classifies parts based on their overhangs
python ./script/parts2fasta.py parts/sd_parts.csv parts/sd_parts.fa

# pre-process ont reads and cutting out TUs (only the sequence between the two primers is kept)
script/cutadapt_bc.sh raw/pTU1b-la-xxx.raw.reads.fastq.gz gaattcgcggccgcttctagag ttctgtggataaccgtattaccgc plasmid/pTU1b-la-xxx-lib.reads.fastq
script/cutadapt_bc.sh raw/pTU2a-ls-xxx.raw.reads.fastq.gz gaattcgcggccgcttctagag ggtgacaccttgcccttttttgccgga plasmid/pTU2a-ls-xxx-lib.reads.fastq

# annotate reads
python ./script/annotate.py -o output/pTU1b-la-xxx-lib/ -a parts/sd_parts.fa -f '{"TU1":[0,1,2,3,4,5]}' plasmid/pTU1b-la-xxx-lib.reads.fastq
python ./script/annotate.py -o output/pTU2a-ls-xxx-lib/ -a parts/sd_parts.fa -f '{"TU1":[0,1,2,4,5],"TU2":[0,1,2,4,5]}' plasmid/pTU2a-ls-xxx-lib.reads.fastq
```
