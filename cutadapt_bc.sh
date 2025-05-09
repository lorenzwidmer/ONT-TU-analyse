TMPD=$(mktemp -d -t "XXXXXX")
CPUS=8

FWD=$2
RVS=$3

FWD_RC=`echo $FWD | tr '[ATUGCYRSWKMBDHNatugcyrswkmbdhn]' '[TAACGRYSWMKVHDNtaacgryswmkvhdn]' | rev`
RVS_RC=`echo $RVS | tr '[ATUGCYRSWKMBDHNatugcyrswkmbdhn]' '[TAACGRYSWMKVHDNtaacgryswmkvhdn]' | rev`

FWD_LEN=${#FWD}
RVS_LEN=${#RVS}

seqkit concat $1 $1 -o $TMPD/duplicated.fastq.gz
cutadapt -g "$FWD;min_overlap=$FWD_LEN...$RVS_RC;min_overlap=$RVS_LEN" -e 0.15 --cores=$CPUS --discard-untrimmed --action=trim --report=minimal -o $TMPD/cutadapt_fwd.fastq.gz $TMPD/duplicated.fastq.gz
cutadapt -g "$RVS;min_overlap=$RVS_LEN...$FWD_RC;min_overlap=$FWD_LEN" -e 0.15 --cores=$CPUS --discard-untrimmed --action=trim --report=minimal -o $TMPD/cutadapt_rvs.fastq.gz $TMPD/duplicated.fastq.gz

seqkit -j $CPUS -t dna -o $TMPD/cutadapt_rvs.rc.fastq.gz seq -r -p $TMPD/cutadapt_rvs.fastq.gz
cat $TMPD/cutadapt_fwd.fastq.gz $TMPD/cutadapt_rvs.rc.fastq.gz > $4

rm -rf $TMPD
