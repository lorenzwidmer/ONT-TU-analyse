TMPD=$(mktemp -d -t "XXXXXX")
CPUS=8

FWD=$1
RVS=$2
OUT=$3

FWD_RC=`echo $FWD | tr '[ATUGCYRSWKMBDHNatugcyrswkmbdhn]' '[TAACGRYSWMKVHDNtaacgryswmkvhdn]' | rev`
RVS_RC=`echo $RVS | tr '[ATUGCYRSWKMBDHNatugcyrswkmbdhn]' '[TAACGRYSWMKVHDNtaacgryswmkvhdn]' | rev`

FWD_LEN=${#FWD}
RVS_LEN=${#RVS}

shift 3
cat $@ > $TMPD/in.fastq.gz

cutadapt -g "$FWD;min_overlap=$FWD_LEN...$RVS_RC;min_overlap=$RVS_LEN" -e 0.15 --revcomp --cores=$CPUS --untrimmed-output $TMPD/ca_full_untrm.fastq.gz --rename='{id}' --action=retain --report=minimal -o $TMPD/ca_full.fastq.gz $TMPD/in.fastq.gz
cutadapt -g "$FWD;min_overlap=$FWD_LEN" -e 0.15 --revcomp --cores=$CPUS --discard-untrimmed --rename='{id}' --action=retain --report=minimal -o $TMPD/ca_fwd.fastq.gz $TMPD/ca_full_untrm.fastq.gz
cutadapt -a "$RVS_RC;min_overlap=$RVS_LEN" -e 0.15 --revcomp --cores=$CPUS --discard-untrimmed --rename='{id}' --action=retain --report=minimal -o $TMPD/ca_rev.fastq.gz $TMPD/ca_full_untrm.fastq.gz

seqkit concat $TMPD/ca_fwd.fastq.gz $TMPD/ca_rev.fastq.gz -o $TMPD/ca_ligated.fastq.gz
zcat $TMPD/ca_fwd.fastq.gz $TMPD/ca_rev.fastq.gz | seqkit rmdup -n -o $TMPD/ca_single.fastq.gz

cat $TMPD/ca_full.fastq.gz $TMPD/ca_ligated.fastq.gz > $OUT

rm -rf $TMPD