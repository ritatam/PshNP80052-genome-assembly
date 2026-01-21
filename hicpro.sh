#!/bin/bash

# bash hicpro.sh -g <ref.fasta> -s <sample_name> --hic1 <read1.fastq.gz> --hic2 <read2.fastq.gz> -o <output_dir>

usage() {
  echo "Usage: $0 -g <ref.fasta> -s <sample_name> --hic1 <read1.fastq.gz> --hic2 <read2.fastq.gz> -o <output_dir>"
  exit 1
}

PARSED_ARGS=$(getopt -o g:s:o: --long hic1:,hic2: -- "$@")
if [[ $? -ne 0 ]]; then usage; fi
eval set -- "$PARSED_ARGS"

while true; do
  case "$1" in
    -g) ref="$2"; shift 2 ;;
    -s) sample_name="$2"; shift 2 ;;
    --hic1) hic1="$2"; shift 2 ;;
    --hic2) hic2="$2"; shift 2 ;;
    -o) outdir="$2"; shift 2 ;;
    --) shift; break ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

if [[ -z "$ref" || -z "$sample_name" || -z "$hic1" || -z "$hic2" || -z "$outdir" ]]; then
  echo "Missing one or more required arguments."
  usage
fi

source /g/data/xf3/miniconda/etc/profile.d/conda.sh
conda activate /g/data/xe2/ht5438/conda_env/hicpro
export PATH=$PATH:/g/data/xe2/ht5438/softwares/HiC-Pro-3.1.0/bin

mkdir -p $outdir
mkdir -p ${outdir}/HiC_input/$sample_name

# check if the input hic files are gzipped. if zipped, unzip them
if [[ $hic1 == *.gz ]]; then
    gunzip $hic1
    hic1=${hic1%.gz}
fi
if [[ $hic2 == *.gz ]]; then
    gunzip $hic2
    hic2=${hic2%.gz}
fi

# softlink the input hic files to the HiC input directory
ln -sr $hic1 ${outdir}/HiC_input/$sample_name/${sample_name}_R1.fastq
ln -sr $hic2 ${outdir}/HiC_input/$sample_name/${sample_name}_R2.fastq

# generate restriction sites annotation file
HICPRO_UTILS=/g/data/xe2/ht5438/softwares/HiC-Pro-3.1.0/bin/utils
$HICPRO_UTILS/digest_genome.py -r ^GATC,G^ANTC,T^TAA,C^TNAG -o ${outdir}/${sample_name}.digested.bed $ref

# generate bowtie2 index
mkdir -p ${outdir}/bowtie2_index
bowtie2-build $ref ${outdir}/bowtie2_index/$(basename "${ref%.fasta}")

# generate chromosome size file
samtools faidx $ref
cut -f1,2 ${ref}.fai > ${outdir}/chrom.sizes

# set up hicpro configuration file
cp /g/data/xe2/ht5438/softwares/HiC-Pro-3.1.0/config-hicpro.txt ${outdir}/config-hicpro.txt
sed -i \
  -e "s|^BOWTIE2_IDX_PATH.*|BOWTIE2_IDX_PATH = ${outdir}/bowtie2_index|" \
  -e "s|^REFERENCE_GENOME.*|REFERENCE_GENOME = $(basename "${ref%.fasta}")|" \
  -e "s|^GENOME_SIZE.*|GENOME_SIZE = ${outdir}/chrom.sizes|" \
  -e "s|^GENOME_FRAGMENT.*|GENOME_FRAGMENT = ${outdir}/${sample_name}.digested.bed|" \
  -e "s|^LIGATION_SITE.*|LIGATION_SITE = GATCGATC,GANTANTC,TTATAA,CTNATNAG|" \
  -e "s|^MIN_MAPQ.*|MIN_MAPQ = 20|" \
  "${outdir}/config-hicpro.txt"

HiC-Pro -i ${outdir}/HiC_input -o ${outdir}/mapq20_output -c ${outdir}/config-hicpro.txt