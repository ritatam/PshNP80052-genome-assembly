## Genome assembly of *Puccinia striiformis* f. sp. *hordei* (*Psh*) isolate NP80052

Rita Tam & Julian Rodriguez-Algaba

**/!\ To keep up with current haplotype naming conventions, we renamed haplotype A to haplotype 1, and haplotype B to haplotype 2, in the manuscript.**

---

### 1. simplex ONT basecalling & read assessment

```bash
# basecalling
model=dorado-0.7.2-linux-x64/models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0
dorado basecaller $model --device cuda:all --emit-fastq pod5/ > PshNP80052.ont.fastq
gzip PshNP80052.ont.fastq

# basic read stats
seqkit stats --all -T PshNP80052.ont.fastq.gz
NanoPlot --fastq PshNP80052.ont.fastq.gz --verbose --N50 -o nanoplot
```

---

### 2. genome size estimation

```bash
jellyfish count -C -m 21 -s 1000000000 PshNP80052.ont.fastq -o PshNP80052.jf
jellyfish histo PshNP80052.jf > PshNP80052.histo
```

then upload k-mer count histogram `PshNP80052.histo` to GenomeScope v2.0 (http://genomescope.org/genomescope2.0/) for visualisation.

---

### 3. genome assembly & clean up

```bash
# genome assembly
hifiasm -t 20 -o PshNP80052 --dual-scaf --telo-m CCCTAA --ont PshNP80052.ont.fastq.gz --h1 PshNP80052.HiC_R1.fastq.gz --h2 PshNP80052.HiC_R2.fastq.gz
awk '/^S/{print ">"$2"\n"$3}' PshNP80052.hic.hap1.p_ctg.gfa > PshNP80052.hic.hap1.p_ctg.fasta
awk '/^S/{print ">"$2"\n"$3}' PshNP80052.hic.hap2.p_ctg.gfa > PshNP80052.hic.hap2.p_ctg.fasta
cat PshNP80052.hic.hap1.p_ctg.fasta PshNP80052.hic.hap2.p_ctg.fasta > PshNP80052.hic.hap12.p_ctg.fasta

# align reads to assembly
minimap2 -t 20 -ax map-ont --secondary=no PshNP80052.hic.hap12.p_ctg.fasta PshNP80052.ont.fastq.gz | samtools sort -@20 -O BAM -o PshNP80052.hifiasm_raw.bam
samtools index PshNP80052.hifiasm_raw.bam

# compute mean per-base coverage across 10kbp windows
samtools faidx PshNP80052.hic.hap12.p_ctg.fasta
cut -f1,2 PshNP80052.hic.hap12.p_ctg.fasta.fai > genome_file.txt
bedtools makewindows -g genome_file.txt -w 10000 > genome_file.window_10kbp.bed
bamtocov --regions genome_file.window_10kbp.bed --report windows-stats.tsv PshNP80052.hifiasm_raw.bam > coverage.bed
```
**assembly blast check**

Contigs/scaffolds in raw hifiasm assembly were split into 1Mbp chunks for parallel blast with snakemake on HPC. Custom scripts can be found in `assembly_blast/`. Coverage depth and all blast results were processed in `assembly_blast/contig_cleanup.ipynb` to filter contigs.

```bash
# split assembly into 1Mbp chunks
python split_genome.py PshNP80052.hic.hap12.p_ctg.fasta PshNP80052_chunks
# edit config.yaml then blast chunks in parallel with snakemake
snakemake -j 1000 --max-jobs-per-second 2 --use-envmodules --keep-going
# concatenate blast results
cat PshNP80052_blast_output/*.blast > all.blast
```

---

### 4. Hi-C phasing & scaffolding

see custom script `juicer_3ddna.sh`. This processes all required input files (e.g. site positions) for Juicer and 3D-DNA. 

**Juicer configuration**

Before using the script, configure Juicer's script `juicer.sh`, such as ligation junction sites. We used a restriction enzyme cocktail (DpnII,HinFI,MseI,DdeI) so need to generate all possible sites on our own. I wrote some python codes that use ligation site sequences (same as for HiC-Pro) to generate all possible combinations for Ns. 

```python
import itertools
re = ['GATCGATC', 'GANTANTC', 'TTATAA', 'CTNATNAG']
result = []
for sequence in re:
    if 'N' in sequence:
        n_indices = [i for i, char in enumerate(sequence) if char == 'N']
        n_count = len(n_indices)
        combinations = list(itertools.product('ATCG', repeat=n_count))
        for comb in combinations:
            new_sequence = sequence
            for i, char in zip(n_indices, comb):
                new_sequence = new_sequence[:i] + char + new_sequence[i+1:]
            result.append(new_sequence)
    else:
        result.append(sequence)
print(result)

"|".join(result)
'GATCGATC|GAATAATC|GAATATTC|GAATACTC|GAATAGTC|GATTAATC|GATTATTC|GATTACTC|GATTAGTC|GACTAATC|GACTATTC|GACTACTC|GACTAGTC|GAGTAATC|GAGTATTC|GAGTACTC|GAGTAGTC|TTATAA|CTAATAAG|CTAATTAG|CTAATCAG|CTAATGAG|CTTATAAG|CTTATTAG|CTTATCAG|CTTATGAG|CTCATAAG|CTCATTAG|CTCATCAG|CTCATGAG|CTGATAAG|CTGATTAG|CTGATCAG|CTGATGAG' 
```

then paste it in `juicer.sh` under ligation junction setting like so:
```bash
...
## Set ligation junction based on restriction enzyme
if [ -z "$ligation" ]
then
    case $site in
        DpnII-HinFI-MseI-DdeI) ligation="'(GATCGATC|GAATAATC|GAATATTC|GAATACTC|GAATAGTC|GATTAATC|GATTATTC|GATTACTC|GATTAGTC|GACTAATC|GACTATTC|GACTACTC|GACTAGTC|GAGTAATC|GAGTATTC|GAGTACTC|GAGTAGTC|TTATAA|CTAATAAG|CTAATTAG|CTAATCAG|CTAATGAG|CTTATAAG|CTTATTAG|CTTATCAG|CTTATGAG|CTCATAAG|CTCATTAG|CTCATCAG|CTCATGAG|CTGATAAG|CTGATTAG|CTGATCAG|CTGATGAG)'" ;;
    ...
...
```
also add recognition sites in Juicer's `generate_site_positions.py`.
```python
patterns = {
    ...
    'DpnII-HinFI-MseI-DdeI': ['GATC', 'GANTC', 'TTAA', 'CTNAG']
    }
```
**Run `juicer_3ddna.sh`**
```bash
# bash juicer_3ddna.sh sample_name restriction_enzyme ref_fasta hic1 hic2 working_dir
bash juicer_3ddna.sh \
    PshNP80052 \
    DpnII-HinFI-MseI-DdeI \
    PshNP80052.hic.hap12.p_ctg.fasta \
    PshNP80052.HiC_R1.fastq.gz  \
    PshNP80052.HiC_R2.fastq.gz \
    .
```

Output files in `3d-dna/PshNP80052.0.assembly` and `3d-dna/PshNP80052.0.hic` were imported into Juicebox for heatmap visualisation, manual scaffolding (only needed for a few chromosomes), haplotype assignment and re-orientation. 

Once happy with the results, run 3d-dna post review to introduce gaps at contig junctions to form new scaffolds. 36 chromosomal scaffolds obtained. 

```bash
run-asm-pipeline-post-review.sh -g 5000 --review PshNP80052.0.reviewed.assembly references/PshNP80052.hic.hap12.p_ctg.fasta aligned/merged_nodups.txt
```

---

### 5. assembly evaluation

**basic assembly statistics**

```bash
# basic assembly stats
seqkit stats --all -T PshNP80052.final.fasta

# count telomeres
FindTelomeres.py PshNP80052.final.fasta

# busco completeness
busco -i PshNP80052.final.fasta -l basidiomycota -m DNA -t 24

# k-mer commpleteness
!!!! to-do !!!!
```

**count within- and cross-haplotype Hi-C contact links**

Generate HiC contact matrices with HiC-Pro. Iutput files prep, configuration (ligation sites set as above; minimum MAPQ threshold set to 20) and HiC-Pro execution are automated in `hicpro.sh`
```bash
mkdir -p hicpro/PshNP80052
bash hicpro.sh -g PshNP80052.final.fasta -s PshNP80052 --hic1 PshNP80052.HiC_R1.fastq.gz --hic2 PshNP80052.HiC_R2.fastq.gz -o hicpro/PshNP80052
```

Count Hi-C contacts from raw matrix for within- and cross-haplotype link statistics and generate a circos plot using codes in https://github.com/ritatam/HiC-Analysis (special thanks to Runpeng Luo who authored the scripts). In the forked repo I only made the haplotype link statistics report more explicit.

```bash
# get hicpro raw matrix (window size 20000) and bed file
cp hicpro/PshNP80052/mapq20_output/hic_results/matrix/PshNP80052/raw/20000/* .

# prepare haplotype chromosome id files
>na.lst
for h in {A,B}; do
    >hap${h}.lst
    for n in {1..18}; do 
        echo ">chr${n}${h}" >> hap${h}.lst
    done
done

# generate report, circos plot and contact histogram across windows
python hic_analysis.py hapA.lst hapB.lst na.lst PshNP80052_20000_abs.bed PshNP80052_20000.matrix . 20000
```
---
### 6. gene annotation

**liftover from reference annotations**

Full dikaryotic gene annotations (hapA + hapB) from reference (Pst104E/AZ2) were lifted onto each PshNP80052 haplotype, one at a time. 
```bash
for n in {1..18}; do echo -e "chr${n}A,chr${n}${hap}\nchr${n}B,chr${n}${hap}" >> chroms_hap${hap}.txt; done
# run liftoff
liftoff -g Pst104E_v3.9.gene_anno.gff3 -chroms chroms_hap${hap}.txt -o output/PshNP80052_hap${hap}.Pst104Ev3.9.liftoff.gff3 -dir output/PshNP80052_hap${hap}.Pst104Ev3.9.liftoff_intermediate_files -p 16 -u output/PshNP80052_hap${hap}.Pst104Ev3.9.unmapped_features.txt PshNP80052_hap${hap}.fasta Pst104Ev3.9.fasta
# keep only models with valid ORFs
agat_sp_filter_feature_by_attribute_value.pl --gff PshNP80052_hap${hap}.Pst104Ev3.9.liftoff.gff3 --attribute valid_ORFs --value 0 --test "=" --type gene -o PshNP80052_hap${hap}.valid_ORFs.gff3 
# re-number the gene models
funannotate util gff-rename -g PshNP80052_hap${hap}.valid_ORFs.gff3 -f PshNP80052_hap${hap}.fasta -o PshNP80052_hap${hap}.valid_ORFs.renum.gff3 --locus_tag PshNP80052 --numbering 1
sed 's/;Alias=.*//g' PshNP80052_hap${hap}.valid_ORFs.renum.gff3 > PshNP80052_hap${hap}.104e-liftoff.gff3
```
Repeated for AZ2. (Note the chromosome assignments were different in AZ2, so `chroms_hap${hap}.txt` was adjusted accordingly for liftoff)

**ab initio predictions**

(Note: I attempted using augustus training parameters generated for evidence-based gene annotations of Pst104E (source [here](https://genome.cshlp.org/content/35/6/1364.full.pdf+html)). However the results missed a lot BUSCOs for some reason; default parameters yield the best results for our Psh so far, so I chose to go with it)

```bash
funannotate predict -i PshNP80052_hap${hap}.masked.fasta -o hap${hap}/funannotate_abinitio  --species "Puccinia striiformis" --isolate PshNP80052_hap${hap} --cpus 24 --optimize_augustus

cat hapA/funannotate_abinitio/predict_results/Puccinia_striiformis_PshNP80052_hapA.gff3 hapB/funannotate_abinitio/predict_results/Puccinia_striiformis_PshNP80052_hapB.gff3 > PshNP80052.abinitio.gff3
```

**Merging of liftoff-derived and ab initio annotations**

See `merge_liftoff_abinitio_annotation.sh`.

**Functional annotation** 

See `functional_annotation.sh`.

---

### 7. karyoplot

!!!to-do!!!

---

### 8. mtDNA assembly

```bash
# concatenate PshNP80052 assembly and AZ2's mtDNA contig for baiting
cat PshNP80052.final.fasta AZ2mtDNA.fasta > PshNP80052_AZ2mtDNA.fasta

# map reads to the concatenated assembly
minimap -ax map-ont --secondary=no PshNP80052_AZ2mtDNA.fasta PshNP80052.ont.fastq.gz | samtools sort -@20 -O BAM -o PshNP80052.PshNP80052_AZ2mtDNA.bam
samtools index PshNP80052.PshNP80052_AZ2mtDNA.bam

# extract reads mapped to AZ2 mtDNA. now these are considered mitochondrial reads
samtools view PshNP80052.PshNP80052_AZ2mtDNA.bam AZ2mtDNA -o PshNP80052.mtDNA.bam
samtools bam2fq PshNP80052.mtDNA.bam > PshNP80052.mtDNA.fastq

# filter mtDNA reads by length
# and subsample down to ~700x to be memory efficient during assembly
seqkit seq --threads 8 --max-len 120000 --min-len 20000 PshNP80052.mtDNA.fastq > PshNP80052.mtDNA.20-120k.fastq
seqkit sample -p 0.3 PshNP80052.mtDNA.20-120k.fastq > PshNP80052.mtDNA.20-120k.subsamp0.3.fastq

# hifiasm
hifiasm --ont --primary -t 24 -f 0 -o hifiasm PshNP80052.mtDNA.20-120k.subsamp0.3.fastq

# extract circular contig
awk '$1 ~ /^>.*c$/ {print; getline; print}' hifiasm.p_ctg.fasta > hifiasm.p_ctg.circular.fasta

# map the mtDNA reads back to the assembled circular contig to inspect alignment
cp hifiasm.p_ctg.circular.fasta PshNP80052.mtDNA_candidate.fasta
minimap2 -ax map-ont PshNP80052.mtDNA_candidate.fasta PshNP80052.mtDNA.20-120k.subsamp0.3.fastq | samtools sort -O BAM -o PshNP80052.mtDNA_candidate.bam
samtools index PshNP80052.mtDNA_candidate.bam
```
