#################################
##### functional annotation #####
#################################

export FUNANNOTATE_DB=gene_annotation/funannotate_db
hap=hapA 
ref=assembly/curated_assembly/PshNP80052_${hap}.final.fasta
numbering=1     ## 17239 for hapB
gff=PshNP80052.liftoff_abinitio_combined.dedup.longest_orf.prot_50aa.renum.utr_clipped.gff3

cd liftoff_abinitio_combined/functional_annotation/01_fix_models
cp ../../preprocessing/$gff .

# fix problematic gene models
funannotate util gff2tbl -g $gff -f $ref > PshNP80052.liftoff_abinitio_merged.${hap}.tbl
funannotate util tbl2gbk --tbl PshNP80052.liftoff_abinitio_merged.${hap}.tbl -f $ref --species "Puccinia striiformis" --isolate NP80052 -o PshNP80052.liftoff_abinitio_merged.${hap}
## --------------------- 
## hapA genes to fix:
## --------------------- 
# PshNP80052_000407       Location: Adjacent intervals in SeqLoc [(lcl|chr1A:1591398-1591885, 1591979-1592028, 1593340-1594001, 1594002-1594310)]
# PshNP80052_000408       Location: Adjacent intervals in SeqLoc [(lcl|chr1A:1592817-1592839, 1593040-1593054, 1593129-1593268, 1593340-1594001, 1594002-1594278, 1594339-1594637, 1594707-1595013, 1595091-1595179)]
# PshNP80052_003766       Coding region extends 1 base(s) past stop codon
# PshNP80052_006171       Coding region extends 2 base(s) past stop codon
# PshNP80052_007341       Coding region extends 1 base(s) past stop codon
# PshNP80052_009661       Coding region extends 1 base(s) past stop codon
# PshNP80052_010464       Location: Adjacent intervals in SeqLoc [(lcl|chr10A:c548599-548522, c548392-548366, c548271-548027, c548026-547278, c547183-546450)]
# PshNP80052_011441       Location: Adjacent intervals in SeqLoc [(lcl|chr11A:c248448-248262, c248140-248036, c247958-247800, c247720-247516, c247515-246261)]
# PshNP80052_012884       Location: Adjacent intervals in SeqLoc [(lcl|chr12A:c2413201-2413184, c2411749-2411731, c2411643-2411425, c2411348-2411028, c2410949-2410787, c2410786-2410659, c2410580-2409922, c2409856-2409087, c2409005-2407214, c2407145-2407080)]
# PshNP80052_014770       Coding region extends 2 base(s) past stop codon
# PshNP80052_015571       Coding region extends 1 base(s) past stop codon
# PshNP80052_015700       Location: Adjacent intervals in SeqLoc [(lcl|chr16A:c1099840-1099574, c1099484-1099057, c1099056-1098015)]
# PshNP80052_015878       Coding region extends 1 base(s) past stop codon

## --------------------- 
## hapB genes to fix:
## --------------------- 
# PshNP80052_018895       Coding region extends 1 base(s) past stop codon
# PshNP80052_021099       Coding region extends 1 base(s) past stop codon
# PshNP80052_024830       Coding region extends 1 base(s) past stop codon
# PshNP80052_027985       Location: Adjacent intervals in SeqLoc [(lcl|chr10B:c564014-563937, c563807-563781, c563686-563442, c563441-562693, c562598-561865)]
# PshNP80052_031249       Location: Adjacent intervals in SeqLoc [(lcl|chr13B:1774445-1774844, 1774908-1775549, 1775624-1775648, 1775746-1775754, 1775823-1775909, 1775910-1775946, 1776037-1776101, 1776177-1776254, 1776339-1776380, 1776444-1776593, 1776660-1776674, 1776793-1776855, 1776923-1777010, 1777097-1777109, 1777183-1777455, 1777552-1777572, 1777670-1777779, 1777872-1777887, 1777959-1778084, 1778165-1778178, 1778274-1778295, 1778387-1778448)]
# PshNP80052_032356       Coding region extends 2 base(s) past stop codon
# PshNP80052_033018       Coding region extends 2 base(s) past stop codon
# PshNP80052_033522       Coding region extends 1 base(s) past stop codon
# PshNP80052_033818       Coding region extends 1 base(s) past stop codon
# PshNP80052_033965       Location: Adjacent intervals in SeqLoc [(lcl|chr17B:c919308-919291, c917850-917832, c917744-917526, c917449-917129, c917050-916888, c916887-916760, c916681-916023, c915957-915188, c915106-913315, c913246-913181)]
# PshNP80052_034499       Coding region extends 1 base(s) past stop codon

# clear functional annotation from previous annotations (product field)
sed -i 's/;product=.*//g' PshNP80052.liftoff_abinitio_merged.${hap}.fixed.gff3
funannotate util gff-rename -g PshNP80052.liftoff_abinitio_merged.${hap}.fixed.gff3 -f $ref -o PshNP80052.liftoff_abinitio_merged.${hap}.fixed.renum.gff3 --locus_tag PshNP80052 --numbering $numbering # this will add back the product field, but all are hypothetical proteins
sed -i 's/;Alias=.*//g' PshNP80052.liftoff_abinitio_merged.${hap}.fixed.renum.gff3
sed -i 's/;$//' PshNP80052.liftoff_abinitio_merged.${hap}.fixed.renum.gff3

# antismash
gff=liftoff_abinitio_combined/functional_annotation/01_fix_models/PshNP80052.liftoff_abinitio_merged.${hap}.fixed.renum.gff3
mkdir -p liftoff_abinitio_combined/functional_annotation/02_antismash
cd liftoff_abinitio_combined/functional_annotation/02_antismash
funannotate util gff2tbl -g $gff -f $ref > PshNP80052.liftoff_abinitio_merged.${hap}.tbl
funannotate util tbl2gbk --tbl PshNP80052.liftoff_abinitio_merged.${hap}.tbl -f $ref --species "Puccinia striiformis" --isolate NP80052 -o PshNP80052.liftoff_abinitio_merged.${hap}
conda activate /home/groups/schwessinger/condaEnvs/antismash
download-antismash-databases  # download antismash databases if not already done
antismash --taxon fungi --output-dir ${hap} --output-basename antismash_${hap} --smcog-trees --cb-knownclusters --asf --cb-subclusters PshNP80052.liftoff_abinitio_merged.${hap}.gbk
antismash_gbk=liftoff_abinitio_combined/functional_annotation/02_antismash/${hap}/antismash_${hap}.gbk

# signalp3 
conda activate /home/groups/schwessinger/condaEnvs/funanntoate
mkdir -p liftoff_abinitio_combined/functional_annotation/03_signalp3
cd liftoff_abinitio_combined/functional_annotation/03_signalp3
ln -sr ../01_fix_models/PshNP80052.liftoff_abinitio_merged.${hap}.fixed.renum.gff3
funannotate gff2prot -g PshNP80052.liftoff_abinitio_merged.${hap}.fixed.renum.gff3 -f $ref --no_stop > PshNP80052.${hap}.proteins.fasta
conda activate /home/groups/schwessinger/condaEnvs/common-tools
export PATH=$PATH:/media/ssd/rita/softwares/signalp-3.0
seqkit split --by-id PshNP80052.${hap}.proteins.fasta --by-id-prefix '' -O ${hap}_flatten
echo running signalp3-nn euk; 
echo 'name	Cmax	pos	?	Ymax	pos	?	Smax	pos	?	Smean	?	D	?' > PshNP80052.${hap}.signalp3.out
find ${hap}_flatten -name "*.fasta" | parallel -j 28 "signalp -t euk -short -m nn {} | sed '1,2d' >> PshNP80052.${hap}.signalp3.out"
conda activate /home/groups/schwessinger/condaEnvs/jupyter
python signalp3_mature_prot_extract.py PshNP80052.${hap}.signalp3.out PshNP80052.${hap}.proteins.fasta PshNP80052.${hap}.signalp3.mature_prot.fasta

# tmhmm
mkdir -p liftoff_abinitio_combined/functional_annotation/04_tmhmm
cd liftoff_abinitio_combined/functional_annotation/04_tmhmm
ln -sr ../03_signalp3/PshNP80052.${hap}.signalp3.mature_prot.fasta
export PATH=$PATH:/media/ssd/rita/softwares/tmhmm-2.0c/bin
tmhmm --workdir . --short PshNP80052.${hap}.signalp3.mature_prot.fasta > PshNP80052.${hap}.tmhmm.out

# phobius on mature protein
mkdir -p liftoff_abinitio_combined/functional_annotation/05_phobius
cd liftoff_abinitio_combined/functional_annotation/05_phobius
ln -sr ../03_signalp3/PshNP80052.${hap}.signalp3.mature_prot.fasta
export PATH=$PATH:/media/ssd/rita/softwares/phobius
phobius.pl -short PshNP80052.${hap}.signalp3.mature_prot.fasta | tail -n+2 > PshNP80052.${hap}.phobius.out 

# filter secretome proteins with transmembrane domains detected by both tmhmm and phobius
mkdir -p liftoff_abinitio_combined/functional_annotation/06_transmembrane_filter
cd liftoff_abinitio_combined/functional_annotation/06_transmembrane_filter
conda activate /home/groups/schwessinger/condaEnvs/jupyter
grep ">" ../03_signalp3/PshNP80052.${hap}.signalp3.mature_prot.fasta | cut -f1 | sed 's/>//' > PshNP80052.${hap}.signalp3.mature_prot.ids
python transmembrane_filter.py --phobius ../05_phobius/PshNP80052.${hap}.phobius.out --tmhmm ../04_tmhmm/PshNP80052.${hap}.tmhmm.out --signalp_union_id PshNP80052.${hap}.signalp3.mature_prot.ids --output PshNP80052.${hap}.secretome_noTM.ids
echo 'name	Cmax	pos	?	Ymax	pos	?	Smax	pos	?	Smean	?	D	?' > PshNP80052.${hap}.signalp3-noTM.out
while read -r id; do
    grep -w $id ../03_signalp3/PshNP80052.${hap}.signalp3.out >> PshNP80052.${hap}.signalp3-noTM.out
done < PshNP80052.${hap}.secretome_noTM.ids

tail -n +2 PshNP80052.${hap}.signalp3-noTM.out > PshNP80052.${hap}.signalp3-noTM.tmp
>PshNP80052.${hap}.signalp3.custom.txt
while read -r id; do
    cleavepos=$(awk -v id="$id" '$1 == id {print $6 - 1; exit}' PshNP80052.${hap}.signalp3-noTM.tmp)
    echo -e "$id\tnote\tSECRETED:SignalP(1-$cleavepos)" >> PshNP80052.${hap}.signalp3.custom.txt
done < PshNP80052.${hap}.secretome_noTM.ids


# iprscan
mkdir -p liftoff_abinitio_combined/functional_annotation/07_iprscan
cd liftoff_abinitio_combined/functional_annotation/07_iprscan
conda activate /home/groups/schwessinger/condaEnvs/funanntoate
ln -sr ../03_signalp3/PshNP80052.${hap}.proteins.fasta
export FUNANNOTATE_DB=gene_annotation/funannotate_db
funannotate iprscan -i PshNP80052.${hap}.proteins.fasta \
	-m local \
	-o PshNP80052.${hap}.iprscan.5.64-96.0.local.xml \
	-c 4 \
	--iprscan_path /media/ssd/rita/softwares/interproscan-5.64-96.0/interproscan.sh


# protein2genome_exonerate
mkdir -p liftoff_abinitio_combined/functional_annotation/08_prot2genome
conda activate /home/groups/schwessinger/condaEnvs/funanntoate
cd liftoff_abinitio_combined/functional_annotation/08_prot2genome
ln -sr ../03_signalp3/PshNP80052.${hap}.proteins.fasta
funannotate util prot2genome \
	-g $ref \
	-p PshNP80052.${hap}.proteins.fasta \
	-o PshNP80052.${hap}.p2g.gff3 \
	--cpus 24 \
	--logfile ${hap}-p2g.log


# phobius on full-length protein
mkdir -p liftoff_abinitio_combined/functional_annotation/09_phobius_fulllen
cd liftoff_abinitio_combined/functional_annotation/09_phobius_fulllen
ln -sr ../03_signalp3/PshNP80052.${hap}.proteins.fasta
export PATH=$PATH:/media/ssd/rita/softwares/phobius
phobius.pl -short PshNP80052.${hap}.proteins.fasta | tail -n+2 > PshNP80052.${hap}.phobius.full-length-prot.out

awk '$2 > 0' PshNP80052.${hap}.phobius.full-length-prot.out > PshNP80052.${hap}.phobius.full-length-prot.TM.out
grep -v -F -f ../06_transmembrane_filter/PshNP80052.${hap}.secretome_noTM.ids PshNP80052.${hap}.phobius.full-length-prot.TM.out > PshNP80052.${hap}.phobius.full-length-prot.TM_noSP.out
>PshNP80052.${hap}.phobius.full-length-prot.TM_noSP.custom.txt
while read -r line; do
    id=$(echo $line | cut -d' ' -f1)
    tm=$(echo $line | cut -d' ' -f2)
    topology=$(echo $line | cut -d' ' -f4)
    echo -e "$id\tnote\tTransmembrane_phobius:$tm($topology)" >> PshNP80052.${hap}.phobius.full-length-prot.TM_noSP.custom.txt
done < PshNP80052.${hap}.phobius.full-length-prot.TM_noSP.out 


# functional annotation
mkdir -p liftoff_abinitio_combined/functional_annotation/10_annotate
mkdir -p liftoff_abinitio_combined/functional_annotation/10_annotate/${hap}
cd liftoff_abinitio_combined/functional_annotation/10_annotate/${hap}
export GENEMARK_PATH=gene_annotation/funannotate_db/gmes_linux_64_4
export PATH=$PATH:gene_annotation/funannotate_db/gmes_linux_64_4
export EGGNOG_DATA_DIR=/media/ssd/rita/softwares/eggnog-mapper-2.1.13/data

funannotate annotate \
    --gff ../../01_fix_models/PshNP80052.liftoff_abinitio_merged.${hap}.fixed.renum.gff3 \
    --fasta $ref \
    --species "Puccinia striiformis" \
    --out . \
    --antismash ../../02_antismash/${hap}/antismash_${hap}.gbk \
    --signalp liftoff_abinitio_combined/functional_annotation/signalp.header \
    --iprscan ../../07_iprscan/PshNP80052.${hap}.iprscan.5.64-96.0.local.xml \
    --phobius ../../09_phobius_fulllen/PshNP80052.${hap}.phobius.full-length-prot.TM_noSP.out\
    --annotations ../../06_transmembrane_filter/PshNP80052.${hap}.signalp3.custom.txt \
    --p2g ../../08_prot2genome/PshNP80052.${hap}.p2g.gff3 \
    --isolate NP80052 \
    --header_length 20 \
    --busco_db basidiomycota \
    --cpus 24


# add tRNA back, but note that all other intermedate files from final funannotate annotate steps will be unusable after this step because numbers will change
export FUNANNOTATE_DB=gene_annotation/funannotate_db
hap=hapB
ref=assembly/curated_assembly/PshNP80052_${hap}.final.fasta
numbering=17810

mkdir -p liftoff_abinitio_combined/functional_annotation/11_add_tRNA
cd liftoff_abinitio_combined/functional_annotation/11_add_tRNA
conda activate /home/groups/schwessinger/condaEnvs/jupyter
python liftoff_abinitio_combined/clip_UTR.py -i ../10_annotate/${hap}/annotate_results/Puccinia_striiformis_NP80052.gff3  -o Puccinia_striiformis_NP80052.funannotate_annotate.${hap}.gff3
grep -v "T2" ../../preprocessing/PshNP80052.liftoff_abinitio_combined.dedup.${hap}.tRNA.gff3 > PshNP80052.liftoff_abinitio_combined.dedup.${hap}.tRNA.gff3
sed -i 's/PshNP80052/tmp/g' PshNP80052.liftoff_abinitio_combined.dedup.${hap}.tRNA.gff3
cat Puccinia_striiformis_NP80052.funannotate_annotate.${hap}.gff3 PshNP80052.liftoff_abinitio_combined.dedup.${hap}.tRNA.gff3 > Puccinia_striiformis_NP80052.funannotate_annotate.with_tRNA.${hap}.gff3
conda activate /home/groups/schwessinger/condaEnvs/funanntoate
funannotate util gff-rename -g Puccinia_striiformis_NP80052.funannotate_annotate.with_tRNA.${hap}.gff3 -f $ref -o Puccinia_striiformis_NP80052.funannotate_annotate.with_tRNA.${hap}.renum.gff3 --locus_tag PshNP80052 --numbering $numbering
sed 's/;Alias=[^;]*//g' Puccinia_striiformis_NP80052.funannotate_annotate.with_tRNA.${hap}.renum.gff3 > Puccinia_striiformis_NP80052.funannotate_annotate.with_tRNA.${hap}.renum.clean.gff3

# finalise gff3
mkdir -p liftoff_abinitio_combined/functional_annotation/12_finalise
cd liftoff_abinitio_combined/functional_annotation/12_finalise
fold $ref > PshNP80052.liftoff_abinitio.final.${hap}.fold.fasta
conda activate /home/groups/schwessinger/condaEnvs/common-tools
cp ../11_add_tRNA/Puccinia_striiformis_NP80052.funannotate_annotate.with_tRNA.${hap}.renum.clean.gff3 PshNP80052.liftoff_abinitio.final.${hap}.gff3
agat_sp_extract_sequences.pl -g PshNP80052.liftoff_abinitio.final.${hap}.gff3 -f PshNP80052.liftoff_abinitio.final.${hap}.fold.fasta -o PshNP80052.liftoff_abinitio.final.${hap}.cds-transcripts.fa -t cds
agat_sp_extract_sequences.pl -g PshNP80052.liftoff_abinitio.final.${hap}.gff3 -f PshNP80052.liftoff_abinitio.final.${hap}.fold.fasta -o PshNP80052.liftoff_abinitio.final.${hap}.proteins.fa -t cds -p
sed -i 's/type=cds/type=prot/g' PshNP80052.liftoff_abinitio.final.${hap}.proteins.fa
conda activate /home/groups/schwessinger/condaEnvs/funanntoate
funannotate gff2tbl -g PshNP80052.liftoff_abinitio.final.${hap}.gff3 -f PshNP80052.liftoff_abinitio.final.${hap}.fold.fasta > PshNP80052.liftoff_abinitio.final.${hap}.tbl
funannotate util tbl2gbk --tbl PshNP80052.liftoff_abinitio.final.${hap}.tbl -f PshNP80052.liftoff_abinitio.final.${hap}.fold.fasta --species "Puccinia striiformis" --isolate NP80052 -o PshNP80052.liftoff_abinitio.final.${hap}
funannotate stats -g PshNP80052.liftoff_abinitio.final.${hap}.gff3 -f PshNP80052.liftoff_abinitio.final.${hap}.fold.fasta -o PshNP80052.liftoff_abinitio.final.${hap}.stats.txt
# get secretome
grep -P "\tmRNA\t" PshNP80052.liftoff_abinitio.final.${hap}.gff3 | grep "SECRETED:SignalP" | cut -f9 | sed -E 's/.*ID=([^;]+)-T[0-9]+.*/\1/' > PshNP80052.liftoff_abinitio.final.${hap}.secretome.tmp
grep -F -f PshNP80052.liftoff_abinitio.final.${hap}.secretome.tmp PshNP80052.liftoff_abinitio.final.${hap}.gff3 > PshNP80052.liftoff_abinitio.final.${hap}.secretome.gff3
rm PshNP80052.liftoff_abinitio.final.${hap}.secretome.tmp

# end of gene annotation.

# combine annotations into full genome set
conda activate /home/groups/schwessinger/condaEnvs/common-tools
cat hapA/PshNP80052.liftoff_abinitio.final.hapA.gff3 hapB/PshNP80052.liftoff_abinitio.final.hapB.gff3 > combined/PshNP80052.liftoff_abinitio.final.gff3
cat hapA/PshNP80052.liftoff_abinitio.final.hapA.proteins.fa hapB/PshNP80052.liftoff_abinitio.final.hapB.proteins.fa > combined/PshNP80052.liftoff_abinitio.final.proteins.fa
cat hapA/PshNP80052.liftoff_abinitio.final.hapA.cds-transcripts.fa hapB/PshNP80052.liftoff_abinitio.final.hapB.cds-transcripts.fa > combined/PshNP80052.liftoff_abinitio.final.cds-transcripts.fa
cat hapA/PshNP80052.liftoff_abinitio.final.hapA.gbk hapB/PshNP80052.liftoff_abinitio.final.hapB.gbk > combined/PshNP80052.liftoff_abinitio.final.gbk
cat hapA/PshNP80052.liftoff_abinitio.final.hapA.tbl hapB/PshNP80052.liftoff_abinitio.final.hapB.tbl > combined/PshNP80052.liftoff_abinitio.final.tbl
cat hapA/PshNP80052.liftoff_abinitio.final.hapA.secretome.gff3 hapB/PshNP80052.liftoff_abinitio.final.hapB.secretome.gff3 > combined/PshNP80052.liftoff_abinitio.final.secretome.gff3 > combined/PshNP80052.liftoff_abinitio.final.secretome.gff3
agat_sp_statistics.pl --gff combined/PshNP80052.liftoff_abinitio.final.gff3 -o combined/PshNP80052.liftoff_abinitio.final.gff3.stats
agat_sp_statistics.pl --gff hapA/PshNP80052.liftoff_abinitio.final.hapA.gff3 -o hapA/PshNP80052.liftoff_abinitio.final.hapA.gff3.stats
agat_sp_statistics.pl --gff hapB/PshNP80052.liftoff_abinitio.final.hapB.gff3 -o hapB/PshNP80052.liftoff_abinitio.final.hapB.gff3.stats

conda activate /home/groups/schwessinger/condaEnvs/busco
busco -i PshNP80052.liftoff_abinitio.final.proteins.fa -l basidiomycota -m protein -c 24
    # ---------------------------------------------------
    # |Results from dataset basidiomycota_odb12          |
    # ---------------------------------------------------
    # |C:94.4%[S:2.2%,D:92.2%],F:1.4%,M:4.2%,n:2409      |
    # |2275    Complete BUSCOs (C)                       |
    # |53    Complete and single-copy BUSCOs (S)         |
    # |2222    Complete and duplicated BUSCOs (D)        |
    # |34    Fragmented BUSCOs (F)                       |
    # |100    Missing BUSCOs (M)                         |
    # |2409    Total BUSCO groups searched               |
    # ---------------------------------------------------