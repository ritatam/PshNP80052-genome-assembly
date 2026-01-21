#############################
##### merge annotations #####
#############################

# merge liftoff and ab initio predictions together
agat_sp_merge_annotations.pl -f PshNP85002.104e-liftoff.gff3 -f PshNP85002.az2-liftoff.gff3 -f PshNP85002.abinitio.gff3 -o PshNP85002.liftoff_abinitio_combined.tmp

# rename and number gene ids
funannotate util gff-rename -g PshNP85002.liftoff_abinitio_combined.tmp -f PshNP85002.final.fasta -o PshNP85002.liftoff_abinitio_combined.gff3 --locus_tag PshNP85002 --numbering 1

# deduplicate cds
python deduplicate_cds_gff.py -i PshNP85002.liftoff_abinitio_combined.gff3 --dedup_cds_gff3_out PshNP85002.liftoff_abinitio_combined.dedup.cds.gff3 --dedup_full_gff3_out PshNP85002.liftoff_abinitio_combined.dedup.gff3

# renumber after deduplication
funannotate util gff-rename -g PshNP85002.liftoff_abinitio_combined.dedup.tmp -f PshNP85002.final.fasta -o PshNP85002.liftoff_abinitio_combined.dedup.gff3 --locus_tag PshNP85002 --numbering 1
sed -i 's/;Alias=.*//g' PshNP85002.liftoff_abinitio_combined.dedup.gff3
agat_sp_extract_sequences.pl -g PshNP85002.liftoff_abinitio_combined.dedup.gff3 -f PshNP85002.final.fasta -o PshNP85002.liftoff_abinitio_combined.dedup.proteins.fa -t cds -p

# keep longest isoform per gene
agat_sp_keep_longest_isoform.pl -gff PshNP85002.liftoff_abinitio_combined.dedup.gff3 -o PshNP85002.liftoff_abinitio_combined.dedup.longest_orf.gff3
agat_sp_extract_sequences.pl -g PshNP85002.liftoff_abinitio_combined.dedup.longest_orf.gff3 -f PshNP85002.final.fasta -o PshNP85002.liftoff_abinitio_combined.dedup.longest_orf.protein.fa -t cds -p

# filter out proteins shorter than 50bp. only 9 genes need to go.
# i can't find ways to do this directly on gff3, so sorted proteins from shortest to longest and manually delete from the top in vim
funannotate gff2prot -g PshNP85002.liftoff_abinitio_combined.dedup.longest_orf.gff3 -f PshNP85002.final.fasta --no_stop > PshNP85002.liftoff_abinitio_combined.dedup.longest_orf.protein.fa
seqkit sort -l PshNP85002.liftoff_abinitio_combined.dedup.longest_orf.protein.fa | less
cp PshNP85002.liftoff_abinitio_combined.dedup.longest_orf.gff3 PshNP85002.liftoff_abinitio_combined.dedup.longest_orf.protein_50aa.gff3
vim PshNP85002.liftoff_abinitio_combined.dedup.longest_orf.protein_50aa.gff3  

# separate filtered annotations into two haplotypes
awk '
$0 ~ /^#/ {
    print > "PshNP85002.liftoff_abinitio_combined.dedup.longest_orf.prot_50aa.hapA.gff3"
    print > "PshNP85002.liftoff_abinitio_combined.dedup.longest_orf.prot_50aa.hapB.gff3"
    next
}
$1 ~ /A$/ { print > "PshNP85002.liftoff_abinitio_combined.dedup.longest_orf.prot_50aa.hapA.gff3" }
$1 ~ /B$/ { print > "PshNP85002.liftoff_abinitio_combined.dedup.longest_orf.prot_50aa.hapB.gff3" }
' PshNP85002.liftoff_abinitio_combined.dedup.longest_orf.protein_50aa.gff3

funannotate util gff-rename -g PshNP85002.liftoff_abinitio_combined.dedup.longest_orf.prot_50aa.hapA.gff3 -f PshNP85002.final.fasta -o PshNP85002.liftoff_abinitio_combined.dedup.longest_orf.prot_50aa.renum.hapA.gff3 --locus_tag PshNP85002 --numbering 1
funannotate util gff-rename -g PshNP85002.liftoff_abinitio_combined.dedup.longest_orf.prot_50aa.hapB.gff3 -f PshNP85002.final.fasta -o PshNP85002.liftoff_abinitio_combined.dedup.longest_orf.prot_50aa.renum.hapB.gff3 --locus_tag PshNP85002 --numbering 17245

sed -i 's/;Alias=.*//g' PshNP85002.liftoff_abinitio_combined.dedup.longest_orf.prot_50aa.renum.hapA.gff3
sed -i 's/;Alias=.*//g' PshNP85002.liftoff_abinitio_combined.dedup.longest_orf.prot_50aa.renum.hapB.gff3

cat PshNP85002.liftoff_abinitio_combined.dedup.longest_orf.prot_50aa.renum.hapA.gff3 PshNP85002.liftoff_abinitio_combined.dedup.longest_orf.prot_50aa.renum.hapB.gff3 > PshNP85002.liftoff_abinitio_combined.dedup.longest_orf.prot_50aa.renum.gff3

# clip UTRs from reference annotations
python clip_UTR.py -i PshNP85002.liftoff_abinitio_combined.dedup.longest_orf.prot_50aa.renum.gff3 -o PshNP85002.liftoff_abinitio_combined.dedup.longest_orf.prot_50aa.renum.utr_clipped.gff3

# extract tRNA from original gff3 (will be added back later)
grep -P "\ttRNA\t" PshNP85002.liftoff_abinitio_combined.dedup.gff3 | awk -F'Parent=' '{print $2}' | cut -d';' -f1 > PshNP85002.liftoff_abinitio_combined.dedup.tRNA.ids
grep -F -f PshNP85002.liftoff_abinitio_combined.dedup.tRNA.ids PshNP85002.liftoff_abinitio_combined.dedup.gff3 > PshNP85002.liftoff_abinitio_combined.dedup.tRNA.gff3
awk '
$0 ~ /^#/ {
    print > "PshNP85002.liftoff_abinitio_combined.dedup.hapA.tRNA.gff3"
    print > "PshNP85002.liftoff_abinitio_combined.dedup.hapB.tRNA.gff3"
    next
}
$1 ~ /A$/ { print > "PshNP85002.liftoff_abinitio_combined.dedup.hapA.tRNA.gff3" }
$1 ~ /B$/ { print > "PshNP85002.liftoff_abinitio_combined.dedup.hapB.tRNA.gff3" }
' PshNP85002.liftoff_abinitio_combined.dedup.tRNA.gff3