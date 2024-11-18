# define dnb4tools is tools/dnbc4tools2.1.3/dnbc4tools
alias dnbc4tools='tools/dnbc4tools2.1.3/dnbc4tools'
dnbc4tools tools mkgtf --action check --ingtf data/reference_data/gencode.vM36.annotation.gtf --output data/reference_data/corrected.gtf
dnbc4tools tools mkgtf --ingtf data/reference_data/corrected.gtf --output data/reference_data/corrected_filter.gtf --type gene_type \
    --include protein_coding
dnbc4tools rna mkref --fasta data/reference_data/GRCm39.primary_assembly.genome.fa \
    --ingtf data/reference_data/corrected_filter.gtf --species Mus_musculus \
    --threads 30 --genomeDir data/reference_data/index
dnbc4tools rna multi --list data/101.raw_data/samplesheet.tsv --genomeDir data/reference_data/index --threads 40