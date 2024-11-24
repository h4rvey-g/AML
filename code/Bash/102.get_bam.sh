# define dnb4tools is tools/dnbc4tools2.1.3/dnbc4tools
alias dnbc4tools='tools/dnbc4tools2.1.3/dnbc4tools'
# unzip /workspaces/AML/data/reference_data/alias/hg38/ensembl_gtf/default/hg38.gtf.gz to data/reference_data/gencode.v47.annotation.gtf
gunzip -c data/reference_data/alias/hg38/gencode_gtf/default/hg38.gtf.gz > data/reference_data/gencode.v47.annotation.gtf
dnbc4tools tools mkgtf --action check --ingtf data/reference_data/gencode.v47.annotation.gtf --output data/reference_data/corrected.gtf
dnbc4tools tools mkgtf --ingtf data/reference_data/corrected.gtf --output data/reference_data/corrected_filter.gtf --type gene_type
dnbc4tools rna mkref --fasta data/reference_data/alias/hg38_primary/fasta/default/hg38_primary.fa \
--ingtf data/reference_data/corrected_filter.gtf --species Homo_sapiens \
--threads 30 --genomeDir data/reference_data/index
dnbc4tools rna multi --list data/101.raw_data/samplesheet.tsv --genomeDir data/reference_data/index --threads 40 \
--outdir data/103.self_workflow
# run bash code/Bash/N1.sh, N2.sh, N4.sh, T1.sh, T2.sh, T4.sh sequentially
for i in 1 2 4; do bash code/Bash/N${i}.sh; done
for i in 1 2 4; do bash code/Bash/T${i}.sh; done

process_bam() {
    bam=$1
    sample=$(basename $(dirname $(dirname "$bam")))
    barcodes="data/103.self_workflow/${sample}/02.count/filter_matrix/barcodes.tsv.gz"
    dnbc4tools tools changetag --inbam "$bam" --outbam "${bam%.bam}.velocyto.bam"
    # rm data/103.self_workflow/${sample}/velocyto content
    rm -rf data/103.self_workflow/${sample}/velocyto/*
    velocyto run -b "$barcodes" -o "data/103.self_workflow/${sample}/velocyto" -m data/reference_data/hg38_rmsk.gtf "${bam%.bam}.velocyto.bam" data/reference_data/corrected_filter.gtf
    # rename the output loom file to the sample name
    find "data/103.self_workflow/${sample}/velocyto" -name "*.loom" -exec mv {} "data/103.self_workflow/${sample}/velocyto/${sample}.loom" \;
}

export -f process_bam
find data/103.self_workflow -name "anno_decon_sorted.bam" | parallel --progress --keep-order --line-buffer process_bam

# test run STARsolo on 