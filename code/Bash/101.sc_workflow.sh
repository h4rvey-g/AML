# # download and unzip https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.annotation.gtf.gz, save to data/reference_data/gencode.vM36.annotation.gtf
# wget -P data/reference_data/ https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.annotation.gtf.gz
# gunzip data/reference_data/gencode.vM36.annotation.gtf.gz
# # download and unzip https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.transcripts.fa.gz
# wget -P data/reference_data/ https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.transcripts.fa.gz
# gunzip data/reference_data/gencode.vM36.transcripts.fa.gz

# nextflow run nf-core/scrnaseq \
#    -profile docker \
#    --input data/101.raw_data/samplesheet.csv \
#    --genome_fasta data/reference_data/gencode.vM36.transcripts.fa \
#    --gtf data/reference_data/gencode.vM36.annotation.gtf \
#    --protocol 10XV2 \
#    --aligner star \
#    --outdir ./data/102.sc_workflow_output \
#    -resume \
#    --save_reference true \
#    --star_feature "Gene Velocyto"
for sample in data/102.expr_data/*; do
    if [ -d "$sample/filter_matrix" ]; then
        mv "$sample/filter_matrix/"* "$sample/"
        rmdir "$sample/filter_matrix"
    fi
done
