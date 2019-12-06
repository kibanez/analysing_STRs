output_folder="/genomes/scratch/kgarikano/GEL_STR/HipSTR/twins_analysis/HipSTR_output/"
hipster_binary="/genomes/scratch/kgarikano/GEL_STR/HipSTR/sw/HipSTR/HipSTR"

# Corresponding fata and regions files for GRCh38
fasta_file="/genomes/resources/genomeref/Illumina/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa"
regions_file="/genomes/scratch/kgarikano/GEL_STR/HipSTR/twins_analysis/hg38.hipstr_reference_forensic.bed"

if [ $# -lt 1 ]
  then
    echo "A list containing BAM files to run through HipSTR is required"
    exit
fi

input_list=$1

cat $input_list | while read line; do

    IFS=?~@~Y,?~@~Y read -ra NAMES <<< "$line"
    lp_id=${NAMES[0]}
    path_to_bam=${NAMES[1]}

    output_vcf=${output_folder}'HipSTR_'${lp_id}'.vcf.gz'
    output_log=${output_folder}'HipSTR_'${lp_id}'.log'
    output_viz=${output_folder}'HipSTR_'${lp_id}'.viz.gz'

    ${hipster_binary} \
    --bams ${lp_id} \
    --fasta ${fasta_file} \
    --regions ${regions_file} \
    --str-vcf  ${output_vcf} \
    --log ${output_log} \
    --viz-out ${output_viz} \
    --lb-tag ID \
    --min-reads 15 \
    --def-stutter-model 