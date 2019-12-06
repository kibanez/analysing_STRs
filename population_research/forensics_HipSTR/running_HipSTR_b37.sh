OUTPUT_FOLDER=/genomes/scratch/kgarikano/GEL_STR/population/HipSTR_output/
HIPSTER_BINARY=/genomes/scratch/kgarikano/GEL_STR/HipSTR/sw/HipSTR_tweak/HipSTR/HipSTR

# Corresponding fasta and forensics regions files for GRCh37
FASTA_FILE=/genomes/resources/genomeref/Illumina/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa
REGIONS_FILE=/genomes/scratch/kgarikano/GEL_STR/population/specs/GRCh37.hipstr_reference_forensic.bed

if [ $# -lt 1 ]
   then
     echo "A list containing <LP_ID, BAM_PATH> is required to run through HipSTR"
     exit
fi

INPUT_LIST=$1

cat ${INPUT_LIST} | while read line; do

 IFS=?~@~Y,?~@~Y read -ra NAMES <<< "$line"
 LP_ID=${NAMES[0]}
 PATH_BAM=${NAMES[1]}

 OUTPUT_VCF=${OUTPUT_FOLDER}'HipSTR_'${LP_ID}'.vcf.gz'
 OUTPUT_LOG=${OUTPUT_FOLDER}'HipSTR_'${LP_ID}'.log'
 OUTPUT_VIZ=${OUTPUT_FOLDER}'HipSTR_'${LP_ID}'.viz.gz'

 TAG_ID='ID'
 MIN_READS=15
 
${HIPSTER_BINARY} \
 --bams ${PATH_BAM} \
 --fasta ${FASTA_FILE} \
 --regions ${REGIONS_FILE} \
 --str-vcf ${OUTPUT_VCF} \
 --log ${OUTPUT_LOG} \
 --viz-out ${OUTPUT_VIZ} \
 --lib-field ${TAG_ID} \
 --min-reads ${MIN_READS} --def-stutter-model

done
