#!/usr/bin/bash --login

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10GB #requirements may be high in case user downloads large images
#SBATCH --partition=work
#SBATCH --time=1-00:00:00
#SBATCH --account=pawsey1142
#SBATCH --mail-type=ALL
#SBATCH --export=NONE

module load nextflow/25.04.6
module load singularity/4.1.0-slurm

SOURCE_PATH=/scratch/pawsey1142/kgagalova/norin-cadenza/wheat-genomes/
TARGET_PATH=/scratch/pawsey1142/kgagalova/norin-cadenza/wheat-genomes/
#EXPR_PATH=/scratch/pawsey1142/kgagalova/norin-cadenza/expression/

export PAWSEY_PROJECT='pawsey1142'

nextflow run . \
    --source_genome     $SOURCE_PATH/Triticum_aestivum_norin61.PGSBv2.1.dna_sm.toplevel.fa \
    --target_genome     $TARGET_PATH/Triticum_aestivum_cadenza.EIv1.1.dna_sm.toplevel.fa \
    --source_gff        $SOURCE_PATH/Triticum_aestivum_norin61.PGSBv2.1.62.gff3 \
    --target_gff        $TARGET_PATH/Triticum_aestivum_cadenza.EIv1.1.62.gff3 \
    --source_proteins   $SOURCE_PATH/Triticum_aestivum_norin61.PGSBv2.1.62.pep.fa \
    --target_proteins   $TARGET_PATH/Triticum_aestivum_cadenza.EIv1.1.62.pep.fa \
    --qtl_bed           qtl_intervals.bed \
    --source_expression RawCountsStar-Norin.tsv \
    --target_expression RawCountsStar-Cadenza.tsv \
    --source_name       Norin \
    --target_name       Cadenza \
    --outdir            ./results \
    -profile pawsey_setonix,singularity \
    --source_annotation Annotations-final-norin.tsv \
    --target_annotation Annotations-final-cadenza.tsv \
    -resume
