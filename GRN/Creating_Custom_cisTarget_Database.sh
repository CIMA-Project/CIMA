git clone https://github.com/aertslab/create_cisTarget_databases
wget https://resources.aertslab.org/cistarget/programs/cbust
chmod a+x cbust

mkdir -p aertslab_motif_colleciton
wget -O aertslab_motif_colleciton/v10nr_clust_public.zip https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/v10nr_clust_public.zip
cd aertslab_motif_colleciton; unzip -q v10nr_clust_public.zip

module load cluster/wice/bigmem
module load BEDTools/2.30.0-GCC-10.3.0

REGION_BED="CIMA/scATAC/Cell_type_l4_Peak.bed"
GENOME_FASTA="CIMA/genomes/hg38.fa"
CHROMSIZES="CIMA/genomes/hg38.chrom.sizes"
SCRIPT_DIR="CIMA/scenicplus/create_cisTarget_databases"

${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh \
        ${GENOME_FASTA} \
        ${CHROMSIZES} \
        ${REGION_BED} \
        hg38.DNBeLab_C4_PBMC.with_1kb_bg_padding.fa \
        1000 \
        yes


ls aertslab_motif_colleciton/v10nr_clust_public/singletons > motifs.txt

OUT_DIR="CIMA/scenicplus"
CBDIR="${OUT_DIR}/aertslab_motif_colleciton/v10nr_clust_public/singletons"
FASTA_FILE="${OUT_DIR}/hg38.DNBeLab_C4_PBMC.with_1kb_bg_padding.fa"
MOTIF_LIST="${OUT_DIR}/motifs.txt"

"${SCRIPT_DIR}/create_cistarget_motif_databases.py" \
    -f ${FASTA_FILE} \
    -M ${CBDIR} \
    -m ${MOTIF_LIST} \
    -o ${OUT_DIR}/${DATABASE_PREFIX} \
    --bgpadding 1000 \
    -t 96