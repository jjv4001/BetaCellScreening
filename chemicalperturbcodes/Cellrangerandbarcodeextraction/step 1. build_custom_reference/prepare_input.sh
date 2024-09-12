
mkdir -p temp/fasta
mkdir -p temp/genes

cat  ../refdata-gex-GRCh38-2020-A/fasta/genome.fa  bfp/bfp.fa  >temp/fasta/genome.fa
cat  ../refdata-gex-GRCh38-2020-A/genes/genes.gtf  bfp/bfp.gtf  >temp/genes/genes.gtf

cellranger=/data/gc-core/sbsuser/software/cellranger-7.1.0/cellranger

${cellranger} mkref --genome=GRCh38-2020-A_bfp \
                 --fasta=temp/fasta/genome.fa \
                 --genes=temp/genes/genes.gtf

rm -rf temp/

