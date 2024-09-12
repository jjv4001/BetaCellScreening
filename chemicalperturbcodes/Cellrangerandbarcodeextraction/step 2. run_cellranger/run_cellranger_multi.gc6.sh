#!/usr/bin/env sh

###cellranger=/data/gc-core/sbsuser/software/cellranger-5.0.0/cellranger
cellranger=/data/gc-core/sbsuser/software/cellranger-7.1.0/cellranger
##cellranger=/data/gc-core/sbsuser/software/cellranger-6.1.1/cellranger

if [ ! -d logs ]
then
	mkdir logs
fi

echo -n "processing..."
${cellranger} multi --id=Chen-XD-14510_230616 \
		--csv=/scratch/seq_data/NovaSeq6000/230616_A00814_0770_AH7H2NDSX7/Analysis.5/CellPlex.csv \
                --localcores=12 \
                --localmem=128 \
                >logs/cellranger_multi.Chen-XD-14510_230616.log 2>&1
echo "ok."

