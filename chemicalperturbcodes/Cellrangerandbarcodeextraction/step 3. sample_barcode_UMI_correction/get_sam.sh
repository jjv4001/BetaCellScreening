
regions=''

for gid in `cut -f 1 /data/gc-core/sbsuser/10X_Refs/refdata-custom-gex-GRCh38-2020-A-jeya-v4/bfp/bfp.gtf`
do
	regions="${regions} ${gid}:1-250"
done

for sid in Jeya_1 Jeya_2 Jeya_3
do
	echo ${sid}
	samtools view Chen-XD-14510_230616/outs/per_sample_outs/${sid}/count/sample_alignments.bam ${regions} >${sid}.sam
	gzip ${sid}.sam
done

echo "Complete!"

