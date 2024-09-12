#python src/extract_barcode_UMIs.py -i raw/test.sam -w raw/barcodes.tsv.gz -b ../design/barcodes.txt -l test.log -o /tmp/standard -t /tmp/standard.barcode.table.tsv.gz -n /tmp/non-standard.barcode.table.tsv.gz

basedir=/data/gc-core/taz2008/scRNAseq/10X_Jeya_14616_230629
workdir=${basedir}/barcode_UMIs_correction
srcdir=${workdir}/src
infodir=${workdir}/info
logdir=${srcdir}/logs

for sid in Jeya_1 Jeya_2 Jeya_3
do
    echo ${sid}
    python ${srcdir}/extract_barcode_UMIs.py -i ${workdir}/raw/${sid}.sam.gz -w ${basedir}/source/${sid}/sample_filtered_feature_bc_matrix/barcodes.tsv.gz -b ${basedir}/design/barcodes.txt -l ${logdir}/extract_barcode_UMIs.${sid}.log -o ${infodir}/${sid}.standard -t ${infodir}/${sid}.standard.table.tsv.gz -n ${infodir}/${sid}.non-standard.table.tsv.gz
    gzip ${infodir}/${sid}.standard.matrix.mtx
done

echo "Complete!"

