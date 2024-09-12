# BetaCellScreening

##FixedRNAprofiling

Tissueseurat.R provides the code for Seurat analysis of the fixed RNA profiling dataset obtained from EndoC-betaH1 grafts subcutaneously transplanted into male and female mice and retrieved after 1 week.

DEseqFLEX.R performs differential expresison analysis of LIP preconditioned grafts transplanted in female versus LIP preconditioned grafts transplanted in male. 

tissues.csv provides the metadata for the assignment of cells to Female Control, Female Treated, Male Control, and Male Treated grafts. 

##chemperturbcodes 

- Cellrangerandbarcodeextraction provides codes to build a custom reference containing lentiviral barcodes, to perform CellRanger analysis, and to extract lentiviral barcode assignments with a maximum of one mismatch.

- Seurat analysis contains codes for the Seurat analysis of the assigned chemical perturb sequencing dataset.

- DEG analysis was used to perform differential expression analysis of small molecule treated samples against untreated control samples.


