Analysis pipeline
====================

## Folder structure

Generally there are following folders:

```
BASE_DIR
|- data
|- enrichedheatmap
|- gviz
|- plots
|- rds
|- rds_methylation
|- scripts
```

-data raw data from Roadmap or UCSC
-enrichedheatmap plots
-gviz plots
-plots general plots
-rds intermediate results
-rds_methylation
-scripts 

## Data

All the data files are stored under `BASE_DIR/data/` folder.

- Sample annotation table (https://docs.google.com/spreadsheets/d/1yikGx4MsO9Ei36b64yOy9Vb6oPC5IBGlFbYEt-N6gOM/edit#gid=15).
  Column B (EID), D (Group), E (Color), R (Type) are used as sample annotations. 

- Methylation data: The methylation rate (http://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/WGBS/FractionalMethylation.tar.gz) and CpG coverage (http://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/WGBS/ReadCoverage.tar.gz) are stored per chromosome. The column names for the methylation rate matrix and CpG coverage matrix are obtained from (http://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/WGBS/EG.mnemonics.name.xls)

- Expression data: Raw counts of all protein coding genes (http://egg2.wustl.edu/roadmap/data/byDataType/rna/expression/57epigenomes.N.pc.gz) and corresponding RPKM values (http://egg2.wustl.edu/roadmap/data/byDataType/rna/expression/57epigenomes.RPKM.pc.gz). The transcriptome annotation is from http://egg2.wustl.edu/roadmap/data/byDataType/rna/annotations/gen10.long.gtf.gz.

- Histone modification data: MACS peak calling results are used: http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/. ChromHMM segmentation based on 15 models http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/all.mnemonics.bedFiles.tgz.

- CpG islands: cpgIslandExt.bed is from [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=512073311_jGgG4W69t5jSJAwKzoyqnGOoCgNB&clade=mammal&org=&db=hg19&hgta_group=regulation&hgta_track=cpgIslandExt&hgta_table=cpgIslandExt&hgta_regionType=genome&position=&hgta_outputType=primaryTable&hgta_outFileName=)

## Preprocess

### filter transcriptome

In the analysis, we only look at protein coding genes. So we first filter out genes and transcripts which are not 
protein coding in the transcriptome annotation file.

```
cd data
sh ../scripts/filter_protein_coding.sh
```

A new file `gen10_long_protein_coding.gtf` is generated in the `data/` folder.

### build txdb

In the analysis, the transcriptome is stored as a `Txdb` object. So here we need to build
it first.

```
Rscript generate_txdb.R
```

A new file `gen10_long_protein_coding.sqlite` is generated in the `data/` folder.

### merge methylation and do smoothing

We observed that some of the samples have too many missing methylation values (represented as -1 in the methylation
files and CpG coverage files). We first calculate NA rate in each sample:

```
cd ../scripts
Rscript all_sample_coverage_NA_rate.R
```

After looking at the plot, we set the cutoff of NA rate to 10% and removed E054 (14%), E070 (%13), E084 (18%) and E085 (41%).

We only look at the methylation on CpG dinucleitide while in the methylation file, methylation
rate is measured on the C on both strand. We calcualte the mean methylation rate of the CpG
weighted by the coverage. CpG sites are removed if in more than 50% of the samples, they have CpG coverage <= 2.
finally we use bsseq to smooth the methylation data.

Following script is executed for each chromosome, e.g.

```
Rscript make_rds_for_methylation_data.R --chr chr1
```

`BASE_DIR/rds_methylation/*_roadmap_merged_bsseq.rds` is generated.

## configuration file
