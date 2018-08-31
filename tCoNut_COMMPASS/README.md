# tCoNuT_COMMPASS
TGen Copy Number Tool version specific for MMRF CoMMpass study

## Summary
tCoNuT is a read depth based comparative copy number tool designed for whole genome, exome and panel NGS data. In addition, tCoNuT pipeline provides scripts for calculating B-allele frequencies.  B-allele frequencies can be used in conjunction with copy number to determine regions of LOH and provide additional evidence of a copy number change. The tool requires a control sample which can be matched or unmatched to the affected tumor sample. The tCoNuT_COMMPASS version of tCoNuT was modified to filter out regions of low-quality to remove an observed ‘waterfall’ effect, add gender aware adjustment to the log2 fold change for chromosome X and Y, and leverage high-coverage exome data to identify high-quality, common SNPs whose positions are used for centering low-coverage long-insert whole-genome data.

The tConut_COMMPASS workflow can be run with the included `ngs_cna2015_WGcenterWithExomesFilt2016_V3.sh` script.  

## System Requirements
The tCoNut_COMMPASS pipeline was developed and run on [CentOS Linux release 7.2.1511](http://vault.centos.org/7.2.1511/) (Core). Most of the scripts are platform independent. Uncompiled MATLAB code (*.m) found in tCoNuT_COMMPASS/pegasusCNA_MMRF
 folder is not platform dependent but would require a license of MATLAB to run. The compiled MATLAB code only requires the MCR.

* 6G of RAM
* 4 processing cores
* [R 3.4.1 or above](https://cran.r-project.org/) with the following packages:
	+ [DNAcopy](https://bioconductor.org/packages/release/bioc/html/DNAcopy.html)
	+ [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
	+ [reshape2](https://cran.r-project.org/web/packages/reshape2/README.html)
	+ [grid](https://stat.ethz.ch/R-manual/R-devel/library/grid/html/grid-package.html)
	+ [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html)
* [perl or above](https://www.perl.org/get.html) with the following modules:
	+ [Statistics::R](https://metacpan.org/pod/release/GMPASSOS/Statistics-R-0.02/lib/Statistics/R.pm)
* [MATLAB Runtime (MCR) v9.0](https://www.mathworks.com/products/compiler/matlab-runtime.html)
* [bedtools 2.26.0 or above](https://bedtools.readthedocs.io/en/latest/)
* [SnpSift 4.2](http://snpeff.sourceforge.net/SnpSift.html)

## Options

The `ngs_cna2015_WGcenterWithExomesFilt2016_V3.sh` workflow options are:

| Option  | Argument  | Required  | Description |
| ------- |:--------- |:---------:|:-------------- |
| -N | string |Yes|Exome normal sample name|
| -T | string |Yes|Exome tumor sample name|
| -b | file |Yes|Padded targets file|
| -n | file |Yes|Genome Normal Dat file|
| -t | file |Yes|Genome Tumor Dat file|
| -h | file |Yes|Exome Normal HS Metrics|
| -H | file |Yes|Exome Tumor HS Metrics|
| -v | file |Yes|Exome Haplotype Caller VCF|
| -m | path |Yes|MCR path|
| -d | file |Yes|dbSNP vcf|
| -s | file |Yes|SnpSift jar file|
| -c | file |Yes|CCDS gtf|
| -f | file |Yes|Copy Number targets file|
| -o | string |Yes|Out file prefix|

## Usage 
To get the arguments and options to the script, run:

`ngs_cna2015_WGcenterWithExomesFilt2016_V3.sh --help`

Prior to running the CoMMpass version of the tCoNuT pipeline, DAT files from the Tumor and Constitutional WGS bams, [Picard HS](https://broadinstitute.github.io/picard/) metrics and a dbSNP annotated [HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php) vcf from the WES bams with the -D option need to be generated. Additionally, the following files will need to be downloaded and provided as arguments to the parameterized workflow.

* [Agilent_SureSelect_V5_plusUTR_hs37d5_GRCh37.74_PaddedTargets_intersect_sorted_padded100.bed](http://tools.tgen.org/Files/CoMMpass_REFfiles/)
* [dbsnp_137.b37.vcf](http://tools.tgen.org/Files/CoMMpass_REFfiles/)
* [Homo_sapiens.GRCh37.74.gtf.hs37d5.EGFRvIII.gtf](http://tools.tgen.org/Files/CoMMpass_REFfiles/)
* [CNA_waterfall_filter_035.txt](http://tools.tgen.org/Files/CoMMpass_REFfiles/)

<b>Step 1 (prior to tCoNuT):</b> Create DAT files using `tgen_CloneCov.pl` for each WGS BAM separately.

```
${tCoNuTdir}/tgen_CloneCov.pl I=${BAMFILE} O=${OUTFILE} M=RG: S=${SAMTOOLS}
```

<b>Step 2 (prior to tCoNuT):</b> Create annotated HC vcf from the WES BAMS.

```
# Step 1 – Haplotype Caller the Match Constitutional and Tumor BAM
java -Djava.io.tmpdir=/scratch/tmp/ -jar -Xmx32g ${GATKPATH}/GenomeAnalysisTK.jar \
	-l INFO \
	-R ${REF} \
	-L ${CHRLIST}/Step${STEP}.list \
	-nct 8 \
	-T HaplotypeCaller \
	${BAMLIST} \
	-D ${KNOWN} \
	-mbq 10 \
	-o ${TRK}_Step${STEP}.HC.vcf 

# Step 2 – Concatenate all VCF's with GATK	
java -cp ${GATKPATH_v}/GenomeAnalysisTK.jar org.broadinstitute.sting.tools.CatVariants \
	-R ${REF} $vcfList \
	-out ${TRK}.HC_All.vcf \
	-assumeSorted
```

<b>Step 3 (prior to tCoNuT):</b> Create Picard HS metrics files.

```
java -Xmx15g -jar ${PICARDPATH}/CalculateHsMetrics.jar \
    REFERENCE_SEQUENCE=${REF} \
    BAIT_INTERVALS=${BAITS} \
    TARGET_INTERVALS=${TARGETS} \
    INPUT=${BAMFILE} \
    OUTPUT=${BAMFILE}.picHSMetrics \
    PER_TARGET_COVERAGE=${BAMFILE}.picStats.HsPerTargetCov \
    TMP_DIR=/scratch/tgenjetstream/tmp/ \
    VALIDATION_STRINGENCY=SILENT
```


<b>Step 4:</b> Run `ngs_cna2015_WGcenterWithExomesFilt2016_V3.sh` with the following options.

```
ngs_cna2015_WGcenterWithExomesFilt2016_V3.sh \
    -N MMRF_1157_1_PB_Whole_C1_KHS5U_L13424 \
    -T MMRF_1157_3_BM_CD138pos_T1_KBS5U_L06509 \
    -b /path/to/Agilent_SureSelect_V5_plusUTR_hs37d5_GRCh37.74_PaddedTargets_intersect_sorted_padded100.bed \
    -n MMRF_1157_3_PB_WBC_C3_KHWGL_L06740.proj.md.jr.bam.clc.cln.dat \
    -t MMRF_1157_3_BM_CD138pos_T1_KHWGL_L06741.proj.md.jr.bam.clc.cln.dat \
    -h MMRF_1157_1_PB_Whole_C1_KHS5U_L13424.proj.md.jr.bam.picHSMetrics \
    -H MMRF_1157_3_BM_CD138pos_T1_KBS5U_L06509.proj.md.jr.bam.picHSMetrics \
    -v MMRF_1157_1_PB_Whole_C1_KHS5U_L13424-MMRF_1157_3_BM_CD138pos_T1_KBS5U_L06509.HC_All.snpEff.vcf \
    -m /path/to/MCR/9.0 \
    -d /path/to/dbsnp_137.b37.vcf \
    -s /path/to/snpEff_4.2_2015-12-05/SnpSift.jar \
    -c /path/to/Homo_sapiens.GRCh37.74.gtf.hs37d5.EGFRvIII.gtf \
    -f /path/to/CNA_waterfall_filter_035.txt \
    -o MMRF_1157_3_PB_WBC_C3_KHWGL_L06740-MMRF_1157_3_BM_CD138pos_T1_KHWGL_L06741

```
