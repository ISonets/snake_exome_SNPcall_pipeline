# snake_exome_SNPcall_pipeline
Snakemake pipeline for SNP and INDEL calling using GATK4 and homemade bash and AWK scripts producing custom output.

## Intro
This pipeline is made in Lab of Mathematical Biology and Bioinformatics of RISBM to to a specific task: produce fast, accurate and annotated with rsid and MAF data SNP calling result tables in ```.txt``` format for specific panels of SNP.  
The biggest difference is that this pipeline also outputs info about SNPs from the panel not only found in each sample, but similar info to SNPs from the panel not found in sample.  
So, pipeline produces 2 files: one with SNPs found in sample, and file with SNPs not found in this sample, and we have info about all SNPs from the panel at once. 

However, at this moment pipeline does **not** produce any plots or visualizations.
 
## Prerequisites
  - Conda([link](https://conda.io/projects/conda/en/latest/index.html)) and mamba ([link](https://github.com/mamba-org/mamba));
  - GATK4 (see [link](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4/) how to install it);
  - hg38 reference genome (see [link](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/) to download it)
  - dbSNP database from NCBI to SNP rsid annotation (see [link](https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz) to FTP download);
  - 1000 Genomes Phase3 VCF file to obtain MAF info (see [link](https://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz)to FTP download);
  - Snakemake(see [link](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) to install it);
  - ```bwa``` and its dependencies (see [link](https://bio-bwa.sourceforge.net/) to download);
  - ```samtools``` and its dependencies (see [link](http://www.htslib.org/) to download).

*Conda* is mandatory for *Snakemake* installation, and we also recommend to install *bwa* and *samtools* in *conda*, however *GATK* installation in *conda* was unstable in our hands, so to achieve maximum stability we recommend to install *GATK* as a standalone version on your machine.

## SNP panel creation
To use the pipeline, you need to get or create 2 files:
1) ```.bed``` file with coordinates of SNP. Its structure should be like this, also according to ```.bed``` format specifications:

| chromosome | start | stop | rsid |
|------------|-------|------|------|
| chr1       | 1000  | 1001 | rs1  |  

etc.

2) 1-column no-header```.txt``` file with rsid of SNPs:

```rs1```  
```rs2```  
```rs3``` 
etc.

## Usage
Before launching this pipeline, you need to modify ```.py``` file, providing absolute paths to :
1) data folder,
2) hg38 location,
3) dbSNP location,
4) 1000 Genomes VCF,
5) ```.bed``` file with SNPs panel coordinates,
6)  ```.txt``` file with rsids of SNPs.

Also you need to index dbSNP ```.vcf.gz``` file and 1000 Genomes ```.vcf.gz``` file with ```tabix``` (part of ```samtools```):  

```tabix -p vcf {vcf.gz}```

And finally you need to create a ```.dict``` (dictionary)  and ```.fai``` (fasta index) files for hg38 reference genome assembly with GATK:  

```gatk CreateSequenceDictionary R={hg38_fastq} O={hg38_dict}```

```samtools faidx {hg38_fastq}``` 

Executing:

```snakemake -s snake_exomes_SNPcall_paired_reads.py -c {cores} ```

## Output format
Output format of this pipeline for SNPs found in this sample looks like this:

| chromosome | coordinate | rsid      | reference_allele | alternative_allele | genotype | depth_ref_allele/depth_alt_allele | total_depth | annotation                                                                                                                                                                                      |
|------------|------------|-----------|------------------|--------------------|----------|-----------------------------------|-------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| chr1       | 67215986   | rs7517847 | T                | G                  | 0/1      | 55/51                             | 106         | dbSNP_151;TSA=SNV;E_Freq;E_Hapmap;E_1000G;E_Cited;E_Phenotype_or_Disease;E_TOPMed;E_gnomAD;MA=G;MAF=0.35603;MAC=1783;AA=T;EAS_AF=0.4236;EUR_AF=0.4205;AMR_AF=0.5663;SAS_AF=0.3313;AFR_AF=0.1634 |

Output format for SNPs not found in this sample is ths same, but has some differences: **genotype** column will always be 0/0 (because only reference allele is found, so no SNP is present).

## Future plans
  - add coverage plots;
  - impore user experience with simpifying installation;
  - migrate from AWK to Python or R.

## Authors
Ignat Sonets, RISBM
