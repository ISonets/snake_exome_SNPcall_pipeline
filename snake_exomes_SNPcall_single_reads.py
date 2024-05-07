#!/usr/bin/env python3
base = "/mnt/disk1/PROJECTS/EXOMS/DATA/external_data/poland_IBD_IonTorrent/ERR1341858/"
suffix_1 = ".fastq"
#suffix_2 = "_2.fastq.gz"

SAMPLES, = glob_wildcards(base +"{sample}"+suffix_1)
want_all = []
want_all.append(expand("res_snake/aln_res/{sample}_aligned_reads.sam",sample=SAMPLES))
want_all.append(expand("res_snake/aln_res/{sample}_dedup_reads.bam",sample=SAMPLES))
want_all.append(expand("res_snake/aln_res/{sample}_dedup_reads_sorted.bam",sample=SAMPLES))
want_all.append(expand("res_snake/vcf/{sample}_raw_variants.vcf",sample=SAMPLES))
want_all.append(expand("res_snake/vcf/{sample}_raw_snps.vcf",sample=SAMPLES))
want_all.append(expand("res_snake/vcf/{sample}_filtered_snps.vcf",sample=SAMPLES))
want_all.append(expand("res_snake/vcf/{sample}_bqsr_snps.vcf",sample=SAMPLES))
want_all.append(expand("res_snake/aln_res/{sample}_recal_reads.bam",sample=SAMPLES))
want_all.append(expand("res_snake/vcf/{sample}_raw_variants_recal.vcf",sample=SAMPLES))
want_all.append(expand("res_snake/vcf/{sample}_raw_snps_recal.vcf",sample=SAMPLES))
want_all.append(expand("res_snake/vcf/{sample}_final_snps.vcf",sample=SAMPLES))
want_all.append(expand("res_snake/vcf/{sample}_final_snps_ann.vcf",sample=SAMPLES))
want_all.append(expand("res_snake/vcf/{sample}_final_snps_ann_af.txt",sample=SAMPLES))
want_all.append(expand("res_snake/vcf/{sample}_unmatched_res.txt",sample=SAMPLES))

rule all:
    input: want_all

rule bwa_aln:  # alignment, specifying header for proper GATK
    input: read1 = base+"{sample}"+suffix_1,
           #read2 = base+"{sample}"+suffix_2
    output: bwa_aln_res = "res_snake/aln_res/{sample}_aligned_reads.sam",
    params: hg38 = "/mnt/disk1/PROJECTS/IMMUNE/DATA/hg38/hg38.fa"
    threads: 32
    shell: """
           bwa mem -Y -t {threads} -R '@RG\\tID:exome_sample\\tLB:sample\\tPL:DNBSEQ\\tSM:exome_sample' {params.hg38} {input.read1} > {output.bwa_aln_res}
    """

rule mark_dupl:  # deduplication
    input: bwa_aln_res = "res_snake/aln_res/{sample}_aligned_reads.sam",
    output: dedup_reads = "res_snake/aln_res/{sample}_dedup_reads.bam",
            dedup_metrics = "res_snake/aln_res/{sample}_dedup_metrics.txt"
    shell: """
           gatk MarkDuplicates -I {input.bwa_aln_res} -M {output.dedup_metrics}	-O {output.dedup_reads}	--ASSUME_SORT_ORDER queryname
    """

rule aln_idx:  # indexing .sam, making .bams, depth calculation, coverage calculation, basic stats
    input: dedup_reads = "res_snake/aln_res/{sample}_dedup_reads.bam",
    output: dedup_reads_sorted = "res_snake/aln_res/{sample}_dedup_reads_sorted.bam",
            dedup_reads_idx = "res_snake/aln_res/{sample}_dedup_reads_sorted.bam.bai",
            dedup_reads_depth = "res_snake/aln_res/{sample}_dedup_reads.depth",
            dedup_reads_cov = "res_snake/aln_res/{sample}_dedup_reads_cov.bedgraph",
            dedup_reads_flagstat="res_snake/aln_res/{sample}_flagstat_res.txt"
    params: hg38 = "/mnt/disk1/PROJECTS/IMMUNE/DATA/hg38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta"
    threads : 32
    shell: """
           samtools sort {input.dedup_reads} -@ {threads} > {output.dedup_reads_sorted}
           samtools index {output.dedup_reads_sorted} -@ {threads} -o {output.dedup_reads_idx}
           samtools depth {output.dedup_reads_sorted} -aa -o {output.dedup_reads_depth}
           samtools flagstat {output.dedup_reads_sorted} > {output.dedup_reads_flagstat}
           bedtools genomecov -ibam {output.dedup_reads_sorted} -g {params.hg38} -bg > {output.dedup_reads_cov}
    """

rule call_1st_round:  # 1st round of SNP/indel calling
    input: dedup_reads_sorted = "res_snake/aln_res/{sample}_dedup_reads_sorted.bam",
    output: raw_vcf = "res_snake/vcf/{sample}_raw_variants.vcf",
    params: hg38 = "/mnt/disk1/PROJECTS/IMMUNE/DATA/hg38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta"
    shell: """
           gatk HaplotypeCaller -R {params.hg38} -I {input.dedup_reads_sorted} -O {output.raw_vcf}
    """

rule extract_SNPs_indels_raw:  # extract only SNP/indels from raw shit
    input: raw_variants = "res_snake/vcf/{sample}_raw_variants.vcf",
    output: raw_snps = "res_snake/vcf/{sample}_raw_snps.vcf",
            raw_indels = "res_snake/vcf/{sample}_raw_indels.vcf"
    params: hg38 = "/mnt/disk1/PROJECTS/IMMUNE/DATA/hg38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta"
    shell: """
           gatk SelectVariants -R {params.hg38} -V {input.raw_variants} -select-type SNP -O {output.raw_snps}
           gatk SelectVariants -R {params.hg38} -V {input.raw_variants} -select-type INDEL -O {output.raw_indels}
    """

rule snps_indels_filtration:  # mild filtering based on GATK/internet manuals
    input: raw_snps = "res_snake/vcf/{sample}_raw_snps.vcf",
           raw_indels = "res_snake/vcf/{sample}_raw_indels.vcf"
    output: filtered_snps = "res_snake/vcf/{sample}_filtered_snps.vcf",
            filtered_indels = "res_snake/vcf/{sample}_filtered_indels.vcf"
    params: hg38 = "/mnt/disk1/PROJECTS/IMMUNE/DATA/hg38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta"
    shell: """
           gatk VariantFiltration -R {params.hg38} -V {input.raw_snps} -O {output.filtered_snps} --filter-name "QD_filter" -filter "QD < 2.0" --filter-name "FS_filter" -filter "FS > 60.0" --filter-name "MQ_filter" -filter "MQ < 40.0" --filter-name "SOR_filter" -filter "SOR > 4.0"
           gatk VariantFiltration -R {params.hg38} -V {input.raw_indels} -O {output.filtered_indels} --filter-name "QD_filter" -filter "QD < 2.0" --filter-name "FS_filter" -filter "FS > 200.0" --filter-name "SOR_filter" -filter "SOR > 10.0"
    """

rule exclude_bad_variants:  # preparing for BQSR
    input: filtered_snps = "res_snake/vcf/{sample}_filtered_snps.vcf",
           filtered_indels = "res_snake/vcf/{sample}_filtered_indels.vcf"
    output: bqsr_snps = "res_snake/vcf/{sample}_bqsr_snps.vcf",
            bqsr_indels = "res_snake/vcf/{sample}_bqsr_indels.vcf"
    shell: """
           gatk SelectVariants 	--exclude-filtered -V {input.filtered_snps} -O {output.bqsr_snps}
           gatk SelectVariants  --exclude-filtered -V {input.filtered_indels} -O {output.bqsr_indels}
    """

rule bqsr_gatk:  # BQSR
    input: bqsr_snps = "res_snake/vcf/{sample}_bqsr_snps.vcf",
           bqsr_indels = "res_snake/vcf/{sample}_bqsr_indels.vcf",
           dedup_reads_sorted = "res_snake/aln_res/{sample}_dedup_reads_sorted.bam"
    output: recal_table = "res_snake/aln_res/{sample}_recal_data.table",
            recal_reads = "res_snake/aln_res/{sample}_recal_reads.bam"
    params: hg38 = "/mnt/disk1/PROJECTS/IMMUNE/DATA/hg38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta"
    shell: """
           gatk BaseRecalibrator -R {params.hg38} -I {input.dedup_reads_sorted} --known-sites {input.bqsr_snps} --known-sites {input.bqsr_indels} -O {output.recal_table}
           gatk ApplyBQSR -R {params.hg38} -I {input.dedup_reads_sorted} -bqsr {output.recal_table} -O {output.recal_reads}
    """

rule call_2nd_round:  # final 2nd round of SNP/indel calling
    input: recal_reads = "res_snake/aln_res/{sample}_recal_reads.bam",
    output: raw_variants_recal = "res_snake/vcf/{sample}_raw_variants_recal.vcf",
    params: hg38 = "/mnt/disk1/PROJECTS/IMMUNE/DATA/hg38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta"
    shell: """
           gatk HaplotypeCaller -R {params.hg38} -I {input.recal_reads} -O {output.raw_variants_recal}
    """

rule extract_SNPs_indels_2nd:  # extract only SNP/indels
    input: raw_variants_recal = "res_snake/vcf/{sample}_raw_variants_recal.vcf",
    output: raw_snps_recal = "res_snake/vcf/{sample}_raw_snps_recal.vcf",
            raw_indels_recal = "res_snake/vcf/{sample}_raw_indels_recal.vcf"
    params: hg38 = "/mnt/disk1/PROJECTS/IMMUNE/DATA/hg38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta"
    shell: """
           gatk SelectVariants -R {params.hg38} -V {input.raw_variants_recal} -select-type SNP -O {output.raw_snps_recal}
           gatk SelectVariants -R {params.hg38} -V {input.raw_variants_recal} -select-type INDEL -O {output.raw_indels_recal}
    """

rule snps_indels_filtration_2nd:  # filtering again
    input: raw_snps_recal = "res_snake/vcf/{sample}_raw_snps_recal.vcf",
           raw_indels_recal = "res_snake/vcf/{sample}_raw_indels_recal.vcf"
    output: final_snps = "res_snake/vcf/{sample}_final_snps.vcf",
            final_indels = "res_snake/vcf/{sample}_final_indels.vcf"
    params: hg38 = "/mnt/disk1/PROJECTS/IMMUNE/DATA/hg38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta"
    shell: """
           gatk VariantFiltration -R {params.hg38} -V {input.raw_snps_recal} -O {output.final_snps} --filter-name "QD_filter" -filter "QD < 2.0" --filter-name "FS_filter" -filter "FS > 60.0" --filter-name "MQ_filter" -filter "MQ < 40.0" --filter-name "SOR_filter" -filter "SOR > 4.0"
           gatk VariantFiltration -R {params.hg38} -V {input.raw_indels_recal} -O {output.final_indels} --filter-name "QD_filter" -filter "QD < 2.0" --filter-name "FS_filter" -filter "FS > 200.0" --filter-name "SOR_filter" -filter "SOR > 10.0"
    """

rule snp_annotate:  # annotate with GATK and dbSNP v151 (2019 version, best consistency of rs_ids and matches rs_ids with 1000 GENOMES db)
    input: final_snps = "res_snake/vcf/{sample}_final_snps.vcf",
    output: final_snps_anno = "res_snake/vcf/{sample}_final_snps_ann.vcf",
    params: hg38 = "/mnt/disk1/PROJECTS/IMMUNE/DATA/hg38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta",
            dbsnp = "/mnt/disk1/DATABASES/dbSNP_ncbi/All_20180418.vcf.gz"
    shell: """
           gatk VariantAnnotator -V {input.final_snps} -R {params.hg38} --dbsnp {params.dbsnp} -O {output.final_snps_anno}
    """

rule snp_annotate_with_af:  #annotate with 1000 GENOMES allele frequency db
    input: final_snps_anno = "res_snake/vcf/{sample}_final_snps_ann.vcf",
    output: final_snps_ready = "res_snake/vcf/{sample}_final_snps_ann_ready.txt",
            final_snps_with_af = "res_snake/vcf/{sample}_final_snps_ann_af.txt"
    params: snps_list = "/mnt/disk1/PROJECTS/EXOMS/RESULTS/panels/rs_txt_lists/crohn_all_nodup.txt",
            genomes1k = "/mnt/disk1/DATABASES/1k_genomes_test/1000GENOMES-phase_3_mod.vcf"
    shell: """
           grep -wf {params.snps_list} {input.final_snps_anno} | cut -f1,2,3,4,5,10 | sed 's/:/\t/g' | cut -f1,2,3,4,5,6,7,8 | awk 'BEGIN{{FS=OFS="\\t"}} {{gsub(/,/, "/", $7)}}1' > {output.final_snps_ready}
           awk '{{if (NR==FNR) {{a[$3]=$8; next}} if ($3 in a) {{print $0, a[$3]}}}}' {params.genomes1k} {output.final_snps_ready} > {output.final_snps_with_af}
    """

rule snp_annotate_unmatched:  # same tactics, but for all SNP from snp list to get full info about all SNPs, even absent in sample.
    input: final_snps_anno = "res_snake/vcf/{sample}_final_snps_ann.vcf",
           dedup_reads_depth = "res_snake/aln_res/{sample}_dedup_reads.depth"
    output: final_snps_unmatched_depth = "res_snake/vcf/{sample}_unmatched_snps_depth.txt",
            final_snps_unmatched_anno = "res_snake/vcf/{sample}_unmatched_res.txt"
    params: snps_list = "/mnt/disk1/PROJECTS/EXOMS/RESULTS/panels/rs_txt_lists/crohn_all_nodup.txt",
            genomes1k = "/mnt/disk1/DATABASES/1k_genomes_test/1000GENOMES-phase_3_mod.vcf",
            snps_coords="/mnt/disk1/PROJECTS/EXOMS/RESULTS/panels/beds/snp_loc_final_0based.bed"
    shell: """
           grep -oFf {params.snps_list} {input.final_snps_anno} |grep -vFf - {params.snps_list} | grep -wf - {params.snps_coords} | cut -f1,2 | grep -wf - {input.dedup_reads_depth} > {output.final_snps_unmatched_depth}
           awk 'NR==FNR{{a[$1,$2]=$3;next}} (($1,$2) in a) {{$6="0/0"; $7=a[$1,$2]; print }}' {output.final_snps_unmatched_depth} {params.genomes1k} | awk '{{print $1,$2,$3,$4,$5,$6,$7,$7,$8}}' | awk '{{$7=$7"/0"; print}}' | sed 's/ /\\t/g' >  {output.final_snps_unmatched_anno}
    """

