ethnicities = ["black", "hispanic", "white"]

rule all:
	input: "best_associations.tsv", "promoters_best_associations.tsv", "summary_table.tsv"

import glob
import gzip

BCF_DIR="source_genotypes"

ALLPWMS="ASCL1_HUMAN.H11MO.0.A,ASCL2_HUMAN.H11MO.0.D,ATF3_HUMAN.H11MO.0.A,ATF4_HUMAN.H11MO.0.A,BACH1_HUMAN.H11MO.0.A,BACH2_HUMAN.H11MO.0.A,BATF3_HUMAN.H11MO.0.B,BATF_HUMAN.H11MO.1.A,BC11A_HUMAN.H11MO.0.A,BCL6_HUMAN.H11MO.0.A,CEBPA_HUMAN.H11MO.0.A,CEBPB_HUMAN.H11MO.0.A,CEBPD_HUMAN.H11MO.0.C,CEBPE_HUMAN.H11MO.0.A,CEBPG_HUMAN.H11MO.0.B,CTCF_HUMAN.H11MO.0.A,CTCFL_HUMAN.H11MO.0.A,DBP_HUMAN.H11MO.0.B,EHF_HUMAN.H11MO.0.B,ELF2_HUMAN.H11MO.0.C,ELF5_HUMAN.H11MO.0.A,EOMES_HUMAN.H11MO.0.D,ERG_HUMAN.H11MO.0.A,ETS1_HUMAN.H11MO.0.A,ETS2_HUMAN.H11MO.0.B,ETV4_HUMAN.H11MO.0.B,ETV6_HUMAN.H11MO.0.D,FIGLA_HUMAN.H11MO.0.D,FLI1_HUMAN.H11MO.1.A,FOSB_HUMAN.H11MO.0.A,FOS_HUMAN.H11MO.0.A,FOSL1_HUMAN.H11MO.0.A,FOSL2_HUMAN.H11MO.0.A,GABPA_HUMAN.H11MO.0.A,GATA1_HUMAN.H11MO.1.A,GATA2_HUMAN.H11MO.1.A,GATA3_HUMAN.H11MO.0.A,GATA4_HUMAN.H11MO.0.A,GATA5_HUMAN.H11MO.0.D,GATA6_HUMAN.H11MO.0.A,GFI1_HUMAN.H11MO.0.C,HEN1_HUMAN.H11MO.0.C,HLF_HUMAN.H11MO.0.C,HTF4_HUMAN.H11MO.0.A,HXC13_HUMAN.H11MO.0.D,ID4_HUMAN.H11MO.0.D,IKZF1_HUMAN.H11MO.0.C,IRF1_HUMAN.H11MO.0.A,IRF2_HUMAN.H11MO.0.A,IRF3_HUMAN.H11MO.0.B,IRF4_HUMAN.H11MO.0.A,IRF5_HUMAN.H11MO.0.D,IRF7_HUMAN.H11MO.0.C,IRF8_HUMAN.H11MO.0.B,IRF9_HUMAN.H11MO.0.C,ITF2_HUMAN.H11MO.0.C,JDP2_HUMAN.H11MO.0.D,JUNB_HUMAN.H11MO.0.A,JUND_HUMAN.H11MO.0.A,JUN_HUMAN.H11MO.0.A,KLF1_HUMAN.H11MO.0.A,KLF4_HUMAN.H11MO.0.A,LYL1_HUMAN.H11MO.0.A,MAFB_HUMAN.H11MO.0.B,MAFF_HUMAN.H11MO.1.B,MAFG_HUMAN.H11MO.1.A,MAFK_HUMAN.H11MO.1.A,MEF2C_HUMAN.H11MO.0.A,MESP1_HUMAN.H11MO.0.D,MITF_HUMAN.H11MO.0.A,MYF6_HUMAN.H11MO.0.C,MYOG_HUMAN.H11MO.0.B,NDF1_HUMAN.H11MO.0.A,NF2L2_HUMAN.H11MO.0.A,NFE2_HUMAN.H11MO.0.A,NFIB_HUMAN.H11MO.0.D,NFIL3_HUMAN.H11MO.0.D,PAX5_HUMAN.H11MO.0.A,PRDM1_HUMAN.H11MO.0.A,RORA_HUMAN.H11MO.0.C,RUNX1_HUMAN.H11MO.0.A,RUNX2_HUMAN.H11MO.0.A,RUNX3_HUMAN.H11MO.0.A,SNAI1_HUMAN.H11MO.0.C,SPI1_HUMAN.H11MO.0.A,SPIB_HUMAN.H11MO.0.A,SPIC_HUMAN.H11MO.0.D,STAT1_HUMAN.H11MO.1.A,STAT2_HUMAN.H11MO.0.A,STAT3_HUMAN.H11MO.0.A,STAT6_HUMAN.H11MO.0.B,TAL1_HUMAN.H11MO.1.A,TBX21_HUMAN.H11MO.0.A,TFAP4_HUMAN.H11MO.0.A,TFE2_HUMAN.H11MO.0.A,TWST1_HUMAN.H11MO.1.A,ZEB1_HUMAN.H11MO.0.A"

BED_FILES="filtered_bed/Bcell-13.bed,filtered_bed/CD4-9.bed,filtered_bed/CD8-10.bed,filtered_bed/CLP-14.bed,filtered_bed/CMP-4.bed,filtered_bed/Erythro-15.bed,filtered_bed/GMP-5.bed,filtered_bed/HSC-1.bed,filtered_bed/LMPP-3.bed,filtered_bed/MEGA1.bed,filtered_bed/MEGA2.bed,filtered_bed/MEP-6.bed,filtered_bed/Mono-7.bed,filtered_bed/MPP-2.bed,filtered_bed/Nkcell-11.bed,filtered_bed/pDC.bed"

chromosomes = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"]


def scan_time(chromosome):
    d = {
        "1": 100,
        "2": 110,
        "3": 71,
        "4": 71,
        "5": 71,
        "6": 72,
        "7": 65,
        "8": 65,
        "9": 60,
        "10": 60,
        "11": 60,
        "12": 60,
        "13": 52,
        "14": 50,
        "15": 47,
        "16": 45,
        "17": 45,
        "18": 30,
        "19": 24,
        "20": 24,
        "21": 13,
        "22": 15,
        "X": 60
    }
    return int(d[chromosome] * 0.5)


localrules: all_samples, african_samples, european_samples, hispanic_samples, reference_genome, ensembl_regulatory_build, ensembl_genes, hocomoco_thresholds, hocomoco_pwm_list, best_associations_per_ethnicity,best_associations, extract_promoters_and_enhancers, conditional_analysis, all


rule all_samples:
	input: "pheno.tab"
	output: "samples"
	shell: "cut -f 1 {input} | grep NWD > {output}"

rule african_samples:
	input: "pheno.tab"
	output: "samples_black"
	shell: "grep Black {input} | cut -f 1 > {output}"

rule european_samples:
	input: "pheno.tab"
	output: "samples_white"
	shell: "grep White {input} | cut -f 1 > {output}"

rule hispanic_samples:
	input: "pheno.tab"
	output: "samples_hispanic"
	shell: "grep Hispanic {input} | cut -f 1 > {output}"

rule reference_genome:
	output: "hg38.fa"
	shell: """
		wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
		gunzip hg38.fa.gz
		rm hg38.fa.gz
		samtools faidx hg38.fa
		"""

rule ensembl_regulatory_build:
	output: promoters = "ensembl_promoters.tsv", enhancers = "ensembl_enhancers.tsv"
	shell: """
		wget ftp://ftp.ensembl.org/pub/release-99/regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.gff.gz
		gunzip homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.gff.gz
		cat homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.gff | awk '{{ if($3 == "promoter") {{print $1,$4,$5}}}}' 'OFS=\t' > {output.promoters}
		cat homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.gff | awk '{{ if($3 == "enhancer") {{print $1,$4,$5}}}}' 'OFS=\t' > {output.enhancer}		
		"""

rule ensembl_genes:
	output: "ensembl_genes.tsv"
	shell: """
		wget ftp://ftp.ensembl.org/pub/release-85/gff3/homo_sapiens/Homo_sapiens.GRCh38.85.gff3.gz
		zcat Homo_sapiens.GRCh38.85.gff3.gz | awk '$3 ~ /gene/ {{ print $1,$4,$5,$3,$9 }}' 'OFS=\t' | grep -v pseudogene > {output}
		"""

rule hocomoco_thresholds:
	output: "thresholds/done"
	shell: """
		wget https://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_thresholds_HUMAN_mono.tar.gz
		tar xzf HOCOMOCOv11_full_thresholds_HUMAN_mono.tar.gz
		rm HOCOMOCOv11_full_thresholds_HUMAN_mono.tar.gz
		touch thresholds/done
		"""

# Genotypes in linkage equilibrium for all populations
rule le_genotypes:
    input: bcf = ancient(BCF_DIR+"/freeze.8.chr{chromosome}.pass_only.phased.bcf"), coordinates = "linkage_equilibrium_coordinates_norsid.tab"
    output: vcf = "genotype_le/chr{chromosome}.vcf.gz"
    threads: 1
    resources: mem=1000, runtime=lambda w: (12*60 if w.chromosome in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10'] else 3*60), threads= 1
    shell:
        """
        mkdir -p genotype_le
        vcftools --bcf {input.bcf} --positions {input.coordinates} --recode --stdout | grep -v '##' | sed 's/^chr//' | bgzip > {output.vcf}
        """

# Genotypes in linkage equilibrium for all populations
rule le_genotypes_one_file:
    input: expand("genotype_le/chr{chromosome}.vcf.gz", chromosome = chromosomes)
    output: "genotype_le/all.vcf.gz"
    threads: 1
    resources: mem=1000, runtime=60, threads= 1
    shell: "cd genotype_le; {{ zcat chr1.vcf.gz; zcat chr2.vcf.gz | grep -v '#'; zcat chr3.vcf.gz | grep -v '#'; zcat chr4.vcf.gz | grep -v '#'; zcat chr5.vcf.gz | grep -v '#'; zcat chr6.vcf.gz | grep -v '#'; zcat chr7.vcf.gz | grep -v '#'; zcat chr8.vcf.gz | grep -v '#'; zcat chr9.vcf.gz | grep -v '#'; zcat chr10.vcf.gz | grep -v '#'; zcat chr11.vcf.gz | grep -v '#'; zcat chr12.vcf.gz | grep -v '#'; zcat chr13.vcf.gz | grep -v '#'; zcat chr14.vcf.gz | grep -v '#'; zcat chr15.vcf.gz | grep -v '#'; zcat chr16.vcf.gz | grep -v '#'; zcat chr17.vcf.gz | grep -v '#'; zcat chr18.vcf.gz | grep -v '#'; zcat chr19.vcf.gz | grep -v '#'; zcat chr20.vcf.gz | grep -v '#'; zcat chr21.vcf.gz | grep -v '#'; zcat chr22.vcf.gz | grep -v '#'; zcat chrX.vcf.gz | grep -v '#'; }} | bgzip > all.vcf.gz"

# Genotypes in linkage equilibrium for each population
rule le_genotypes_pop:
    input: vcf = ancient("genotype_le/all.vcf.gz"), samples = "samples_{ethnicity}"
    output: vcf = "{ethnicity}/le.vcf.gz", tbi = "{ethnicity}/le.vcf.gz.tbi"
    threads: 1
    resources: mem=1000, runtime=3*60, threads= 1
    shell:
        """
        vcftools --gzvcf genotype_le/all.vcf.gz --recode --recode-INFO-all --keep {input.samples} --stdout | gzip -c > {output.vcf} &&
        tabix -f -p vcf {output.vcf}
        """

# Kinship matrix for this population ####################
rule pop_kinship_matrix:
   input: vcf = ancient("{ethnicity}/le.vcf.gz"), tbi = ancient("{ethnicity}/le.vcf.gz.tbi"), samples = ancient("samples_{ethnicity}")
   output: "{ethnicity}/kinship.kinf"
   threads: 4
   resources: mem=lambda w: (30000 if w.pop not in ['white'] else 60000), runtime=lambda w: (3*60 if w.pop not in ['white'] else 9*60), threads= 4
   shell:
        """
		(cat "#FAM_ID,IND_ID,FAT_ID,MOT_ID"; (cat {input.samples} | awk '{{print "",$1,"",""}}' 'OFS=\t')) | tr ',' '\t' > {input.samples}.ped
        epacts make-kin --vcf {input.vcf} --ped {input.samples}.ped --min-maf 0.01 --min-callrate 0.95 --out {output}
        make -f {output}.Makefile -j 4
        """

rule hocomoco_pwm_list:
	output: "HOCOMOCOv11_full_pwms_HUMAN_mono.txt"
	shell: "wget https://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_pwms_HUMAN_mono.txt"

rule scan:
	input: all_samples=ancient("samples"), bcf=ancient("%s/freeze.8.chr{chromosome}.pass_only.phased.bcf" % BCF_DIR), reference_genome=ancient("hg38.fa"), thresholds=ancient("thresholds/done"), pwms=ancient("HOCOMOCOv11_full_pwms_HUMAN_mono.txt")
	output: "chr{chromosome}.vcf.gz"
	resources: mem=1000, threads=2, runtime=lambda wildcards: (60*scan_time(wildcards.chromosome))
	threads: 2
	shell: """
		find-tfbs --threads 2 -c chr{wildcards.chromosome} -o {output} --min_maf 5 --samples {input.all_samples} -r {input.reference_genome} -i {input.bcf} -n %s -b %s -p {input.pwms} --pwm_threshold_directory thresholds --pwm_threshold 0.0001
		""" % (ALLPWMS, BED_FILES)

rule ethnicity_vcf:
	input: vcf="chr{chromosome}.vcf.gz", samples="samples_{ethnicity}"
	output: vcf = "{ethnicity}_chr{chromosome}.vcf.gz", tbi = "{ethnicity}_chr{chromosome}.vcf.gz.tbi", positions = "{ethnicity}_chr{chromosome}.positions"
	resources: mem=2000, threads=1, runtime=lambda w: (3*60 if (w.ethnicity == 'white' and w.chromosome in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10']) else 60)
	threads: 1
	shell: """
		vcftools --gzvcf {input.vcf} --keep {input.samples} --recode --stdout | grep -v '##' | awk '{{if (NR!=1) {{$2 = NR-1; $7="PASS";}} print;}}' OFS='\t' | bgzip > {output.vcf}
		tabix -f -p vcf {output.vcf}
		zcat {output.vcf} | cut -f 1,2 > {output.positions}
		"""

rule ethnicity_vcf_rbc:
	input: vcf="{ethnicity}_chr{chromosome}.vcf.gz", samples="samples_{ethnicity}"
	output: vcf = "{ethnicity}_chr{chromosome}_rbc.vcf.gz", tbi = "{ethnicity}_chr{chromosome}_rbc.vcf.gz.tbi", positions = "{ethnicity}_chr{chromosome}_rbc.positions"
	resources: mem=1000, threads=1, runtime=60
	threads: 1
	shell: """
		zcat {input.vcf} | grep -E '(CHROM|Erythro-15|CMP-4|MEP-6|MPP-2|HSC-1)' | bgzip > {output.vcf}
		tabix -f -p vcf {output.vcf}
		zcat {output.vcf} | cut -f 1,2 > {output.positions}
		"""

rule ethnicity_vcf_plt:
	input: vcf="{ethnicity}_chr{chromosome}.vcf.gz", samples="samples_{ethnicity}"
	output: vcf = "{ethnicity}_chr{chromosome}_plt.vcf.gz", tbi = "{ethnicity}_chr{chromosome}_plt.vcf.gz.tbi", positions = "{ethnicity}_chr{chromosome}_plt.positions"
	resources: mem=1000, threads=1, runtime=60
	threads: 1
	shell: """
		zcat {input.vcf} | grep -E '(CHROM|MEGA1|MEGA2|MEP-6|CMP-4|MPP-2|HSC-1)' | bgzip > {output.vcf}
		tabix -f -p vcf {output.vcf}
		zcat {output.vcf} | cut -f 1,2 > {output.positions}
		"""

rule ethnicity_vcf_ben:
	input: vcf="{ethnicity}_chr{chromosome}.vcf.gz", samples="samples_{ethnicity}"
	output: vcf = "{ethnicity}_chr{chromosome}_ben.vcf.gz", tbi = "{ethnicity}_chr{chromosome}_ben.vcf.gz.tbi", positions = "{ethnicity}_chr{chromosome}_ben.positions"
	resources: mem=1000, threads=1, runtime=60
	threads: 1
	shell: """
		zcat {input.vcf} | grep -E '(CHROM|GMP-5|CMP-4|LMPP-3|MPP-2|HSC-1)' | bgzip > {output.vcf}
		tabix -f -p vcf {output.vcf}
		zcat {output.vcf} | cut -f 1,2 > {output.positions}
		"""

rule ethnicity_vcf_wbc:
	input: vcf="{ethnicity}_chr{chromosome}.vcf.gz", samples="samples_{ethnicity}"
	output: vcf = "{ethnicity}_chr{chromosome}_wbc.vcf.gz", tbi = "{ethnicity}_chr{chromosome}_wbc.vcf.gz.tbi", positions = "{ethnicity}_chr{chromosome}_wbc.positions"
	resources: mem=1000, threads=1, runtime=60
	threads: 1
	shell: """
		zcat {input.vcf} | grep -E '(CHROM|Bcell-13|CD4-9|Nkcell-11|Mono-7|pDC|CD8-10|GMP-5|CLP-14|MPP-2|LMPP-3|HSC-1|CMP-4|mDC)' | bgzip > {output.vcf}
		tabix -f -p vcf {output.vcf}
		zcat {output.vcf} | cut -f 1,2 > {output.positions}
		"""

rule ethnicity_vcf_lymph:
	input: vcf="{ethnicity}_chr{chromosome}.vcf.gz", samples="samples_{ethnicity}"
	output: vcf = "{ethnicity}_chr{chromosome}_lymph.vcf.gz", tbi = "{ethnicity}_chr{chromosome}_lymph.vcf.gz.tbi", positions = "{ethnicity}_chr{chromosome}_lymph.positions"
	resources: mem=1000, threads=1, runtime=60
	threads: 1
	shell: """
		zcat {input.vcf} | grep -E '(CHROM|CD4-9|CD8-10|Nkcell-11|CLP-14|LMPP-3|MPP-2|HSC-1)' | bgzip > {output.vcf}
		tabix -f -p vcf {output.vcf}
		zcat {output.vcf} | cut -f 1,2 > {output.positions}
		"""

rule all_vcf:
	input: vcf=expand("{{ethnicity}}_chr{chromosome}.vcf.gz", chromosome=chromosomes), 
	       vcf_rbc=expand("{{ethnicity}}_chr{chromosome}_rbc.vcf.gz", chromosome=chromosomes),
		   vcf_plt=expand("{{ethnicity}}_chr{chromosome}_plt.vcf.gz", chromosome=chromosomes),
		   vcf_ben=expand("{{ethnicity}}_chr{chromosome}_ben.vcf.gz", chromosome=chromosomes),
		   vcf_wbc=expand("{{ethnicity}}_chr{chromosome}_wbc.vcf.gz", chromosome=chromosomes),
		   vcf_lymph=expand("{{ethnicity}}_chr{chromosome}_lymph.vcf.gz", chromosome=chromosomes),
	output: "{ethnicity}_all_vcf_ready"


def vcf_suffix(trait):
	if trait in ["MCH", "MCV", "MCHC", "HCT", "HGB", "RBC", "RDW"]:
		return "_rbc"
	elif trait in ["MPV", "PLT"]:
		return "_plt"
	elif trait in ["BASO", "EOSIN", "NEUTRO"]:
		return "_ben"
	elif trait in ["WBC"]:
		return "_wbc"
	elif trait in ["LYMPH"]:
		return "_lymph"
	else:
		return ""

rule ethnicity_epacts:
	input: vcf=expand("{{ethnicity}}_chr{chromosome}.vcf.gz", chromosome=chromosomes), 
	       vcf_rbc=expand("{{ethnicity}}_chr{chromosome}_rbc.vcf.gz", chromosome=chromosomes),
		   vcf_plt=expand("{{ethnicity}}_chr{chromosome}_plt.vcf.gz", chromosome=chromosomes),
		   vcf_ben=expand("{{ethnicity}}_chr{chromosome}_ben.vcf.gz", chromosome=chromosomes),
		   vcf_wbc=expand("{{ethnicity}}_chr{chromosome}_wbc.vcf.gz", chromosome=chromosomes),
		   vcf_lymph=expand("{{ethnicity}}_chr{chromosome}_lymph.vcf.gz", chromosome=chromosomes),
		   ped = "residuals/{ethnicity}/{trait}.ped",
		   kin = "{ethnicity}/kinship.kinf"
	output: "epacts/{ethnicity}/EPACTS_{trait}.epacts.gz"
	params: prefix = "epacts/{ethnicity}/EPACTS_{trait}", window = lambda wildcards: 10000, threads = lambda wildcards: (1 if wildcards.pop == "white" else 6), input_vcf_suffix = lambda wildcards: vcf_suffix(wildcards.trait)
	resources: mem=lambda wildcards: (15000 if wildcards.pop == "white" else 5000), threads=lambda wildcards: (6 if wildcards.pop == "white" else 6), runtime=lambda wildcards: (40*60 if wildcards.pop == "white" else 15*60)
	threads: lambda wildcards: (6 if wildcards.pop == "white" else 6)
	shell: """
		mkdir -p epacts/{wildcards.ethnicity}/
		find epacts/{wildcards.ethnicity}/ -name "EPACTS_{wildcards.trait}*" -delete
		epacts single --sepchr --kinf {input.kin} --min-maf 0 --min-mac 0 --min-callrate 0 --vcf {wildcards.ethnicity}_chr1{params.input_vcf_suffix}.vcf.gz --ped {input.ped} --pheno PHENO --field DS --unit {params.window} --test q.emmax --out {params.prefix} --cov PC1 --cov PC2 --cov PC3 --cov PC4 --cov PC5 --cov PC6 --cov PC7 --cov PC8 --cov PC9 --cov PC10
		python optimize_epacts_makefile.py {params.prefix}.Makefile {wildcards.ethnicity}_chr1{params.input_vcf_suffix}.positions {params.window}
		make -f {params.prefix}.Makefile.optimized -j {resources.threads}
		#epacts-plot --in {output} --ignore-maf --out epacts/{wildcards.ethnicity}/EPACTS_{wildcards.trait}
		#find epacts/{wildcards.ethnicity}/ -name "*{wildcards.trait}*.epacts" -delete
		#rm {params.prefix}.Makefile {params.prefix}.Makefile.optimized {params.prefix}.cov {params.prefix}.phe {params.prefix}.ind #{params.prefix}.epacts.conf
		touch {output}
		"""

# Find the associations with lowest p-values, in each ethnicity
rule best_associations_per_ethnicity:
	input: expand("epacts/{{ethnicity}}/EPACTS_{trait}.epacts.gz", trait=["HGB", "HCT", "MPV", "PLT", "WBC", "BASO", "EOSIN", "LYMPH", "MONO", "NEUTRO", "MCH", "MCHC", "MCV", "RDW", "RBC"])
	output: "{ethnicity}/best_associations.tsv"
	shell: "python extract_best_associations.py phenotypes_per_cell_type.tab transcription_factors_per_cell_type.tab {output} {wildcards.ethnicity}"


# Find the associations with lowest p-values, in all ethnicities
rule best_associations:
	input: black = "black/best_associations.tsv", hispanic = "hispanic/best_associations.tsv", white = "white/best_associations.tsv"
	output: "best_associations.tsv"
	shell: """ (head -1 {input.white}; grep -v "POP" {input.white}; grep -v "POP" {input.black}; grep -v "POP" {input.hispanic};) | sed -e /^$/d | sort -r -k 6,6g > {output}"""

rule extract_promoters_and_enhancers:
	input: associations = "best_associations.tsv", ensembl_promoters = "ensembl_promoters.tsv", ensembl_enhancers = "ensembl_enhancers.tsv", ensembl_genes = "ensembl_genes.tsv"
	output: "promoters_enhancers_best_associations.tsv"
	shell: "python extract_associations_that_overlap_promoters_enhancers.py {input.associations} {input.ensembl_promoters} {input.ensembl_enhancers} {input.ensembl_genes} {output}"


rule conditional_analysis:
	input: "promoters_enhancers_best_associations.tsv"
	output: "conditional_analysis.done", "summary_table.tsv"
        resources: mem=3500, runtime=10*60, threads= 10
        threads: 10
	shell: """
		module load nixpkgs/16.09  gcc/7.3.0 r/3.4.4
		module load python/3.8.2
		module load bcftools

		python3 conditional_analysis_promoters.py 10 %s
		touch conditional_analysis.done""" % BCF_DIR
