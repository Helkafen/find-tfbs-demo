import os
import os.path
import sys
import re
import gzip
from multiprocessing import Pool

def add_covars(ped_file, additional_covar_file, modified_ped_file):
    #print("# %s %s -> %s" % (ped_file, additional_covar_file, modified_ped_file))
    covars = {}
    covar_names = []
    f = open(additional_covar_file, "rt")
    individuals = f.readline().strip().split("\t")[4:]
    for individual in individuals:
        covars[individual] = []
    covar_names = []
    for line in f.readlines():
        if len(line) < 3:
            continue
        values = [x.strip() for x in line.split("\t")][4:]
        for (individual,value) in zip(individuals,values):
            covars[individual].append(value)
        covar_names.append("%s_%s_%s_%s" % (line.split("\t")[0],line.split("\t")[1],line.split("\t")[2],line.split("\t")[3]))
    f.close()

    f = open(ped_file, "rt")
    o = open(modified_ped_file, "wt")
    o.write(f.readline().strip()+ "\t" + "\t".join(covar_names))
    for l in f:
        if len(l) < 3:
            continue
        individual = l.split("\t")[1]
        if individual in covars:
            l = "\n" + l.strip() + "\t" + "\t".join(covars[individual])
            o.write(l)
        #else:
        #    print("Not in BCF: %s" % individual)
    o.close()
    f.close()

    #print("# Done %s %s -> %s" % (ped_file, additional_covar_file, modified_ped_file))
    return covar_names

class ConditionalTest:
    def __init__(self, pop, trait, cell, region, pwm, originalpvalue, originalbeta):
        self.gene = None
        self.cond_variants = None
        self.pop = pop
        self.trait = trait
        self.cell = cell
        self.region = region
        self.pwm = pwm
        self.chromosome = self.region.split(":")[0]
        self.dir = "conditional_analysis/%s" % self.name()
        if pop == "European":
            self.pop_lower = "white"
        elif pop == "African":
            self.pop_lower = "black"
        elif pop == "Hispanic":
            self.pop_lower = "hispanic"
        else:
            raise
        self.originalpvalue = originalpvalue
        self.originalbeta = originalbeta
    def set_cond_variants(self, cond_variants):
        self.cond_variants = list(sorted(set(cond_variants),key=lambda x:(x.chromosome, x.position)))
    def name(self):
        return "test_%s_%s_%s_%s_%s" % (self.gene, self.pop, self.trait, self.pwm, self.region.replace(":","_").replace("-","to"))
    def __repr__(self):
        return "Test"
    def __eq__(self, other):
        return self.name() == other.name() #(self.gene == other.gene) and (self.pop == other.pop) and (self.trait == other.trait) and (self.cell == other.cell) and self

class Variant:
    def __init__(self, chromosome, position): # , ref, alt):
        self.chromosome = chromosome
        self.position = position
        #self.ref = ref
        #self.alt = alt
    def __repr__(self):
        return "%s_%s_%s_%s" % (self.chromosome, self.position) #, self.ref, self.alt)
    def __hash__(self):
        return hash(self.position)
    def __eq__(self, other):
        return ((self.chromosome, self.position) == (other.chromosome, other.position))

def parse_variant(s):
    c,p,r,a = s.replace(":"," ").replace("_"," ").split(" ")
    return Variant(c,int(p),r,a)

def execute(command, verbose = True):
    if verbose:
        print(command)
    for l in os.popen(command).readlines():
        if verbose:
            print(l)

cells = ['Bcell-13', 'CD4-9', 'CD8-10', 'CLP-14', 'CMP-4', 'Erythro-15', 'GMP-5', 'HSC-1', 'LMPP-3', 'MEGA1', 'MEGA2', 'MEP-6', 'MPP-2', 'Mono-7', 'Nkcell-11', 'mDC', 'pDC']

def close_enough(x,y, tolerance = 0):
    (chromosome, position) = x
    (c,s,e) = y
    if not (c == chromosome):
        return False
    return min([abs(x-position) for x in [s,e]]) <= tolerance


def overlap(x,y):
    (c,s,e) = x
    (chromosome, start, end) = y
    if not (c == chromosome):
        return False
    if start <= e and start >= s:
        return True
    if end >= s and end <= e:
        return True
    return False

def parse_bed_line(l):
    c,s,e,q = l.strip().split("\t")
    return (c, int(s), int(e))

def parse_coordinates(l):
    c,rest = l.strip().split(":")
    s,e = rest.split("-")
    return (c, int(s), int(e))

def extract_counts(scan_file):
    fields = os.popen("zcat %s | grep -v CHROM | cut -f 8,10-" % scan_file).read().strip().split("\t")
    distinct_counts = [int(x) for x in fields[0].split(";")[0].replace("COUNTS=","").split(",")]
    min_count = min(distinct_counts)
    max_count = max(distinct_counts)
    counts = [min_count + round(float(x.split(":")[1])/2 * (max_count-min_count)) for x in fields[1:]]
    counts_str = ""
    for i in range(min_count,max_count+1):
        if i in counts:
            counts_str=counts_str+"%s:%s/" % (i, counts.count(i))
    return counts_str.strip("/")

def process_test(test):
    result = Result(test, False, False, False, True, True, None, None)

    print("Processing %s" % test.name())
    vcf = "%s_%s.vcf.gz" % (test.pop_lower, test.chromosome)
    original_epacts_gz = "epacts/%s/EPACTS_%s.epacts.gz" % (test.pop_lower, test.trait)
    #original_line = os.popen("zcat %s | grep '%s'" % (original_epacts_gz, test.region.split(":")[1])).read().strip() # TODO slow

    # Paths of the files created during the conditional analysis
    execute("mkdir -p %s" % test.dir, verbose = False)
    additional_covar_file = "%s/covar.vcf" % test.dir
    log = "%s/epacts.log" % test.dir
    new_top5000 = "%s/conditional.epacts.top5000" % test.dir
    non_conditional_top5000 = "%s/non_conditional.epacts.top5000" % test.dir
    summary_file = "%s/summary" % test.dir
    non_conditional_makefile = "%s/non_conditional.Makefile" % test.dir
    makefile = "%s/conditional.Makefile" % test.dir
    ped_file = "residuals/%s/%s.ped" % (test.pop_lower, test.trait)
    modified_ped_file = "%s/conditional.ped" % test.dir
    region_bed = "%s/region.bed" % test.dir
    sample_file = "samples_%s" % test.pop_lower
    filtered_sample_file = "%s/samples" % test.dir
    scan_file = "%s/scan.vcf.gz" % test.dir
    # Removed --min-maf 5 because this limit is applied to the whole population. Here we only look at one ethnicity
    scan_command =                 "find-tfbs --threads 1 -c %s -o %s.temp --min_maf 0 --samples %s -r hg38.fa -i %s/freeze.8.%s.pass_only.phased.bcf -n %s -b %s -p HOCOMOCOv11_full_pwms_HUMAN_mono.txt --pwm_threshold_directory thresholds --pwm_threshold 0.0001 --verbose > %s/scan.log" % (test.chromosome, scan_file, filtered_sample_file, BCF_DIR, test.chromosome, test.pwm, region_bed, test.dir)

    # Create bed files
    if not os.path.isfile(region_bed):
        print("Creating %s" % region_bed)
        (chromosome, rest) = test.region.split(":")
        (start, end) = rest.split("-")
        o = open(region_bed, "wt")
        o.write("%s\t%s\t%s\t" % (chromosome, start, end))
        o.close()

    bcf_sample_file = "bcf_header"
    if not os.path.isfile(bcf_sample_file):
        print("Create %s" % bcf_sample_file)
        os.popen("bcftools view --h BCF_DIR/freeze.8.chr16.pass_only.phased.bcf | grep -v '##' | cut -f 10 | tr '\t' '\n' > %s" % (BCF_DIR, bcf_sample_file))
    
    if not os.path.isfile(filtered_sample_file):
        print("Create %s" % filtered_sample_file)
        bcf_samples = set(open(bcf_sample_file, "rt").read().strip().split("\t"))
        ethnicity_samples = open(sample_file, "rt").read().strip().split("\n")
        f = open(filtered_sample_file, "wt")
        f.write("\n".join([x for x in ethnicity_samples if x in bcf_samples]))
        f.close()
    
    samples = open(sample_file, "rt").read().strip().split("\t")

    can_do_conditional_analysis = False
    print("Reached 1")
    # Create covariable phenotype file, with additional covariates (the conditional variants)
    if not test.cond_variants:
        result.conditional_didnt_run = "No known variants for conditional analysis"
        print(result.conditional_didnt_run)
        return result
    print("Reached 2")
    if not os.path.isfile(additional_covar_file+".OK"):
        for (i, cond_variant) in enumerate(test.cond_variants):
            if i == 0:
                prefilter = ""
                end = ">"
            else:
                prefilter = "| grep -v CHROM"   
                end = ">>"
            execute("bcftools view -f '%%LINE\\n' -r %s:%s -S %s %s/freeze.8.%s.pass_only.phased.bcf | grep -v '##' %s | awk '{if ($1 == \"#CHROM\" || (length($4) == 1 && length($5) == 1)) {print;}}' | cut -f 1,2,4,5,10- | sed 's/0|0/0/g' | sed 's/0|1/1/g' | sed 's/1|0/1/g' | sed 's/1|1/2/g' %s %s" % (cond_variant.chromosome, cond_variant.position, filtered_sample_file, BCF_DIR, cond_variant.chromosome, prefilter, end, additional_covar_file))
        execute("touch %s.OK" % additional_covar_file)
    print("Reached 3")
    if not os.path.isfile(modified_ped_file):
        number_of_lines_in_additional_covar_file = int(os.popen("grep -v CHROM %s | wc -l" % additional_covar_file).read().strip().split(" ")[0])
        if number_of_lines_in_additional_covar_file > 1:
            can_do_conditional_analysis = True
            print("Create %s" % modified_ped_file)
            covar_names = add_covars(ped_file, additional_covar_file, modified_ped_file)
        else:
            result.conditional_didnt_run = "No known variants in TOPMed for conditional analysis"
            print(result.conditional_didnt_run)
            return result
    else:
        can_do_conditional_analysis = True
        covar_names = add_covars(ped_file, additional_covar_file, modified_ped_file)
    print("Reached 4")
    number_of_lines_in_additional_covar_file = int(os.popen("grep -v CHROM %s | wc -l" % additional_covar_file).read().strip().split(" ")[0])
    if not (number_of_lines_in_additional_covar_file > 1):
        result.conditional_didnt_run = "No known variants in TOPMed for conditional analysis"
        print(result.conditional_didnt_run)
        return result
        
    # Extract the variants from the ATAC-seq peak
    print("Reached 5 %s" % scan_file)
    if not os.path.isfile(scan_file):
        print("Create %s" % scan_file)
        execute(scan_command)
        execute("zcat %s.temp | grep 'CHROM\|%s' | bgzip > %s" % (scan_file, test.region.split(":")[1], scan_file))
        execute("rm %s.temp" % scan_file)
        execute("tabix -f -p vcf %s" % scan_file)

    makefile_linear = "%s/conditional_linear.Makefile" % test.dir
    if not os.path.isfile(makefile_linear):
        print("Create %s" % makefile_linear)
        cov_str = " --cov PC1 --cov PC2 --cov PC3 --cov PC4 --cov PC5 --cov PC6 --cov PC7 --cov PC8 --cov PC9 --cov PC10 " + " ".join("--cov %s" % x for x in covar_names)
        execute("epacts single --vcf %s --region %s:1-100 --ped %s --field DS --pheno PHENO %s --min-maf 0 --min-mac 0 --min-callrate 0 --test q.linear --out %s &> %s" % (scan_file, test.chromosome.replace("chr",""), modified_ped_file, cov_str, makefile_linear.replace(".Makefile",""), log))
        execute("sed -i 's/epacts-plot/epacts-plot --ignore-maf/g' %s &>> %s" % (makefile_linear, log))

    if not os.path.isfile(makefile):
        print("Create %s" % makefile)
        cov_str = " --cov PC1 --cov PC2 --cov PC3 --cov PC4 --cov PC5 --cov PC6 --cov PC7 --cov PC8 --cov PC9 --cov PC10 " + " ".join("--cov %s" % x for x in covar_names)
        execute("epacts single --vcf %s --region %s:1-100 --ped %s --field DS --pheno PHENO %s --min-maf 0 --min-mac 0 --min-callrate 0 --test q.emmax --kinf %s/kinship.kinf --out %s &> %s" % (scan_file, test.chromosome.replace("chr",""), modified_ped_file, cov_str, test.pop_lower, makefile.replace(".Makefile",""), log))
        execute("sed -i 's/epacts-plot/epacts-plot --ignore-maf/g' %s &>> %s" % (makefile, log))

    conditional_linear_epacts_file = makefile_linear.replace(".Makefile", ".epacts.gz")
    if not os.path.isfile(conditional_linear_epacts_file):
        execute("make -f %s &>> %s" % (makefile_linear, log))

    linear_pvalue_file = makefile_linear.replace(".Makefile", ".pvalue")
    emmax_pvalue_file = makefile.replace(".Makefile", ".pvalue")

    if not os.path.isfile(linear_pvalue_file):
        lines = [l.strip() for l in gzip.open(conditional_linear_epacts_file, "rt").readlines() if l][1:]
        assert(len(lines) == 1)
        linear_pvalue = float(lines[0].split("\t")[-5])
        print(linear_pvalue)
        f = open(linear_pvalue_file, "wt")
        f.write("%s" % linear_pvalue)
        f.close()

    linear_pvalue = float(open(linear_pvalue_file, "rt").read())
    result.linear_pvalue = linear_pvalue
    emmax_pvalue = None
    if linear_pvalue <= 0.05:
        conditional_epacts_file = makefile.replace(".Makefile", ".epacts.gz")
        if not os.path.isfile(conditional_epacts_file):
            print("EMMAX needs to run %s" % test.pop_lower)
            execute("make -f %s &>> %s" % (makefile, log))
        
        if not os.path.isfile(emmax_pvalue_file):
            lines = [l.strip() for l in gzip.open(conditional_epacts_file, "rt").readlines() if l][1:]
            assert(len(lines) == 1)
            emmax_pvalue = float(lines[0].split("\t")[-4])
            f = open(emmax_pvalue_file, "wt")
            f.write("%s" % emmax_pvalue)
            f.close()

    if os.path.isfile(emmax_pvalue_file):
        emmax_pvalue = float(open(emmax_pvalue_file, "rt").read())
        result.emmax_pvalue = emmax_pvalue


    return result

    if not os.path.isfile(summary_file):
        f = open(summary_file, "wt")
        if test.cond_variants and (number_of_lines_in_additional_covar_file >= 1):
            f.write("Region: %s\n\n" % (test.region,))
            f.write("Original from %s:\n" % original_epacts_gz)
            f.write("%s\n\n" % original_line)
            f.write("After correcting for sentinel variants:\n")
            f.write(os.popen("grep -E '(%s)|CHROM' %s" % (test.region.split(":")[1], new_top5000)).read())
            f.write("Without correcting for sentinel variants:\n")
            f.write(os.popen("grep -E '(%s)|CHROM' %s" % (test.region.split(":")[1], non_conditional_top5000)).read())
            f.write("\n")

            f.write("Conditional variants at that are in TOPMed or not:\n")
            pf = open(modified_ped_file, "rt")
            header = pf.readline().strip().split("\t")
            for cond_variant in test.cond_variants:
                try:
                    index = header.index(str(cond_variant))
                except:
                    f.write("%s not in TOPMed\n" % cond_variant)
                    continue
                counts = [l.split("\t")[index].strip() for l in open(modified_ped_file, "rt").readlines()[1:] if l]
                f.write("%s\t0:%s 1:%s 2:%s\n" % (cond_variant, counts.count("0"), counts.count("1"), counts.count("2")))

            buffer = os.popen("grep '^score' %s/scan.log | cut -d' ' -f 8,10 | sort | uniq -c" % test.dir).read()
            beginnings_of_match = list(sorted(list(set([int(x.split(" ")[1]) for x in re.findall("\d+ \d+", buffer)]))))

            f.write("\n\nBeginning of matches and how many unique haplotypes have these matches:\n")
            f.write(buffer)

            f.write("\nRaw scanner data:\n")
            f.write(os.popen("zcat %s | cut -f 1-13" % scan_file).read())

            p_value = original_line.split("\t")[8]
            conditional_p_value = os.popen("zcat %s/conditional.epacts.gz | grep %s | cut -f 9" % (test.dir, test.region.split(":")[1])).read().strip()
            counts = extract_counts(scan_file)
            f.write("\n#TABLE %s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (test.gene, test.pop, test.trait, test.pwm, test.region, counts, "N/A", p_value, conditional_p_value))
            
            f.write("\nVariants right after these matches:\n")
            for beginning in beginnings_of_match:
                f.write("#AUX Match at position %s\n" % beginning)
                for line in os.popen("bcftools view -f '%%LINE\\n' -r %s:%s-%s -S %s %s/freeze.8.%s.pass_only.phased.bcf | grep -v '##'" % (test.chromosome, beginning, beginning+30, filtered_sample_file, BCF_DIR, test.chromosome)).readlines():
                    if "CHROM" in line:
                        continue
                    fields = line.strip().split("\t")
                    f.write("#AUX %s\t%s\t%s\t%s\t0:%s 1:%s 2:%s\n" % (fields[0], fields[1], fields[3], fields[4], fields[10:].count("0|0"), fields[10:].count("1|0")+fields[10:].count("0|1"), fields[10:].count("1|1") * 2))
        else:
            p_value = original_line.split("\t")[8]
            counts = extract_counts(scan_file)
            f.write("#TABLE %s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (test.gene, test.pop, test.trait, test.pwm, test.region, counts, "N/A", p_value, "N/A"))
            # Gene Ethnicity Trait TF peak_coords Counts SNP_frequency pvalue conditional_pvalue

        f.close()


def is_duffy(region):
    chromosome, start, end = parse_coordinates(region)
    return (chromosome == "chr1")
    #center = (start + end) / 2
    #return (chromosome == "chr1") and (center >= 140000000) and (center <= 180000000)

def is_hemoglobin_chr16(region):
    chromosome, start, end = parse_coordinates(region)
    center = (start + end) / 2
    return ((chromosome == "chr16") and (center <= 400000))

def is_hla(region):
    chromosome, start, end = parse_coordinates(region)
    center = (start + end) / 2
    return (chromosome == "chr6") and (center >= 26000000) and (center <= 33000000)

class Result:
    def __init__(self, test, duffy, hemoglobin_chr16, hla, is_promoter, is_enhancer, linear_pvalue, emmax_pvalue):
        self.test = test
        self.duffy = duffy
        self.hemoglobin_chr16 = hemoglobin_chr16
        self.hla = hla
        self.is_promoter = is_promoter
        self.is_enhancer = is_enhancer
        self.linear_pvalue = linear_pvalue
        self.emmax_pvalue = emmax_pvalue
        self.conditional_didnt_run = None # Reason why conditional analysis was not run
    def to_line(self):
        test = self.test
        comment = ""
        if self.duffy:
            comment = "Duffy/DARC locus"
        elif self.hemoglobin_chr16:
            comment = "Hemoglobin A/B locus"
        elif self.hla:
            comment = "HLA complex"
        elif float(test.originalpvalue) >= 1e-9:
            comment = "Not study-wide significant"
        if self.conditional_didnt_run:
            if comment:
                comment = comment + ". " + self.conditional_didnt_run
            else:
                comment = self.conditional_didnt_run
        #if self.emmax_pvalue:
        #    if comment:
        #        comment = comment + ". " + "Gene: " + self.test.gene.replace("_",",")
        #    else:
        #        comment = "Gene: " + self.test.gene.replace("_",",")
        regulatory = ""
        genes = ""
        if (self.is_promoter and self.test.gene and len(self.test.gene.replace("_",",")) > 2):
            regulatory = "Promoter"
            genes = self.test.gene.replace("_",",") if self.test.gene and len(self.test.gene.replace("_",",")) > 2 else ""
        elif self.is_enhancer:
            regulatory = "Enhancer"
            genes = self.test.gene.replace("_",",") if self.test.gene and len(self.test.gene.replace("_",",")) > 2 else ""
        #"Promoter" if (self.is_promoter and self.test.gene and len(self.test.gene.replace("_",",")) > 2) else "Not a promoter"

        return ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (test.pop, test.trait, test.cell, test.region, test.pwm, test.originalpvalue, test.originalbeta, regulatory, genes, self.linear_pvalue if self.linear_pvalue else "NA", self.emmax_pvalue if self.emmax_pvalue else "NA", comment))


def main():
    results = []
    associated_regions_that_overlap_promoters_enhancers = dict([(l.split(",")[0], (l.split(",")[1], l.split(",")[2])) for l in open("promoters_enhancers_best_associations.tsv", "rt").read().split("\n") if l])
    bcx2_variants_coordinates = [Variant(c,s) for (c,s,e) in [parse_bed_line(l) for l in open("BCX2_GWAS_sentinel.bed", "rt") if l and "track name" not in l and "N/A" not in l] if e-s == 1] # We ignore the 3 BCX variants that are not SNPs, as they are not in the vicinity of any of our regions
    tests = []
    for l in [x for x in open("best_associations.tsv", "rt").read().split("\n") if x][1:]:
        fields = l.strip().split("\t")
        pop,trait,cell,region,pwm,pvalue,beta,sebeta = fields
        duffy = is_duffy(region) and pop in ['African', 'Hispanic'] and trait in ["WBC", "BASO", "EOSIN", "LYMPH", "MONO", "NEUTRO"] # This region is known to affect white blood cell phenotypes in African populations
        hemoglobin_chr16 = is_hemoglobin_chr16(region) and pop in ['African', 'Hispanic'] and trait in ["HGB", "HCT", "MCH", "MCHC", "MCV", "RDW", "RBC"]          # This region is known to affect red blood cell phenotypes
        hla = is_hla(region)

        test = ConditionalTest(pop, trait, cell, region, pwm, pvalue, beta)
        if region in associated_regions_that_overlap_promoters_enhancers:
            # We the region overlaps a promoter, we want to run a conditional test over the BCX2 hematopoietic variants that are in the vicinity
            # Whether or not there's any variant to run a conditional test, we want to know the frequency of the variants that create/disrupt the TFBS in this region. So we need to run find-tfbs in verbose mode and parse the logs
            test.set_cond_variants(set([bcx2_variant for bcx2_variant in bcx2_variants_coordinates if close_enough((bcx2_variant.chromosome, bcx2_variant.position), parse_coordinates(region), tolerance = 1000000)]))
            test.gene, regulatory = associated_regions_that_overlap_promoters_enhancers[region]
            is_promoter = (regulatory == "Promoter")
            is_enhancer = (regulatory == "Enhancer")
            if is_promoter:
                if test in tests:
                    pass
                elif duffy or hemoglobin_chr16 or hla:
                    results.append(Result(test, duffy, hemoglobin_chr16, hla, is_promoter, is_enhancer, None, None)) # No conditional analysis for these regions
                elif float(pvalue) < 1e-9:
                    print("Loaded conditional test definition: %s" % test.name())
                    if test.name() == "test_RENBP_African_RBC_CTCFL_HUMAN.H11MO.0.A_chrX_153945909to153946614":
                        tests.append(test)
                else:
                    results.append(Result(test, duffy, hemoglobin_chr16, hla, is_promoter, is_enhancer, None, None))
            else:
                results.append(Result(test, duffy, hemoglobin_chr16, hla, is_promoter, is_enhancer, None, None))
        else:
            results.append(Result(test, duffy, hemoglobin_chr16, hla, False, False, None, None))

    print("Number of tests: %s" % len(tests))

    thread_number = int(sys.argv[1])
    #with Pool(thread_number) as p:
    #    results.extend(p.map(process_test, tests))
    for t in tests:
        print(process_test(t))

    f = open("summary_table.tsv", "wt")
    header = "POPULATION	TRAIT	CELL	COORDS	TF	PVALUE_EMMAX	BETA_EMMAX	REGULATORY_FEATURE	GENE	CONDITIONAL_LINEAR_PVALUE	CONDITIONAL_EMMAX_PVALUE	COMMENT"
    f.write("\n".join([header] + [r.to_line() for r in results]))
    f.close()

    

if __name__ == '__main__':
    BCF_DIR = sys.argv[2]
    main()
