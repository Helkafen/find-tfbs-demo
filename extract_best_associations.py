import glob
import gzip
import sys

phenotypes_per_cell_type_filename = sys.argv[1]
transcription_factors_per_cell_type_filename = sys.argv[2]
output_filename = sys.argv[3]  # Best results
ethnicity = sys.argv[4]

# type: {trait, (cell, pwm)}
def load_tests():
    # Build the list of association tests that need to be run
    traits_per_cell = {}
    for (c,t) in [(l.split(",")[0], l.split(",")[1].strip()) for l in open(phenotypes_per_cell_type_filename).readlines()[1:]]:
        x = traits_per_cell.get(c,set([]))
        x.add(t)
        traits_per_cell[c] = x
    # type: [(trait, (cell, pwm))]
    tests_to_run = [(trait, (l.split("\t")[0], l.split("\t")[1].strip())) for l in open(transcription_factors_per_cell_type_filename).readlines()[2:] if l for trait in list(traits_per_cell[l.split("\t")[0]])]
    d = {}
    for trait, cell_pwm in tests_to_run:
        if trait not in d:
            d[trait] = set()
        d[trait].add(cell_pwm)
    return d

def upper(ethnicity):
    if ethnicity == "white":
        ethnicity = "European"
    elif ethnicity == "black":
        ethnicity = "African"
    elif ethnicity == "hispanic":
        ethnicity = "Hispanic"
    else:
        raise
    return ethnicity

header = False
output = open(output_filename,"wt")
file_list = glob.glob("epacts/%s/EPACTS_*.epacts.gz" % ethnicity)
best_associations = set()
tests_to_run = load_tests()
for i,file in enumerate(file_list):
    trait = file.split("/")[-1].replace("EPACTS_","").replace(".epacts.gz","")
    ethnicity = file.split("/")[1]
    
    print("Extract from file %s (trait=%s)" % (file, trait))
    f = gzip.open(file, "rt")
    for line in f:
        if not line:
            continue
        if line.startswith("#"):
            if not header:
                output.write("POPULATION\tTRAIT\tCELL\tCOORDS\tTF\tPVALUE\tBETA\tSEBETA\n")
                header = True
        else:
            fields = line.split("\t")
            marker_id = fields[3]
            cell = marker_id.split(",")[0].replace(".bed","").split("_")[-1]
            pwm = marker_id.split(",")[1]
            coords = "chr"+fields[0] + ":" + marker_id.split(",")[-1]
            if (cell,pwm) not in tests_to_run[trait]:
                continue
            try:
                p_value = fields[10]
                beta = fields[11]
                sebeta = fields[12]
            except IndexError:
                print(file)
                print(fields)
                raise IndexError
            if p_value == "NA":
                continue
            if p_value == "-nan":
                continue
            if p_value == "nan":
                continue
            if float(p_value) < 0.00001:
                best_associations.add("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (upper(ethnicity), trait, cell, coords, pwm, p_value, beta, sebeta))
    f.close()
best_associations = list(sorted(best_associations, reverse = False, key = lambda x: float(x.split("\t")[5])))
output.write("\n".join(best_associations))
output.close()
