import sys

makefile = sys.argv[1]
positions_file = sys.argv[2]
step = int(sys.argv[3])

chromosomes = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"]

print("Optimization: Fool EPACTS by removing some temporary files that would be empty anyway")

def load_max_position_by_chromosome(filename):
    d = {}
    for l in open(filename, "rt"):
        x = l.split("\t")
        c = x[0]
        if x[1].startswith("POS"):
            continue
        p = int(x[1].strip())
        if c not in d:
                d[c] = 0
        d[c] = max(p, d[c])
    return d

def load_positions_by_chromosome(filename):
    d = {}
    for chromosome in chromosomes:
        filename_chr = filename.replace("chr1", "chr" + chromosome)
        for l in open(filename_chr, "rt"):
            x = l.split("\t")
            c = x[0]
            if x[1].startswith("POS"):
                continue
            p = int(x[1].strip())
            if c not in d:
                d[c] = set()
            d[c].add(round_step(p))
    return d

def round_step(x):
    y = x - 1 # We need 0-based positions
    return y - (y % step)

max_positions = load_positions_by_chromosome(positions_file)

f = open(makefile, "rt")
o = open(makefile + ".optimized", "wt")
for l in f:
    # Shorten list of required files. We don't need to generate empty .epacts files
    if ".epacts.OK: " in l:
        fields = l.strip("\n").split(" ")
        good_fields = []
        for field in fields:
            if field.endswith(".epacts"):
                (c,s,e) = field.split(".")[-4:-1]
                if c not in max_positions:
                    pass # Not a chromosome we're considering
                elif round_step(int(s)) in max_positions[c]:
                    good_fields.append(field)
                else:
                    pass # We don't need to produce this temporary epacts file. It would be empty
            else:
                good_fields.append(field)
        o.write(" ".join(good_fields) + "\n")
    # Generate the .epacts.gz with a different method. This is useful when the list of .epacts files is too large
    elif "epacts-cat" in l:
        target = l.split(" ")[-1].strip()
        pattern = target.replace(".epacts.gz", ".*.epacts") #/home/seb/dell/tfbs/tfbs_replicable/epacts/black/EPACTS_HGB.epacts.gz
        o.write("""\tls  %s | sort -t'.' -k3,3n -k4,4n -r | xargs cat | awk  'BEGIN {C=0} {if ($$1 == "#CHROM" && C==0) {C = 1; print;} else if ($$1 != "#CHROM") {print;}}' | bgzip -c > %s""" % (pattern, target) + "\n")
    elif "tabix -f -pbed" in l:
        pass
    elif "epacts-plot" in l:
        pass
    elif "rm -f" in l:
        pass
    else:
        o.write(l.strip("\n") + "\n")
o.close()
f.close()
