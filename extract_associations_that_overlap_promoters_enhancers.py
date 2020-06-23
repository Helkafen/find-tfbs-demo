import sys

associations_filename = sys.argv[1]
ensembl_promoters_filename = sys.argv[2]
ensembl_enhancers_filename = sys.argv[3]
ensembl_genes_filename = sys.argv[4]
output_filename = sys.argv[5]

def overlaps(x,y):
    (c,s,e) = x
    (chromosome, start, end) = y
    if not (c == chromosome):
        return False
    if start <= e and start >= s:
        return True
    if end >= s and end <= e:
        return True
    return False

def close_enough(x,y):
    (c,s,e) = x
    (chromosome, start, end) = y
    if not (c == chromosome):
        return False
    center1 = (s+e) / 2
    center2 = (start+end) / 2
    return (abs(center1 - center2) < 1000000)

def parse_coordinates(s):
    chromosome,rest = s.split(":")
    start,end = [int(x) for x in rest.split("-")]
    return (chromosome, start, end)

def parse_coordinates_tab(s):
    (chrom, s, e) = s.strip().split("\t")[:3]
    return ("chr"+chrom, int(s), int(e))

f = open(associations_filename, "rt")
g = open(ensembl_promoters_filename, "rt")
h = open(ensembl_genes_filename, "rt")
e = open(ensembl_enhancers_filename, "rt")
o = open(output_filename, "wt")
ensembl_promoters = set([parse_coordinates_tab(x) for x in g.readlines()])
ensembl_enhancers = set([parse_coordinates_tab(x) for x in e.readlines()])
ensembl_genes = set([(parse_coordinates_tab(x), x.split("\t")[3], x.split("\t")[4]) for x in h.readlines()])
header = f.readline()
for l in f:
    coordinates = parse_coordinates(l.split("\t")[3])
    overlaps_promoter = False
    overlaps_enhancer = False
    for ensembl_promoter in ensembl_promoters:
        if overlaps(ensembl_promoter, coordinates):
            overlaps_promoter = True
            break
    for ensembl_enhancer in ensembl_enhancers:
        if overlaps(ensembl_enhancer, coordinates):
            overlaps_enhancer = True
            break
    if overlaps_promoter:
        all_gene_overlaps = []
        for (gene_coords,type,description) in ensembl_genes:
            if overlaps(gene_coords, ensembl_promoter):
                name_fields = [x.replace("Name=","") for x in description.split(";") if x.startswith("Name=")]
                assert(len(name_fields) == 1)
                gene_name = name_fields[0]
                all_gene_overlaps.append(gene_name)
        o.write("%s:%s-%s,%s,Promoter\n" % (coordinates[0], coordinates[1], coordinates[2], "_".join(all_gene_overlaps)))
    elif overlaps_enhancer:
        all_gene_nearby = []
        for (gene_coords,type,description) in ensembl_genes:
            if close_enough(gene_coords, ensembl_enhancer):
                name_fields = [x.replace("Name=","") for x in description.split(";") if x.startswith("Name=")]
                assert(len(name_fields) == 1)
                gene_name = name_fields[0]
                all_gene_nearby.append(gene_name)
        o.write("%s:%s-%s,%s,Enhancer\n" % (coordinates[0], coordinates[1], coordinates[2], "_".join(all_gene_nearby)))

o.close()
f.close()
g.close()
e.close()