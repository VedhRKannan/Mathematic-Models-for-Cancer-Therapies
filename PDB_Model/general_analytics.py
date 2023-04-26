with open('1nql_biopython.tsv', 'r') as f:
    lines = f.readlines()

for line in lines:
    columns = line.strip().split('\t')
    residue = str(columns[3])
    if residue not in ["General", "Glycine", "Pre-Pro", "Proline"]:
        print(residue)
