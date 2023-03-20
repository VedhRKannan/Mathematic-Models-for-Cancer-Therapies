input_file = "1nql_biopython.tsv"
output_file = "new.tsv"

with open(input_file, "r") as f_in, open(output_file, "w") as f_out:
    for line in f_in:
        cols = line.strip().split("\t")
        # Extract amino acid type from column 1
        amino_acid = cols[0].split(":")[2]
        # Append amino acid type as a new column
        new_line = "\t".join(cols + [amino_acid]) + "\n"
        f_out.write(new_line)
