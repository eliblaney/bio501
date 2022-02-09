# Import the regular expressions package
import re

# Open the messed-up file
grape_promoters = open('mangled_promoters.txt')

# These two variables keep track of the gene with the highest
# number of As and Ts. We keep updating it until we have the
# single gene with the highest number of all.
gene_most_at = ""
most_at = -1
# Iterate through all of the mangled promoter sequences
for line in grape_promoters:
    # We can separate the gene name and DNA sequence into two variables
    gene, seq = line.rstrip().split(",")
    # This condensed line will remove everything except As and Ts, making
    # the DNA sequence uppercase in the process. This takes care of two problems
    # at once: removing invalid base characters and counting the number of As
    # and Ts. By removing everything that's not an A or T and taking the
    # length of the string, we accomplish both tasks in one short operation.
    at = len(re.sub(r"[^AT]", "", seq.upper()))
    # Here, we update our previous variables if we find a gene with a higher
    # number of As and Ts than anything we've found before.
    if at > most_at:
        gene_most_at = gene
        most_at = at

# Print out the gene name and the total number of As and Ts that we found.
print(gene_most_at, most_at)

# We don't need the file anymore, close it.
grape_promoters.close()