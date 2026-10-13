from Bio import SeqIO
from argnorm.lib import homolog containment  # or argNorm's main normalization function
# Note: Import the exact argNorm module function based on your input tool/database

# 1. Use Biopython to read a sequence file or parse results
fasta_file = "contigs.fasta"
for record in SeqIO.parse(fasta_file, "fasta"):
    # Suppose you extracted or predicted a gene name / tool output from sequence analysis
    raw_gene_name = "blaOXA-1"  # Example raw name from an annotation step
    source_tool = "amrfinderplus"  # e.g., amrfinderplus, resfinder, deeparg
    
    # 2. Normalize using argNorm
    # argNorm normalizes gene names/categories to the Antibiotic Resistance Ontology (ARO)
    from argnorm.normalization import normalize
    
    normalized_result = normalize(raw_gene_name, catalog=source_tool)
    
    print(f"Original: {raw_gene_name} ({source_tool}) -> ARO: {normalized_result}")
