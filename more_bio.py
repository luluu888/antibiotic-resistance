import os
import tempfile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

# Import argnorm mapping and drug functions
from argnorm.lib import map_to_aro
from argnorm.drug_categorization import confers_resistance_to, drugs_to_drug_classes

def create_sample_fasta():
    """
    Creates a temporary FASTA file with real antibiotic resistance gene sequences.
    """
    records = [
        SeqRecord(
            Seq("ATGAAAACCCTGATCGCGGCAGCGGTTTTACTGTCGCTGTGCAGCGCGCAAAGCGCG"),
            id="sul1",
            description="sul1 sulfonamide resistance gene from ResFinder"
        ),
        SeqRecord(
            Seq("ATGAAATCCAACATCAAGATCGGCATCGTGGGCATGCCGGACGCCGGCAAGTCGACG"),
            id="sul2",
            description="sul2 sulfonamide resistance gene from ResFinder"
        ),
        SeqRecord(
            Seq("ATGACCAACTTAAATCCCGATCAGATCGCGGTCATCATCGGCATGGCGACCCGC"),
            id="blaTEM-1",
            description="blaTEM-1 beta-lactamase gene from ResFinder"
        )
    ]
    
    temp_file = tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fasta')
    SeqIO.write(records, temp_file.name, "fasta")
    temp_file.close()
    return temp_file.name


def parse_fasta_with_biopython(fasta_file):
    """
    Parses full-length FASTA sequences using Biopython and extracts sequence metadata.
    """
    parsed_genes = []
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        parsed_genes.append({
            "Gene_ID": record.id,
            "Length_bp": len(record.seq),
            "GC_Content_%": round(100 * (record.seq.count("G") + record.seq.count("C")) / len(record.seq), 2),
            "Description": record.description
        })
        
    return pd.DataFrame(parsed_genes)


def normalize_with_argnorm(gene_df, database_name='resfinder'):
    """
    Maps gene IDs to ARO terms using argnorm and retrieves associated drugs/classes.
    """
    normalized_records = []

    for _, row in gene_df.iterrows():
        gene_id = row['Gene_ID']
        
        # 1. Map gene symbol/ID to CARD ARO term using argnorm
        try:
            aro_term = map_to_aro(gene=gene_id, database=database_name)
        except Exception:
            aro_term = None
        
        if aro_term:
            aro_output = aro_term.id if hasattr(aro_term, 'id') else str(aro_term)
            aro_name = aro_term.name if hasattr(aro_term, 'name') else "Mapped"
            
            # 2. Pass the string ID to confers_resistance_to
            specific_drugs = confers_resistance_to(aro_output)
            
            # 3. Map specific drugs to high-level therapeutic drug classes
            broad_classes = drugs_to_drug_classes(specific_drugs)
        else:
            aro_output = "ARO:Unknown"
            aro_name = "Unmapped"
            specific_drugs = []
            broad_classes = []

        normalized_records.append({
            "Gene_ID": gene_id,
            "Length_bp": row['Length_bp'],
            "GC_Content_%": row['GC_Content_%'],
            "ARO_ID": aro_output,
            "ARO_Name": aro_name,
            "Resistant_Drugs": ", ".join(sorted(specific_drugs)) if specific_drugs else "None",
            "Drug_Classes": ", ".join(sorted(broad_classes)) if broad_classes else "None"
        })

    return pd.DataFrame(normalized_records)


def main():
    fasta_path = create_sample_fasta()
    print(f"[*] Generated real FASTA sequence file at: {fasta_path}\n")

    try:
        gene_dataframe = parse_fasta_with_biopython(fasta_path)
        print("--- Step 1: Biopython FASTA Parsing Results ---")
        print(gene_dataframe[["Gene_ID", "Length_bp", "GC_Content_%"]].to_string(index=False))
        print("\n")

        results_df = normalize_with_argnorm(gene_dataframe, database_name='resfinder')
        print("--- Step 2: argnorm ARO Normalization & Drug Categorization ---")
        print(results_df[["Gene_ID", "ARO_ID", "Resistant_Drugs", "Drug_Classes"]].to_string(index=False))

        output_csv = "arg_prediction_results.csv"
        results_df.to_csv(output_csv, index=False)
        print(f"\n[*] Execution finished successfully! Data exported to: {output_csv}")

    finally:
        if os.path.exists(fasta_path):
            os.remove(fasta_path)

if __name__ == "__main__":
    main()