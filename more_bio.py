import os
import tempfile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd


from argnorm.lib import map_to_aro
from argnorm.drug_categorization import confers_resistance_to, drugs_to_drug_classes

def create_sample_fasta():

    records = [
        SeqRecord(
            Seq("ATGAAAACCCTGATCGCGGCAGCGGTTTTACTGTCGCTGTGCAGCGCGCAAAGCGCG"),
            id="sul1_2_U12338",  
            description="sul1 sulfonamide resistance gene"
        ),
        SeqRecord(
            Seq("ATGAAATCCAACATCAAGATCGGCATCGTGGGCATGCCGGACGCCGGCAAGTCGACG"),
            id="sul2_1_AF542061",  
            description="sul2 sulfonamide resistance gene"
        )
    ]
    
    temp_file = tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fasta')
    SeqIO.write(records, temp_file.name, "fasta")
    temp_file.close()
    return temp_file.name


def parse_fasta_with_biopython(fasta_file):

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
   
    normalized_records = []

    for _, row in gene_df.iterrows():
        gene_id = row['Gene_ID']
        aro_output = "ARO:Unknown"
        aro_name = "Unmapped"
        specific_drugs = []
        broad_classes = []
        
        try:
            # 1. Map full database identifier to CARD ARO term object
            aro_term = map_to_aro(gene=gene_id, database=database_name)
            
            if aro_term:
                aro_output = aro_term.id if hasattr(aro_term, 'id') else str(aro_term)
                aro_name = aro_term.name if hasattr(aro_term, 'name') else "Mapped"
                
                # 2. Get specific resistant drugs
                specific_drugs = confers_resistance_to(aro_output)
                
                # 3. Map specific drugs to high-level therapeutic drug classes
                broad_classes = drugs_to_drug_classes(specific_drugs)
                
        except Exception:
            pass

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
    print(f"[*] Generated FASTA sequence file at: {fasta_path}\n")

    try:
       
        gene_dataframe = parse_fasta_with_biopython(fasta_path)
        print("--- Step 1: Biopython FASTA Parsing Results ---")
        print(gene_dataframe[["Gene_ID", "Length_bp", "GC_Content_%"]].to_string(index=False))
        print("\n")

       
        results_df = normalize_with_argnorm(gene_dataframe, database_name='resfinder')
        print("--- Step 2: argnorm ARO Normalization & Drug Categorization ---")
        print(results_df[["Gene_ID", "ARO_ID", "Resistant_Drugs", "Drug_Classes"]].to_string(index=False))

        # Step 3: Export results to CSV
        output_csv = "arg_prediction_results.csv"
        results_df.to_csv(output_csv, index=False)
        print(f"\n[*] Execution finished successfully! Data exported to: {output_csv}")

    finally:
        # Clean up temporary file
        if os.path.exists(fasta_path):
            os.remove(fasta_path)

if __name__ == "__main__":
    main()