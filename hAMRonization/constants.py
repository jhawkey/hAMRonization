# genetic_variation type associated with detected resistance
NUCLEOTIDE_VARIANT = "nucleotide_variant_detected"
AMINO_ACID_VARIANT = "protein_variant_detected"
GENE_PRESENCE = "gene_presence_detected"
INACTIVATING_VARIANT = "inactivating_mutation_detected"
INACTIVATING_VARIANT_PARTIAL = "inactivating_mutation_partial"
PROMOTER_VARIANT = "promoter_variant_detected"

# conversion from single letter amino acid to three letter code, in line with HGVS
aa_conversion = {'G': 'Gly', 'A': 'Ala', 'S': 'Ser', 'P': 'Pro', 'T': 'Thr', 'C': 'Cys', 'V': 'Val', 'L': 'Leu', 'I': 'Ile', 
                 'M': 'Met', 'N': 'Asn', 'Q': 'Gln', 'K': 'Lys', 'R': 'Arg', 'H': 'His', 'D': 'Asp', 'E': 'Glu', 'W': 'Trp', 
                 'Y': 'Tyr', 'F': 'Phe', '*': 'Ter', 'STOP': 'Ter'}