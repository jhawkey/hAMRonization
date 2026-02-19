#!/usr/bin/env python

import csv
import warnings
import re
from .Interfaces import hAMRonizedResultIterator
from hAMRonization.constants import (
    NUCLEOTIDE_VARIANT,
    AMINO_ACID_VARIANT,
    GENE_PRESENCE,
    INACTIVATING_VARIANT,
    PROMOTER_VARIANT,
    INACTIVATING_VARIANT_PARTIAL,
    aa_conversion
)

required_metadata = [
    "analysis_software_version",
    "reference_database_version",
    "input_file_name",
]

# Alternative metadata when using Name column
required_metadata_use_name = [
    "analysis_software_version",
    "reference_database_version",
]


class AmrFinderPlusIterator(hAMRonizedResultIterator):

    # these aliases are to support newer versions of the AMRFinderPlus output, which have different column names for some fields
    column_aliases = {
        "Element type": "Type",
        "Element subtype": "Subtype",
        "Gene symbol": "Element symbol",
        "Sequence name": "Element name",
        "Protein identifier": "Protein id",
        "% Coverage of reference sequence": "% Coverage of reference",
        "% Identity to reference sequence": "% Identity to reference",
        "Accession of closest sequence": "Closest reference accession",
        "Name of closest sequence": "Closest reference name",
        "HMM id": "HMM accession"
    }

    nuc_field_map = {
        "Protein id": None,
        "Contig id": "input_sequence_id",
        "Start": "input_gene_start",
        "Stop": "input_gene_stop",
        "Strand": "strand_orientation",
        "Element symbol": "gene_symbol",
        "Element name": "gene_name",
        "Scope": None,
        "Type": None,
        "Subtype": None,
        "Class": "drug_class",
        "Subclass": "antimicrobial_agent",
        "Method": None,
        "Target length": "input_gene_length",
        "Reference sequence length": "reference_gene_length",
        "% Coverage of reference": "coverage_percentage",
        "% Identity to reference": "sequence_identity",
        "Alignment length": None,
        "Closest reference accession": "reference_accession",
        "Closest reference name": None,
        "HMM accession": None,
        "HMM description": None,
        "Hierarchy node": "gene_reference_id",
        # Fields we compute below (not in TSV)
        "amino_acid_mutation": "amino_acid_mutation",
        "nucleotide_mutation": "nucleotide_mutation",
        "genetic_variation_type": "genetic_variation_type",
    }

    # AMP outputs the same column set for nuc and prot detections,
    # with Start and Stop always in nt units; however target and
    # reference length are reported in AA for proteins.
    prot_field_map = nuc_field_map.copy()
    prot_field_map.update({
        "Target length": "input_protein_length",
        "Reference sequence length": "reference_protein_length"
    })

    def __init__(self, source, metadata):
        metadata["analysis_software_name"] = "amrfinderplus"
        metadata["reference_database_name"] = "NCBI Reference Gene Database"
        self.metadata = metadata
        self.use_name_column = metadata.get("input_file_name") == "useName"

        # We pass None for the field_map as it differs depending on
        # whether we return a nucleotide or protein variant detection.
        # TODO: refactor field_map out of super's constructor, and make
        # it a parameter on super's hARMonize().
        super().__init__(source, None, self.metadata)

    def parse(self, handle):
        """
        Read each and return it
        """
        skipped_truncated = 0
        reader = csv.DictReader(handle, delimiter="\t")
        
        # Check if we're using Name column and validate
        if self.use_name_column:
            if "Name" not in reader.fieldnames:
                raise ValueError(
                    "input_file_name was set to 'useName' but the input file does not contain a 'Name' column. "
                    "Please provide a valid input_file_name value or use an AMRFinderPlus file that includes the 'Name' column."
                )
        
        for result in reader:

            # Replace NA value with None for consistency
            for field, value in result.items():
                if value == "NA":
                    result[field] = None

            self._apply_column_aliases(result)

            # Skip reported virulence genes
            if result.get("Type") == "VIRULENCE":
                continue

            # "POINT" indicates mutational resistance
            # amrfinderplus has no special fields but the mutation itself is
            # appended to the symbol name so we want to split this
            result['amino_acid_mutation'] = None
            result['nucleotide_mutation'] = None
            result['genetic_variation_type'] = GENE_PRESENCE

            if result.get("Subtype").startswith("POINT"):
                # get the mutation type if there is one
                self._parse_mutation(result)

            # no longer skipping internal stops, we want to
            # report these as inactivating mutaations
            if "INTERNAL_STOP" in result.get('Method', ''):
                #skipped_truncated += 1
                result['genetic_variation_type'] = INACTIVATING_VARIANT
            # now we just have to work out what to do with partial hits...
            # we want to keep these as inactivated variants, but we also
            # need some way to record them as partial hits, as this will be
            # important for interpretation
            if result['Method'].startswith("PARTIAL"):
                result['genetic_variation_type'] = INACTIVATING_VARIANT_PARTIAL

            # Determine the field_map to use depending on the method used
            # The following seems to cover all bases with a minimum of fuss
            have_prot = result.get("Protein id") is not None
            method = result.get("Method") or ""
            if method.endswith('P') or method.endswith('X'):
                field_map = self.prot_field_map
            elif method.endswith('N'):
                field_map = self.nuc_field_map
            elif method in ['COMPLETE', 'HMM']:
                field_map = self.prot_field_map if have_prot else self.nuc_field_map
            else:
                warnings.warn(f"Assuming unknown method {method} implies a protein detection"
                              f" in {self.metadata['input_file_name']}")
                field_map = self.prot_field_map

            # This uses the "override hack" that should perhaps be cleaned up
            
            # Handle input_file_name based on whether we're using Name column
            if self.use_name_column:
                # Use the Name column value for input_file_name
                self.metadata["input_file_name"] = result.get("Name", "unknown")
            
            yield self.hAMRonize(result, self.metadata, field_map)

        if skipped_truncated > 0:
            warnings.warn(f"Skipping {skipped_truncated} records with INTERNAL_STOP "
                          f"from {self.metadata['input_file_name']}")

    def _apply_column_aliases(self, result):
        """Apply the column aliases so we can support different version of AMRFinderPlus output with different column names"""
        for alias, canonical in self.column_aliases.items():
            if canonical not in result or result[canonical] is None:
                if alias in result and result[alias] is not None:
                    result[canonical] = result[alias]
    
    def _parse_mutation(self, result):

        value_in_symbol = result.get("Element symbol")
        gene_symbol, mutation = value_in_symbol.rsplit("_", 1)

        # this means it is a protein mutation
        if result['Method'] in ["POINTX", "POINTP"]:
            # variation type is protein variant, unless subtype is POINT_DISRUPT
            # which means the variation type is inactivation mutation
            if result['Subtype'] == "POINT_DISRUPT":
                result['genetic_variation_type'] = INACTIVATING_VARIANT
            else:
                result['genetic_variation_type'] = AMINO_ACID_VARIANT
            # extract the relevant parts of the mutation
            pattern = re.compile(r"(\D+)(\d+)(\D+)")
            ref, pos, alt = pattern.match(mutation).groups()
            # convert the single letter AA code to the 3 letter code
            # note that we need to determine if we've got a simple substitution of ref to alt
            # or do we have a deletion or an insertion, or a frameshift?
            # gyrA_S83L -> p.Ser83Leu, this is a substitution
            # penA_D346DD -> p.345_346insAsp, this is an insertion
            # ompK35_E42RfsTer47 -> p.Glu42ArgfsTer47, this is an inactivating frameshift
            # ompK35_Y36Ter -> p.Tyr36Ter, this is an inactivating frameshift

            # okay so if there are two or more characters in alt, and no Ter, then we have an insertion
            if len(alt) > 1 and alt != 'STOP' and 'Ter' not in alt:
                # then we have an insertion, and the inserted AAs is the second character onwards of alt
                # need to get all the inserted AAs converted
                alt_string = ''
                for char in alt[1:]:
                    alt_string += aa_conversion.get(char)
                # our coordinates are the original pos, and original - 1
                pos_coords = str(int(pos) - 1) + "_" + str(pos)
                result['amino_acid_mutation'] = f"p.{pos_coords}ins{alt_string}"
                return result
            if len(alt) > 1 and 'Ter' in alt:
                # then we have a frameshift leading to a stop codon
                # if the start of alt is Ter, then it's simply ref pos *
                if alt.startswith('Ter'):
                    result['amino_acid_mutation'] = f"p.{aa_conversion.get(ref)}{pos}Ter"
                    return result
                else:
                    # otherwise we need to get the first aa and convert it
                    alt = aa_conversion.get(alt[0])
                    # then grab the number after Ter
                    fs_pos = re.search(r'Ter(\d+)', mutation).group(1)
                    result['amino_acid_mutation'] = f"p.{aa_conversion.get(ref)}{pos}{alt}fsTer{fs_pos}"
                    return result
            # this is the option if we have a simple conversion of one aa to another
            else:
                ref = aa_conversion.get(ref)
                alt = aa_conversion.get(alt)
                result['amino_acid_mutation'] = f"p.{ref}{pos}{alt}"
                return result
        
        elif result['Method'] == "POINTN":
            # we need to extract the relevant parts, this will be different because we may have promoter mutations
            pattern = re.compile(r'^([A-Za-z]+)(-?\d+)([A-Za-z]+)$')
            ref, pos, alt = pattern.match(mutation).groups()
            if '-' in pos:
                result['genetic_variation_type'] = PROMOTER_VARIANT
            else:
                result['genetic_variation_type'] = NUCLEOTIDE_VARIANT
            # if there's a '-' in the position, then this is a promoter mutation
            # if 'del' is in mutation, then we need to convert to c.PosNTdel.
            if alt == 'del':
                result['nucleotide_mutation'] = f"c.{pos}{ref}del"
                return result
            # otherwise it's more like 23S_G2032T -> c.2032G>T, with a - if it's in the promoter.
            else:
                result['nucleotide_mutation'] = f"c.{pos}{ref}>{alt}"
                return result
