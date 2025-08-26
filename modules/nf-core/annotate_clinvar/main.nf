process ANNOTATE_CLINVAR {
    tag "$meta.id"
    label 'process_low'

    container "docker.io/kmhzamir/dragena:v1.4"
    
    input:
    tuple val(meta), path(merged_tab)  // Concatenated tab file from previous step
    path(tsv)                           

    output:
    tuple val(meta), path("${meta.id}_annotated_clinvar.tab.gz")   , emit: annotated_clinvar

    script:
    """
    python3 <<CODE
import pandas as pd

# Load ClinVar file
clinvar = pd.read_csv("${tsv}", sep="\\t")

# Keep only desired columns + RSID
clinvar = clinvar[["RSID", "AlternateAlleleVCF", "Type", "Name", "ClinicalSignificance"]]

# Load VEP file (tab.gz)
vep = pd.read_csv("${merged_tab}", sep="\\t", compression="gzip")

# Standardize column name for merging
vep.rename(columns={"#Uploaded_variation": "RSID"}, inplace=True)
vep.rename(columns={"ALT": "AlternateAlleleVCF"}, inplace=True)


# Merge on RSID
merged = vep.merge(clinvar, on=["RSID", "AlternateAlleleVCF"], how="left")

# Save as tab-separated file
merged.to_csv("${meta.id}_annotated_clinvar.tab", sep="\\t", index=False)
CODE

# Compress using bgzip for Tabix compatibility
bgzip -c ${meta.id}_annotated_clinvar.tab > ${meta.id}_annotated_clinvar.tab.gz
    """
}
