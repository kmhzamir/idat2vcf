process MERGE_VCF_COLUMNS {
    tag "$meta.id"
    label 'process_low'

    container "docker.io/kmhzamir/dragena:v1.4"
    
    input:
    tuple val(meta), path(merged_tab)  // Concatenated tab file from previous step
    tuple val(meta), path(vcf), path(txt), path(tbi)                           // Raw VCF file

    output:
    tuple val(meta), path("${meta.id}_cleaned.tab.gz")   , emit: cleaned_output

    script:
    """
    python3 <<CODE
import gzip
import pandas as pd

# Define file paths
vcf_file_path = '${vcf}'
annotated_file_path = '${merged_tab}'
output_file_path = '${meta.id}_merged_with_vcf_data.tab.gz'

# Step 1: Get sample name from VCF header
with gzip.open(vcf_file_path, 'rt') as f:
    for line in f:
        if line.startswith('#CHROM'):
            vcf_columns = line.strip().split('\\t')
            sample_name = vcf_columns[9]
            break

# Step 2: Load VCF variant data (exclude ID from output, keep RSID for merging)
vcf_data = []
with gzip.open(vcf_file_path, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\\t')
        vcf_data.append({
            '#Uploaded_variation': fields[2],  # rsID used for merging
            'REF': fields[3],
            'ALT': fields[4],
            'FORMAT': fields[8],
            sample_name: fields[9]
        })

vcf_df = pd.DataFrame(vcf_data)

# Step 3: Read annotated tab file (header is first row, no ##)
annotated_df = pd.read_csv(
    annotated_file_path,
    sep='\\t',
    compression='gzip'
)

# Step 4: Merge on RSID and drop rows that didn’t match if needed (optional)
merged_df = pd.merge(
    annotated_df,
    vcf_df,
    on='#Uploaded_variation',
    how='left'  # use 'inner' if you want to keep only matching RSIDs
)

# Step 5: Write to compressed output file
merged_df.to_csv(output_file_path, sep='\\t', index=False, compression='gzip')
print("Merged file saved to:", output_file_path)
CODE

# Optional: Clean extra # lines if needed
zcat ${meta.id}_merged_with_vcf_data.tab.gz | awk '!/^##/ && (!/^#/ || !seen++)' | bgzip > ${meta.id}_cleaned.tab.gz
    """
}