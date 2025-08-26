process POSTPROCESS_VEP {
    tag "$meta.id"
    label 'process_low'

    container "docker.io/kmhzamir/dragena:v1.4"
    
    input:
    tuple val(meta), path(merged_tab)  // Concatenated tab file from previous step

    output:
    tuple val(meta), path("${meta.id}_cleaned.tab.gz")   , emit: cleaned_output

    script:
"""
python3 <<'CODE'
import pandas as pd
import gzip
from io import StringIO

def merge_vep_annotations(input_gz_file, temp_output_file):
    with gzip.open(input_gz_file, 'rt') as f:
        header_lines = []
        col_header = ""
        data_lines = []

        for line in f:
            if line.startswith("##"):
                continue  # Skip lines starting with ##
            elif line.startswith("#"):
                col_header = line.lstrip("#").strip()
            else:
                data_lines.append(line)

    if not col_header:
        raise ValueError("Column header line not found in the input file.")

    # Combine all into one string to load into DataFrame
    full_data = col_header + "\\n" + "".join(data_lines)
    df = pd.read_csv(StringIO(full_data), sep="\\t", dtype=str).fillna("-")

    # Merge logic
    group_cols = ["Uploaded_variation", "Location"]
    grouped = df.groupby(group_cols, sort=False)
    merged_rows = []

    for keys, group in grouped:
        merged_row = dict(zip(group_cols, keys))
        for col in df.columns:
            if col in group_cols:
                continue
            values = group[col].unique()
            values = [v for v in values if v != "-"]
            if not values:
                merged_value = "-"
            elif len(values) == 1:
                merged_value = values[0]
            else:
                merged_value = "|".join(sorted(set(values)))
            merged_row[col] = merged_value
        merged_rows.append(merged_row)

    # Write plain .tab file (with only one header line)
    with open(temp_output_file, "w") as out:
        out.write("#" + col_header + "\\n")
        merged_df = pd.DataFrame(merged_rows)
        merged_df.to_csv(out, sep="\\t", index=False, header=False)

# Run it
merge_vep_annotations("${merged_tab}", "${meta.id}_merged_output.tab")
CODE

# Compress using bgzip (ready for tabix)
bgzip -c ${meta.id}_merged_output.tab > ${meta.id}_cleaned.tab.gz
"""
}