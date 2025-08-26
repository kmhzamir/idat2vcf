process ANNOTATE_GT {
    tag "$meta.id"
    label 'process_low'

    container "docker.io/kmhzamir/dragena:v1.4"

    input:
    tuple val(meta), path(cleaned_tab)  // Assume this is a .tab.gz or .tab file

    output:
    tuple val(meta), path("${meta.id}_annotated_gt.tab.gz"), emit: gt_output

    script:
    """
    gunzip -c ${cleaned_tab} > temp_input.tab  # decompress if gzipped

    python3 - <<EOF
import csv

input_file = "temp_input.tab"
output_file = "${meta.id}_annotated_gt.tab"

with open(input_file, "r") as fin, open(output_file, "w", newline="") as fout:
    reader = csv.reader(fin, delimiter="\t")
    writer = csv.writer(fout, delimiter="\t")

    header = next(reader)
    ref_index = header.index("REF")
    alt_index = header.index("AlternateAlleleVCF")
    format_index = header.index("FORMAT")
    sample_column_index = format_index + 1  # first sample after FORMAT

    # Add both "GT" and "GT_Allele" after FORMAT
    new_header = header[:format_index] + ["GT", "GT_Allele"] + header[format_index:]
    writer.writerow(new_header)

    for row in reader:
        ref = row[ref_index]
        alt = row[alt_index]
        format_fields = row[format_index].split(":")
        sample_data = row[sample_column_index].split(":")

        # Default values
        gt_value = "./."
        allele_gt = "./."

        # Extract GT (numeric)
        try:
            gt_index = format_fields.index("GT")
            gt_value = sample_data[gt_index]  # e.g. "0/0"
        except (ValueError, IndexError):
            pass

        # Convert GT numeric to allele form
        alleles = [ref] + alt.split(",")
        if "/" in gt_value:
            try:
                a1, a2 = map(int, gt_value.split("/"))
                allele_gt = f"{alleles[a1]}/{alleles[a2]}"
            except (ValueError, IndexError):
                pass

        # Insert GT and GT_Allele
        new_row = row[:format_index] + [gt_value, allele_gt] + row[format_index:]
        writer.writerow(new_row)

EOF

    # Compress output using bgzip
    bgzip -f ${meta.id}_annotated_gt.tab
    """
}
