process ANNOTATE_PHENOTYPE {
    tag "$meta.id"
    label 'process_low'

    container "docker.io/kmhzamir/dragena:v1.4"
    
    input:
    tuple val(meta), path(merged_tab)  // Concatenated tab file from previous step
    path(xml)                           

    output:
    tuple val(meta), path("${meta.id}_annotated_clinvar_phenotype.tab.gz")   , emit: annotated_clinvar_phenotype

    script:
    """
    python3 <<CODE
import csv
import gzip
import xml.etree.ElementTree as ET

# === CONFIG ===
tab_input = "${merged_tab}"             
tab_output = "${meta.id}_annotated_clinvar_phenotype.tab"     
clinvar_xml = "${xml}"        
# ==============

# Step 1: Parse ClinVar XML and index annotations by RSID
def build_clinvar_rsid_dict(xml_path):
    rsid_dict = {}

    context = ET.iterparse(xml_path, events=("end",))
    for event, elem in context:
        if elem.tag.endswith("VariationArchive"):
            rsids = [xref.attrib.get("ID") for xref in elem.findall(".//XRef") if xref.attrib.get("DB") == "dbSNP"]
            for rsid in rsids:
                annotations = []
                for rcv in elem.findall(".//RCVAccession"):
                    classification = rcv.findtext(".//RCVClassifications/GermlineClassification/Description")
                    conditions = [c.text for c in rcv.findall(".//ClassifiedCondition") if c.text]
                    for condition in conditions:
                        if classification and condition:
                            annotations.append(f"{classification}:{condition}")
                if annotations:
                    rsid_dict[rsid] = " | ".join(annotations)
            elem.clear()
    return rsid_dict

# Step 2: Read .tab.gz, match RSID, and write output
def annotate_tab_file(input_path, output_path, rsid_annotations):
    with gzip.open(input_path, 'rt') as infile, open(output_path, 'w', newline='') as outfile:
        reader = csv.DictReader(infile, delimiter='\t')
        fieldnames = reader.fieldnames + ["Clinvar Annotations"]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        for row in reader:
            raw_rsid = row['RSID']
            rsid = raw_rsid.replace('rs', '')
            annotation = rsid_annotations.get(rsid, "-")
            row["Clinvar Annotations"] = annotation
            writer.writerow(row)

# === RUN ===
rsid_annotations = build_clinvar_rsid_dict(clinvar_xml)
annotate_tab_file(tab_input, tab_output, rsid_annotations)
CODE

# Compress using bgzip for Tabix compatibility
bgzip -c ${meta.id}_annotated_clinvar_phenotype.tab > ${meta.id}_annotated_clinvar_phenotype.tab.gz
    """
}
