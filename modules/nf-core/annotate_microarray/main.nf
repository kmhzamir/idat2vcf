process ANNOTATE_MICROARRAY {
    tag "$meta.id"
    label 'process_low'

    container "docker.io/kmhzamir/dragena:v1.4"

    input:
    tuple val(meta), path(vcf), path(txt), path(tbi)
    tuple val(meta), path(gt_tab)

    output:
    tuple val(meta), path("${meta.id}_final.tab.gz"), emit: final_output

    script:
    """
    python3 <<CODE
import pandas as pd
import subprocess

def parse_gsgt_file(filepath):
    gsgt_data = {}
    with open(filepath, 'r') as file:
        in_data_section = False
        for line in file:
            line = line.strip()
            if line == "[Data]":
                in_data_section = True
                continue
            if in_data_section:
                if not line or line.startswith("SNP Name"):
                    continue
                parts = line.split('\\t')
                if len(parts) >= 4:
                    rsid, _, allele1, allele2 = parts[:4]
                    gsgt_data[rsid] = f"{allele1}/{allele2}"
    return gsgt_data

def annotate_tab_file(tab_gz_path, gsgt_dict, output_tab_path):
    df = pd.read_csv(tab_gz_path, sep='\\t', compression='gzip')
    df['Microarray Result'] = df['RSID'].map(gsgt_dict).fillna('NA')
    cols = list(df.columns)
    if 'FORMAT' in cols:
        cols.remove('Microarray Result')  # avoid duplication
        insert_idx = cols.index('FORMAT')
        cols.insert(insert_idx, 'Microarray Result')
        df = df[cols]
    df.to_csv(output_tab_path, sep='\\t', index=False)

def compress_with_bgzip(input_tab_path):
    subprocess.run(['bgzip', '-f', input_tab_path], check=True)

if __name__ == "__main__":
    gsgt_path = '${txt}'
    tab_gz_path = '${gt_tab}'
    output_tab_path = '${meta.id}_final.tab'

    gsgt_dict = parse_gsgt_file(gsgt_path)
    annotate_tab_file(tab_gz_path, gsgt_dict, output_tab_path)
    compress_with_bgzip(output_tab_path)
CODE
    """
}
