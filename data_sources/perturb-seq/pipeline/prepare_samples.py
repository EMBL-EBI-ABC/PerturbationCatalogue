import pandas as pd
import os
import sys
import re
import json

def main(metadata_path, fastq_dir, output_csv):
    df = pd.read_csv(metadata_path, sep='\t')
    
    # Extract Sample ID and Modality from library_name
    # Pattern: jurkat_<modality>_<sample_id>_L<lane>
    # Example: jurkat_mRNA_8_1_L003 -> modality=mRNA, sample_id=8_1
    def parse_library(name):
        match = re.search(r'jurkat_(mRNA|sgRNA)_([\d_]+)_L\d+', name)
        if match:
            return match.group(1), match.group(2)
        return None, None

    df[['modality', 'sample_id']] = df['library_name'].apply(lambda x: pd.Series(parse_library(x)))
    df = df.dropna(subset=['sample_id'])

    # Find present files
    files = os.listdir(fastq_dir)
    present_srrs = set()
    for f in files:
        if f.endswith('.fastq.gz'):
            srr = f.split('_')[0]
            present_srrs.add(srr)

    df = df[df['run_accession'].isin(present_srrs)]

    # Group by sample_id
    samples = {}
    for _, row in df.iterrows():
        sid = row['sample_id']
        mod = row['modality']
        srr = row['run_accession']
        
        if sid not in samples:
            samples[sid] = {'mRNA': [], 'sgRNA': []}
        
        samples[sid][mod].append(srr)

    # Convert to CSV for Nextflow
    # Format: sample_id,mRNA_srrs,sgRNA_srrs
    # mRNA_srrs is semicolon-separated list of SRRs
    rows = []
    for sid, mods in samples.items():
        if mods['mRNA'] or mods['sgRNA']:
            rows.append({
                'sample_id': sid,
                'mRNA_srrs': ';'.join(sorted(set(mods['mRNA']))),
                'sgRNA_srrs': ';'.join(sorted(set(mods['sgRNA'])))
            })

    pd.DataFrame(rows).to_csv(output_csv, index=False)
    print(f"Generated sample sheet with {len(rows)} samples to {output_csv}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python prepare_samples.py <metadata_tsv> <fastq_dir> <output_csv>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
