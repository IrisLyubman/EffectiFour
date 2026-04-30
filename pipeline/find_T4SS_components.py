import argparse
import os
import subprocess
import pandas as pd
from Bio import SeqIO

RUN_WITH_CONDA = False
E_VALUE_CUT_OFF = 1e-10
QUERY_COVERAGE_PERCENTAGE_CUT_OFF = 0.3
MMSEQS_OUTPUT_FORMAT = 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,qlen,qcov,tstart,tend,evalue,bits'

COLUMN_NAMES = [
    'T4SS_protein', 'bacterial_protein', 'identity_percent', 'alignment_length',
    'mismatch', 'gapopen', 'query_start', 'query_end', 'query_length',
    'query_coverage_percentage', 'target_start', 'target_end', 'e_value', 'bit_score'
]

HEADERS = HEADERS = ["T4SS Protein", "T4SS Hit Annotation", "Bacterial Protein ID"]


def run_mmseqs(query_fasta, output_mmseqs, bacterial_proteome, tmp_mmseqs, threads=1):
    mmseqs_command = (
        f"mmseqs easy-search {query_fasta} {bacterial_proteome} {output_mmseqs} "
        f"{tmp_mmseqs} --format-output {MMSEQS_OUTPUT_FORMAT} --threads {threads}"
    )

    if RUN_WITH_CONDA:
        conda_activate_command = ". ~/miniconda3/etc/profile.d/conda.sh; conda activate test; "
        mmseqs_command = conda_activate_command + mmseqs_command

    subprocess.run(mmseqs_command, shell=True, check=True)


def get_best_hits_dictionary(mmseqs_results_file):
    """
    Returns:
        {T4SS_protein: (bacterial_protein, bit_score)}
    """
    best_hits_dict = {}

    if not os.path.exists(mmseqs_results_file) or os.path.getsize(mmseqs_results_file) == 0:
        return best_hits_dict

    df = pd.read_csv(
        mmseqs_results_file,
        sep='\t',
        names=COLUMN_NAMES,
        header=None,
        dtype={'bacterial_protein': str}
    )

    filtered_df = df[
        (df['e_value'] < E_VALUE_CUT_OFF) &
        (df['query_coverage_percentage'] > QUERY_COVERAGE_PERCENTAGE_CUT_OFF)
    ]

    for _, row in filtered_df.iterrows():
        T4SS_protein = row['T4SS_protein']
        bacterial_protein = row['bacterial_protein']
        bit_score = row['bit_score']

        if T4SS_protein not in best_hits_dict:
            best_hits_dict[T4SS_protein] = (bacterial_protein, bit_score)
        elif bit_score > best_hits_dict[T4SS_protein][1]:
            best_hits_dict[T4SS_protein] = (bacterial_protein, bit_score)

    return best_hits_dict


def build_full_results_dict(T4SS_fasta, best_hits_dict):
    """
    Returns:
        {T4SS_protein: (annotation, bacterial_protein_or_None)}
    """
    full_results_dict = {}

    for rec in SeqIO.parse(T4SS_fasta, "fasta"):
        protein_id = rec.id

        # remove the accession from the full FASTA header and "MULTISPECIES:" if existst
        annotation = rec.description.replace(protein_id, "", 1).strip().replace("MULTISPECIES:", "").strip()

        if "[" in annotation and "]" in annotation:
            annotation = annotation.split("[")[0].strip()

        if protein_id in best_hits_dict:
            bacterial_protein = best_hits_dict[protein_id][0]
        else:
            bacterial_protein = None

        full_results_dict[protein_id] = (annotation, bacterial_protein)

    return full_results_dict


def write_dict_to_output_file(full_results_dict, output_file):
    data = []

    for T4SS_protein, (annotation, bacterial_protein) in full_results_dict.items():
        data.append([T4SS_protein, annotation, bacterial_protein])

    df = pd.DataFrame(data, columns=HEADERS)
    df.to_csv(output_file, index=False)


def main():
    parser = argparse.ArgumentParser(
        description='Run mmseqs for a Type IV protein set against a bacterial proteome and output a CSV of best hits.'
    )
    parser.add_argument('working_directory', type=str, help='Path to the working directory')
    parser.add_argument('bacterial_proteome', type=str, help='Name of the bacterial proteome file')
    parser.add_argument('T4SS_data', type=str, help='Directory containing T4SS fasta')
    parser.add_argument('--threads', type=int, default=1, help='Number of mmseqs threads')
    args = parser.parse_args()

    working_directory = args.working_directory
    bacterial_proteome = os.path.join(working_directory, args.bacterial_proteome)
    T4SS_dir = os.path.join(working_directory, "T4SS_data")
    T4SS_fasta = os.path.join(T4SS_dir, "T4SS.faa")

    tmp_mmseqs = os.path.join(working_directory, "tmp_mmseqs")
    output_mmseqs = os.path.join(working_directory, "output_mmseqs")
    output_file = os.path.join(working_directory, "T4SS.csv")

    os.makedirs(tmp_mmseqs, exist_ok=True)

    run_mmseqs(T4SS_fasta, output_mmseqs, bacterial_proteome, tmp_mmseqs, threads=args.threads)
    best_hits_dict = get_best_hits_dictionary(output_mmseqs)
    full_results_dict = build_full_results_dict(T4SS_fasta, best_hits_dict)
    write_dict_to_output_file(full_results_dict, output_file)


if __name__ == "__main__":
    main()