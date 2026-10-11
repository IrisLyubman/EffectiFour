import fasta_parser
from sys import argv
import os
from Bio import SeqIO, SeqRecord
import subprocess
import csv
import shutil
import time
import sys
sys.path.append('/effectidor/type4b/script/auxiliaries')
from auxiliaries import fail
import effectidor_CONSTANTS as CONSTS
import logging
import warnings

logging.captureWarnings(True)

working_directory = os.path.abspath(argv[1])
identity_cutoff = argv[2]  # default: '50'
coverage_cutoff = argv[3]  # default: '60'
error_path = os.path.join(working_directory, 'error_OGs.txt')
scripts_dir = CONSTS.EFFECTIDOR_AUXILIARIES
os.chdir(working_directory)

# Prepare all genomes for Microbializer according to genome_sources.csv
genome_sources_f = os.path.join(working_directory, 'genome_sources.csv')
microbializer_input_dir = os.path.join(working_directory, 'genomes_for_Microbializer')

if not os.path.exists(genome_sources_f):
    fail(
        f'genome_sources.csv was not found in {working_directory}',
        error_path
    )

# Remove previous input directory so old genomes from a previous run
# cannot accidentally remain in the Microbializer input.
if os.path.exists(microbializer_input_dir):
    shutil.rmtree(microbializer_input_dir)

os.makedirs(microbializer_input_dir)

seen_gcfs = set()
seen_orfs_paths = set()

with open(genome_sources_f, newline='') as sources_f:
    reader = csv.DictReader(sources_f)

    for row in reader:
        genome = row['genome'].strip()
        gcf = row.get('gcf', '').strip().upper()
        orfs_path = row.get('orfs_path', '').strip()

        if not genome:
            continue

        if not orfs_path:
            fail(
                f'No ORFs path was provided for genome {genome}',
                error_path
            )

        # Support both absolute and relative paths
        if not os.path.isabs(orfs_path) and not os.path.exists(orfs_path):
            orfs_path = os.path.join(working_directory, orfs_path)
        orfs_path = os.path.abspath(orfs_path)

        if not os.path.exists(orfs_path):
            fail(f'ORFs file was not found for genome {genome}: {orfs_path}', error_path)

        # genome_sources.csvshould already be deduplicated in main_T4Es.py, but keep this protection here as well.
        if gcf:
            if gcf in seen_gcfs:
                continue
            seen_gcfs.add(gcf)
        else:
            # For genomes without a GCF accession, prevent the exact same
            # ORFs file from being added twice.
            if orfs_path in seen_orfs_paths:
                continue

        seen_orfs_paths.add(orfs_path)

        locus_dic = fasta_parser.parse_ORFs(orfs_path)

        recs = [
            SeqRecord.SeqRecord(
                locus_dic[locus],
                id=locus,
                description=''
            )
            for locus in locus_dic
        ]

        output_f = os.path.join(
            microbializer_input_dir,
            f'{genome}.fasta'
        )

        SeqIO.write(recs, output_f, 'fasta')


# Run Microbializer block Josef replace --contigs_dir with --input
cmds_file_content = f'''source {CONSTS.MICROBIALIZER_SCRIPTS}/miniconda3/etc/profile.d/conda.sh!@#\
conda activate {CONSTS.MICROBIALIZER_SCRIPTS}/miniconda3/envs/microbializer!@#\
export PATH=$CONDA_PREFIX/bin:$PATH!@#!@#\
python {CONSTS.MICROBIALIZER_SCRIPTS}/pipeline/main.py --input \
{working_directory}/genomes_for_Microbializer/ --output_dir output_OGs --add_orphan_genes_to_ogs true --only_calc_ogs true \
--inputs_fasta_type orfs --bypass_number_of_genomes_limit true --identity_cutoff {identity_cutoff} --coverage_cutoff \
{coverage_cutoff} --account_name pupko-users_v2 -q pupko-pool\tMicrobializer_for_Effectidor'''
with open('search_OGs.cmds', 'w') as cmds_f:
    cmds_f.write(cmds_file_content)
cmd = f'python {os.path.join(scripts_dir, "q_submitter_power.py")} search_OGs.cmds {working_directory} --memory 4'
subprocess.call(cmd, shell=True)

Microbializer_output_f = 'M1CR0B1AL1Z3R_output_OGs/05a_orthogroups/orthogroups.csv'
while not os.path.exists(Microbializer_output_f):
    if os.path.exists('M1CR0B1AL1Z3R_output_OGs/error.txt'):
        with open('M1CR0B1AL1Z3R_output_OGs/error.txt') as error:
            error_msg = error.read()
        fail(error_msg, error_path)
    # while the job hasn't finished
    time.sleep(60)

# move Microbializer output to the working directory for future use
shutil.copy(Microbializer_output_f, working_directory)

# subprocess.call("rm -r M1CR0B1AL1Z3R_output_OGs",shell=True)
# subprocess.call("rm -r output_OGs",shell=True)

with open('orthogroups.csv') as ortho_table:
    ortho_reader = csv.reader(ortho_table)
    header = next(ortho_reader)
    with open(os.path.join(working_directory, 'clean_orthologs_table.csv'), 'w', newline='') as out_F:
        writer = csv.writer(out_F)
        writer.writerow(header)
        for row in ortho_reader:
            for i in range(1, len(header)):
                row[i] = row[i].replace(f'{header[i]}:', '')
            writer.writerow(row)
