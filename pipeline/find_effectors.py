from Bio import SeqIO
from sys import argv
import os
import csv
from protein_mmseq import protein_mmseqs_all_vs_all

effectors = argv[1]
bacterial_proteome = argv[2]
effectors_prots = argv[3]
minimal_coverage = float(argv[4])

if not os.path.exists('blast_outputs'):
    os.makedirs('blast_outputs')

effector_annotations = {}
for rec in SeqIO.parse(effectors, 'fasta'):
    annotation = rec.description[len(rec.id):].strip()
    effector_annotations[rec.id] = annotation

in_f = open(effectors)
content = in_f.read()
replaced_pipes = content.replace('|', '^')  # MMSeq changes the IDs when pipes are available in the rec.id,
# so I replace it now with a character that is not expected to be found in the ID, to replace it back while parsing
# the MMSeq result.
in_f.close()
with open(effectors, 'w') as effectors_f:
    effectors_f.write(replaced_pipes)
effectors_recs = SeqIO.parse(effectors, 'fasta')
effectors_lens = {effector_rec.id: len(effector_rec.seq) for effector_rec in effectors_recs}


def parse_blast_out(blast_out, e_val=10 ** (-10), min_identity=0.7, min_coverage=50.0):
    blast_out_dic = {}
    with open(blast_out, 'r') as blast_out_f:
        best_hits = {}
        for row in blast_out_f:
            row = row.split('\t')
            prot_id = row[0]
            hit = row[1]
            e_value = float(row[-2])
            identity = float(row[2])
            coverage = (float(row[3]) / effectors_lens[prot_id]) * 100
            prot_id = prot_id.replace('^', '|')

            if e_value <= e_val and identity >= min_identity and coverage >= min_coverage:
                if prot_id not in blast_out_dic:
                    blast_out_dic[prot_id] = {hit}
                else:
                    blast_out_dic[prot_id].add(hit)

                score = float(row[-1]) # MMseqs bitscore
                if hit not in best_hits or score > best_hits[hit]['score']:
                    best_hits[hit] = {'score': score, 'effector': prot_id}

    with open('blast_outputs/effectors_hits.csv', 'w') as out_f:
        for hit, info in best_hits.items():
            out_f.write(f"{hit},{info['score']}\n")

    with open('blast_outputs/effectors_hits_with_annotations.csv', 'w', newline='') as out_f:
        writer = csv.writer(out_f)
        writer.writerow(['locus', 'database_effector', 'database_annotation', 'bitscore'])

        for hit, info in best_hits.items():
            effector_id = info['effector']
            writer.writerow([hit, effector_id, effector_annotations.get(effector_id, ''), info['score']])

    return blast_out_dic

if not os.path.exists('tmp_mmseq'):
    os.makedirs('tmp_mmseq')
protein_mmseqs_all_vs_all(effectors, bacterial_proteome, 'blast_outputs/effectorsDB.blast', 'tmp_mmseq')

k = 0.6
effectors_homologs = set()
homologs = list(parse_blast_out('blast_outputs/effectorsDB.blast', min_identity=k, min_coverage=minimal_coverage).values())
for h in homologs:
    effectors_homologs = set.union(effectors_homologs, h)

recs = SeqIO.parse(bacterial_proteome, 'fasta')
effectors_recs = []
for rec in recs:
    if rec.id in effectors_homologs:
        effectors_recs.append(rec)
SeqIO.write(effectors_recs, effectors_prots, 'fasta')
