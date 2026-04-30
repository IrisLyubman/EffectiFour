from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import subprocess
import math
import re
import csv
import os
import fasta_parser
from itertools import combinations


#%%

def parse_gff(gff_f, locus_dic):
    CDS_l = []
    RNA = []
    with open(gff_f) as in_f:
        for line in in_f:
            if line.startswith('#'): #header
                continue
            
            line_l = line.strip().split('\t')
            if len(line_l) > 2:
                if line_l[2] == 'CDS':
                    is_locus_tag = False
                    features_l = line_l[-1].split(';')
                    for feature in features_l:
                        if feature.startswith('locus'):
                            locus_tag = feature.split('=')[1]
                            CDS_l.append(locus_tag)
                            is_locus_tag = True
                            break
                    if not is_locus_tag:
                        for feature in features_l:
                            if feature.startswith('ID=CDS:'):
                                locus_tag = feature.split(':')[1]
                                CDS_l.append(locus_tag)
                                is_locus_tag = True
                                break
                    if not is_locus_tag:
                        for locus in locus_dic:
                            if f'{locus};' in line_l[-1]:
                                locus_tag = locus
                                CDS_l.append(locus_tag)
                                break
                else:
                    features_l = line_l[-1].split(';')
                    is_locus_tag = False
                    for feature in features_l:
                        if feature.startswith('locus'):
                            locus_tag = feature.split('=')[1]
                            RNA.append(locus_tag)
                            is_locus_tag = True
                            break
                    if not is_locus_tag:
                        for feature in features_l:
                            if feature.startswith('ID=transcript:'):
                                locus_tag = feature.split(':')[1]
                                RNA.append(locus_tag)
                                is_locus_tag = True
                                break
                    if not is_locus_tag:
                        for locus in locus_dic:
                            if f'{locus};' in line_l[-1]:
                                locus_tag = locus
                                RNA.append(locus_tag)
                                break
    RNA_set = set(RNA).difference(set(CDS_l))
    return set(CDS_l), RNA_set

def parse_gff_to_CDS_loc(gff_f, locus_dic):
    regions = []
    circulars = []
    linears = []
    locus_area_d = {}
    with open(gff_f) as in_f:
        for line in in_f:
            if line.startswith('#'): #header
                line_l = line.split()
                if 'sequence-region' in line_l[0]:
                    regions.append(line_l[1])
                    locus_area_d[line_l[1]] = {}            
            line_l = line.strip().split('\t')
            if len(line_l) > 2:
                if line_l[2] == 'region' or line_l[2] == "chromosome" or line_l[2] == "plasmid":
                    if 'Is_circular=true' in line_l[-1]:
                        circulars.append(line_l[0])
                    else:
                        linears.append(line_l[0])
                if line_l[2] == 'CDS':
                    region = line_l[0]
                    area = [(int(line_l[3])-1,int(line_l[4]))]
                    if line_l[6] == '+':
                        area.append(1)
                    elif line_l[6] == '-':
                        area.append(-1)
                    features_l = line_l[-1].split(';')
                    is_locus_tag = False
                    for feature in features_l:
                        if feature.startswith('locus'):
                            is_locus_tag = True
                            locus_tag = feature.split('=')[1]
                            if region not in locus_area_d:
                                locus_area_d[region] = {}
                            locus_area_d[region][locus_tag] = area
                            break
                    if not is_locus_tag:
                        for feature in features_l:
                            if feature.startswith('ID=CDS:'):
                                locus_tag = feature.split(':')[1]
                                if region not in locus_area_d:
                                    locus_area_d[region] = {}
                                locus_area_d[region][locus_tag] = area
                                is_locus_tag = True
                                break
                    if not is_locus_tag:
                        for locus in locus_dic:
                            if f'{locus};' in line_l[-1]:
                                locus_tag = locus
                                if region not in locus_area_d:
                                    locus_area_d[region] = {}
                                locus_area_d[region][locus_tag]=area
                                break
    return locus_area_d,circulars

def get_promoters(gene_area_dic, circular_contigs, genome_f, promoter_length=150):
    promoters_d = {}
    genome = SeqIO.to_dict(SeqIO.parse(genome_f, 'fasta'))
    for region_dic in gene_area_dic:
        for gene in gene_area_dic[region_dic]:
            if gene_area_dic[region_dic][gene][1] == 1:  # forward
                start, end = gene_area_dic[region_dic][gene][0][0], gene_area_dic[region_dic][gene][0][1]
                if start >= promoter_length:
                    promoter = genome[region_dic].seq[start-promoter_length:start]
                else:
                    if region_dic in circular_contigs:
                        to_add = promoter_length-start
                        promoter = genome[region_dic].seq[-to_add:]+genome[region_dic].seq[:start]
                    else:
                        promoter = genome[region_dic].seq[:start]
            else: #-1 reverse
                start, end = gene_area_dic[region_dic][gene][0][0], gene_area_dic[region_dic][gene][0][1]
                if len(genome[region_dic].seq)-end >= promoter_length:
                    promoter = genome[region_dic].seq[end:end+promoter_length].reverse_complement()
                else:
                    if region_dic in circular_contigs:
                        to_add = promoter_length - (len(genome[region_dic].seq)-end)
                        promoter = genome[region_dic].seq[end:]+genome[region_dic].seq[:to_add]
                        promoter = promoter.reverse_complement()
                    else:
                        promoter = genome[region_dic].seq[end:].reverse_complement()
            promoters_d[gene] = promoter
    return promoters_d


def create_box_mismatches(box_list, num_mismatches=2):
    motifs = []
    for mismatch_positions in combinations(range(len(box_list)), num_mismatches):
        li = list(box_list)
        for pos in mismatch_positions:
            li[pos] = '[ATGC]'
        motifs.append(''.join(li))

    return '|'.join(motifs)

exs_box = '[AT]{3}[AC][AT]{2}[AC]{3}C[GT][GTA]CC[GA]A[AT][ATC][CT]TG[GA][TC]A'
exs_box_l = ['[AT]']*3+['[AC]']+['[AT]']*2+['[AC]']*3+['C', '[GT]', '[GTA]']+['C']*2+['[GA]', 'A', '[AT]', '[ATC]', '[CT]', 'T', 'G', '[GA]', '[TC]', 'A']
exs_box_mismatch = create_box_mismatches(exs_box_l)

#CpxR = cpxr, will change later
cpxr = 'GTAAA[ATGC]{6}G[AT]AAA'
# 'GTAAA[ATGC][AT]{2}[ATGC]{3}G[AT]AAA' less permissive CpxR
cpxr_l = ['G','T','A','A','A'] + ['[ATGC]']*6 + ['G','[AT]','A','A','A']
cpxr_mismatch = create_box_mismatches(cpxr_l, num_mismatches=1)

pmra = 'CTTAA[TG][AG][TC]T[ATGC]{2}[ACG]TTAA[TG][ATG][ATC][TA]'
pmra_l = (list("CTTAA") + ['[TG]', '[AG]', '[TC]'] + list("T") + ['[ATGC]']*2 + ['[ACG]'] + list("TTAA") + ['[TG]', '[ATG]', '[ATC]', '[TA]'])
pmra_mismatch = create_box_mismatches(pmra_l)

csra = 'A?[ATGC]GGA[ATGC]{1,30}A?[ATGC]GGA'
csra_l = (['A?'] + ['[ATGC]', 'G', 'G', 'A'] + ['[ATGC]{1,30}'] + ['A?'] + ['[ATGC]', 'G', 'G', 'A'])
csra_mismatch = create_box_mismatches(csra_l)  

tts_box = 'GTCAG[TCG]T[TCAG]{4}G[AT][AC]AG[CGT][TAC][ATCG]{3}[CTG]{2}[ATCG]{4}A'
tts_box_l = ['G', 'T', 'C', 'A', 'G', '[TGC]', 'T']+['[ATGC]']*4+['G', '[AT]', '[AC]', 'A', 'G', '[CGT]', '[TAC]']+['[ATGC]']*3+['[CTG]']*2+['[ATGC]']*4+['A']
tts_box_mismatch = create_box_mismatches(tts_box_l)


def existence_upstream_to_AUG(locus, pattern, promoter_dic,promoter_length=150):
    promoter = promoter_dic[locus][-promoter_length:]
    if re.search(pattern, str(promoter), re.I):
        return 1
    else:
        return 0
                
    
def main(ORFs_file, working_directory, gff_f, genome_f, PmrA=False, CpxR=False, CsrA=False, exs=False, tts=False):
    os.chdir(working_directory)
    locus_dic = fasta_parser.parse_ORFs(ORFs_file)
    locus_area_d, circulars = parse_gff_to_CDS_loc(gff_f, locus_dic)
    if CsrA:
        promoters_d = get_promoters(locus_area_d, circulars, genome_f, promoter_length=50)
    else:
        promoters_d = get_promoters(locus_area_d, circulars, genome_f)
        
    with open(f'regulatory_features.csv', 'w', newline='') as f:
        csv_writer = csv.writer(f)
        header = ['locus']
        if PmrA:
            header += ['PmrA']
        if CpxR:
            header += ['CpxR']
        if CsrA:
            header += ['CsrA']
        if exs:
            header += ['exs_box']
        if tts:
            header += ['tts_box']
        csv_writer.writerow(header)
        for locus in locus_dic:
            l = [locus]
            if PmrA:
                l.append(existence_upstream_to_AUG(locus, pmra_mismatch, promoters_d))
                #l.append(existence_upstream_to_AUG(locus,pmra_mismatch))
            if CpxR:
                l.append(existence_upstream_to_AUG(locus, cpxr_mismatch, promoters_d, promoter_length=100))
                #l.append(existence_upstream_to_AUG(locus,cpxr_mismatch))
            if CsrA:
                l.append(existence_upstream_to_AUG(locus, csra, promoters_d))
                #l.append(existence_upstream_to_AUG(locus,csra_mismatch))
            if exs:
                l.append(existence_upstream_to_AUG(locus, exs_box, promoters_d))
                #l.append(existence_upstream_to_AUG(locus,exs_box_mismatch))
            if tts:
                l.append(existence_upstream_to_AUG(locus, tts_box, promoters_d))
                #l.append(existence_upstream_to_AUG(locus,tts_box_mismatch))
                
            csv_writer.writerow(l)
    endfile = open('regulatory_features.done', 'w')
    endfile.close()

if __name__ == '__main__':

        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('input_ORFs_path',
                            help='A path to a DNA ORFs file.',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
        parser.add_argument('output_dir_path',
                            help='A path to a folder in which the output files will be created.',
                            type=lambda path: path.rstrip('/'))
        parser.add_argument('gff_path',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'),
                            help='A path to GFF file.')
        parser.add_argument('genome_path',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'),
                            help='A path to fasta file with full genome records.')
        parser.add_argument('--PmrA', help='look for PmrA in promoters', action='store_true')
        parser.add_argument('--CpxR', help='look for CpxR in promoters', action='store_true')
        parser.add_argument('--CsrA', help='look for CsrA in promoters', action='store_true')
        parser.add_argument('--exs', help='look for exs-box in promoters', action='store_true')
        parser.add_argument('--tts',help='look for tts-box in promoters', action='store_true')
        
        args = parser.parse_args()
        ORFs_file = args.input_ORFs_path
        working_directory = args.output_dir_path
        gff_f = args.gff_path
        genome_f = args.genome_path
        PmrA_flag = args.PmrA
        CpxR_flag = args.CpxR
        CsrA_flag = args.CsrA
        exs_flag = args.exs
        tts_flag = args.tts
        main(ORFs_file, working_directory, gff_f, genome_f, PmrA=PmrA_flag, CpxR=CpxR_flag, CsrA=CsrA_flag, exs=exs_flag,
             tts=tts_flag)
