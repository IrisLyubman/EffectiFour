from sys import argv
import os
import csv
import pandas as pd
import math
import numpy as np
from functools import reduce

working_directory = argv[1]
os.chdir(working_directory)

Effectidor_features_d = os.path.join(working_directory, 'Effectidor_runs')


def get_ortho_dict(wd, ortho_f='clean_orthologs_table.csv'):
    '''from the output file of Microbializer, present in the working directory,
    this function creates a dictionary that maps in every genome the loci to OGs'''
    genomes_orthogroup_dict = {}
    with open(os.path.join(wd, ortho_f)) as ortho_table:
        ortho_reader = csv.reader(ortho_table)
        header = next(ortho_reader)
        for row in ortho_reader:
            OG = row[0]
            for i in range(1, len(header)):
                genome = header[i]
                if genome not in genomes_orthogroup_dict:
                    genomes_orthogroup_dict[genome] = {}
                for locus in row[i].split(';'):
                    genomes_orthogroup_dict[genome][locus] = OG
    return genomes_orthogroup_dict


genomes_orthogroup_dict = get_ortho_dict(working_directory)
flatten_ortho_dict = {key: genomes_orthogroup_dict[dic_name][key] for dic_name in genomes_orthogroup_dict
                      for key in genomes_orthogroup_dict[dic_name]}


# prepare the orthologs' table for final merge and remove from it the T4SS components:
if os.path.exists(Effectidor_features_d):
    pseudogenes = []
    genomes = os.listdir(Effectidor_features_d)
    for genome in genomes:
        pseudogenes_path = os.path.join(Effectidor_features_d, genome, 'pseudogenes.txt')
        with open(pseudogenes_path) as pseudo_f:
            content = pseudo_f.read()
            if content != '':
                pseudogenes.extend(content.split('\n'))
    T4SS_tables = []
    for genome in genomes:
        t4ss_path = os.path.join(Effectidor_features_d, genome, 'T4SS.csv')
        df = pd.read_csv(t4ss_path, dtype={'Bacterial Protein ID': str})
        df.insert(0, 'Genome', genome)
        T4SS_tables.append(df)

    merged_T4SS_df = pd.concat(T4SS_tables, ignore_index=True)
    merged_T4SS_df.to_csv('T4SS.csv', index=False)

else:
    with open('pseudogenes.txt') as pseudo_f:
        pseudogenes = pseudo_f.read().split('\n')
    merged_T4SS_df = pd.read_csv('T4SS.csv', index_col='T4SS Protein', dtype={'Bacterial Protein ID': str})
with open(os.path.join(working_directory, 'pseudogenes.txt'), 'w') as pseudo_f:
    pseudo_f.write('\n'.join(pseudogenes))

T4SS_loci_to_remove = []
if os.path.exists(Effectidor_features_d):
    for genome in genomes:
        t4ss_path = os.path.join(Effectidor_features_d, genome, 'T4SS.csv')
        df = pd.read_csv(t4ss_path, dtype={'Bacterial Protein ID': str})
        if 'Bacterial Protein ID' in df.columns:
            ids = df['Bacterial Protein ID'].dropna()
            ids = ids[ids.str.strip() != '']
            T4SS_loci_to_remove.extend(ids.tolist())
else:
    df = pd.read_csv('T4SS.csv', dtype={'Bacterial Protein ID': str})
    if 'Bacterial Protein ID' in df.columns:
        ids = df['Bacterial Protein ID'].dropna()
        ids = ids[ids.str.strip() != '']
        T4SS_loci_to_remove.extend(ids.tolist())
        
if len(T4SS_loci_to_remove) > 0:
    flattened_T4SS_loci_no_nan = T4SS_loci_to_remove
    T4SS_OGs_to_remove = set([flatten_ortho_dict[locus] for locus in flattened_T4SS_loci_no_nan])
    OGs_table = pd.read_csv('clean_orthologs_table.csv', index_col='OG_name', dtype=str)
    OGs_no_T4SS = OGs_table.drop(index=list(T4SS_OGs_to_remove), errors="ignore")
else:
    OGs_no_T4SS = pd.read_csv('clean_orthologs_table.csv', index_col='OG_name', dtype=str)
OGs_no_T4SS.to_csv(os.path.join(working_directory, 'clean_orthologs_table_noT4SS.csv'))

genomes_orthogroup_dict = get_ortho_dict(working_directory, ortho_f='clean_orthologs_table_noT4SS.csv')

with open(os.path.join(working_directory, 'clean_orthologs_table.csv')) as in_f:
    with open('clean_orthologs_table_with_pseudo.csv', 'w') as out_f:
        header = next(in_f)
        header = header.replace('OG_name,', 'OG,')
        out_f.write(header)
        for row in in_f:
            to_remove = []
            for pseudo in pseudogenes:
                if f'{pseudo},' in row or f'{pseudo};' in row or f'{pseudo}\n' in row or f'{pseudo}\r\n' in row:
                    row = row.replace(f'{pseudo},', f'{pseudo}(pseudogene),')
                    row = row.replace(f'{pseudo};', f'{pseudo}(pseudogene);')
                    row = row.replace(f'{pseudo}\n', f'{pseudo}(pseudogene)\n')
                    row = row.replace(f'{pseudo}\r\n', f'{pseudo}(pseudogene)\r\n')
                    to_remove.append(pseudo)
            out_f.write(row)
            for pseudo in to_remove:
                pseudogenes.remove(pseudo)


def combine_all_genomes_data(out_f_path, in_f_name):
    with open(out_f_path, 'w', newline='') as out_path:
        writer = csv.writer(out_path)
        if os.path.exists(Effectidor_features_d):
            genomes = os.listdir(Effectidor_features_d)
            first_genome = genomes[0]
            with open(os.path.join(Effectidor_features_d, first_genome, f'{in_f_name}.csv')) as in_f:
                reader = csv.reader(in_f)
                header = ['OG'] + next(reader)
                writer.writerow(header)
                for row in reader:
                    if row[0] in genomes_orthogroup_dict[first_genome]:
                        row = [genomes_orthogroup_dict[first_genome][row[0]]] + row
                        writer.writerow(row)
            for genome in genomes[1:]:
                with open(os.path.join(Effectidor_features_d, genome, f'{in_f_name}.csv')) as in_f:
                    reader = csv.reader(in_f)
                    next(reader)
                    for row in reader:
                        if row[0] in genomes_orthogroup_dict[genome]:
                            row = [genomes_orthogroup_dict[genome][row[0]]] + row
                            writer.writerow(row)
        else:
            with open(f'{in_f_name}.csv') as in_f:
                reader = csv.reader(in_f)
                header = ['OG'] + next(reader)
                writer.writerow(header)
                for row in reader:
                    if row[0] in genomes_orthogroup_dict['genome_ORFs']:
                        row = [genomes_orthogroup_dict['genome_ORFs'][row[0]]] + row
                        writer.writerow(row)


combine_all_genomes_data('full_data.csv', 'features')
combine_all_genomes_data('full_OGs_annotations.csv', 'annotations')
combine_all_genomes_data('full_effector_homologs.csv', 'closest_effector_homologs')


def groupbyMode(in_f, out_f, target, new_col_name):
    df = pd.read_csv(in_f)
    try:
        g_df = df.groupby(['OG'])[target].agg(pd.Series.mode)
    except ValueError:
        g_df = df.groupby(['OG'], dropna=False)[target].apply(lambda x:None)
    g_df = g_df.apply(lambda x: ', '.join(x) if type(x) == np.ndarray else x)
    g_df.rename(new_col_name, inplace=True)
    g_df.to_csv(out_f)


groupbyMode('full_OGs_annotations.csv', 'OGs_annotations.csv', 'annotation', 'Annotation(s)')
groupbyMode('full_effector_homologs.csv', 'OG_effector_homologs.csv', 'Effector_ID', 'Effector_homolog(s)')


def label(iterable_arg):
    '''define a label of an OG, based on the labels of its members'''
    if any(iterable_arg == 'effector'):
        return 'effector'
    elif any(iterable_arg == 'no'):
        return 'no'
    else:
        return '?'   
        

df = pd.read_csv('full_data.csv')
features = list(df.columns[2:-1])
# defining manipulation per feature in the transformation to OGs
median = []
if 'distance_from_closest_effector' in features:
    median.append('distance_from_closest_effector')
if 'distance_to_mobile_genetic_element' in features:
    median.append('distance_to_mobile_genetic_element')
maximum_features = ['T4_signal', 'PmrA', 'cpxr', 'csra']
maximum = [f for f in features if any(k in f for k in maximum_features)]
for feature in maximum:
    features.remove(feature)
for feature in median:
    features.remove(feature)

#%% transformation to OGs
df.sort_values(by=features, inplace=True)
grouped = df.groupby('OG').agg({**{feature: ['mean'] for feature in features},
                                **{feature: ['median'] for feature in median},
                                **{feature: ['max'] for feature in maximum},
                                **{'is_effector': label}})

grouped.sort_values(by=list(grouped.columns[:-1]), inplace=True)
#%% adding the "similarity_to_effectors_vs_non_effectors" features for the OGs
aa_freqs = [feature for feature in features if 'full_protein' in feature]
effectors = grouped[grouped.is_effector.label == 'effector'][aa_freqs]
non_effectors = grouped[grouped.is_effector.label == 'no'][aa_freqs]
effectors_mean = effectors.mean()
non_effectors_mean = non_effectors.mean()


def aa_profile(OG):
    OG_aa_freqs = grouped.loc[OG][aa_freqs]
    if OG in effectors.index:
        effector_1 = effectors.drop(index=OG)
        effector_1_mean = effector_1.mean()
        dis_to_eff = math.sqrt(sum([(a-b)**2 for a, b in zip(OG_aa_freqs.values, effector_1_mean.values)]))
    else:
        dis_to_eff = math.sqrt(sum([(a-b)**2 for a, b in zip(OG_aa_freqs.values, effectors_mean.values)]))
    if OG in non_effectors.index:
        non_effector_1 = non_effectors.drop(index=OG)
        non_effector_1_mean = non_effector_1.mean()
        dis_to_non_eff = math.sqrt(sum([(a-b)**2 for a, b in zip(OG_aa_freqs.values, non_effector_1_mean.values)]))
    else:
        dis_to_non_eff = math.sqrt(sum([(a-b)**2 for a, b in zip(OG_aa_freqs.values, non_effectors_mean.values)]))
    
    return dis_to_eff - dis_to_non_eff


OGs = list(grouped.index)
aa_profiles_dis = [aa_profile(OG) for OG in OGs]
grouped.insert(len(grouped.columns)-1, 'similarity_to_effectors_vs_non_effectors', aa_profiles_dis)
 
#%% flatten columns and indices to a simple features table
updated_features = grouped.reset_index()
updated_features.columns = ['_'.join(col) for col in updated_features.columns.values]
updated_features.columns = [col.replace('label', '_').strip('_') for col in updated_features.columns]
updated_features.sort_values(by=list(updated_features.columns[1:]), inplace=True)
updated_features.to_csv('OGs_features.csv', index=False)

f = open('merge_OGs.done', 'w')
f.close()