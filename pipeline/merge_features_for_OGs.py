from sys import argv
import os
import csv
import pandas as pd
import math
import numpy as np
from functools import reduce

working_directory = os.path.abspath(argv[1])
os.chdir(working_directory)

genome_sources_f = os.path.join(working_directory, 'genome_sources.csv')

if not os.path.exists(genome_sources_f):
    raise FileNotFoundError(
        f'genome_sources.csv was not found in {working_directory}'
    )

genome_sources = pd.read_csv(genome_sources_f, dtype=str).fillna('')

genome_sources['requested'] = (
    genome_sources['requested']
    .str.strip()
    .str.lower()
    .eq('true')
)

requested_genomes = genome_sources.loc[genome_sources['requested'], 'genome'].tolist()


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
                    locus = locus.strip()
                    if locus:
                        genomes_orthogroup_dict[genome][locus] = OG
    return genomes_orthogroup_dict


genomes_orthogroup_dict = get_ortho_dict(working_directory)

# Collect pseudogenes and T4SS data from ALL genomes participating in the run.
# All genomes are used for determining which OGs are T4SS OGs.
# Only requested genomes are written to the user-facing T4SS.csv.

pseudogenes_by_genome = {}
requested_T4SS_tables = []
T4SS_OGs_to_remove = set()

for _, source in genome_sources.iterrows():

    genome = source['genome']
    source_dir = source['source_dir']
    requested = source['requested']

    # Pseudogenes
    pseudogenes_path = os.path.join(source_dir, 'pseudogenes.txt')

    if not os.path.exists(pseudogenes_path):
        raise FileNotFoundError(
            f'pseudogenes.txt was not found for genome {genome}: '
            f'{pseudogenes_path}'
        )

    with open(pseudogenes_path) as pseudo_f:
        pseudogenes = [
            p.strip()
            for p in pseudo_f.read().splitlines()
            if p.strip()
        ]

    pseudogenes_by_genome[genome] = set(pseudogenes)

    # T4SS
    t4ss_path = os.path.join(source_dir, 'T4SS.csv')

    if not os.path.exists(t4ss_path):
        raise FileNotFoundError(
            f'T4SS.csv was not found for genome {genome}: {t4ss_path}'
        )

    df = pd.read_csv(
        t4ss_path,
        dtype={'Bacterial Locus Tag': str}
    )

    # User-facing T4SS table should contain only requested genomes
    if requested:
        requested_df = df.copy()
        requested_df.insert(0, 'Genome', genome)
        requested_T4SS_tables.append(requested_df)

    # But T4SS OG removal is based on ALL genomes
    if 'Bacterial Locus Tag' in df.columns:

        ids = df['Bacterial Locus Tag'].dropna().astype(str)
        ids = ids[ids.str.strip() != '']

        for locus in ids:

            locus = locus.strip()

            OG = genomes_orthogroup_dict.get(
                genome, {}
            ).get(locus)

            if OG is not None:
                T4SS_OGs_to_remove.add(OG)


# Write T4SS.csv only for requested genomes
if requested_T4SS_tables:

    merged_T4SS_df = pd.concat(
        requested_T4SS_tables,
        ignore_index=True
    )

else:
    merged_T4SS_df = pd.DataFrame()

merged_T4SS_df.to_csv(
    os.path.join(working_directory, 'T4SS.csv'),
    index=False
)


# Write pseudogenes.txt only for requested genomes
requested_pseudogenes = []

for genome in requested_genomes:
    requested_pseudogenes.extend(
        sorted(pseudogenes_by_genome.get(genome, set()))
    )

with open(
    os.path.join(working_directory, 'pseudogenes.txt'),
    'w'
) as pseudo_f:
    pseudo_f.write('\n'.join(requested_pseudogenes))


# Remove OGs containing T4SS components
OGs_table = pd.read_csv(
    'clean_orthologs_table.csv',
    index_col='OG_name',
    dtype=str
)

OGs_no_T4SS = OGs_table.drop(
    index=list(T4SS_OGs_to_remove),
    errors='ignore'
)

OGs_no_T4SS.to_csv(
    os.path.join(
        working_directory,
        'clean_orthologs_table_noT4SS.csv'
    )
)

genomes_orthogroup_dict = get_ortho_dict(working_directory, ortho_f='clean_orthologs_table_noT4SS.csv')

# OGs that belong to at least one requested genome
requested_OGs = set()

for genome in requested_genomes:

    if genome not in genomes_orthogroup_dict:
        raise RuntimeError(
            f'Requested genome {genome} was not found in the '
            f'Microbializer orthogroups table.'
        )

    requested_OGs.update(genomes_orthogroup_dict[genome].values())

with open(os.path.join(working_directory, 'requested_OGs.txt'), 'w') as out_f:
    for OG in sorted(requested_OGs):
        out_f.write(f'{OG}\n')

ortho_df = pd.read_csv(os.path.join(working_directory, 'clean_orthologs_table.csv'), dtype=str).fillna('')
ortho_df = ortho_df[['OG_name'] + requested_genomes].copy()

for genome in requested_genomes:
    genome_pseudogenes = pseudogenes_by_genome.get(genome, set())

    for pseudo in genome_pseudogenes:
        ortho_df[genome] = ortho_df[genome].apply(
            lambda x: ';'.join([f'{locus}(pseudogene)' if locus == pseudo else locus for locus in x.split(';')])
        )

ortho_df.rename(columns={'OG_name': 'OG'}, inplace=True)
ortho_df.to_csv('clean_orthologs_table_with_pseudo.csv', index=False)


def combine_all_genomes_data(out_f_path, in_f_name, requested_only=False):
    sources = genome_sources[genome_sources['requested']] if requested_only else genome_sources
    dfs = []

    for _, source in sources.iterrows():
        genome = source['genome']
        source_dir = source['source_dir']
        input_f = os.path.join(source_dir, f'{in_f_name}.csv')

        if not os.path.exists(input_f):
            raise FileNotFoundError(f'{in_f_name}.csv was not found for genome {genome}: {input_f}')

        if genome not in genomes_orthogroup_dict:
            raise RuntimeError(f'Genome {genome} was not found in the orthogroups table.')

        df = pd.read_csv(input_f)
        locus_col = df.columns[0]
        df = df[df[locus_col].isin(genomes_orthogroup_dict[genome])].copy()
        df.insert(0, 'OG', df[locus_col].map(genomes_orthogroup_dict[genome]))
        dfs.append(df)

    combined_df = pd.concat(dfs, ignore_index=True, sort=False)

    if 'is_effector' in combined_df.columns:
        cols = [col for col in combined_df.columns if col != 'is_effector'] + ['is_effector']
        combined_df = combined_df[cols]

    combined_df.to_csv(out_f_path, index=False)


# Training features: ALL genomes
combine_all_genomes_data('full_data.csv', 'features', requested_only=False)

# User-facing annotation information: requested genomes only
combine_all_genomes_data('full_OGs_annotations.csv', 'annotations', requested_only=True)
combine_all_genomes_data('full_effector_homologs.csv', 'closest_effector_homologs', requested_only=True)


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
maximum_features = ['T4_signal', 'PmrA', 'CpxR', 'CsrA', 'gnlp_cosine_similarity']
maximum = [f for f in features if any(k in f for k in maximum_features)]
minimum_features = [] # comment if euclidian distance is not relevant or add 'gnlp_euclidean_distance' when want to use
minimum = [f for f in features if any(k in f for k in minimum_features)] # comment if euclidian distance is not relevant

for feature in maximum:
    features.remove(feature)
for feature in median:
    features.remove(feature)
for feature in minimum: # comment if euclidian distance is not relevant
    features.remove(feature)

#%% transformation to OGs
df.sort_values(by=features, inplace=True)
grouped = df.groupby('OG').agg({**{feature: ['mean'] for feature in features},
                                **{feature: ['median'] for feature in median},
                                **{feature: ['max'] for feature in maximum},
                                **{feature: ['min'] for feature in minimum},
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

# Complete feature table used for training - all pan-genome + requested OGs
updated_features.to_csv('OGs_features_training.csv', index=False)
# User-facing feature table - only OGs present in requested genomes
requested_features = updated_features[updated_features['OG'].isin(requested_OGs)].copy()
requested_features.to_csv('OGs_features.csv', index=False)
f = open('merge_OGs.done', 'w')
f.close()
