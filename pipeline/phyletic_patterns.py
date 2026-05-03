import pandas as pd
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from sys import argv
import re

working_directory = argv[1]

def create_phyletic_patterns(in_file, phyletic_csv, phyletic_text, gene_order, start_col=4, delimiter=',', T4Es=False):
    df = pd.read_csv(in_file, sep=delimiter)
    header = df.columns
    species = header[start_col:]

    df.set_index(header[0], inplace=True)
    df.dropna(how='all', subset=species, inplace=True)
    df.fillna('0', inplace=True)
    df[species] = np.where(df[species] != '0', '1', '0')

    if T4Es:
        predicted = df[(df['is_effector'] == '?') & (df['score'] > 0.5)][species]
        positives = df[df.is_effector == 'yes'][species]
        positives, predicted = positives.T, predicted.T
        phyletic_df = positives.merge(predicted, left_index=True, right_index=True)
    else:
        phyletic_df = df.T

    phyletic_df.to_csv(phyletic_csv)

    phyletic_sequence = phyletic_df.apply(lambda row: ''.join(row.values.astype(str)), axis=1)

    with open(phyletic_text, 'w') as out_phyletic:
        for i in range(len(phyletic_sequence)):
            out_phyletic.write(f'>{phyletic_sequence.index[i]}\n')
            out_phyletic.write(f'{phyletic_sequence.iloc[i]}\n')

    with open(gene_order, 'w') as out_gene_order:
        out_gene_order.write(','.join(phyletic_df.columns))

def simplify_T4SS_component(annotation):
    if pd.isna(annotation):
        return annotation

    annotation = str(annotation)

    component_match = re.search(
        r'((?:Icm[A-Za-z0-9]+|Dot[A-Za-z0-9]+|LvgA)(?:/(?:Icm[A-Za-z0-9]+|Dot[A-Za-z0-9]+|LvgA))*)',
        annotation
    )

    if component_match:
        return component_match.group(1)

    return annotation

def create_T4SS_phyletic_patterns(in_file, phyletic_csv, phyletic_text, gene_order):
    df = pd.read_csv(in_file)

    df = df.dropna(subset=['Bacterial Protein ID'])
    df = df[df['Bacterial Protein ID'].astype(str).str.strip() != '']

    df['T4SS component'] = df['T4SS Hit Annotation'].apply(simplify_T4SS_component)
    df['present'] = 1

    phyletic_df = df.pivot_table(
        index='Genome',
        columns='T4SS component',
        values='present',
        aggfunc='max',
        fill_value=0
    )

    phyletic_df.to_csv(phyletic_csv)

    phyletic_sequence = phyletic_df.apply(lambda row: ''.join(row.values.astype(str)), axis=1)

    with open(phyletic_text, 'w') as out_phyletic:
        for i in range(len(phyletic_sequence)):
            out_phyletic.write(f'>{phyletic_sequence.index[i]}\n')
            out_phyletic.write(f'{phyletic_sequence.iloc[i]}\n')

    with open(gene_order, 'w') as out_gene_order:
        out_gene_order.write(','.join(phyletic_df.columns))

def create_presence_absence_map(
        phyletic_csv,
        x_label,
        out_fig_path,
        colors=['silver', 'lightseagreen'],
        base_font_size=9,
        col_cluster=True,
        width_per_col=0.08,
        height_per_row=0.35,
        min_width=12,
        max_width=30,
        min_height=6,
        max_height=16):

    df = pd.read_csv(phyletic_csv)
    df.set_index(df.columns[0], inplace=True)

    n_rows, n_cols = df.shape
    fig_width = max(n_cols * width_per_col, min_width)
    fig_height = max(n_rows * height_per_row, min_height)

    if max_width is not None:
        fig_width = min(fig_width, max_width)
    if max_height is not None:
        fig_height = min(fig_height, max_height)

    cmap = LinearSegmentedColormap.from_list('Custom', colors, len(colors))

    clustermap = sns.clustermap(
        df,
        cmap=cmap,
        figsize=(fig_width, fig_height),
        tree_kws={"linewidths": 0.},
        col_cluster=col_cluster,
        row_cluster=False,
        vmin=0,
        vmax=1,
        cbar_pos=(0.03, 0.75, 0.025, 0.12)
    )

    clustermap.ax_heatmap.xaxis.set_label_position('top')
    clustermap.ax_heatmap.xaxis.tick_top()

    plt.setp(
        clustermap.ax_heatmap.xaxis.get_majorticklabels(),
        rotation=60,
        ha='left',
        fontsize=base_font_size
    )

    plt.setp(
        clustermap.ax_heatmap.yaxis.get_majorticklabels(),
        fontsize=base_font_size
    )

    clustermap.ax_heatmap.set_xlabel(x_label, fontsize=base_font_size + 2, labelpad=12)
    clustermap.ax_heatmap.set_ylabel('Genome', fontsize=base_font_size + 2, labelpad=10)

    colorbar = clustermap.ax_heatmap.collections[0].colorbar
    colorbar.set_ticks([0.25, 0.75])
    colorbar.set_ticklabels(['Absent', 'Present'])
    colorbar.ax.tick_params(labelsize=base_font_size)

    clustermap.fig.savefig(
        out_fig_path,
        dpi=180,
        bbox_inches='tight',
        pad_inches=0.15
    )

    plt.close(clustermap.fig)

# T4Es

in_f_T4Es = os.path.join(
    working_directory,
    'out_learning',
    'consensus_predictions_with_annotations_and_ortho_table.csv'
)

phyletic_csv_T4Es = os.path.join(working_directory, 'PresenceAbsence_T4Es.csv')
phyletic_text_T4Es = os.path.join(working_directory, 'Effectors_phyletic_pattern.txt')
gene_order_T4Es = os.path.join(working_directory, 'T4Es_order_in_PhyleticPattern.txt')

create_phyletic_patterns(
    in_f_T4Es,
    phyletic_csv_T4Es,
    phyletic_text_T4Es,
    gene_order_T4Es,
    T4Es=True
)

create_presence_absence_map(
    phyletic_csv_T4Es,
    'T4E (ortholog group)',
    os.path.join(working_directory, 'T4Es_presence_absence.png'),
    base_font_size=8,
    width_per_col=0.12,
    height_per_row=0.55,
    min_width=18,
    max_width=45,
    min_height=10,
    max_height=28
)

# T4SS

in_file_T4SS = os.path.join(working_directory, 'T4SS.csv')
phyletic_csv_T4SS = os.path.join(working_directory, 'PresenceAbsence_T4SS.csv')
phyletic_text_T4SS = os.path.join(working_directory, 'T4SS_phyletic_pattern.txt')
gene_order_T4SS = os.path.join(working_directory, 'T4SS_order_in_PhyleticPattern.txt')

create_T4SS_phyletic_patterns(
    in_file_T4SS,
    phyletic_csv_T4SS,
    phyletic_text_T4SS,
    gene_order_T4SS
)

create_presence_absence_map(
    phyletic_csv_T4SS,
    'T4SS component',
    os.path.join(working_directory, 'T4SS_presence_absence.png'),
    col_cluster=False
)

# chaperones

in_f_chaperones = os.path.join(working_directory, 'chaperones.csv')

if os.path.exists(in_f_chaperones):
    df = pd.read_csv(in_f_chaperones)
    n_rows, n_cols = df.shape

    if n_rows > 1 and n_cols > 1:
        phyletic_csv_chaperones = os.path.join(working_directory, 'PresenceAbsence_chaperones.csv')
        phyletic_text_chaperones = os.path.join(working_directory, 'chaperones_phyletic_pattern.txt')
        gene_order_chaperones = os.path.join(working_directory, 'chaperones_order_in_PhyleticPattern.txt')

        create_phyletic_patterns(
            in_f_chaperones,
            phyletic_csv_chaperones,
            phyletic_text_chaperones,
            gene_order_chaperones,
            start_col=1
        )

        create_presence_absence_map(
            phyletic_csv_chaperones,
            'chaperone',
            os.path.join(working_directory, 'chaperones_presence_absence.png')
        )