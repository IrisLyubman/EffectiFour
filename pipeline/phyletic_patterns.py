import pandas as pd
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm  # NEW: better for category colors
from matplotlib.patches import Patch  # NEW: for custom legend
from sys import argv
import re

working_directory = argv[1]
features_f = os.path.join(working_directory, 'features_filled_nan.csv')

def create_phyletic_patterns(in_file, phyletic_csv, phyletic_text, gene_order,
                             start_col=4, delimiter=',', T4Es=False):
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

        positives = positives.T
        predicted = predicted.T

        # safer than merge, keeps both positive and predicted columns
        phyletic_df = pd.concat([positives, predicted], axis=1)

    else:
        phyletic_df = df.T

    phyletic_df.to_csv(phyletic_csv)

    phyletic_sequence = phyletic_df.apply(
        lambda row: ''.join(row.values.astype(str)),
        axis=1
    )

    with open(phyletic_text, 'w') as out_phyletic:
        for i in range(len(phyletic_sequence)):
            out_phyletic.write(f'>{phyletic_sequence.index[i]}\n')
            out_phyletic.write(f'{phyletic_sequence.iloc[i]}\n')

    with open(gene_order, 'w') as out_gene_order:
        out_gene_order.write(','.join(phyletic_df.columns))


# NEW: adds regulatory-element categories to the T4E presence/absence table
def add_regulation_to_T4E_phyletic_csv(phyletic_csv, features_file):
    """
    Converts T4E heatmap values:
    0 = absent
    1 = present T4E
    2 = T4E + PmrA
    3 = T4E + CpxR
    4 = T4E + CsrA
    5 = T4E + multiple regulatory elements
    """

    phyletic_df = pd.read_csv(phyletic_csv)
    phyletic_df.set_index(phyletic_df.columns[0], inplace=True)

    features = pd.read_csv(features_file)
    id_col = features.columns[0]

    regulatory_elements = ['PmrA', 'CpxR', 'CsrA']

    reg_cols = {
        element: [
            col for col in features.columns
            if col == element or col.startswith(f'{element}_')
        ]
        for element in regulatory_elements
    }

    reg_by_og = {}

    for _, row in features.iterrows():
        og = row[id_col]
        present_regs = []

        for reg_name, cols in reg_cols.items():
            for col in cols:
                if pd.notna(row[col]) and float(row[col]) == 1:
                    present_regs.append(reg_name)
                    break

        reg_by_og[og] = present_regs

    for og in phyletic_df.columns:
        if og not in reg_by_og:
            continue

        regs = reg_by_og[og]

        if len(regs) == 0:
            continue
        elif len(regs) > 1:
            value = 5
        elif regs[0] == 'PmrA':
            value = 2
        elif regs[0] == 'CpxR':
            value = 3
        elif regs[0] == 'CsrA':
            value = 4
        else:
            continue

        phyletic_df.loc[phyletic_df[og].astype(int) == 1, og] = value

    phyletic_df.to_csv(phyletic_csv)


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

    df = df.dropna(subset=['Bacterial Locus Tag'])
    df = df[df['Bacterial Locus Tag'].astype(str).str.strip() != '']

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
        colors=['#E6E6E6', '#4C78A8'],
        legend_labels=None,  # NEW: optional category legend
        base_font_size=9,
        show_x_labels=True,
        width_per_col=0.35,
        height_per_row=0.22,
        min_width=14,
        max_width=40,
        min_height=8,
        max_height=28):

    df = pd.read_csv(phyletic_csv)
    df.set_index(df.columns[0], inplace=True)
    df = df.astype(int)

    n_rows, n_cols = df.shape

    fig_width = min(max(n_cols * width_per_col, min_width), max_width)
    fig_height = min(max(n_rows * height_per_row, min_height), max_height)

    # NEW: categorical color handling
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(
        boundaries=np.arange(len(colors) + 1) - 0.5,
        ncolors=len(colors)
    )

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    sns.heatmap(
        df,
        cmap=cmap,
        norm=norm,
        linewidths=0,
        cbar=(legend_labels is None),
        cbar_kws={
            "shrink": 0.15,
            "aspect": 10,
            "pad": 0.02
        },
        ax=ax
    )

    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()

    ax.set_xlabel(
        x_label,
        fontsize=base_font_size + 2,
        labelpad=10
    )

    ax.set_ylabel(
        'Genome',
        fontsize=base_font_size + 2
    )

    if show_x_labels:
        ax.set_xticks(np.arange(len(df.columns)) + 0.5)
        ax.set_xticklabels(
            df.columns,
            rotation=90,
            fontsize=base_font_size
        )
    else:
        ax.set_xticks([])
        ax.set_xticklabels([])

    ax.set_yticks(np.arange(len(df.index)) + 0.5)
    ax.set_yticklabels(
        df.index,
        rotation=0,
        fontsize=max(4, base_font_size - 2)
    )

    # NEW: custom legend for T4E regulation categories
    if legend_labels is not None:
        handles = [
            Patch(facecolor=colors[i], edgecolor='none', label=legend_labels[i])
            for i in range(len(legend_labels))
        ]

        ax.legend(
            handles=handles,
            title='Category',
            bbox_to_anchor=(1.01, 1),
            loc='upper left',
            borderaxespad=0,
            fontsize=max(6, base_font_size - 1),
            title_fontsize=base_font_size
        )

    else:
        colorbar = ax.collections[0].colorbar
        colorbar.set_ticks([0, 1])
        colorbar.set_ticklabels(['Absent', 'Present'])

    plt.tight_layout()

    plt.savefig(
        out_fig_path,
        dpi=300,
        bbox_inches='tight',
        pad_inches=0.05
    )

    plt.close()


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

# NEW: add PmrA/CpxR/CsrA regulation information into the T4E heatmap values
add_regulation_to_T4E_phyletic_csv(
    phyletic_csv_T4Es,
    features_f
)

create_presence_absence_map(
    phyletic_csv_T4Es,
    'T4E (ortholog group)',
    os.path.join(working_directory, 'T4Es_presence_absence.png'),
    colors=[
        '#F2F2F2',  # absent
        '#4C78A8',  # present T4E only
        '#F58518',  # PmrA
        '#54A24B',  # CpxR
        '#B279A2',  # CsrA
        '#E45756'   # multiple regulatory elements
    ],
    legend_labels=[
        'Absent',
        'T4E',
        'T4E + PmrA',
        'T4E + CpxR',
        'T4E + CsrA',
        'T4E + multiple'
    ],
    show_x_labels=False,
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
