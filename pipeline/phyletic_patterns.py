import pandas as pd
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from sys import argv

working_directory = argv[1]


def create_phyletic_patterns(in_file, phyletic_csv, phyletic_text, gene_order, start_col=4, delimiter=',',
                             T4Es=False):
    '''
    in_file: path to csv file with the data to transform to phyletic pattern, where the columns are genomes, and the
    rows are genes/OGs
    phyletic_csv: path to a csv file that will be created with the phyletic pattern. it will be used later for the
    presence/absense map
    phyletic_text: a path to the output phyletic pattern  that will be created, in txt format
    gene_order: a path to an output txt file that will hold a list of the gene order in the phyletic pattern
    start_col: position of the first column in the input file that represents a genome
    delimiter: delimiter to use to read the dataframe (, for csv and \t for tsv)
    T4Es: whether the input is T4Es prediction (this input requires dropping columns that do not exist elsewhere)
    '''
    df = pd.read_csv(in_file, sep=delimiter)
    header = df.columns
    species = header[start_col:]
    df.set_index(header[0], inplace=True)
    df.dropna(how='all', subset=species, inplace=True)
    df.fillna('0', inplace=True)
    df[species] = np.where(df[species] != '0', '1', '0')  # put 1 wherever there's a locus
    if T4Es:
        predicted = df[(df['is_effector'] == '?') & (df['score'] > 0.5)][species]
        positives = df[df.is_effector == 'yes'][species]
        positives, predicted = positives.T, predicted.T
        phyletic_df = positives.merge(predicted, left_index=True, right_index=True)
    else:  # T4SS
        phyletic_df = df.T
    phyletic_df.to_csv(phyletic_csv)
    phyletic_sequence = phyletic_df.apply(lambda row: ''.join(row.values.astype(str)), axis=1)
    with open(phyletic_text, 'w') as out_phyletic:
        for i in range(len(phyletic_sequence)):
            out_phyletic.write(f'>{phyletic_sequence.index[i]}\n')
            out_phyletic.write(f'{phyletic_sequence.iloc[i]}\n')
    with open(gene_order, 'w') as out_gene_order:
        out_gene_order.write(','.join(phyletic_df.columns))


# generating clustered presence/absence map
# working_directory = r'C:\Users\TalPNB2\Downloads'


def create_presence_absence_map(phyletic_csv, x_label, out_fig_path, colors=['silver', 'lightseagreen'],
                                base_font_size=8, col_cluster=True):
    df = pd.read_csv(phyletic_csv)
    df.set_index(df.columns[0], inplace=True)
    n_rows, n_cols = df.shape
    fig_width = max(n_cols * 0.18, 14)
    fig_height = max(n_rows * 0.22, 4)

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
        cbar_pos=(0.02, 0.75, 0.03, 0.15)
    )
    clustermap.ax_heatmap.xaxis.set_label_position('top')
    clustermap.ax_heatmap.xaxis.tick_top()
    if n_cols > 120:
        clustermap.ax_heatmap.set_xticklabels([])
    else:
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
    clustermap.ax_heatmap.set_xlabel(x_label, fontsize=base_font_size + 2, labelpad=10)
    clustermap.ax_heatmap.set_ylabel('Genome', fontsize=base_font_size + 2, labelpad=10)
    colorbar = clustermap.ax_heatmap.collections[0].colorbar
    colorbar.set_ticks([0.25, 0.75])
    colorbar.set_ticklabels(['Absent', 'Present'])
    colorbar.ax.tick_params(labelsize=base_font_size)
    clustermap.fig.savefig(
        out_fig_path,
        dpi=300,
        bbox_inches='tight',
        pad_inches=0.1
    )

    plt.close(clustermap.fig)


# T4Es:

in_f_T4Es = os.path.join(working_directory, 'out_learning',
                         'consensus_predictions_with_annotations_and_ortho_table.csv')
phyletic_csv_T4Es = os.path.join(working_directory, 'PresenceAbsence_T4Es.csv')
phyletic_text_T4Es = os.path.join(working_directory, 'Effectors_phyletic_pattern.txt')
gene_order_T4Es = os.path.join(working_directory, 'T4Es_order_in_PhyleticPattern.txt')
create_phyletic_patterns(in_f_T4Es, phyletic_csv_T4Es, phyletic_text_T4Es, gene_order_T4Es, T4Es=True)
create_presence_absence_map(phyletic_csv_T4Es, 'T4E (ortholog group)', os.path.join(working_directory,
                                                                                    'T4Es_presence_absence.png'))

# T4SS:

in_file_T4SS = os.path.join(working_directory, 'T4SS.csv')
phyletic_csv_T4SS = os.path.join(working_directory, 'PresenceAbsence_T4SS.csv')
phyletic_text_T4SS = os.path.join(working_directory, 'T4SS_phyletic_pattern.txt')
gene_order_T4SS = os.path.join(working_directory, 'T4SS_order_in_PhyleticPattern.txt')
create_phyletic_patterns(in_file_T4SS, phyletic_csv_T4SS, phyletic_text_T4SS, gene_order_T4SS, start_col=1)
create_presence_absence_map(phyletic_csv_T4SS, 'T4SS component', os.path.join(working_directory,
                                                                              'T4SS_presence_absence.png'),
                            col_cluster=False)

# chaperones:
in_f_chaperones = os.path.join(working_directory, 'chaperones.csv')

if os.path.exists(in_f_chaperones):
    df = pd.read_csv(in_f_chaperones)
    n_rows, n_cols = df.shape
    if n_rows > 1 and n_cols > 1:
        phyletic_csv_chaperones = os.path.join(working_directory, 'PresenceAbsence_chaperones.csv')
        phyletic_text_chaperones = os.path.join(working_directory, 'chaperones_phyletic_pattern.txt')
        gene_order_chaperones = os.path.join(working_directory, 'chaperones_order_in_PhyleticPattern.txt')
        create_phyletic_patterns(in_f_chaperones, phyletic_csv_chaperones, phyletic_text_chaperones, gene_order_chaperones, start_col=1)
        create_presence_absence_map(phyletic_csv_chaperones, 'chaperone', os.path.join(working_directory, 'chaperones_presence_absence.png'))