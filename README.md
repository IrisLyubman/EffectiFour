# EffectiFour

**A web server for the prediction of type IVB secretion system effectors**

EffectiFour is a machine-learning-based pipeline for identifying bacterial proteins translocated through the Icm/Dot type IVB secretion system (T4BSS).

EffectiFour supports the simultaneous analysis of multiple bacterial genomes, enabling the identification and comparison of candidate effector repertoires across orthologous groups (OGs).

## Web Server

https://effectidor.tau.ac.il/type4b/

## Main Features

- **Pan-genome analysis:** Identification of orthologous groups across multiple bacterial genomes.
- **Machine-learning predictions:** Integration of protein sequence properties, homology, genomic organization, and regulatory information.
- **C-terminal secretion-signal prediction:** Incorporation of ProtT5-derived protein embeddings and a trained classifier.
- **Genomic-context analysis:** GeNLP-based embeddings to identify similarities between the genomic neighborhoods of candidate proteins and known effectors.
- **Regulatory features:** Identification of regulatory motifs associated with T4BSS.
- **Icm/Dot system detection:** Identification of secretion machinery components.

## Input

EffectiFour requires protein-coding sequences in FASTA format, with one file per analyzed genome.

Optional inputs include:

- Complete genomic sequences (FASTA)
- Genome annotation files (GFF)
- Host proteomes
- Proteomes of related bacteria lacking the corresponding secretion system
- Additional known type IVB effectors

## Output

EffectiFour generates ranked candidate effector predictions, including prediction scores, functional annotations, and orthologous-group assignments.

Additional outputs include presence–absence profiles of predicted effectors across genomes, regulatory-feature information, and detected Icm/Dot secretion-system components.

## External Dependencies

EffectiFour uses external tools and models, including MMseqs2, Microbializer 2, ProtT5, and GeNLP. Running the source code locally requires configuring the relevant software environments, reference datasets, and computational resources.

The ProtT5-based secretion-signal classifier developed for EffectiFour is maintained in a separate repository:

https://github.com/IrisLyubman/T4Es_secretion_signal_prediction

## Citation

Lyubman I., Schers A., Segal G., and Pupko T.

*EffectiFour: a web server for the prediction of type IVB secretion system effectors.*

Manuscript in preparation.
