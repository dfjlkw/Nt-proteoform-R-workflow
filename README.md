# Nt-proteoform-R-workflow
This repository contains data and code of the R workflow to select N-terminal proteoform pairs for interactome analysis published in:

**[N-terminal proteoforms may engage in different protein complexes](https://pubmed.ncbi.nlm.nih.gov/37316325/)**

Annelies Bogaert 1,2 , Daria Fijalkowska 1,2, An Staes 1,2, Tessa Van de Steene 1,2, Marnik Vuylsteke 3, Charlotte Stadler 4 , Sven Eyckerman 1,2, Kerstin Spirohn 5,6,7, Tong Hao 5,6,7,  Michael A. Calderwoo d5,6,7 and Kris Gevaert 1,2,*

1 VIB Center for Medical Biotechnology, VIB, Ghent, B9052, Belgium

2 Department of Biomolecular Medicine, Ghent University, Ghent, B9052, Belgium

3 Gnomixx, Melle, Belgium

4 Department of Protein Science, KTH Royal Institute of Technology and Science for Life Laboratories, Stockholm, Sweden

5 Center for Cancer Systems Biology (CCSB), Dana-Farber Cancer Institute, Boston, USA

6 Department of Genetics, Blavatnik Institute, Harvard Medical School, Boston, USA

7 Department of Cancer Biology, Dana-Farber Cancer Institute, Boston, USA

\* To whom correspondence should be addressed. Tel: +32 (0) 9 224 98 35; Email: kris.gevaert@vib-ugent.be



## Description of analysis steps
### Selection of N-terminal proteoform pairs for interactome analysis
Nt-peptides selected using the KNIME workflow outlined in [Bogaert, A., et al., Limited Evidence for Protein Products of Noncoding Transcripts in the HEK293T Cellular Cytosol. Mol Cell Proteomics, 2022](https://pubmed.ncbi.nlm.nih.gov/35788065/) and Supplementary Materials and Methods were further annotated and curated in R version 4.1.0. 

**Part I** – annotation of genomic features. Human gene and transcript annotations were downloaded from Ensembl BioMart Archive Release 98 (September 2019) and proteoform accessions were linked to gene identifiers, names, biotypes and transcript support levels (TSL). Subsequently, Nt-peptides were mapped to genomic positions and a peptide BED file was created. Furthermore, we report the exact Nt-proteoform exon rank, start codon, frame and distance to the corresponding aTIS. 

**Part II** – annotation and prediction of protein sequence features. UniProt annotations of human canonical proteins were downloaded in January 2020, including the following features: signal peptide, transit peptide, propeptide, region, motif, coiled coil, compositional bias, repeat and zinc finger. To determine sequence features annotated in UniProt that Nt-proteoforms had lost, Nt-peptides were exactly matched to canonical human UniProt proteins using dbtoolkit (version 4.2.5). Nt-peptides not exactly mapping to any UniProt protein, alignment of the entire Nt-proteoform to one UniProt reference sequence was performed. From these alignments, sequence ranges that were lost compared to the UniProt reference protein were extracted. We also determined if Nt-truncated proteoforms could be derived from N-terminal processing (removed signal, transit or pro- peptide) considering a ± one amino acid margin of error. For Nt-peptides without an exact match to an UniProt protein, we scanned for short linear motifs in the proteoform and their selected UniProt reference using Eukaryotic Linear Motif (ELM) resource API. TopFIND was used to determine if the Nt-peptides could have been derived from post-translational processing. Additionally, Nt-proteoforms that could derive from N-terminal dipeptidase cleavage are marked; see Table S4 column “dipeptidase”. Finally, we performed a BLAST analysis against the human UniProt or the human non-redundant proteins (NCBI, July 2020). Hits with over 80% of protein sequence identity over 50% of the length are reported.

**Part III** – survey of complementary data(bases). Curated, experimentally determined protein interactions were downloaded from BioGRID version 3.5.182, human genetic phenotypes and disorders available from OMIM were matched by gene identifier using biomaRt version 2.50.3, release 98. Prior experimental evidence for Nt-proteoform expression in human (primary) cells reported by Van Damme et al. 2014 was included for matching N-termini. Cytosolic expression associated was confirmed by three independent sources: our own cytosolic proteomics data, gene ontology GOSlim annotation and cytosolic subcellular localization determined by immunostaining from The Human Protein Atlas.

Genes with multiple proteoforms were classified into four categories: 1. annotated + alternative TIS, 2. multiple alternative TIS, 3. multiple annotated TIS and 4. one TIS, where canonical UniProt and Ensembl aTIS were considered annotated TIS. Subsequently, we calculated a TIS score for each alternative Nt-proteoform (see 12082022-C13-FullMergedResults-R_output.txt). If the gene score was >1 it is included in (12082022-MOST-INTERESTING-C13-FullMergedResults-R_output.txt). For more details, please consult the publication.
