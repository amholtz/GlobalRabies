# About

This repository contains scripts and data that were used to analyze the global spread of canine-rabies virus. The paper can be found here: LINK

Despite the rapid growth in viral sequencing, statistical methods face challenges in handling historical viral endemic diseases with vast amounts of underutilized partial sequence data. We propose a phylogenetic pipeline harnessing full genomes as well as partial sequences to investigate historical pathogen spread between countries. Its application to Rabies virus (RABV) yields precise dating and confident estimates of geographic dispersal. By using full genomes and partial sequences, we limit both geographic and genetic bias which often hinder studies that focus on specific genes. Our pipeline reveals an emergence of the present canine-mediated RABV between 1301 and 1401 and estimates regional introductions over a 700-year period. This geographic reconstruction enables us to locate episodes of human-mediated introductions of RABV around the globe and examine the role that European colonization played in its spread. Our work enables phylogeographic reconstruction on large and genetically diverse datasets for many viral pathogens.



---
# Work Flow Description

**Data can either be downloaded from a link or found in [data folder](https://github.com/amholtz/GlobalRabies/tree/main/data)**

### Set up
1. Download fasta and metadata file from [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Lyssavirus%20rabies,%20taxid:11292)
2. Initial Data cleaning was done using [clean_rabv.R](https://github.com/amholtz/GlobalRabies/blob/main/R/clean_RABV.R) to remove sequences with data integrity issues (see methods)

### Sequence alignment
1.  Global alignment by [MAFFT(v7.505)](https://doi.org/10.1093/nar/gkf436)
```
mafft --reorder --keeplength --compactmapout --maxambiguous 0.05 --addfragments fragments --auto allRABV.fasta > allRABV_aligned.fasta
```
2.  NC_001542, reference genome was cut at the positions in the table below and sequences were categorized into separate gene specific fasta files from the position cut offs which represent the start codon to the end of the coding region of the gene (mRNA). Sequences were added to the fasta gene files they belonged (for example, WGS will appear in each of the 5 gene FA files since they contain each of the five genes)  ([clean_rabv.R](https://github.com/amholtz/GlobalRabies/blob/main/R/clean_RABV.R)). The reference gene (from above step) was added to these files.

  | Gene      | Position Start | Position End |
  |-----------|----------------|--------------|
  | N protein | 71             | 1423         |
  | P protein | 1514           | 2407         |
  | M protein | 2496           | 3104         |
  | G protein | 3318           | 4892         |
  | L protein | 5418           | 11846        |

  ```
Rscript --vanilla clean_rabv.R --meta ../data/meta_full_exclusion_clade_simple.tab --aln ../data/with_keeplength_RABV.fasta --out_wgs_text ../data/sequence_alignments/gene_specific_analysis/wgs.txt --out_n ../data/sequence_alignments/gene_specific_analysis/n.txt --out_p ../data/sequence_alignments/gene_specific_analysis/p.txt --out_m ../data/sequence_alignments/gene_specific_analysis/m.txt --out_g ../data/sequence_alignments/gene_specific_analysis/g.txt --out_l ../data/sequence_alignments/gene_specific_analysis/l.txt
  ```

3.  Original alignment is then subsected according to genetic categorization by [Goalign(v0.3.5)](https://github.com/evolbioinfo/goalign) (G gene example)

  ```
  goalign subset -i allRABV.fasta -f G.txt --unaligned -o Ggene.fasta

  ```

4.  Each gene is then aligned independently by [MAFFT(v7.505)](https://doi.org/10.1093/nar/gkf436) according to the cut reference sequence
Ex: G gene from reference + all RABV sequences that were classified as G gene example with G gene)
```
mafft --reorder --keeplength --compactmapout --maxambiguous 0.05 --addfragments fragments --auto Ggene.fasta > Ggene_aln.fa
```

5. [Goalign(v0.3.5)](https://github.com/evolbioinfo/goalign) concat was then used to concatenate the aligned sequences back together (without noncoding regions)
```
goalign concat -i N_gene_aln.fa P_gene_aln.fa M_gene_aln.fa G_gene_aln.fa L_gene_aln.fa -o concat_seq_genes.fa
```

![Alt text](https://github.com/amholtz/GlobalRabies/blob/main/concatenation_genes.png)

### Phylogenetic Tree Reconstruction & Dating
1. FastTree reconstruction on all sequences - [iTol link to tree result](https://itol.embl.de/tree/15799174202126551652369486#)
```
~/FastTreeMP -gtr -gamma -nt concat_seq_genes.fasta > genewise_aln_RABV.nwk
```
2.  Canine Cluster Subsected- Sequences IDs under cluster defining node were identified in itol and saved as a text file (canine_ids.txt)
```
Gotree prune -f canine_ids.txt -r -i genewise_aln_RABV.nwk -o genewise_aln_RABV_canine.nwk
```
3.  Canine Tree Dating - Evolutionary rate from tree of just canine-WGS
```
lsd2 -i wgs_TempestRooted1327_canine.nwk -d canine_lsd2_dates.txt -o wgs_TempestRooted1327_canine_CI.result -e 3 -s 10860 -f 1000
```
  **Evolutionary rate: 0.000199876 [0.000194762; 0.000221632]**

4. Rate used on entire canine-RABV Tree
```
lsd2 -i TempestRooted1327_canine_10199Taxa_collapsed.nwk-d fullCanine_lsd2.tab -o TempestRooted1327_canine_10199_WGSRate.result-e 5 -s 10860 -f 1000 -w rate.txt
```

### Purification and Diversifying Selection

```
$mpirun -np 6 HYPHYMPI meme --alignment NT_macse_out.fa --tree troupin_tree.nwk
$mpirun -np 6 HYPHYMPI absrel --alignment NT_macse_out.fa --tree troupin_tree.nwk
$mpirun -np 6 HYPHYMPI fel --alignment NT_macse_out.fa --tree troupin_tree.nwk --branches Internal --ci Yes

```


### Phylogeography


#### PastML Country Level (Full Tree)
```
Final_Full_10044_FinalACR_WGSrate_OutRemove.LSD2.nwk -d meta_RABV_cleaned_clade_gene.tab -c Country --prediction_method MPPA --root_date 1356.74 --html_compressed HTML_compressed_canine_MPPA_nexus_100.html --upload_to_itol -o canine_MPPA_nexus_pastML --parameters --tip_size_threshold 100
```
**[iTol Tree with ACR Estimation annotations](https://itol.embl.de/tree/1579917420235811657296942#)** & **[PastML Visualization- ACR Country Results](https://github.com/amholtz/GlobalRabies/tree/main/data/ACR_Results/Country/Full_Tree)**

#### Subsampling by country
Subsampling was done using the custom script [py_subsampling.py](https://github.com/amholtz/GlobalRabies/blob/main/python/py_subsampling.py)


#### Tree Reconstruction, dating, and comparison of subsampled tree (Subsample 5 example)

1.  IQTREE2 Reconstruction
```
iqtree2 -s subsampled_5000_5.fa -st DNA -nt 8 -alrt 0 -m GTR+I+G4 -B 1000 -p gene_partition.txt
```
2.  LSD2 Dating
```
lsd2 -i TempEstRooted_subsampled_5000_5_OutRem.treefile -d fullCanine_lsd2.tab -s 10860 -o sub1_CI -f 1000 -e 3 -w rate.txt
```
3. Triplet Distance Calculation with full-tree with custom script [triplet_distance.R](https://github.com/amholtz/GlobalRabies/blob/main/R/triplet_distance.R)

#### PastML Country Level (Subsample 5 example)
```
named.tree_Subsample5000_5_TempestRooted1327_WGSRate_OutRem.date.nexus -d meta_RABV_cleaned_clade_gene.tab -c Country --prediction_method MPPA --root_date 1365 --html_compressed HTML_compressed_canine_5000_5_MPPA_nexus_100.html --upload_to_itol -o canine_5000_5_subsample_MPPA_nexus_pastML --tip_size_threshold 100
```
**[iTol Tree with ACR Estimation annotations](https://itol.embl.de/tree/15799174109116831658497579)** & **[PastML Visualization- ACR Country Results](https://github.com/amholtz/GlobalRabies/blob/main/data/ACR_Results/Country/Sub5)**

#### Country Level ACR Consensus Tree
Consensus tree was created by comparing country estimations for each node across all subsamples and the full canine tree following custom script [ACR_Sub_comparsion.R](https://github.com/amholtz/GlobalRabies/blob/main/R/ACR_Sub_comparison.R)

###### PastML Consensus ACR

```
pastml -t named.tree_dated.Fullpruned_subsampledTips.nwk --prediction_method COPY --root_date 1356.74 -o canine_AdaptedFullSUBSTATE_pastML -d inSUB_Full_Sub_states_2Col.tab --upload_to_itol --columns Country_Full agg --html_compressed FullSubAdapted_2Col_30.html --tip_size_threshold 30 --colours manual_colours.character_Country_Full.tab colours.character_agg.tab
```
**[PastML Visualization- ACR Country Consensus Results](https://github.com/amholtz/GlobalRabies/tree/main/data/ACR_Results/Country/Consensus_Tree)**


#### PastML Regional Level
```
named.tree_dated.Fullpruned_subsampledTips.nwk -d meta_full_exclusion_clade_simple.tab -c region23 --prediction_method MPPA --root_date 1356.74 --html_compressed pastml_compressed_visualisation_region23.html --upload_to_itol -o pastml_region23 --tip_size_threshold 100
```
**[PastML Visualization- ACR Regional Results](https://github.com/amholtz/GlobalRabies/tree/main/data/ACR_Results/Region)**

#### PastML Colony Level
```
named.tree_dated.Fullpruned_subsampledTips.nwk -d meta_full_exclusion_clade_simple.tab -c colony --prediction_method MPPA --root_date 1356.74 --html_compressed pastml_compressed_visualisation_colony.html --upload_to_itol -o pastml_colony --tip_size_threshold 100
```
**[PastML Visualization- ACR Colony Results](https://github.com/amholtz/GlobalRabies/tree/main/data/ACR_Results/Colony)**


#### PastML Clade Level
```
named.tree_dated.Fullpruned_subsampledTips.nwk -d meta_full_exclusion_clade_simple.tab -c clade --prediction_method MPPA --root_date 1356.74 --html_compressed pastml_compressed_visualisation_clade.html --upload_to_itol -o pastml_clade --tip_size_threshold 100
```
**[PastML Visualization- ACR Clade Results](https://github.com/amholtz/GlobalRabies/tree/main/data/ACR_Results/Clade)**


### Human-Mediated Introductions
Human-mediated transmission were inferred using the custom script [introduction_events_FullTree.R](https://github.com/amholtz/GlobalRabies/blob/main/R/introduction_events_FullTree.R)

## Interactive Visulizations and Tables can be viewed **[here](https://amholtz.github.io/GlobalRabies/)**

# Get in touch

To report any issues, please [open an isssue](https://github.com/amholtz/GlobalRabies/issues).
