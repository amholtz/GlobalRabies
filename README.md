# About

This repository contains scripts and data that were used to analyze the global spread of canine-rabies virus. The paper can be found here: LINK

Despite the rapid growth in viral sequencing, statistical methods face challenges in handling historical viral endemic diseases with vast amounts of underutilized partial sequence data. We propose a phylogenetic pipeline harnessing full genomes as well as partial sequences to investigate historical pathogen spread between countries. Its application to Rabies virus (RABV) yields precise dating and confident estimates of geographic dispersal. By using full genomes and partial sequences, we limit both geographic and genetic bias which often hinder studies that focus on specific genes. Our pipeline reveals an emergence of the present canine-mediated RABV between 1301 and 1401 and estimates regional introductions over a 700-year period. This geographic reconstruction enables us to locate episodes of human-mediated introductions of RABV around the globe and examine the role that European colonization played in its spread. Our work enables phylogeographic reconstruction on large and genetically diverse datasets for many viral pathogens.



---
# Introduction

The following work was performed with [MAFFT(v7.505)](https://doi.org/10.1093/nar/gkf436), [FastTree(v2.1.11)](https://doi.org/10.1371/journal.pone.0009490), [Goalign(v0.3.5)](https://github.com/evolbioinfo/goalign), [Gotree(v0.4.4)](https://github.com/evolbioinfo/gotree), [TempEst(v1.5.3)](https://doi.org/10.1093/ve/vew007), [HyPhy(v2.5.40)](https://github.com/veg/hyphy), [IQTREE2(v2.2.2.2)](10.1093/molbev/msaa015), [LSD2(v1.8.8)](https://doi.org/10.1093/sysbio/syv068), [PastML(v.1.9.34)](10.1093/molbev/msz131), and [iTol](https://itol.embl.de/tree/1579917420235811657296942#). In addition custom scripts in R and Python were used, which can be found in [R](https://github.com/amholtz/GlobalRabies/tree/main/R) and [Python](https://github.com/amholtz/GlobalRabies/tree/main/python) folders. R version 4.2.1 was used with the following packages, dplyr, tidyverse, ggplot2, plotly, treeio, phangorn, cepiigeogist, countrycode, reshape, data.table, DT, optparse, lubridate, seqinR, readr, taxize, rworldmap, googleVis,rgdal, scales, wesanderson, ape, and Quartet. Python version 3.8 was used with the following packages numpy, pandas, random, pastml.tree, and collections   

The intermediate data files can be found in the [data folder](https://github.com/amholtz/GlobalRabies/tree/main/data). To reproduce the analyses from scratch follow the instructions below.


### Set up

Data can either be downloaded in the [data folder](https://github.com/amholtz/GlobalRabies/tree/main/data). Large alignment files can be downloaded [here](https://www.dropbox.com/scl/fo/nnaz349rkwqew3qsf2fhp/h?dl=0&rlkey=bd0b65ql5ewy8szn29j4wndvu)

Initial sequence dataset was downloaded from [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Lyssavirus%20rabies,%20taxid:11292) by  searching for taxid:11292. Download [fasta](https://www.dropbox.com/s/qflcjo106h2tnw6/allRABV.fasta?dl=0) and [metadata file](https://github.com/amholtz/GlobalRabies/blob/main/data/metadata_edited.tab) and save in data folder as allRABV.fasta


### Sequence alignment
1.  Global alignment by [MAFFT(v7.505)](https://doi.org/10.1093/nar/gkf436)
  ```
  mafft --reorder --keeplength --maxambiguous 0.05 --addfragments ../data/allRABV.fasta --auto ../data/rabv_reference_1988.fasta  > with_keeplength_RABV.fasta
  ```
2. Host species were categorized by family and order to simplify by custom script [species_host_table.R](https://github.com/amholtz/GlobalRabies/blob/main/R/species_host_table.R)
  ```
  Rscript --vanilla species_host_table.R --meta ../data/metadata.tab --table ../data/species_host_table.csv
  ```

3.  The custom script, ([clean_rabv.R](https://github.com/amholtz/GlobalRabies/blob/main/R/clean_RABV.R)), organized sequences by subgenomic region. NC_001542, the reference genome was cut at the positions in the table below [(partition_RABVGenes.txt)](https://github.com/amholtz/GlobalRabies/blob/main/data/sequence_alignments/gene_specific_analysis/partition_RABVGenes.txt) which represent start and stop codons for each gene. As an example, sequences categorized as G gene, contain more than 200 nucleotides between start and stop codons and were saved as a new line in a text file. A quality check was conducted to remove sequences that were (1) missing date and country information, (2) older than 1972, (3) identified as vaccine or laboratory strains, (4) with coding regions shorter than 200 nucleotides. As a result, 14,752 sequences were retained for this study.

  | Gene      | Position Start | Position End |
  |-----------|----------------|--------------|
  | N protein | 71             | 1423         |
  | P protein | 1514           | 2407         |
  | M protein | 2496           | 3104         |
  | G protein | 3318           | 4892         |
  | L protein | 5418           | 11846        |

  ```
  Rscript --vanilla clean_RABV.R --meta ../data/metadata.tab --aln ../data/with_keeplength_RABV.fasta --host_table ../data/species_host_table.csv  --out_n_text ../data/sequence_alignments/gene_specific_analysis/n.txt --out_p_text ../data/sequence_alignments/gene_specific_analysis/p.txt --out_m_text ../data/sequence_alignments/gene_specific_analysis/m.txt --out_g_text ../data/sequence_alignments/gene_specific_analysis/g.txt --out_l_text ../data/sequence_alignments/gene_specific_analysis/l.txt --out_wgs_text ../data/sequence_alignments/gene_specific_analysis/wgs.txt --meta_out ../data/metadata_edited.tab
  ```


4. The sequence alignment was  split into 4 different files, representing the coding regions of the 4 genes defined in [partition_RABVGenes.txt](https://github.com/amholtz/GlobalRabies/blob/main/data/sequence_alignments/gene_specific_analysis/partition_RABVGenes.txt).
  ```
  goalign split -i ../data/with_keeplength_RABV.fasta --partition ../data/sequence_alignments/gene_specific_analysis/partition_RABVGenes.txt --out-prefix ../data/sequence_alignments/gene_specific_analysis/cutalign_
  ```

5.  Each gene-specific alignment is then subsected for sequences grouped in step 3 using [Goalign(v0.3.5)](https://github.com/evolbioinfo/goalign) along with WGS  (G gene example)

  ```
  cat ../data/sequence_alignments/gene_specific_analysis/G.txt ../data/sequence_alignments/gene_specific_analysis/wgs.txt > ../data/sequence_alignments/gene_specific_analysis/G_wgs.txt

  goalign subset -i ../data/sequence_alignments/gene_specific_analysis/cutalign_G.fa -f ../data/sequence_alignments/gene_specific_analysis/G_wgs.txt --unaligned -o ../data/sequence_alignments/gene_specific_analysis/cutalign_G_only.fa
  ```

5.  Each gene is then aligned independently by [MAFFT(v7.505)](https://doi.org/10.1093/nar/gkf436) according to the cut reference sequence (Ex: G gene from reference + all RABV sequences that were classified as G gene)
```
mafft --reorder --keeplength --maxambiguous 0.05 --addfragments ../data/sequence_alignments/gene_specific_analysis/cutalign_G_only.fa --auto ../data/sequence_alignments/gene_specific_analysis/Ggene_Ref.fa > ../data/sequence_alignments/gene_specific_analysis/Ggene_aln.fa
```

6. [Goalign(v0.3.5)](https://github.com/evolbioinfo/goalign) concat was then used to concatenate the aligned sequences back together (without noncoding regions)
```
goalign concat -i ../data/sequence_alignments/gene_specific_analysis/Ngene_aln.fa ../data/sequence_alignments/gene_specific_analysis/Pgene_aln.fa ../data/sequence_alignments/gene_specific_analysis/Mgene_aln.fa ../data/sequence_alignments/gene_specific_analysis/Ggene_aln.fa ../data/sequence_alignments/gene_specific_analysis/Lgene_aln.fa -o ../data/concat_seq_genes.fasta
```
Note: Download [concat_seq_genes](https://www.dropbox.com/s/517xr3gl38ysi43/concat_seq_genes.fasta?dl=0)

![Alt text](https://github.com/amholtz/GlobalRabies/blob/main/concatenation_genes.png)

### Phylogenetic Tree Reconstruction & Dating
1. A global phylogenetic tree was reconstructed using [FastTree(v2.1.11)](https://doi.org/10.1371/journal.pone.0009490)  on all sequences - [iTol link to tree result](https://itol.embl.de/tree/15799174202126551652369486#)
```
~/FastTreeMP -gtr -gamma -nt ../data/concat_seq_genes.fasta > ../data/genewise_aln_RABV.nwk
```
2.  Canine Cluster Subsected- Sequences IDs under cluster defining node were identified (n=10209) in [iTOL](https://doi.org/10.1093/nar/gkab301) and saved as a text file [(canine_ids.txt)](https://github.com/amholtz/GlobalRabies/blob/main/data/canine_ids.txt)

3. A subset alignment of canine-mediate sequences (n=10209) was selected from the original alignment [(concat_seq_genes.fasta)](https://www.dropbox.com/s/517xr3gl38ysi43/concat_seq_genes.fasta?dl=0) using [Goalign(v0.3.5)](https://github.com/evolbioinfo/goalign)
```
goalign subset -i ../data/concat_seq_genes.fasta -f ../data/canine_ids.txt -o ../data/canine.fa
```

3. A phylogenetic tree of all canine-mediated sequences ([canine.fa](https://www.dropbox.com/s/xji8hmmzykae27x/canine.fa?dl=0)) was reconstructed using [FastTree(v2.1.11)](https://doi.org/10.1371/journal.pone.0009490)
```
~/FastTreeMP -gtr -gamma -nt ../data/canine.fa > ../data/RABV_canine10209.nwk
```

4. Rooting was accomplished via [TempEst(v1.5.3)](https://doi.org/10.1093/ve/vew007) by residual-mean square function. Date file ([Tempest_fullCanine.tab](https://github.com/amholtz/GlobalRabies/blob/main/data/Tempest_fullCanine.tab)) is adapted to [Tempest format](https://beast.community/tempest_tutorial)  from our [metadata file](https://github.com/amholtz/GlobalRabies/blob/main/data/metadata_edited.tab)

  **Input**: [RABV_canine10209.nwk](https://github.com/amholtz/GlobalRabies/blob/main/data/RABV_canine10209.nwk), [Tempest_fullCanine.tab](https://github.com/amholtz/GlobalRabies/blob/main/data/Tempest_fullCanine.tab)

  **Output**: [TempestRooted_RABV_canine.nwk](https://github.com/amholtz/GlobalRabies/blob/main/data/TempestRooted_RABV_canine.nwk)

4.  Canine Tree Dating - Evolutionary rate from tree of just canine-WGS
###### Pruning Canine Tree for just WGS by [Gotree(v0.4.4)](https://github.com/evolbioinfo/gotree)
```
gotree prune -i ../data/TempestRooted_RABV_canine.nwk -f ../data/sequence_alignments/gene_specific_analysis/wgs.txt -r -o ../data/wgs_TempestRooted1327_canine.nwk
```
###### Evolutionary Rate [(rate.txt)](https://github.com/amholtz/GlobalRabies/blob/main/data/rate.txt) from WGS Pruned Tree estimated by [LSD2(v1.8.8)](https://doi.org/10.1093/sysbio/syv068)
```
lsd2 -i ../data/wgs_TempestRooted1327_canine.nwk -d ../data/canine_lsd2_dates.txt -o ../data/wgs_TempestRooted1327_canine_CI.result -e 3 -s 10860 -f 1000
```
Evolutionary rate: 0.000199876 [0.000194762; 0.000221632]

5. Rate applied on entire canine-RABV Tree with [LSD2(v1.8.8)](https://doi.org/10.1093/sysbio/syv068)
```
lsd2 -i ../data/TempestRooted_RABV_canine.nwk -d ../data/canine_lsd2_dates.txt -o ../data/TempestRooted1327_WGSRate_OutRem.date.nwk -e 5 -s 10860 -f 1000 -w ../data/rate.txt
```

### Purification and Diversifying Selection was performed using [HyPhy(v2.5.40)](https://github.com/veg/hyphy)

```
$mpirun -np 6 HYPHYMPI meme --alignment ../data/selection_model/NT_macse_out.fa --tree ../data/selection_model/troupin_tree.nwk
$mpirun -np 6 HYPHYMPI absrel --alignment ../data/selection_model/NT_macse_out.fa --tree ../data/selection_model/troupin_tree.nwk
$mpirun -np 6 HYPHYMPI fel --alignment ../data/selection_model/NT_macse_out.fa --tree ../data/selection_model/troupin_tree.nwk --branches Internal --ci Yes

```


### Phylogeography


#### Ancestral Character Reconstruction on Country Level (Full Tree) was estimated with [PastML(v.1.9.34)](10.1093/molbev/msz131)
```
pastml -t ../data/TempestRooted1327_WGSRate_OutRem.date.nwk -d metadata_edited.tab -c Country --prediction_method MPPA --root_date 1356.74 --html_compressed HTML_compressed_canine_MPPA_nexus_100.html --upload_to_itol -o canine_MPPA_nexus_pastML --tip_size_threshold 100
```
**[iTol Tree with ACR Estimation annotations](https://itol.embl.de/tree/1579917420235811657296942#)** & **[PastML Visualization- ACR Country Results](https://github.com/amholtz/GlobalRabies/tree/main/data/ACR_Results/Country/Full_Tree)**

#### Subsampling by country

1.  Smart-Subsampling by country was performed using the custom script [py_subsampling.py](https://github.com/amholtz/GlobalRabies/blob/main/python/py_subsampling.py)
```
python3 py_subsampling.py --input_tree ../data/TempestRooted1327_WGSRate_OutRem.date.nwk --input_locs ../data/metadata_edited.tab --size 5500 --output_ids subsampled_5500_1 subsampled_5500_2 subsampled_5500_3 subsampled_5500_4 subsampled_5500_5
```


#### Tree Reconstruction, dating, and comparison of subsampled tree (Example: Subsample 5)

1. Subsection of alignment for subsample 5 sequence IDs by [Goalign(v0.3.5)](https://github.com/evolbioinfo/goalign)

  ```
  goalign subset -i ../data/concat_seq_genes.fasta -f ../data/subsampled_5500_5.txt -o ../data/subsample5.fa

  ```
2.  Phylogenetic Reconstruction by [IQTREE2(v2.2.2.2)](10.1093/molbev/msaa015) GTR+I+G4 and partioning (Example: Subsample 5)
```
iqtree2 -s ../data/subsample5.fa -st DNA -nt 8 -alrt 0 -m GTR+I+G4 -B 1000 -p ../data/gene_partition.txt -pre ../data/subsample5.fa
```

3. Rooting via [TempEst(v1.5.3)](https://doi.org/10.1093/ve/vew007) residual-mean square function to infer the best-fitting root)

  **Input**:[subsample5.fa.treefile](https://github.com/amholtz/GlobalRabies/blob/main/data/subsample5.fa.treefile), [Tempest_fullCanine.tab](https://github.com/amholtz/GlobalRabies/blob/main/data/Tempest_fullCanine.tab)

  **Output**:[TempEstRooted_subsampled_5000_5.fa.treefile](https://github.com/amholtz/GlobalRabies/blob/main/data/TempEstRooted_subsampled_5000_5.fa.treefile)

4.  Rate from WGS [(rate.txt)](https://github.com/amholtz/GlobalRabies/blob/main/data/rate.txt) applied on Subsample 5 Tree with [LSD2(v1.8.8)](https://doi.org/10.1093/sysbio/syv068)
```
lsd2 -i ../data/TempEstRooted_subsampled_5000_5.fa.treefile -d ../data/fullCanine_lsd2.tab -s 10860 -o sub5_CI -f 1000 -e 3 -w rate.txt
```


#### Comparing Subsampled tree reconstructions by Triplet Distance Calculations (Example: Subsample 5)

1. Prune the full-tree so the tips match the tips contained in the subsample
```
gotree prune -i ../data/TempestRooted_RABV_canine.nwk -f ../data/subsampled_5500_5.txt -r -o ../data/full_tree_Sub5Tips.nwk
```
2.  Resolve multifurcations with [Gotree(v0.4.4)](https://github.com/evolbioinfo/gotree)
```
gotree resolve -i ../data/full_tree_Sub5Tips.nwk -o ../data/full_tree_Sub5Tips.nwk_bi.nwk  
```

3.  Calculate Triplet Distance with full-tree by custom script [triplet_distance.R](https://github.com/amholtz/GlobalRabies/blob/main/R/triplet_distance.R)
```
Rscript --vanilla triplet_distance.R --subtree ../data/TempEstRooted_subsampled_5000_5.fa.treefile --fulltree ../data/full_tree_Sub5Tips_bi.nwk --tips 5484 -o ../data/subsample5_tripletdistance.csv
```


#### Ancestral Character Reconstruction on Country Level (Subsample 5) was estimated with [PastML(v.1.9.34)](10.1093/molbev/msz131)
(working directory ~data/ACR_Results/Country/Sub5/)
```
pastml -t ../data/sub5_CI.nwk -d ../data/metadata_edited.tab -c Country --prediction_method MPPA --root_date 1365 --html_compressed HTML_compressed_canine_5000_5_MPPA_nexus_100.html --upload_to_itol -o canine_5000_5_subsample_MPPA_nexus_pastML --tip_size_threshold 100
```
**[iTol Tree with ACR Estimation annotations](https://itol.embl.de/tree/15799174109116831658497579)** & **[PastML Visualization- ACR Country Results](https://github.com/amholtz/GlobalRabies/blob/main/data/ACR_Results/Country/Sub5)**

#### Consensus tree was created from the aggregation of ACR results on the country level
Consensus tree was created by comparing country estimations for each node across all subsamples and the full canine tree following custom script [ACR_Sub_comparsion.R](https://github.com/amholtz/GlobalRabies/blob/main/R/ACR_Sub_comparison.R)

```
Rscript --vanilla ACR_Sub_comparison.R --meta ../data/metadata_edited.tab --path ../data/ACR_Results/Country/ --prefix consensus
```
### Ancestral Character Reconstructions by [PastML(v.1.9.34)](10.1093/molbev/msz131) on consensus tree

###### 1. Country Level (Consensus) was estimated with [PastML(v.1.9.34)](10.1093/molbev/msz131)

**[PastML Visualization- ACR Country Consensus Results](https://github.com/amholtz/GlobalRabies/tree/main/data/ACR_Results/Country/Consensus_Tree)**
```
pastml -t ../data/ACR_Results/Country/Consensus_Tree/consensus_tree6096.nwkconsensus_tree6096.nwk --prediction_method COPY --root_date 1356.74 -o canine_AdaptedFullSUBSTATE_pastML -d ../data/ACR_Results/Country/Consensus_Tree/consensus_inSUB_Full_Sub_states_2Col.tab --upload_to_itol --columns Country_Full agg --html_compressed ../data/ACR_Results/Country/Consensus_Tree/FullSubAdapted_2Col_50.html --tip_size_threshold 50 --colours manual_colours.character_Country_Full.tab colours.character_agg.tab
```



###### 2. Ancestral Character Reconstruction on Regional Level (Consensus) was estimated with [PastML(v.1.9.34)](10.1093/molbev/msz131)
**[PastML Visualization- ACR Regional Results](https://github.com/amholtz/GlobalRabies/tree/main/data/ACR_Results/Region)**
```
pastml -t ../data/ACR_Results/Country/Consensus_Tree/consensus_tree6096.nwk -d metadata_edited.tab -c region23 --prediction_method MPPA --root_date 1356.74 --html_compressed ../data/ACR_Results/Region/pastml_compressed_visualisation_region23.html --upload_to_itol -o pastml_region23 --tip_size_threshold 100
```


###### 3. Ancestral Character Reconstruction on Colony Level (Consensus) was estimated with [PastML(v.1.9.34)](10.1093/molbev/msz131)
**[PastML Visualization- ACR Colony Results](https://github.com/amholtz/GlobalRabies/tree/main/data/ACR_Results/Colony)**

```
pastml -t ../data/ACR_Results/Country/Consensus_Tree/consensus_tree6096.nwk -d metadata_edited.tab -c colony --prediction_method MPPA --root_date 1356.74 --html_compressed ../data/ACR_Results/Colony/pastml_compressed_visualisation_colony.html --upload_to_itol -o pastml_colony --tip_size_threshold 100
```

###### 4. Ancestral Character Reconstruction on Clade Level (Concensus) was estimated with [PastML(v.1.9.34)](10.1093/molbev/msz131)
**[PastML Visualization- ACR Clade Results](https://github.com/amholtz/GlobalRabies/tree/main/data/ACR_Results/Clade)**

```
pastml -t ../data/ACR_Results/Country/Consensus_Tree/consensus_tree6096.nwk -d metadata_edited.tab -c clade --prediction_method MPPA --root_date 1356.74 --html_compressed ../data/ACR_Results/Clade/pastml_compressed_visualisation_clade.html --upload_to_itol -o pastml_clade --tip_size_threshold 100
```


### Human-Mediated Introductions
Human-mediated transmission were inferred using the custom script [introduction_events_FullTree.R](https://github.com/amholtz/GlobalRabies/blob/main/R/introduction_events_FullTree.R)

```
Rscript --vanilla introduction_events_FullTree.R --meta ../data/metadata_edited.tab --tree ../data/ACR_Results/Country/Full_Tree/named.tree_TempestRooted1327_WGSRate_OutRem.date.nexus --annotations ../data/ACR_Results/Country/Full_Tree/combined_ancestral_states.tab --probabilities ../data/ACR_Results/Country/Full_Tree/marginal_probabilities.character_Country.model_F81.tab --newick ../data/ACR_Results/Country/Full_Tree/named.tree_TempestRooted1327_WGSRate_OutRem.date.nwk --prefix fulltree
```

### Final Metadata file with exclusion and inclusion information
```
Rscript --vanilla metadata_exclusioncriteria.R --meta ../data/metadata_edited.tab --full_tree ../data/genewise_aln_RABV.nwk --canine_orig ../data/RABV_canine10209.nwk --canine_lsd ../data/TempestRooted1327_WGSRate_OutRem.date.nwk --canine_cons ../data/ACR_Results/Country/Consensus_Tree/consensus_tree6096.nwk --meta_out ../data/metadata_edited_exclusion.tab
```


## Interactive Visualizations and Tables can be viewed **[here](https://amholtz.github.io/GlobalRabies/)**

# Get in touch

To report any issues, please [open an isssue](https://github.com/amholtz/GlobalRabies/issues).
