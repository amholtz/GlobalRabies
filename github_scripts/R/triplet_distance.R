####################################################
##                           Triplet Distances
## ####   Calculate Triplet Distances between subsampled tree 
## ####   by IQTREE + Gene Partitioning 
## ####   and 
## ####   Full Tree from FastTree reconstruction without partitioning
## author: Andrew Holtz
## creation date: 2022/03/16
###################################################

# FullTree has been adapted to only contain tips in the subsample and is converted
# to a bifurcating tree

sub1 <- TripletDistance(
  "full_tree_Sub1Tips_bi.nwk", "TempEstRooted_subsampled_5000_1_oR.treefile")
sub2 <- TripletDistance(
  "full_tree_Sub2Tips_bi.nwk", "TempEstRooted_subsampled_5000_2_oR.treefile")
sub3 <- TripletDistance(
  "full_tree_Sub3Tips_bi.nwk", "TempEstRooted_subsampled_5000_3_oR.treefile")
sub4 <- TripletDistance(
  "full_tree_Sub4Tips_bi.nwk", "TempEstRooted_subsampled_5000_4_oR.treefile")
sub5 <- TripletDistance(
  "full_tree_Sub5Tips_bi.nwk", "TempEstRooted_subsampled_5000_5_oR.treefile")

# Binomial coefficient for 5484 tips: 27472834684

per_norm_diff_1 <- sub1/27472834684
per_norm_diff_2 <- sub1/27472834684
per_norm_diff_3 <- sub1/27472834684
per_norm_diff_4 <- sub1/27472834684
per_norm_diff_5 <- sub1/27472834684


