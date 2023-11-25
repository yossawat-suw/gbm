
# Choose group to import either control or gbm
radiation <- c("control", "radiated")
runs <- c("run1", "run2")

# Parameter
# cell input to be analyse

object <- "all"

sigs <- c("6metamodules", "4metamodules_new_50", "3metamodules_new_50", "4metamodules_new_50_mes")

sig <- sigs[2]

# For reference based
merges <- c("6metamodules", "4_merge_metamodules", "4_merge_metamodules_3celltypes", "4_merge_metamodules_mes")

merge <- merges[2]
# pick which celltype to be analyse

# all_celltypes <- c("AClike", "MESlike", "NPClike", "OPClike")
# chosen_celltypes <- all_celltypes[c(2)]
#annotation

# For tool to be analized in annotation
anno.picks <- c(by_type = 1, by_consistency = 2, only_diverge = 3, marker_based = 4)
anno.pick <- anno.picks[3]



