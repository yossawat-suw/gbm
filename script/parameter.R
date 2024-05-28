
## common

object <- "all"



#marker-based
sigs <- c( "4metamodules_new_50_mes",
           "2metamodules_NPC_OPC", 
           "3metamodules_new_50", 
           "4metamodules_new_50",
           "npcopc_ac_mes",
           "npcopc_acmes",
           "4pathways",
           "neurodev_bulk")


# For reference based
merges <- c("4_merge_metamodules_mes",
            "2_merge_metamodules", 
            "4_merge_metamodules_3celltypes",
            "4_merge_metamodules",
            "npcopc_ac_mes",
            "npcopc_acmes",
            "4pathways",
            "neurodev")


all_cellstates <- c("MESlike",
                   "AClike", 
                   "NPClike", 
                   "OPClike",
                   "NPC_OPClike",
                   "AC_MESlike",
                   "GPM","MTC","NEU","PPR",
                   "Astro","Mesenchymal","Neuronal","Oligo","Progenitor","Unassigned")

chosing <- list(c(1),c(2,5),c(2,3,4),c(1,2,3,4),c(1,2,5),c(5,6),c(7,8,9,10),c(11,12,13,14,15))




## common

object <- "all"

run_eachs <- c(TRUE,FALSE)
runs <- c("each","whole")


## consensus
unknowns <- c("","_no_unknown")


#choosing tools
radiations <- c("control","radiated")

radiation <- c("control","radiated")


cellstate_colors <- c(
  "OPClike" = "#1f77b4",  # Replace 'cellstate1' with actual cellstate names
  "unknown" = "grey",
  "NPClike" = "#2ca02c",
  "AClike" = "#d62728",
  "MESlike" = "#9467bd"
  # Add more colors as needed for each cellstate
)
