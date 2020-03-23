library(readr)
library(tidyr)
library(ape)
library(tidytree)
library(ggplot2)
library(ggtree)
# reconstruct the tree from the node_to_children output file ------------------
file_base <- snakemake@params[['output_dir']]
node_to_children <- read_csv(snakemake@input[['node_to_children']])
# transform so that there is a one-to-one mapping of node to child
node_to_children <- pivot_longer(node_to_children, cols = -node_id, 
                                 names_to = "original", values_to = "child")
node_to_children <- filter(node_to_children, child != -1)

# parent, node, branch_length, label
tree_tibble <- data.frame(parent = node_to_children$node_id, node = node_to_children$child, 
                          branch_length = 1, label = node_to_children$child)


# add metadata ------------------------------------------------------------
# using other csv files, add metadata to the tree_tibble, including hash 
# identity, taxonomic assignment

leaves_to_hash <- read_csv(snakemake@input[['leaves_to_hash']])
hashes_to_fragment <- read_csv(snakemake@input[["hash_to_frag"]])
# read in and parse lineage of hashes_to_tax
hashes_to_tax <- read_csv(snakemake@input[['hash_to_tax']])
hashes_to_tax <- separate(data = hashes_to_tax, col = "lineage", sep = ";",
                          into = c("superkingdom", "phylum", "class", "order", 
                                   "family", "genus", "species", "strain"))

# "label" currently holds node info
tree_tibble <- full_join(tree_tibble, leaves_to_hash, by = c("label" = "leaf_id"))
tree_tibble <- full_join(tree_tibble, hashes_to_fragment, by = "hashval")
tree_tibble <- full_join(tree_tibble, hashes_to_tax, by = "hashval")

# add group clade based on lca lineage. 
# group by family to test
tree_tibble$group <- tree_tibble$phylum
tree_tibble$group <- ifelse(tree_tibble$group == "", NA, tree_tibble$group)
#print(tree_tibble$group)

# plot the tree -----------------------------------------------------------

# convert the tibble to phylo and info objects
tree <- as.phylo(tree_tibble) 
tree_info <- tree_tibble %>%
  select(-parent, -branch_length, -label)

# create a plot
p <- ggtree(tree, layout="circular") 

# add color to plot
svg(filename = snakemake@output[['svg']], width = 7, height = 7)
p %<+%
  tree_info +
  geom_tippoint(aes(color=group)) +
  labs(title = snakemake@wildcards[['f']])
dev.off()

pdf(file = snakemake@output[['pdf']], width = 7, height = 7)
p %<+%
  tree_info +
  geom_tippoint(aes(color=group)) +
  labs(title = snakemake@wildcards[['f']])
dev.off()

