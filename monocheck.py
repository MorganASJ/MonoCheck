import argparse
from taxonomygen import fetch_taxonomy
from ete3 import Tree

# Root tree based on taxa listed in file
def root_tree(tree, outgroup_file):
    with open(outgroup_file) as f:
        outgroup_taxa = f.read().strip().split("\n")

    print(outgroup_taxa) # Works

    # if set of taxa, find commmon ancestor node, else root on leaf node
    if len(outgroup_taxa) > 1:
        tree.set_outgroup(tree.get_common_ancestor(outgroup_taxa))
    else:
        tree.set_outgroup(outgroup_taxa[0])

    # return transformed tree object
    return tree

# Entry point
def main():

    # Parse arguments python monocheck [tree file] -group [taxonomic group] -outgroup [outgroup file OR outgroup taxon]
    parser = argparse.ArgumentParser(description="Check monophyly of a group in a phylogenetic tree.")
    parser.add_argument("tree", help="Newick file containing the tree")
    parser.add_argument("-group", required=True, help="Group to check for monophyly")
    # parser.add_argument("-family", required=False, help="Group to check for monophyly")
    # parser.add_argument("-order", required=False, help="Group to check for monophyly")
    # parser.add_argument("-genus", required=False, help="Group to check for monophyly")
    # parser.add_argument("-species", required=False, help="Group to check for monophyly")
    # parser.add_argument("-phylum", required=False, help="Group to check for monophyly")
    parser.add_argument("-outgroup", required=True, help="File containing outgroup taxa (one per line)")
    args = parser.parse_args()

    # Load tree
    tree = Tree(args.tree, format=1)

    print(tree) # Works

    # Root the tree
    tree = root_tree(tree, args.outgroup)

    print(tree) # Works

    # Load or generate taxonomy
    taxonomy_file = "taxonomy.tax"
    taxonomy_path = fetch_taxonomy([leaf.name for leaf in tree.get_leaves()], taxonomy_file)

    return 0


if __name__ == "__main__":
    main()