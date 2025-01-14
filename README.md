# MonoCheck

**MonoCheck**  is a Python utility script that generates a taxonomy reference file for taxa within a phylogenetic tree and determines whether one or more specified clades are monophyletic.

This tool is particularly useful for researchers studying phylogenetics who need to quickly evaluate the monophyly of groups across large numbers of trees or datasets.

### Usage

MonoCheck is run from the terminal and can be used to quickly determine monophyly or produce a csv report for more complex requests.

    $ monocheck.py [-h] -clade CLADE [CLADE ...] -outgroup OUTGROUP -email EMAIL [-taxonomy TAXONOMY] [-out OUT] [--showtree | --no-showtree] [--resettaxonomy | --no-resettaxonomy] tree
    
**Positional Arguments:**

 - **tree**: Path to a Newick file OR directory containing Newick files.
	Trees input must be in the standard Newick format, and tips must be prefixed with the taxon name i.e.: Genus_species_* . Information following the taxon name will be ignored.

**Optional Arguments:**

 - **-clade**: Specify one or more taxonomic clades to check for monophyly.

 - **-outgroup**: File containing outgroup tips (one taxon per line).

 - **-email**: Your email address for querying NCBI Taxonomy.

 - **-taxonomy**: Manually point the program to an existing taxonomy file. (Optional)

 - **-out**: Path to save the output as a CSV file. (Optional)

 - **--showtree**: Display the rooted tree before checking monophyly. (Optional)

 - **--resettaxonomy**: Force regeneration of the taxonomy file by querying NCBI. (Optional)

# Requirements

MonoCheck requires **BioPython** for the Entrez libraries, **Pandas** for data frame management, and **ete3** for phylogenetic tree data structures and functions.

These packages can be installed by running the following:

    $ pip install BioPython Pandas ete3

Or through creating a conda environment using the supplied environment.yml file:

    $ conda env create -f environment.yml

Activate the conda environment for MonoCheck:

    $ conda activate monocheck

# Using MonoCheck

Researchers studying phylogenetics regularly produce a large number of trees and it can be difficult to inspect each one independently. You may for instance have a set of over 1000 gene trees and need to know what proportion of these trees support a particular evolutionary hypothesis. 

Using MonoCheck we can quickly determine if a tree, or set of trees contain a given clade. This can be useful to see the proportion of trees supporting a particular hypothesis under different gene selections, or phylogenetic model parameters.

MonoCheck works through association with a taxonomy file. You can provide a taxonomy file with the -taxonomy flag, or let MonoCheck generate one by searching NCBI Taxonomy. It is important to note that NCBI Taxonomy as a resource may not necessarily be up to date, and you may instead need to manually change the taxonomy if the clade you wish to test is not recognised within the Taxonomy database.

To test a single tree for the presence or absence of a clade we need to have the tree file and an outgroup file.

### Creating an outgroup file

In order to root the tree we need an outgroup file. This should be a normal text file containing the tip labels of the outgroup. For example:

    Priapulus_caudatus_sc
	Limulus_polyphemus_sc

These will be my outgroup taxa for the test tree in the main directory `test_tree.tre`

### Running MonoCheck on a single tree

With the outgroup generated, we can run MonoCheck with the test tree. The tree will be rooted on the taxa above and then we will be checking if Insecta is monophyletic.

    $ python3 monocheck.py test_tree.tre -clade Insecta -outgroup outgroup.txt -email myemail@university.ac.uk

Since there is no taxonomy file in the directory MonoCheck will search for the taxonomy records for each taxa in the tree and save the taxonomy file to the working directory. This can take a few seconds and up to a few minutes for larger trees. If a taxonomy file is already present in the directory MonoCheck will use the one already there. This means for larger trees you will only need to compute the taxonomy once.

Output from this command:

    Retrieving taxonomy information from NCBI Taxonomy database took 11.2 seconds.

    ----- Checking Monophyly of clades -----
    test_tree.tre>>> Is Insecta monophyletic? Yes
    ----------------------------------------

If we provide the `--showtree` flag we will get a print out of the rooted tree:


             /-Tachypleus_tridentatus_sc
          /-|
       /-|   \-Limulus_polyphemus_sc
      |  |
      |   \-Priapulus_caudatus_sc
      |
    --|      /-Daphnia_pulex_sc
      |   /-|
      |  |  |   /-Tribolium_castaneum_sc
      |  |   \-|
      |  |     |   /-Drosophila_melanogaster_sc
       \-|      \-|
         |         \-Aedes_aegypti_sc
         |
         |   /-Paramacrobiotus_metropolitanus_sc
          \-|
             \-Caenorhabditis_elegans_sc


 ### Running MonoCheck on multiple trees

In the directory trees there are a set of orthofinder gene trees for the same taxa as above. If I wanted to check the monophyly of clades within all these trees manually it could take significant time.

To run MonoCheck on a folder simply type the folder name instead of the tree. Unless given a direct path MonoCheck will check for the outgroup file and taxonomy file in the Working and Tree directory.

Example:

    python3 monocheck.py trees -clade Arthropoda Insecta Diptera -outgroup trees/outgroup.txt -email myemail@university.ac.uk

after running this command we will get output telling us the state of each of the three clades in each tree:

    trees/OG0000138_tree.txt>>> Arthropoda | polyphyletic | 3
    trees/OG0000138_tree.txt>>> Insecta | monophyletic | *
    trees/OG0000138_tree.txt>>> Diptera | polyphyletic | 1

The first column before the arrows is the tree being assessed. Then we have the clade followed by its state. The number following this is the number of taxa that prevent this clade being monophyletic. An Asterix means that the clade is monophyletic.

With so many trees viewing this information in the terminal is not very useful. Instead we can send it to an output file with the -out option. This will create a csv file containing the above information ready for exploration in R or Pandas.

|Tree|Clade|Status|NumberDisrupting|DisruptingTaxa|
|---|---|---|---|---|
|trees/OG0000206_tree.txt|Arthropoda|paraphyletic|2|"['Caenorhabditis_elegans', 'Paramacrobiotus_metropolitanus']"|
|trees/OG0000206_tree.txt|Insecta|monophyletic|2|"['Caenorhabditis_elegans', 'Paramacrobiotus_metropolitanus']"|
|trees/OG0000206_tree.txt|Diptera|paraphyletic|2|"['Caenorhabditis_elegans', 'Paramacrobiotus_metropolitanus']"|
|trees/OG0000215_tree.txt|Arthropoda|paraphyletic|2|"['Caenorhabditis_elegans', 'Paramacrobiotus_metropolitanus']"|
|trees/OG0000215_tree.txt|Insecta|monophyletic|0|[]|
|trees/OG0000215_tree.txt|Diptera|monophyletic|0|[]|
|...|...|...|...|...|

A small script included in the repo; analyse.py, shows me monophyly stats when the tree_stats.csv is input:

    Clade: Arthropoda, Monophyletic: 14, Total: 46, Percentage: 30.43%
    Clade: Insecta, Monophyletic: 38, Total: 46, Percentage: 82.61%
    Clade: Diptera, Monophyletic: 41, Total: 46, Percentage: 89.13%

### Error Checking

I have built error recovery into the program. If you are missing an input, or there is an issue with file structure the program should tell you. If there are any issues please feel free to contact me.

### Library References:
_Jaime Huerta-Cepas, Francois Serra and Peer Bork._ ETE 3: Reconstruction, analysis and visualization of phylogenomic data. **Mol Biol Evol 2016;**  [doi: 10.1093/molbev/msw046](http://mbe.oxfordjournals.org/content/early/2016/03/21/molbev.msw046 "link to citation reference")

Peter J. A. Cock, Tiago Antao, Jeffrey T. Chang, Brad A. Chapman, Cymon J. Cox, Andrew Dalke, Iddo Friedberg, Thomas Hamelryck, Frank Kauff, Bartek Wilczynski, Michiel J. L. de Hoon, Biopython: freely available Python tools for computational molecular biology and bioinformatics, _Bioinformatics_, Volume 25, Issue 11, June 2009, Pages 1422–1423, [https://doi.org/10.1093/bioinformatics/btp163](https://doi.org/10.1093/bioinformatics/btp163)

McKinney, W., & others. (2010). Data structures for statistical computing in python. In _Proceedings of the 9th Python in Science Conference_ (Vol. 445, pp. 51–56).
