#!/bin/bash

echo "----------------- Running through whole program example -----------------"

echo "Running MonoCheck on the first single tree with command:"

echo "$ python3 monocheck.py test_tree.tre -clade Insecta -outgroup outgroup.txt -email myemail@university.ac.uk"

python3 monocheck.py test_tree.tre -clade Insecta -outgroup outgroup.txt -email myemail@university.ac.uk > monocheck_output_test_tree.log

echo "Output saved in monocheck_output_test_tree.txt"

echo "Running MonoCheck on all trees in the trees directory with command:"

echo "$ python3 monocheck.py trees -clade Arthropoda Insecta Diptera -outgroup trees/outgroup.txt -email myemail@university.ac.uk"

python3 monocheck.py trees -clade Arthropoda Insecta Diptera -outgroup trees/outgroup.txt -email myemail@university.ac.uk -out tree_stats.csv > monocheck_output_trees.log

echo "Output saved in monocheck_output_trees.txt"

echo "Running analyse.py on the output of the previous command: "

python3 analyse.py tree_stats.csv

echo "Tests complete"

echo "----------------- End of run -----------------"

#EOF