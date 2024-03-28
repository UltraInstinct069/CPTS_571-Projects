#!/bin/bash
echo "s1.fas running"
python suffix_tree_main.py s1.fas English_alphabet.txt
echo "s2.fas running"
python suffix_tree_main.py s2.fas English_alphabet.txt
echo "colorblind_human_gene.fasta running"
python suffix_tree_main.py colorblind_human_gene.fasta DNA_alphabet.txt
echo "colorblind_mouse_gene.fasta running"
python suffix_tree_main.py colorblind_mouse_gene.fasta DNA_alphabet.txt
echo "Slyco.fas running"
python suffix_tree_main.py Slyco.fas DNA_alphabet.txt
echo "chr12.fas running"
python suffix_tree_main.py chr12.fas DNA_alphabet.txt
echo "Covid_Wuhan.fasta running"
python suffix_tree_main.py Covid_Wuhan.fasta DNA_alphabet.txt
echo "Covid_Brazil.fasta running"
python suffix_tree_main.py Covid_Brazil.fasta DNA_alphabet.txt
echo "Covid_India.fasta running"
python suffix_tree_main.py Covid_India.fasta DNA_alphabet.txt
echo "Covid_Australia.fasta running"
python suffix_tree_main.py Covid_Australia.fasta DNA_alphabet.txt
echo "Covid_USA-CA4.fasta running"
python suffix_tree_main.py Covid_USA-CA4.fasta DNA_alphabet.txt
