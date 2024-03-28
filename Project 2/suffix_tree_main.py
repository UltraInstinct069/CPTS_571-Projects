
from suffix_tree_header import *
import time
import sys
import suffix_tree_API as sf
from node_header import *

if __name__ == "__main__":
    args= list(sys.argv)
    fasta_file=args[1]
    alphabet_file=args[2]

    # loading inputs
    seq=sf.load_sequences(fasta_file)
    sf.load_alphabet(alphabet_file)
    output_file='output_report.txt'
    output_file_bwt=fasta_file.replace('.','_')+'_BWT.txt'
    
    # printing input file information
    file_ptr=open(output_file, "a")
    sf.print_input_info(file_ptr,fasta_file,alphabet_file)
        
    # building the tree and printing tree statistics 
    start = time.time()
    tree = sf.suffix_tree()    
    end = time.time()
     
    sf.print_tree_report(file_ptr,end-start)
    file_ptr.close()
    
    # builiding bwt
    bwt_file=open(output_file_bwt, "w")
    sf.bwt(tree,bwt_file)
    bwt_file.close()
    print('Suffix tree construction complete : ',fasta_file)