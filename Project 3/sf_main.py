import time
import random
import sys
from Bio import SeqIO
import gc
import configparser
import cell
import alignment as align
import gen_suffix_tree as gst
import suffix_node as sn
import input_loader as loader

def task1(file_ptr,input, alphabet,total_colors,seq_alphas,file_names):
    file_ptr.write(f"\nTask 1) Detecting DNA fingerprints: ")
    sf_build_st = time.time()
    gen_suff_tree = gst.SuffixTree()
    gen_suff_tree.SuffixTree(alphabet, input, total_colors)
    sf_build_ed = time.time()
    file_ptr.write(f'\nGenralized suffix tree construction time: {sf_build_ed-sf_build_st} sec')
    try:
        gen_suff_tree.color_node(gen_suff_tree.root)
        print('Coloring Complete!')
    except Exception as e:
        print('Coloring: ',e)

    try:
        print('seq_alphas: ',seq_alphas)
        print('total_colors: ',total_colors)
        fp_st = time.time()
        gen_suff_tree.find_alpha(gen_suff_tree.root, seq_alphas)
        fp_ed = time.time() 
            
        for fi in range(gen_suff_tree.total_colors):
            gen_suff_tree.print_fingerprints(file_ptr,fi, seq_alphas,file_names[fi])
        
        print('find_alpha Complete')    
        file_ptr.write(f'\nFingerprint generation time: {fp_ed-fp_st} sec')

    except Exception as e:
        print('nFingerprint: ',e)            

def task2(file_ptr, seq_list,file_names):
    file_ptr.write(f"\n\nTask 2) Computing similarity matrix: ")

    num_seq = len(seq_list)
    similarity_matrix = [[0.0 for j in range(num_seq)] for i in range(num_seq)]

    sim_mat_st = time.time()
    for i in range(num_seq):
        similarity_matrix[i][i] = len(seq_list[i])

    for i in range(num_seq):
        for j in range(i + 1, num_seq):
            file_ptr.write(f"\nProcessing {file_names[i].replace('input_files/','')}  &  {file_names[j].replace('input_files/','')}")
            score_a = 0 # prefixes alignment score
            score_b = 0 # lcs score
            score_c = 0 # suffixed alignment score
            input_strings = seq_list[i] + dollar_sym + seq_list[j] + dollar_sym
    
            gen_t_cons_st = time.time()
            gen_suf_tree_lcs= gst.SuffixTree()
            gen_suf_tree_lcs.SuffixTree(alphabet, input_strings, 2)
            gen_t_cons_ed = time.time()

            gen_t_cons_t = gen_t_cons_ed-gen_t_cons_st
            #file_ptr.write("Genralized suffix tree construction time for ",file_names[i].replace('input_files/','') ,' & ',file_names[j].replace('input_files/',''),' : ', gen_t_cons_t , " sec.")
            file_ptr.write(f"\nGenralized suffix tree construction time for {file_names[i].replace('input_files/','') } & {file_names[j].replace('input_files/','')} : { gen_t_cons_t }  sec.")

            gen_suf_tree_lcs.color_node(gen_suf_tree_lcs.root)
            lcs = None
            lcs_start_s1 = 0
            lcs_start_s2 = 0
            lcs_end_ind = 0
            lcs,lcs_start_s1,lcs_start_s2,lcs_end_ind = gen_suf_tree_lcs.find_lcs(gen_suf_tree_lcs.root, lcs, lcs_start_s1, lcs_start_s2, lcs_end_ind, 2)

            score_b = len(seq_list[i][lcs_start_s1:lcs_start_s1+lcs_end_ind])
            file_ptr.write(f"\nLongest Common Sequence :  {score_b}  characters")

            lcs_start_s2 -= len(seq_list[i]) + 1

            align_time_st = time.time()
            if lcs_start_s1 == 0 or lcs_start_s2 == 0:
                score_a = 0
            else:
                s_i_rev = seq_list[i][0:lcs_start_s1][::-1]
                s_j_rev = seq_list[j][0:lcs_start_s2][::-1]

                rev_i_j = align.alignement("parameters.config")
                rev_i_j.alignement(s_i_rev, s_j_rev)
                rev_i_j.calculate_scores()
                score_a = rev_i_j.align_max_score
                rev_i_j.dp_table = [[]]
                del rev_i_j
                gc.collect()

            if lcs_start_s1 + lcs_end_ind == len(seq_list[i]) or lcs_start_s2 + lcs_end_ind == len(seq_list[j]):
                score_c = 0
            else:
                s_i_fwd = seq_list[i][lcs_start_s1+lcs_end_ind:len(seq_list[i])]
                s_j_fwd = seq_list[j][lcs_start_s2+lcs_end_ind:len(seq_list[j])]

                fwd_i_j = align.alignement("parameters.config")
                fwd_i_j.alignement(s_i_fwd, s_j_fwd)
                fwd_i_j.calculate_scores()
                score_c = fwd_i_j.align_max_score
                fwd_i_j.dp_table = [[]]
                del fwd_i_j
                gc.collect()

            align_time_ed = time.time()

            similarity_matrix[i][j] = score_a + score_b + score_c
            similarity_matrix[j][i] = score_a + score_b + score_c

            Align_time = align_time_ed-align_time_st
            # file_ptr.write("Total Alignment time for ",file_names[i].replace('input_files/','') ,' & ',file_names[j].replace('input_files/',''),' : ', Align_time , " Sec.")
            file_ptr.write(f"\nTotal Alignment time for {file_names[i].replace('input_files/','')} & {file_names[j].replace('input_files/','')} : {Align_time}  Sec.")
            del gen_suf_tree_lcs
            gc.collect()

    sim_mat_ed = time.time()
    file_ptr.write(f"\n\n Similarity Matrix generation Time: {sim_mat_ed-sim_mat_st}")
    file_ptr.write(f"\n\n{str(similarity_matrix)}")
    

if __name__ == "__main__":
    file_ptr=open('output_file.txt', "w")
    dollar_sym = '$'
    concat_str = ""
    input_alphabet = ""
    

    input_alphabet = 'input_files/DNA_alphabet.txt'
    alphabet = loader.load_alphabet(input_alphabet)
    file_names=['input_files/Covid_Wuhan.fasta','input_files/Covid_USA-CA4.fasta'
                ]
    file_names=['input_files/Covid_Wuhan.fasta','input_files/Covid_USA-CA4.fasta',
                'input_files/Covid_Australia.fasta','input_files/Covid_India.fasta',
                'input_files/Covid_Brazil.fasta','input_files/SARS_2017_MK062179.fasta',
                'input_files/SARS_2003_GU553363.fasta','input_files/MERS_2014_USA_KP223131.fasta',
                'input_files/MERS_2014_KY581694.fasta','input_files/MERS_2012_KF600620.fasta'
                ]
    seq_alphas = []
    seq_list = []
    for num_seq in range(len(file_names)):
        seq = loader.load_sequences(file_names[num_seq])
        seq_list.append(seq)
        concat_str += seq
        concat_str += dollar_sym
        print(float('inf'))

        node = sn.sf_node(-10, None, None, 0, 0, 5, 100000, "")
        seq_alphas.append(node)

    total_colors = len(file_names)

    strrep = str(concat_str) 
    try:
        #task 1
        t_task1_st=time.time()
        task1(file_ptr,strrep, alphabet,total_colors,seq_alphas,file_names)
        t_task1_ed=time.time()
        file_ptr.write(f"\n\nTotal Time for Finger Print task: {t_task1_ed-t_task1_st}")
        
        
        #task 2
        t_task2_st=time.time()
        task2(file_ptr,seq_list,file_names)
        t_task2_ed=time.time()
        file_ptr.write(f"\n\nTotal Time for Alignment task: {t_task2_ed-t_task2_st}")
    except Exception as e:
        print(e)
    
    file_ptr.close()
    











