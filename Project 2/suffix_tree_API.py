from numpy import * 
from Bio import SeqIO
from suffix_tree_header import *
from node_header import *
import sys

# loading sequences
def load_sequences(filename):
    sequences= list(SeqIO.parse(filename, "fasta"))
    global seq
    seq=sequences[0].seq

    seq=seq+'$'
    global sequence_len
    sequence_len=len(seq)
    return seq
    
# loading alphabets
def load_alphabet(file):
    f = open(file, "r")
    alphabet=f.read()
    alphabet = alphabet.replace(' ','')
    alphabet='$'+alphabet
    global alphabet_len
    alphabet_len = len(alphabet)

# building the suffix tree
def suffix_tree():
    global sequence_len
    # creating the root
    root = create_node(0, None, 0, 0, 0)
    root.suffix_link = root
    leaf = root
    for i in range(sequence_len):
        #inserting suffixes
        leaf = insert_node(i, root, leaf)
        if leaf is None:
            return None
    return root

# creating node
def create_node(id, parent, suffix_head, suffix_tail, stringDepth):
    node = Node()
    
    if suffix_head > suffix_tail:
        print("Invalid value appears fro head and tail ", id)
        exit(-1)
    node.id = id
    node.suffix_link = None
    node.parent = node if parent == None else parent
    node.suffix_head = suffix_head
    node.suffix_tail = suffix_tail
    node.num_child = 0
    node.childs = [None]*alphabet_len 
    node.depth = stringDepth
    
    return node


def insert_child(parent, child):
    global seq
    parent.childs[parent.num_child] = child
    parent.num_child += 1
    child.parent = parent

    for i in range(parent.num_child - 1):
        for j in range(parent.num_child - 1 - i):
            if parent.childs[j] == None or parent.childs[j + 1]== None:
                continue
            if seq[parent.childs[j].suffix_head] > seq[parent.childs[j + 1].suffix_head]:
                temp = parent.childs[j + 1]
                parent.childs[j + 1] = parent.childs[j]
                parent.childs[j] = temp
    return

def find_child(n, suffix, i):
    global seq
    cur_node = None
    
    for i in range(n.num_child):
        if n.num_child<=0:
            break
        cur_node = n.childs[i]
        if cur_node== None:
            continue        
        if seq[cur_node.suffix_head] == seq[suffix]:
            return i, cur_node
        
    return i, None

def create_internal_node(current, head, tail):
    global seq, internal_nodes, string_depth, max_depth,leafs,max_depth_node,sequence_len
    
    ch_num = 0
    ch_num,_= find_child(current.parent, head, ch_num)
    for i in range(current.suffix_head, tail + 1):
        if seq[i] != seq[head + (i - current.suffix_head)]:
            new_internal_node = create_node(sequence_len + internal_nodes + 1, current.parent, current.suffix_head, i - 1,(current.parent.depth + i) - (current.suffix_head))
            internal_nodes += 1
            string_depth += new_internal_node.depth

            if new_internal_node.depth > max_depth:
                max_depth = new_internal_node.depth
                max_depth_node = new_internal_node

            insert_child(new_internal_node, current)
            new_internal_node.parent.childs[ch_num] = new_internal_node
            current.suffix_head = new_internal_node.suffix_tail + 1
            newLeaf = create_node(leafs, new_internal_node, head + new_internal_node.suffix_tail - new_internal_node.suffix_head + 1, tail, tail - (head + new_internal_node.suffix_tail - new_internal_node.suffix_head) + new_internal_node.depth)
            leafs += 1
            insert_child(new_internal_node, newLeaf)
            return newLeaf

    print("Error in create_internal_node")
    exit(1)
    return current

def hop(n, head, tail):
    global seq
    num_child = 0
    i = 0
    num_child,node_child = find_child(n, head, num_child)
    if node_child == None:
        return n

    min_n = min((tail - head), (node_child.suffix_tail - node_child.suffix_head)) +1
    for i in range(min_n):
        if seq[head + i] != seq[node_child.suffix_head + i]:
            return n
    
    return hop(node_child, head + min_n, tail)

def find_path(v, head):
    global leafs,sequence_len
    child_num=0
    tail = sequence_len - 1
    hop_dis = hop(v, head, sequence_len - 1)
    num_hops = hop_dis.depth - v.depth
    child_num, child= find_child(hop_dis, head + num_hops,child_num)
    if child == None:
        child = create_node(leafs, hop_dis, head + num_hops, tail, hop_dis.depth + (tail - head) + 1)
        try:
            insert_child(hop_dis, child)
        except Exception as ex:
            print('Error in insert_child in find_path: ',ex)       
        leafs += 1
    else:
        try:
            child = create_internal_node(child, head + num_hops, tail)
        except Exception as ex:
            print('Error in create_internal_node in find_path: ',ex)
    return child

def node_hops (v_prime,u, beta_head, beta_end, suffix):
    global sequence_len
    r=0
    child_num=0
    beta_len = beta_end - beta_head + 1
    curr_node = v_prime
    e = None
    i = None
    v = None
    try:
        while r <= beta_len:
            child_num,e = find_child(curr_node, beta_head + r, child_num)
            if e:
                edge_len = e.suffix_tail - e.suffix_head + 1
                if edge_len + r > beta_len:
                    i = create_internal_node(e, suffix + curr_node.depth, sequence_len - 1)
                    v = i.parent
                    u.suffix_link = v
                    return i
                elif edge_len + r == beta_len:
                    v = e
                    u.suffix_link = v
                    k = u.depth
                    i = find_path(v, suffix + k - 1)
                    return i
                else:
                    r += edge_len
                    curr_node = e
                    continue
    except Exception as ex:
        print('Error in Node hops: ',ex)
    
    return i


# node u : suffix leaf i-1's parent
# node v : suffix link SL(u)
# node u': parent of u (if exists)
# node v': SL(u')  - may or may not be the parent of v

def insert_node(i, root, leaf):
    u = leaf.parent
    Case = -1

    if u.suffix_link is not None:
        if u != root:           
            Case = 0 # case IA
        else:                   
            Case = 1 # case IB
    elif u.parent != root:      
        Case = 2 # case IIA
    else:                       
        Case = 3 # case IIB

    if Case == 0: # case IA
        alpha_len = u.depth
        v = u.suffix_link
        return find_path(v, i + alpha_len - 1)
    elif Case == 1: # case IB
        return find_path(u, i)
    elif Case == 2:  # case IIA
        u_prime = u.parent
        beta_head = u.suffix_head
        beta_tail = u.suffix_tail
        if u_prime == None:
            v_prime=None
        else:
            v_prime = u_prime.suffix_link
        return node_hops(v_prime, u, beta_head, beta_tail, i)
    elif Case == 3: # case IIB
        u_prime = u.parent
        beta_head = u.suffix_head
        beta_tail = u.suffix_tail
        beta_prime_head = beta_head + 1
        if beta_tail == beta_head: 
            u.suffix_link = u_prime            
            return find_path(u_prime, i)
        else:
            return node_hops(u_prime, u, beta_prime_head, beta_tail, i)
    else:
        print("ERROR: Couldn't insert")
        exit(1)
    return 0  


def dfs(node):
    global dfs
    dfs += 1
    
    if dfs % 10 == 0: # Limiting tree depth in case of larger tree
        print("Depth: {}".format(node.depth))
    else:
        print("Depth: {}\t".format(node.depth), end="")
    for i in range(node.num_child):        
        dfs(node.children[i])
    
    return 0


def display_children(node):
    print(f"Showing children of Node ID: {node.id}:")
    for i in range(node.num_child):
        print(f"Child: {i}, id: {node.childs[i].id} ; ", end="")
    print()
    return 0

# building BWT
def bwt(node, ptr):
    global sequence_len
    if node== None:
        return
    if node.num_child == 0:
        value = node.id - 1
        ptr.write(f"{seq[value - 1] if value > 0 else seq[sequence_len - 1]}\r")
    else:
        for child in node.childs:
            bwt(child, ptr)

def print_tree_report(file_ptr,total_time):
    global internal_nodes,leafs,string_depth,max_depth
    
    file_ptr.write("\nTree Statistics:\n")
    file_ptr.write("----------------\n")
    file_ptr.write(f"Total time to generate tree: {(total_time)*1000} MSec\n")
    file_ptr.write(f"Implment constant: {sys.getsizeof(Node())} bytes\n")
    file_ptr.write(f"Total Tree Size: {sys.getsizeof(Node())*(internal_nodes+leafs)} bytes\n")
    file_ptr.write(f"Internal Nodes: {internal_nodes +1}\n")    
    file_ptr.write(f"Total Number of Nodes: {leafs-1 }\n")    
    file_ptr.write(f"Total Nodes: {internal_nodes + leafs}\n")    
    file_ptr.write(f"Average string-depth of an internal node: {string_depth / (internal_nodes+1)}\n")
    file_ptr.write(f"The string-depth of the deepest internal node: {max_depth}\n")    
    
    # printing Exact match report
    exact_matching_repeat_report(file_ptr)  

def exact_matching_repeat_report(file_ptr):
    global max_depth_node,max_depth
    n = max_depth_node.childs[0]
    file_ptr.write("\nExact Match Report:\n")
    file_ptr.write("---------------------\n")
    file_ptr.write("Exact Matching Repeat: ")
    # getting the exact match repeat from the deepest internal node
    for i in range(n.id - 1, max_depth_node.depth + n.id - 1):
        file_ptr.write(seq[i])
    
    file_ptr.write("\n")
    file_ptr.write(f"Length of longest exact matching repeat: {max_depth}\n")
    
    output_str='Exact Matching Repeat positions: '
    for i in range(max_depth_node.num_child):        
        output_str+=f"{max_depth_node.childs[i].id}, "

    file_ptr.write(output_str[:-2])    

    file_ptr.write("\n")
    file_ptr.write("______________________________________________________________________\n")
    file_ptr.write("______________________________________________________________________\n")

def print_input_info(file_ptr,fasta_file,alphabet_file):
    file_ptr.write(f"Input Information: \n")
    file_ptr.write("---------------------\n")
    file_ptr.write(f"Input: {fasta_file}\n")
    file_ptr.write(f"Input Characters: {len(seq)-1}\n")
    file_ptr.write(f"Alphabet: {alphabet_file}\n")
    


    
    