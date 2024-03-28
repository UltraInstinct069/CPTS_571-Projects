# node structure
from suffix_tree_header import *
class Node:
    def __init__(self):
        self.id = 0
        self.depth = 0
        self.num_child = 0
        self.suffix_head = 0
        self.suffix_tail = 0
        self.parent = None
        self.suffix_link = None
        self.childs = [None]*alphabet_len