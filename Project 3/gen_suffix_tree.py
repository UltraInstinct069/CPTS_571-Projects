import suffix_node as sn

class SuffixTree:
    def __init__(self):
        self.root = None
        self.join_string = ""
        self.alpha = ""
        self.alpha_map = {}
        self.in_id = 0
        self.total_colors = 0

    def SuffixTree(self, alpha, Seq, colors):
        alpha = '$' + alpha
        self.alpha_map = {'$': 0}
        for i in range(1, len(alpha)):
            self.alpha_map[alpha[i]] = i
        self.total_colors = colors
        self.join_string = Seq
        self.alpha = alpha
        self.in_id = len(Seq)
        self.build_gen_suffix_tree()
    
    def build_gen_suffix_tree(self):
        self.root = sn.sf_node(self.in_id, None, None, 0, 0, len(self.alpha), 0, self.join_string)
        self.in_id+=1
        self.root.parent = self.root
        self.root.color = self.total_colors
        for i in range(len(self.join_string)):
            self.find_path(self.root, i, i)
    
        return self.root
    
    def find_path(self,u, string_index, id):
        child_index = self.alpha_map[self.join_string[string_index]]
        if u.children[child_index] == None:
            string_depth_l = u.string_depth + self.join_string.find('$', string_index) - string_index + 1
            cur_node = sn.sf_node(id, u, None, string_index, self.join_string.find('$', string_index),
                    len(self.alpha), string_depth_l, self.join_string)
            u.children[child_index] = cur_node
        else:
            v = u.children[child_index]
            for i in range(v.parent_edge_label[0], v.parent_edge_label[1] + 1):
                f_df = i - v.parent_edge_label[0]
                n_df = string_index + f_df
                a_df = len(self.join_string) - string_index
                a_df_offset = len(self.join_string)
                if self.join_string[n_df] != self.join_string[i]:
                    internal_node = sn.sf_node(self.in_id, u, None,
                        v.parent_edge_label[0], i - 1, len(self.alpha), u.string_depth + f_df, self.join_string)
                    self.in_id+=1
                    if id > self.join_string.find('$', 0):
                        a_df -= a_df_offset
                        a_df += self.join_string.find('$', id)+1
                    leafNode = sn.sf_node(id, internal_node, None,
                        n_df, self.join_string.find('$', string_index), len(self.alpha), u.string_depth + a_df, self.join_string)
                    leafInd = self.alpha_map[self.join_string[n_df]]
                    vInd = self.alpha_map[self.join_string[i]]

                    internal_node.children[leafInd] = leafNode
                    internal_node.children[vInd] = v
                    v.parent_edge_label[0] = i
                    v.parent = internal_node
                    u.children[child_index] = internal_node
                    return leafNode
            if v.id < len(self.join_string):
                v.color = self.total_colors
                v.id = id
                return v
            else:
                return self.find_path(v, (string_index + (v.parent_edge_label[1] - v.parent_edge_label[0] + 1)), id)

    def color_node(self,Node):
        for i in range(len(self.alpha)):
            if Node.children[i] != None:
                self.color_node(Node.children[i])
        if Node.id > len(self.join_string):
            child_colors = -1
            for i in range(len(self.alpha)):
                if Node.children[i] != None:
                    if child_colors == -1:
                        child_colors = Node.children[i].color
                    elif child_colors != Node.children[i].color:
                        child_colors = self.total_colors
            Node.color = child_colors

    def find_alpha(self,Node, a_alpha):
        if Node.color == self.total_colors:
            for i in range(1, len(self.alpha)):
                if Node.children[i] is not None and Node.children[i].color != self.total_colors:
                    if a_alpha[Node.children[i].color].string_depth > Node.string_depth:
                        a_alpha[Node.children[i].color] = Node
        for i in range(1, len(self.alpha)):
            if Node.children[i] is not None:
                self.find_alpha(Node.children[i], a_alpha)
    def print_fingerprints(self,file_ptr,color, alphas,file_name):
        c = 0
        if alphas[color].id == -10:
            file_ptr.write(f"\n{file_name.replace('input_files/','') } 's fingerprint: NONE")
            return
        for i in range(1, len(self.alpha)):
            if alphas[color].children[i] != None and alphas[color].children[i].color == color:
                c = alphas[color].children[i].parent_edge_label[0]
                break
        alpha = self.join_string[c - alphas[color].string_depth: 1 +c] 
        file_ptr.write(f"\n{file_name.replace('input_files/','') } 's fingerprint: {alpha}")
        file_ptr.write(f"\nFingerprint Length:  {str(len(alpha))}")

        
    def find_lcs(self,Node, result_node, start_seq1, start_seq2, lcs_len, num_colors):
        ci = 0
        if Node.id > len(self.join_string):
            if Node.color == num_colors:
                if lcs_len < Node.string_depth:
                    lcs_len = Node.string_depth
                    while ci < len(self.alpha):
                        if Node.children[ci] != None and Node.children[ci].color == 0:
                            start_seq1 = Node.children[ci].parent_edge_label[0] - lcs_len
                        if Node.children[ci] != None and Node.children[ci].color == 1:
                            start_seq2 = Node.children[ci].parent_edge_label[0] - lcs_len
                        ci += 1
                    result_node = Node
        for i in range(len(self.alpha)):
            if Node.children[i] != None:
                result_node,start_seq1,start_seq2,lcs_len=self.find_lcs(Node.children[i], result_node, start_seq1, start_seq2, lcs_len, num_colors)

        return result_node,start_seq1,start_seq2,lcs_len