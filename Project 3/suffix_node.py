class sf_node:
    def __init__(self, id, parent, suffix_link, parent_edge_label, parent_edge_label_end, alphabet_size, string_depth, full_string):
        self.id = id
        self.parent = parent
        self.suffix_link = suffix_link
        self.parent_edge_label = [parent_edge_label, parent_edge_label_end]
        self.children = [None] * alphabet_size
        self.string_depth = string_depth
        if id < len(full_string):
            color_end = -1
            color = -1
            #print('color_end: ',color_end)
            while id > color_end:
                #print('color_end: ', color_end)
                color_end = full_string.find('$', color_end + 1)
                color += 1
            self.color = color
        else:
            self.color = -1
