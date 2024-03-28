import configparser
import cell as dp_c
import gc

class alignement:
    def __init__(self,filename):
        self.dp_table = []

        self.align_max_score = 0
        self.local_align_indices = [0, 0]
        config_file = configparser.ConfigParser()
        config_file.read(filename)
        
        self.h=int(config_file['params']['h'])
        self.g=int(config_file['params']['g'])
        self.match=int(config_file['params']['match'])
        self.mismatch=int(config_file['params']['mismatch'])

    def alignement(self, s1, s2):
        self.s1 = s1
        self.s2 = s2
        self.len_s1 = len(self.s1)
        self.len_s2 = len(self.s2)
        # print('S1 length: ',self.len_s1)
        # print('S2 length: ',self.len_s2 )
        temp_table = [[dp_c.cell(float('-inf'), float('-inf'), float('-inf')) for j in range(self.len_s2 + 1)] for i in range(self.len_s1 + 1)]

        temp_table[0][0] = dp_c.cell(float('-inf'), float('-inf'), 0)
        for i in range(1, self.len_s1 + 1):
            temp_table[i][0] = dp_c.cell(float('-inf'), (self.h + i * self.g), float('-inf'))
        for j in range(1, self.len_s2 + 1):
            temp_table[0][j] = dp_c.cell((self.h + j * self.g), float('-inf'), float('-inf'))

        self.dp_table = temp_table
        del temp_table
        gc.collect()

    def actions(self, i, d, s):
        if s >= max([i, d]):
            return 'D'
        elif i >= max([s, d]):
            return 'I'
        else:
            return 'D'
        
    def calculate_scores(self):
        # print('S1 length: ',self.len_s1)
        # print('S2 length: ',self.len_s2)
        for i in range(1, self.len_s1 + 1):
            for j in range(1, self.len_s2 + 1):
                align_dir = ''
                self.dp_table[i][j].score_subs = max(self.dp_table[i - 1][j - 1].score_d, self.dp_table[i - 1][j - 1].score_i, self.dp_table[i - 1][j - 1].score_subs)
                if self.s1[i - 1] == self.s2[j - 1]:
                    self.dp_table[i][j].score_subs += self.match
                else:
                    self.dp_table[i][j].score_subs += self.mismatch
                
                self.dp_table[i][j].score_d = max(self.dp_table[i - 1][j].score_d + self.g, self.dp_table[i - 1][j].score_i + self.g + self.h, self.dp_table[i - 1][j].score_subs + self.g + self.h)
                
                self.dp_table[i][j].score_i = max(self.dp_table[i][j - 1].score_i + self.g, self.dp_table[i][j - 1].score_d + self.g + self.h, self.dp_table[i][j - 1].score_subs + self.g + self.h)
                
                self.dp_table[i][j].max_score()
                
                if self.dp_table[i][j].score >= self.align_max_score:
                    self.align_max_score = self.dp_table[i][j].score
                    self.local_align_indices[0] = i
                    self.local_align_indices[1] = j
                
                align_dir = self.actions(self.dp_table[i][j].score_i, self.dp_table[i][j].score_d, self.dp_table[i][j].score_subs)
                if align_dir == 'D':
                    self.dp_table[i][j].m_score = self.dp_table[i - 1][j].m_score
                elif align_dir == 'I':
                    self.dp_table[i][j].m_score = self.dp_table[i][j - 1].m_score
                elif align_dir == 'S':
                    if self.s1[i - 1] == self.s2[j - 1]:
                        self.dp_table[i][j].m_score = self.dp_table[i - 1][j - 1].m_score + self.match