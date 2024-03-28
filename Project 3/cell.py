
class cell:
    def __init__(self, i, d, s):
        self.score_i = i # insertion Score
        self.score_d = d # deletion Score
        self.score_subs = s # substituion score
        self.score = max(self.score_d, self.score_i, self.score_subs)
        self.m_score=0

    def max_score(self):
        self.score = max(self.score_d, self.score_i, self.score_subs)