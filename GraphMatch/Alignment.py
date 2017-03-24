__author__ = 'laurencoombe'

class Alignment:

    def __init__(self, score, correspondingVertices, indelVertices, q_edges, r_edges):
        self.score = score
        self.correspVertices = correspondingVertices
        self.indelVertices = indelVertices
        self.query_edges = q_edges
        self.ref_edges = r_edges

