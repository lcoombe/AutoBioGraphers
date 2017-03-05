__author__ = 'laurencoombe'

class Corr_Vertex:

    def __init__(self, name, score, queryVertex):
        self.name = name
        self.score = float(score)
        self.queryVertex = queryVertex

    def __str__(self):
        return self.name + ", " + str(self.score) + ", " + str(self.queryVertex)