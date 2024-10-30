from multiprocessing import Pool
from tqdm import tqdm

"""
Graphlet counting code isn't compatible with standard svmlight format
"""

class SVML_Reformatter:
    def __init__(self):
        self.graphlet_map = {}
        self.hash_set = set()
        self.n_graphlets = 0
    def translate_tup(self,tup):
        hash_val, count = tup.split(":")
        if hash_val in self.hash_set:
            return self.graphlet_map[hash_val],count
        self.hash_set.add(hash_val)
        self.graphlet_map[hash_val] = self.n_graphlets
        v = self.n_graphlets
        self.n_graphlets += 1
        return v,count
    def sort_line(self,parts):
        tups = sorted(parts[1:-1],key=lambda tup: tup[0])
        out = " ".join([parts[0],*[":".join([str(ti) for ti in tup]) for tup in tups],parts[-1]])
        return out
    def format_file(self,fp):
        line_lists = []
        with open(fp) as f:
            rawlines = f.readlines()
        for line in tqdm(rawlines,leave=False):
            parts = line.strip().split(" ")
            tups = [self.translate_tup(tup) for tup in parts[1:-1]]
            line_lists.append([parts[0],*tups, parts[-1]])
        outlines = []
        for l in tqdm(line_lists,leave=False):
            outlines.append(self.sort_line(l))
        return outlines
