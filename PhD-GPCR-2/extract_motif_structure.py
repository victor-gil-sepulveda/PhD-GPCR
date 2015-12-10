'''
Created on Dec 10, 2015

@author: victor
'''
from prody.proteins.pdbfile import parsePDB, writePDB
import sys
from histogram import parse_motifs
import urllib

if __name__ == '__main__':
    pdb = parsePDB(sys.argv[1])
    motifs = parse_motifs(sys.argv[2])
    prot_name = sys.argv[3]
    for motif in motifs[prot_name]:
        motif_struct = pdb.select("resid %d to %d"%motifs[prot_name][motif])
        writePDB("%s.pdb"%urllib.quote(motif, safe=""),motif_struct)
        print motif, sorted(list(set(motif_struct.getResnums())))
