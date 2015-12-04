'''
Created on Dec 3, 2015

@author: victor
'''
import glob
import os
import errno
from optparse import OptionParser
from prody.proteins.pdbfile import parsePDB, writePDB

def create_directory(directory_path, ensure_writability = False):
    """
    Creates a directory (with subdirs) if it doesn't exist.
    
    @param directory_path: the path of the directory and subdirectories to be created. 
    """
    if ensure_writability:
        if not os.access(os.path.dirname(directory_path), os.W_OK):
            return False
    try:
        os.makedirs(directory_path)
        return True
    except OSError, e:
        if e.errno != errno.EEXIST:
            raise
    return False

if __name__ == '__main__':
    """
    Given a folder or list of folders into a file, it creates the files needed to represent the clusters in a results
    folder. The pdbs have only the ligand. 
    """
    parser = OptionParser()
    parser.add_option("-i", dest="input")
    parser.add_option("-r", dest="results")
    parser.add_option("-o", dest="output")
    (options, args) = parser.parse_args()

    if options.input is None:
        parser.error('Input file is needed. Line format: "protein_name drug_name clusters_folder"')
    
    if options.output is None:
        parser.error('Output file is needed. The new path for clusters will be stored there.')
    
    if options.results is None:
        parser.error('You must define the results folder')
    
    output = open(options.output,"w")
    for line in open(options.input):
        protein, drug, folder = line.strip().split()
        files = glob.glob(os.path.join(folder, "cluster_*.pdb"))
        results_folder = os.path.join(options.results, drug, protein)
        create_directory(results_folder)
        output.write("%s %s %s\n"%(protein, drug, results_folder))
        for i, filename in enumerate(sorted(files)):
            pdb = parsePDB(filename)
            if i == 0:
                # Extract first frame
                prot = pdb.select("protein")
                writePDB(os.path.join(results_folder, "%s.pdb"%protein), prot, csets=[0])
            # Extract ligands
            ligands = pdb.select("resname %s"%drug)
            writePDB(os.path.join(results_folder, "ligand_%s"%os.path.basename(filename)), ligands)
    output.close()
    