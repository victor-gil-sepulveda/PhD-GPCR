import sys
import tabulate
import os
from plopmetrics import processDir, filterRecords, genSingleTraj, genMetricsFile, genMetrics
from tools import load_dic_in_json, save_dic_in_json, get_best_clustering,\
create_dir, define_tableau20_cm
from pyproct.clustering.clustering import Clustering
import numpy
from prody.proteins.pdbfile import parsePDB
from prody.measure.contacts import Contacts
import matplotlib.pyplot as pl
import itertools

if __name__ == "__main__":

    RESULTS_PATH = os.path.join("Results", sys.argv[1], sys.argv[2])
    create_dir(RESULTS_PATH)

    #--------------------------------
    # Perform the filtering
    #--------------------------------
    PDB_FILE = os.path.join(sys.argv[1], "%s.pdb"%sys.argv[2])
    FILTERED_PDB_FILE = os.path.join(RESULTS_PATH,"%s.filtered.pdb"%(sys.argv[2]))
    METRICS_FILE = os.path.join(RESULTS_PATH,"%s.metrics.dat"%(sys.argv[2]))
    records = processDir(sys.argv[1], sys.argv[2])
    selection = filterRecords("'L1  Binding Ene' < -226 and 'L1  Binding Ene' > -424 and 'L1(15.360.555.4)' < 6.5 and 'L1(15.360.555.4)' > 1.5", records)
    genSingleTraj(FILTERED_PDB_FILE, records, selection)
#    genMetricsFile(METRICS_FILE, ["L1(15.360.555.4)","L1  Binding Ene"], selection)
    metrics = genMetrics(["L1(15.360.555.4)","L1  Binding Ene"], selection).T
  
    #--------------------------------  
    # Prepare the clustering for this guy
    #--------------------------------
    ## Load template and modify its contents for this case
    CLUSTERING_PATH = os.path.join(RESULTS_PATH,"%s_%s_clustering"%(sys.argv[1], sys.argv[2]))
    SCRIPT_PATH = os.path.join(RESULTS_PATH,"clustering.json")
    OUT_FILE = os.path.join(RESULTS_PATH, "clustering.out")
    script = load_dic_in_json(sys.argv[3])
    script["global"]["workspace"]["base"] = CLUSTERING_PATH
    script["data"]["files"].append(FILTERED_PDB_FILE)
    save_dic_in_json(script, SCRIPT_PATH)
    os.system("python -m pyproct.main %s > %s"%(SCRIPT_PATH, OUT_FILE))
#     os.system("python -m pyproct.main %s"%(SCRIPT_PATH))
    best_clustering = Clustering.from_dic(get_best_clustering(CLUSTERING_PATH)["clustering"])
     
    #--------------------------------
    # Now calculate the values
    #--------------------------------
    results = {}
    for cluster in best_clustering.clusters:
        energies = metrics[1][cluster.all_elements]
        distances = metrics[0][cluster.all_elements]
        results[cluster.id] = {}
        results[cluster.id]["max_energy"] = numpy.max(energies)
        results[cluster.id]["min_energy"] = numpy.min(energies)
        results[cluster.id]["mean_energy"] = numpy.mean(energies)
        results[cluster.id]["mean_distance"] = numpy.mean(distances)
        results[cluster.id]["population"] = len(cluster.all_elements)
        
    #--------------------------------
    # Tabulate!
    #--------------------------------
    TABLE_PATH = os.path.join(RESULTS_PATH,"table.rst")
    items = ["max_energy", "min_energy", "mean_energy", "mean_distance", "population"]
    headers = ["cluster id"]; headers.extend(items)
    table = []
    for cluster_id in results:
        row = [cluster_id]
        for item in items:
            row.append(results[cluster_id][item])
        table.append(row)
    open(TABLE_PATH,"w").write(tabulate.tabulate(table, headers, tablefmt="rst"))
    
    #--------------------------------
    # Get Residues < 2A of the ligand
    #--------------------------------
    CLUSTERS_PATH = os.path.join(CLUSTERING_PATH, "clusters")
    CONTACTS_FILE = os.path.join(RESULTS_PATH,"residue_contacts.txt")
    cluster_pdbs = os.listdir(CLUSTERS_PATH)
    residues_file = open(CONTACTS_FILE,"w")
    residues_per_id = {}
    for cluster_pdb in cluster_pdbs:
        if ".pdb" in cluster_pdb:
            cluster_id,_ = os.path.splitext(os.path.basename(cluster_pdb))
            struct = parsePDB(os.path.join(CLUSTERS_PATH,cluster_pdb))
            protein = struct.select("protein")
            ligand = struct.select("chain C")
            residues = []
            for i in range(protein.numCoordsets()):
                protein.setCoords(protein.getCoordsets()[i])
                ligand.setCoords(ligand.getCoordsets()[i])
                prot_contacts = Contacts(protein)
                contacts =  prot_contacts.select(2, ligand)
                
                # Get residue names
                for atom in iter(contacts):
                    residues.append("%s:%s:%s"%(atom.getResindex(),atom.getResname(),atom.getChid()))
            #residues_file.write("%s %s \n"%(cluster_id, " ".join(set(residues))))
            residues_file.write("%s %s \n"%(cluster_id, " ".join(residues)))
            residues_per_id[cluster_id] = residues
    residues_file.close()
    
    #--------------------------------
    # Plot distribution of the residues
    #--------------------------------
    labels = residues_per_id.keys()
    all_res_tags = []
    for label in labels:
        all_res_tags.extend(residues_per_id[label])
    all_res_tags = list(set(all_res_tags)); all_res_tags.sort()
    res_tag_id = dict(zip(all_res_tags, range(len(all_res_tags))))
    
    # bar data
    data = []
    for label in labels:
        row = [0]*len(all_res_tags)
        for tag in residues_per_id[label]:
            row[res_tag_id[tag]] += 1
        data.append(row)
    data = numpy.array(data)
    print data

    # plot data    
    fig, ax = pl.subplots()
    tab20cm = itertools.cycle(define_tableau20_cm())
    base = numpy.array([0]*len(all_res_tags)) 
    for i, row in enumerate(data):
        pl.bar(range(len(all_res_tags)), 
                  row, 
                  label = labels[i],
                  bottom = base,
                  color = next(tab20cm))
        base += row
    ax.set_xticklabels(all_res_tags ,rotation=45)
    
    pl.legend()
    pl.savefig(os.path.join(RESULTS_PATH,"histogram.svg"))
    
    
    