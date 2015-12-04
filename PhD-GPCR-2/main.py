import os
import tabulate
from plopmetrics import processDir, filterRecords, genSingleTrajFast, genMetricsFile, genMetrics
from tools import load_dic_in_json, save_dic_in_json, get_best_clustering,\
create_dir
from pyproct.clustering.clustering import Clustering
import numpy
from prody.proteins.pdbfile import parsePDB
from prody.measure.contacts import Contacts
from histogram import get_labels, filter_less_contacts_than,\
    plot_histogram, get_num_contacts_per_residue
from optparse import OptionParser

if __name__ == "__main__":

    # argv
    # 1 - directory where things are 
    # 2 - name of the pdb
    # 3 - template script for pyProCT
    # Changed to optparse
    
    parser = OptionParser()
    parser.add_option("--drug", dest="drug")
    parser.add_option("--protein", dest="protein")
    parser.add_option("--results", dest="results")
    parser.add_option("--template", dest="template")
    parser.add_option("--allosteric", action = "store_true", default= False, dest="do_allosteric")
    parser.add_option("--plot", action = "store_true", default= False, dest="do_plots")
    (options, args) = parser.parse_args()

    if options.drug is None or options.protein is None or options.template is None or options.results is None:  
        parser.error('Base directory, protein name, results directory and template are mandatory arguments.')

    RESULTS_PATH = os.path.join(options.results, options.drug, options.protein)
    create_dir(RESULTS_PATH)

    #--------------------------------
    # Perform the filtering
    #--------------------------------
    FILTERED_PDB_FILE = os.path.join(RESULTS_PATH,"%s.filtered.pdb"%(options.protein))
    METRICS_FILE = os.path.join(RESULTS_PATH,"%s.metrics.dat"%(options.protein))
    records = processDir(options.drug, options.protein)
    if not options.allosteric:
        selection = filterRecords("'L1  Binding Ene' < -226 and 'L1  Binding Ene' > -424 and 'L1(15.360.555.4)' < 6.5 and 'L1(15.360.555.4)' > 1.5", records)
        genSingleTrajFast(FILTERED_PDB_FILE, records, selection)
        genMetricsFile(METRICS_FILE, ["L1(15.360.555.4)","L1  Binding Ene"], selection)
        metrics = genMetrics(["L1(15.360.555.4)","L1  Binding Ene"], selection).T
    else: # range 6A - 14A
        selection = filterRecords("'L1  Binding Ene' < -226 and 'L1  Binding Ene' > -424 and 'L1(24.954.352.7)' < 14.1 and 'L1(24.954.352.7)' > 5.9", records)
        genSingleTrajFast(FILTERED_PDB_FILE, records, selection)
        genMetricsFile(METRICS_FILE, ["L1(24.954.352.7)","L1  Binding Ene"], selection)
        metrics = genMetrics(["L1(24.954.352.7)","L1  Binding Ene"], selection).T

    metrics = numpy.loadtxt(METRICS_FILE).T
    
    #--------------------------------  
    # Prepare the clustering for this guy
    #--------------------------------
    ## Load template and modify its contents for this case
    CLUSTERING_PATH = os.path.join(RESULTS_PATH,"%s_%s_clustering"%(options.drug, options.protein))
    MAX_CLUSTERS = 10
    SCRIPT_PATH = os.path.join(RESULTS_PATH,"clustering.json")
    OUT_FILE = os.path.join(RESULTS_PATH, "clustering.out")
    script = load_dic_in_json(options.template)
    script["global"]["workspace"]["base"] = CLUSTERING_PATH
    script["data"]["files"].append(FILTERED_PDB_FILE)
    script["clustering"]["evaluation"]["maximum_clusters"] = MAX_CLUSTERS
    save_dic_in_json(script, SCRIPT_PATH)
    os.system("python -m pyproct.main %s > %s"%(SCRIPT_PATH, OUT_FILE))
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
    CONTACTS_DISTANCE = 2.0

    cluster_pdbs = os.listdir(CLUSTERS_PATH)
    residues_file = open(CONTACTS_FILE,"w")
    contacts_per_cluster = {}
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
                
                contacts =  prot_contacts.select(CONTACTS_DISTANCE, ligand)
                
                if contacts is not None:
                    # Get residue names
                    for atom in iter(contacts):
                        residues.append("%s:%s:%s"%(atom.getResnum(), atom.getResname(), atom.getChid()))
            
            #residues_file.write("%s %s \n"%(cluster_id, " ".join(set(residues))))
            residues_file.write("%s %s \n"%(cluster_id, " ".join(residues)))
            contacts_per_cluster[cluster_id] = residues
    residues_file.close()
    
    if options.do_plots:
        #--------------------------------
        # Plot distribution of the residues
        #--------------------------------
        contacts_per_residue = get_num_contacts_per_residue(contacts_per_cluster)
        
        contact_residue_labels =  get_labels (contacts_per_cluster)
        
        # A normal plot
        target = os.path.join(RESULTS_PATH, "histogram.svg")
        plot_histogram(contacts_per_cluster, contact_residue_labels, target, False)
        
        # A plot averaging
        target = os.path.join(RESULTS_PATH, "histogram_a.svg")
        plot_histogram(contacts_per_cluster, contact_residue_labels, target, True)
    
        
        # A plot filtering
        filtered_contact_residue_labels =  filter_less_contacts_than(2000, contact_residue_labels, contacts_per_residue)
        target = os.path.join(RESULTS_PATH, "histogram_f.svg")
        plot_histogram(contacts_per_cluster, filtered_contact_residue_labels, target, False)
    
        # A plot filtering + averaging
        filtered_contact_residue_labels =  filter_less_contacts_than(2000, contact_residue_labels, contacts_per_residue)
        target = os.path.join(RESULTS_PATH, "histogram_fa.svg")
        plot_histogram(contacts_per_cluster, filtered_contact_residue_labels, target, True)
