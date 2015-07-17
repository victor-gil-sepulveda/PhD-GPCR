import json
from pyproct.tools.commonTools import convert_to_utf8
import numpy
import os
import math


PYPROCT = "/home/victor/workspaces/Python/pyProClust/pyproct/main.py"

def load_dic_in_json(filename):
    return convert_to_utf8(json.loads("".join(open(filename,"r").readlines())))

def save_dic_in_json(this, there):
    open(there, "w").write(json.dumps(this, sort_keys=False, indent=4, separators=(',', ': ')))

def create_dir(_dir):
    if not os.path.exists(_dir):
        os.makedirs(_dir)

def use_pyproct(directory, script):
    create_dir(directory)
    os.system("python -m pyproct.main %s "%(script))

def get_best_clustering(clustering_directory):
    results = load_dic_in_json(os.path.join(clustering_directory,"results","results.json"))
    best_clustering_id = results['best_clustering']
    return results['selected'][best_clustering_id]

def score_cluster(cluster, matrix, metrics):
    all_rmsd_distances = [metrics[element][0] for element in cluster.all_elements]
    all_binding_energies = [metrics[element][1] for element in cluster.all_elements]
    return numpy.mean(all_rmsd_distances),numpy.mean(all_binding_energies)

def find_most_negative_be_half(scores):
    # Retrieve all scores and do the mean
    be_scores_mean = numpy.mean([x[0][1] for x in scores])
    filtered_scores = []
    for score in scores:
        rmsd,be = score[0]
        if be <= be_scores_mean:
            filtered_scores.append(score)
    return filtered_scores

def find_most_negative_rmsd_half(scores):
    # Retrieve all scores and do the mean
    rmsd_scores_mean = numpy.mean([x[0][0] for x in scores])
    filtered_scores = []
    for score in scores:
        rmsd,be = score[0]
        if rmsd <= rmsd_scores_mean:
            filtered_scores.append(score)
    return filtered_scores

def find_most_negative_cluster(scores):
    filtered = scores
    while len(filtered)>1:
        filtered = find_most_negative_be_half(filtered)
        filtered = find_most_negative_rmsd_half(filtered)
    return filtered[0]

def find_5_clusters_with_less_energy(scores):
    # Mean energy will be the second element of each score
    energies = [(en,cluster) for ((dist, en),cluster) in scores]  # @UnusedVariable
    energies.sort()
    clusters = []
    for i in range(5):
        clusters.append(energies[i])
    return clusters

def normalize_metrics(metrics):
    normalized_metrics = []
    metricsT = metrics.T
    for metric_array in metricsT:
        m_max = numpy.max(metric_array)
        m_min = numpy.min(metric_array)
        normalized_metrics.append( (metric_array - m_min) / (m_max - m_min))
    return numpy.array(normalized_metrics).T

def distance(first,second):
    r = first - second
    return math.sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])

def define_tableau20_cm ():
    """
    http://www.randalolson.com/2014/06/28/how-to-make-beautiful-data-visualizations-in-python-with-matplotlib/
    """
    # These are the "Tableau 20" colors as RGB.  
    tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
                 (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
                 (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
                 (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
                 (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]  
      
    # Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.  
    for i in range(len(tableau20)):  
        r, g, b = tableau20[i]  
        tableau20[i] = (r / 255., g / 255., b / 255.) 
    
    return tableau20
