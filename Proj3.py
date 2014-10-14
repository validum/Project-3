# -*- coding: utf-8 -*-
"""
Created on Fri Oct 03 21:31:31 2014

@author: Eric

Template for Project 3
Student will implement four functions:

slow_closest_pairs(cluster_list)
fast_closest_pair(cluster_list) - implement fast_helper()
hierarchical_clustering(cluster_list, num_clusters)
kmeans_clustering(cluster_list, num_clusters, num_iterations)

where cluster_list is a list of clusters in the plane
"""

import math
#import copy
import Provided_Code_Proj3 as alg_cluster
#import alg_cluster

def pair_distance(cluster_list, idx1, idx2):
    """
    Helper function to compute Euclidean distance between two clusters
    in cluster_list with indices idx1 and idx2
    
    Returns tuple (dist, idx1, idx2) with idx1 < idx2 where dist is distance between
    cluster_list[idx1] and cluster_list[idx2]
    """
    return (cluster_list[idx1].distance(cluster_list[idx2]), min(idx1, idx2), max(idx1, idx2))


def slow_closest_pairs(cluster_list):
    """
    Compute the set of closest pairs of cluster in list of clusters
    using O(n^2) all pairs algorithm
    
    Returns the set of all tuples of the form (dist, idx1, idx2) 
    where the cluster_list[idx1] and cluster_list[idx2] have minimum distance dist.   
    
    """
    shortest_dist = float("inf")
    set_pairs = set([])
    for idx1 in range(len(cluster_list)):
        for idx2 in range(len(cluster_list)):
            if idx2 != idx1:
                pair = pair_distance(cluster_list, idx2, idx1)
                if pair[0] == shortest_dist:
                    set_pairs.add(pair)
                elif pair[0] < shortest_dist:
                    shortest_dist = pair[0]
                    set_pairs.clear()
                    set_pairs.add(pair)
    return set_pairs


def fast_closest_pair(cluster_list):
    """
    Compute a closest pair of clusters in cluster_list
    using O(n log(n)) divide and conquer algorithm
    
    Returns a tuple (distance, idx1, idx2) with idx1 < idx 2 where
    cluster_list[idx1] and cluster_list[idx2]
    have the smallest distance dist of any pair of clusters
    """
        
    def fast_helper(cluster_list, horiz_order, vert_order):
        """
        Divide and conquer method for computing distance between closest pair of points
        Running time is O(n * log(n))
        
        horiz_order and vert_order are lists of indices for clusters
        ordered horizontally and vertically
        
        Returns a tuple (distance, idx1, idx2) with idx1 < idx 2 where
        cluster_list[idx1] and cluster_list[idx2]
        have the smallest distance dist of any pair of clusters
    
        """
        
        # base case
        if len(horiz_order) <= 3:
            dist = (float("inf"), -1, -1)
            for clust1 in horiz_order:
                for clust2 in horiz_order:
                    if clust1 != clust2:
                        dist = min(dist, pair_distance(cluster_list, clust1, clust2))
            return dist
        
        # divide
        else:
            half_num_clust = int(math.ceil(len(horiz_order)/2))
            mid_horiz = 0.5*(cluster_list[horiz_order[half_num_clust-1]].horiz_center() + cluster_list[horiz_order[half_num_clust]].horiz_center())
            horiz_order_left = horiz_order[:half_num_clust]
            horiz_order_right = horiz_order[half_num_clust:]
            vert_order_left = []
            vert_order_right = []
            horiz_left_set = set(horiz_order_left)
            for cluster_indx in vert_order:
                 if cluster_indx in horiz_left_set:
                     vert_order_left.append(cluster_indx)
                 else:
                     vert_order_right.append(cluster_indx)
            min_of_pair = min(fast_helper(cluster_list, horiz_order_left, vert_order_left), fast_helper(cluster_list, horiz_order_right, vert_order_right))

        # conquer
            vlist = []
            for clusterv in vert_order:
                if abs(cluster_list[clusterv].horiz_center() - mid_horiz) <= min_of_pair[0]:
                    vlist.append(clusterv)
        
            for indx1 in range(len(vlist)-1):
                for indx2 in range(indx1+1, min(indx1+3, len(vlist))):
                    min_of_pair = min(min_of_pair, pair_distance(cluster_list, vlist[indx1], vlist[indx2]))

        return min_of_pair
            
    # compute list of indices for the clusters ordered in the horizontal direction
    hcoord_and_index = [(cluster_list[idx].horiz_center(), idx) 
                        for idx in range(len(cluster_list))]    
    hcoord_and_index.sort()
    horiz_order = [hcoord_and_index[idx][1] for idx in range(len(hcoord_and_index))]
     
    # compute list of indices for the clusters ordered in vertical direction
    vcoord_and_index = [(cluster_list[idx].vert_center(), idx) 
                        for idx in range(len(cluster_list))]    
    vcoord_and_index.sort()
    vert_order = [vcoord_and_index[idx][1] for idx in range(len(vcoord_and_index))]

    # compute answer recursively
    answer = fast_helper(cluster_list, horiz_order, vert_order) 
    return (answer[0], min(answer[1:]), max(answer[1:]))

    

def hierarchical_clustering(cluster_list, num_clusters):
    """
    Compute a hierarchical clustering of a set of clusters
    Note: the function mutates cluster_list
    
    Input: List of clusters, number of clusters
    Output: List of clusters whose length is num_clusters
    """
    while len(cluster_list) > num_clusters:
        closest_pair = fast_closest_pair(cluster_list)
        cluster_list[closest_pair[1]].merge_clusters(cluster_list[closest_pair[2]])
        cluster_list.pop(closest_pair[2])
    return cluster_list



    
def kmeans_clustering(cluster_list, num_clusters, num_iterations):
    """
    Compute the k-means clustering of a set of clusters
    
    Input: List of clusters, number of clusters, number of iterations
    Output: List of clusters whose length is num_clusters
    """
    
    # initialize k-means clusters to be initial clusters with largest populations
    cluster_list_sorted = sorted(cluster_list, key = lambda x: x.total_population(), reverse=True)
    num_cluster_list = [alg_cluster.Cluster(set([]), cluster_list_sorted[x].horiz_center(), cluster_list_sorted[x].vert_center(), 0, 0) for x in range(num_clusters)]
    num_cluster_list_copy = [alg_cluster.Cluster(set([]), cluster_list[x].horiz_center(), cluster_list[x].vert_center(), 0, 0) for x in range(num_clusters)]
    for iteration in range(num_iterations):
        for idx in range(len(num_cluster_list)):
            num_cluster_list_copy[idx] = num_cluster_list[idx].copy()
        
        for cluster in cluster_list:
            dist = float("inf")
            kindex = 0
            for ncluster in num_cluster_list:
                if ncluster.distance(cluster) < dist:
                    dist = ncluster.distance(cluster)
                    kindex = num_cluster_list.index(ncluster)
            num_cluster_list_copy[kindex].merge_clusters(cluster)
        
        num_cluster_list = [alg_cluster.Cluster(set([]), x.horiz_center(), x.vert_center(), 0, 0) for x in num_cluster_list_copy]

    return num_cluster_list_copy




#print fast_closest_pair([alg_cluster.Cluster(set([]), 0, 0, 1, 0), alg_cluster.Cluster(set([]), 1, 0, 1, 0), alg_cluster.Cluster(set([]), 2, 0, 1, 0), alg_cluster.Cluster(set([]), 3, 0, 1, 0), alg_cluster.Cluster(set([]), 4, 0, 1, 0), alg_cluster.Cluster(set([]), 5, 0, 1, 0), alg_cluster.Cluster(set([]), 6, 0, 1, 0), alg_cluster.Cluster(set([]), 7, 0, 1, 0), alg_cluster.Cluster(set([]), 8, 0, 1, 0), alg_cluster.Cluster(set([]), 9, 0, 1, 0), alg_cluster.Cluster(set([]), 10, 0, 1, 0), alg_cluster.Cluster(set([]), 11, 0, 1, 0), alg_cluster.Cluster(set([]), 12, 0, 1, 0), alg_cluster.Cluster(set([]), 13, 0, 1, 0), alg_cluster.Cluster(set([]), 14, 0, 1, 0), alg_cluster.Cluster(set([]), 15, 0, 1, 0), alg_cluster.Cluster(set([]), 16, 0, 1, 0), alg_cluster.Cluster(set([]), 17, 0, 1, 0), alg_cluster.Cluster(set([]), 18, 0, 1, 0), alg_cluster.Cluster(set([]), 19, 0, 1, 0)])

#print slow_closest_pairs([alg_cluster.Cluster(set(['01101']), 720.281573781, 440.436162917, 223510, 5.7e-05), alg_cluster.Cluster(set(['01117']), 709.193528999, 417.394467797, 143293, 5.6e-05), alg_cluster.Cluster(set(['01073']), 704.191210749, 411.014665198, 662047, 7.3e-05), alg_cluster.Cluster(set(['01015']), 723.907941153, 403.837487318, 112249, 5.6e-05), alg_cluster.Cluster(set(['01113']), 740.385154867, 436.939588695, 49756, 5.6e-05), alg_cluster.Cluster(set(['04013']), 214.128077618, 396.893960776, 3072149, 6.8e-05), alg_cluster.Cluster(set(['06025']), 156.397958859, 393.161127277, 142361, 5.6e-05), alg_cluster.Cluster(set(['06065']), 146.410389633, 374.21707964, 1545387, 6.1e-05), alg_cluster.Cluster(set(['06073']), 129.2075529, 387.064888184, 2813833, 6.6e-05), alg_cluster.Cluster(set(['06059']), 113.997715586, 368.503452566, 2846289, 9.8e-05), alg_cluster.Cluster(set(['06037']), 105.369854549, 359.050126004, 9519338, 0.00011), alg_cluster.Cluster(set(['06111']), 93.4973310868, 344.590570899, 753197, 5.8e-05), alg_cluster.Cluster(set(['06083']), 76.0382837186, 340.420376302, 399347, 6.4e-05), alg_cluster.Cluster(set(['06029']), 103.787886113, 326.006585349, 661645, 9.7e-05), alg_cluster.Cluster(set(['06071']), 148.402461892, 350.061039619, 1709434, 7.7e-05), alg_cluster.Cluster(set(['06107']), 108.085024898, 306.351832438, 368021, 5.8e-05), alg_cluster.Cluster(set(['06039']), 97.2145136451, 278.975077449, 123109, 6e-05), alg_cluster.Cluster(set(['06019']), 95.6093812211, 290.162708843, 799407, 6.4e-05), alg_cluster.Cluster(set(['06081']), 52.6171444847, 262.707477827, 707161, 5.6e-05), alg_cluster.Cluster(set(['06001']), 61.782098866, 259.312457296, 1443741, 7e-05), alg_cluster.Cluster(set(['06085']), 63.1509653633, 270.516712105, 1682585, 6.3e-05), alg_cluster.Cluster(set(['06067']), 74.3547338322, 245.49501455, 1223499, 6.1e-05), alg_cluster.Cluster(set(['06075']), 52.7404001225, 254.517429395, 776733, 8.4e-05), alg_cluster.Cluster(set(['06113']), 68.2602083189, 236.862609218, 168660, 5.9e-05), alg_cluster.Cluster(set(['06101']), 74.2003718491, 229.646592975, 78930, 5.6e-05), alg_cluster.Cluster(set(['06021']), 65.2043358182, 213.245337355, 26453, 6.9e-05), alg_cluster.Cluster(set(['06089']), 77.359494209, 188.945068958, 163256, 5.7e-05), alg_cluster.Cluster(set(['08005']), 380.281283151, 270.268826873, 487967, 5.9e-05), alg_cluster.Cluster(set(['08001']), 379.950978294, 265.078784954, 363857, 6.6e-05), alg_cluster.Cluster(set(['08031']), 371.038986573, 266.847932979, 554636, 7.9e-05), alg_cluster.Cluster(set(['09003']), 925.917212741, 177.152290276, 857183, 5.7e-05), alg_cluster.Cluster(set(['12073']), 762.463896365, 477.365342219, 239452, 6.1e-05), alg_cluster.Cluster(set(['13313']), 737.308367745, 378.040993858, 83525, 5.6e-05), alg_cluster.Cluster(set(['13215']), 745.265661102, 430.987078939, 186291, 5.9e-05), alg_cluster.Cluster(set(['13067']), 747.238620236, 397.293799252, 607751, 6.4e-05), alg_cluster.Cluster(set(['13063']), 752.853876848, 406.722877803, 236517, 6.6e-05), alg_cluster.Cluster(set(['13151']), 756.589546538, 407.288873768, 119341, 5.6e-05), alg_cluster.Cluster(set(['13089']), 754.465443436, 400.059456026, 665865, 6.8e-05), alg_cluster.Cluster(set(['13121']), 750.160287596, 399.907752014, 816006, 7e-05), alg_cluster.Cluster(set(['13135']), 758.038826857, 395.110327675, 588448, 6.3e-05)])

#print fast_closest_pair([alg_cluster.Cluster(set(['01101']), 720.281573781, 440.436162917, 223510, 5.7e-05), alg_cluster.Cluster(set(['01117']), 709.193528999, 417.394467797, 143293, 5.6e-05), alg_cluster.Cluster(set(['01073']), 704.191210749, 411.014665198, 662047, 7.3e-05), alg_cluster.Cluster(set(['01015']), 723.907941153, 403.837487318, 112249, 5.6e-05), alg_cluster.Cluster(set(['01113']), 740.385154867, 436.939588695, 49756, 5.6e-05), alg_cluster.Cluster(set(['04013']), 214.128077618, 396.893960776, 3072149, 6.8e-05), alg_cluster.Cluster(set(['06025']), 156.397958859, 393.161127277, 142361, 5.6e-05), alg_cluster.Cluster(set(['06065']), 146.410389633, 374.21707964, 1545387, 6.1e-05), alg_cluster.Cluster(set(['06073']), 129.2075529, 387.064888184, 2813833, 6.6e-05), alg_cluster.Cluster(set(['06059']), 113.997715586, 368.503452566, 2846289, 9.8e-05), alg_cluster.Cluster(set(['06037']), 105.369854549, 359.050126004, 9519338, 0.00011), alg_cluster.Cluster(set(['06111']), 93.4973310868, 344.590570899, 753197, 5.8e-05), alg_cluster.Cluster(set(['06083']), 76.0382837186, 340.420376302, 399347, 6.4e-05), alg_cluster.Cluster(set(['06029']), 103.787886113, 326.006585349, 661645, 9.7e-05), alg_cluster.Cluster(set(['06071']), 148.402461892, 350.061039619, 1709434, 7.7e-05), alg_cluster.Cluster(set(['06107']), 108.085024898, 306.351832438, 368021, 5.8e-05), alg_cluster.Cluster(set(['06039']), 97.2145136451, 278.975077449, 123109, 6e-05), alg_cluster.Cluster(set(['06019']), 95.6093812211, 290.162708843, 799407, 6.4e-05), alg_cluster.Cluster(set(['06081']), 52.6171444847, 262.707477827, 707161, 5.6e-05), alg_cluster.Cluster(set(['06001']), 61.782098866, 259.312457296, 1443741, 7e-05), alg_cluster.Cluster(set(['06085']), 63.1509653633, 270.516712105, 1682585, 6.3e-05), alg_cluster.Cluster(set(['06067']), 74.3547338322, 245.49501455, 1223499, 6.1e-05), alg_cluster.Cluster(set(['06075']), 52.7404001225, 254.517429395, 776733, 8.4e-05), alg_cluster.Cluster(set(['06113']), 68.2602083189, 236.862609218, 168660, 5.9e-05), alg_cluster.Cluster(set(['06101']), 74.2003718491, 229.646592975, 78930, 5.6e-05), alg_cluster.Cluster(set(['06021']), 65.2043358182, 213.245337355, 26453, 6.9e-05), alg_cluster.Cluster(set(['06089']), 77.359494209, 188.945068958, 163256, 5.7e-05), alg_cluster.Cluster(set(['08005']), 380.281283151, 270.268826873, 487967, 5.9e-05), alg_cluster.Cluster(set(['08001']), 379.950978294, 265.078784954, 363857, 6.6e-05), alg_cluster.Cluster(set(['08031']), 371.038986573, 266.847932979, 554636, 7.9e-05), alg_cluster.Cluster(set(['09003']), 925.917212741, 177.152290276, 857183, 5.7e-05), alg_cluster.Cluster(set(['12073']), 762.463896365, 477.365342219, 239452, 6.1e-05), alg_cluster.Cluster(set(['13313']), 737.308367745, 378.040993858, 83525, 5.6e-05), alg_cluster.Cluster(set(['13215']), 745.265661102, 430.987078939, 186291, 5.9e-05), alg_cluster.Cluster(set(['13067']), 747.238620236, 397.293799252, 607751, 6.4e-05), alg_cluster.Cluster(set(['13063']), 752.853876848, 406.722877803, 236517, 6.6e-05), alg_cluster.Cluster(set(['13151']), 756.589546538, 407.288873768, 119341, 5.6e-05), alg_cluster.Cluster(set(['13089']), 754.465443436, 400.059456026, 665865, 6.8e-05), alg_cluster.Cluster(set(['13121']), 750.160287596, 399.907752014, 816006, 7e-05), alg_cluster.Cluster(set(['13135']), 758.038826857, 395.110327675, 588448, 6.3e-05)])


#print slow_closest_pairs([alg_cluster.Cluster(set([]), 214.128077618, 396.893960776, 1, 0), alg_cluster.Cluster(set([]), 117.98714746, 355.203359224, 1, 0), alg_cluster.Cluster(set([]), 63.6146242678, 255.772864302, 1, 0), alg_cluster.Cluster(set([]), 376.551142131, 267.577115777, 1, 0), alg_cluster.Cluster(set([]), 754.616924589, 430.694788671, 1, 0), alg_cluster.Cluster(set([]), 572.136841573, 151.345524697, 1, 0), alg_cluster.Cluster(set([]), 646.807657539, 479.18344992, 1, 0), alg_cluster.Cluster(set([]), 575.468951442, 250.939543294, 1, 0), alg_cluster.Cluster(set([]), 912.229343552, 204.714030654, 1, 0), alg_cluster.Cluster(set([]), 297.479700714, 166.095561196, 1, 0), alg_cluster.Cluster(set([]), 544.129272622, 500.592855456, 1, 0), alg_cluster.Cluster(set([]), 746.701109685, 334.01735233, 1, 0), alg_cluster.Cluster(set([]), 389.120681395, 359.676072945, 1, 0), alg_cluster.Cluster(set([]), 125.27486023, 39.1497730391, 1, 0), alg_cluster.Cluster(set([]), 658.593243544, 197.564417673, 1, 0), alg_cluster.Cluster(set([]), 766.421393637, 209.960342327, 1, 0), alg_cluster.Cluster(set([]), 501.979819355, 329.536983654, 1, 0)])

#print fast_closest_pair([alg_cluster.Cluster(set([]), 214.128077618, 396.893960776, 1, 0), alg_cluster.Cluster(set([]), 117.98714746, 355.203359224, 1, 0), alg_cluster.Cluster(set([]), 63.6146242678, 255.772864302, 1, 0), alg_cluster.Cluster(set([]), 376.551142131, 267.577115777, 1, 0), alg_cluster.Cluster(set([]), 754.616924589, 430.694788671, 1, 0), alg_cluster.Cluster(set([]), 572.136841573, 151.345524697, 1, 0), alg_cluster.Cluster(set([]), 646.807657539, 479.18344992, 1, 0), alg_cluster.Cluster(set([]), 575.468951442, 250.939543294, 1, 0), alg_cluster.Cluster(set([]), 912.229343552, 204.714030654, 1, 0), alg_cluster.Cluster(set([]), 297.479700714, 166.095561196, 1, 0), alg_cluster.Cluster(set([]), 544.129272622, 500.592855456, 1, 0), alg_cluster.Cluster(set([]), 746.701109685, 334.01735233, 1, 0), alg_cluster.Cluster(set([]), 389.120681395, 359.676072945, 1, 0), alg_cluster.Cluster(set([]), 125.27486023, 39.1497730391, 1, 0), alg_cluster.Cluster(set([]), 658.593243544, 197.564417673, 1, 0), alg_cluster.Cluster(set([]), 766.421393637, 209.960342327, 1, 0), alg_cluster.Cluster(set([]), 501.979819355, 329.536983654, 1, 0)])

