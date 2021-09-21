#Collection of functions to use with to collated images to more easily view and work with clusters
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt


#Subfunction used in some of the more complex functions below
def dist_func(a1, a2):
    tot = 0
    for i in range(len(a1)):
        tot += (a1[i]-a2[i])**2
    return tot**0.5


#Functions to use directly on collated hduls
def hdultocluster(hduls, name="CAT", tablename="OBJ"):
    """
    Takes list of collated HDULs and returns a list of predicted objects.
    Each entry contains an ndarray of sources associated with a specific object.
    ndarrays contain source data taken from the HDU used to cluster sources.
    In cases where this function takes a long time to run, it is recommended to pickle the clusters
    to avoid having to rerun this function.
    Arguments:
        hduls -- list of fits HDULs
        name -- name of HDU for each HDUL containing catalog of sources and used to collate sources
        tablename -- name of TableHDU appended to HDUL by collate function
    """
    data = [hdul[name].data for hdul in hduls]
    datatype = data[0].dtype
    clusters = [np.zeros(len(hduls), dtype=datatype) for i in range(len(hduls[0][tablename].data))]
    for i in range(len(hduls)):
        indexes = hduls[i][tablename].data
        for j in range(len(indexes)):
            index = indexes[j][0]
            if index != -1:
                for k in range(len(datatype)):
                    clusters[j][i][k] = data[i][index][k]
    return clusters


def cluster_ratio(hduls, name="CAT", tablename="OBJ"):
    """
    Prints number of clustered sources out of total sources as well as the percent of clustered sources.
    Arguments:
        hduls -- list of fits HDULSs
        name -- name of HDU for each HDUL containing catalog of sources that were collated
        tablename -- name of TableHDU appended by collate
    """
    points = 0
    clustered_points = 0
    for hdul in hduls:
        points += len(hdul[name].data)
        noise = 0
        for point in hdul[tablename].data:
            if point[0] == -1:
                noise += 1
        clustered_points += len(hdul[tablename].data) - noise
    print("{} out of {} points clustered".format(clustered_points, points), clustered_points/points)



#Functions to be used on either individual clusters or the whole list of clusters returned by hdultocluster()
def cluster_mean(cluster, dim='flux'):
    """
    Returns the mean value of a particular dimension for a cluster.
    Arguments:
        cluster -- singular entry(cluster) from hdultocluster()
        dim -- dimension for which mean will be calculated and returned
    """
    sources = len(cluster[dim]) - list(cluster[dim]).count(0)
    return sum(cluster[dim])/sources


def cluster_std(cluster, dim='flux'):
    """
    Returns the standard deviation of a particular dimension for a cluster.
    Arguments:
        cluster -- singular entry(cluster) from hdultocluster()
        dim -- dimension for which standard deviation will be calculated and returned
    """
    values = [value for value in cluster[dim] if value != 0]
    return np.std(values)


def cluster_sm_ratio(cluster, dim='flux'):
    """
    Returns the ratio of mean to standard deviation for a dimension of a cluster.
    Can be used to sort a list of clusters and find stars of unusually high flux variance.
    Arguments:
        cluster -- singular entry(cluster) from hdultocluster()
        dim -- dimension for which mean and standard deviation will be calculated and returned
    """
    return cluster_std(cluster, dim)/cluster_mean(cluster, dim)


def cluster_search(clusters, x, y):
    """
    Returns closest cluster to a specified x,y coordinate.
    Arguments:
        clusters -- list of clusters returned by hdultocluster()
        x -- x coordinate
        y -- y coordinate
    """
    means = [[cluster_mean(cluster, 'x'), cluster_mean(cluster, 'y')] for cluster in clusters]
    dist = 10000000
    smallest = None
    for i in range(len(means)):
        if dist_func(means[i], [x,y]) < dist:
            dist = dist_func(means[i], [x,y])
            smallest = i
    return clusters[smallest]


def cluster_plot(clusters):
    """
    Plots that shows each cluster's source count on the y-axis and each cluster's mean flux value on the x-axis.
    Useful alongside cluster_ratio() for evaluating clustering quality.
    Arguments:
        clusters -- list of clusters returned by hdultocluster()
    """
    x = [cluster_mean(cluster) for cluster in clusters]
    y = [len(cluster)-list(cluster['flux']).count(0) for cluster in clusters]
    plt.plot(x,y,'go')
    plt.xscale('log')
    plt.xlabel('flux')
    plt.ylabel('record sources count')
    plt.show()

