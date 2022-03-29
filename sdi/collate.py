import click
from . import _cli as cli
import sys
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import dbscan, OPTICS
from astropy.io import fits

#Subfunctions utilized in clustering function
def _norm(array):
    array -= min(array)
    array *= 1/max(array)
    return array

def dist_func(a1, a2):
    tot = 0
    for i in range(len(a1)):
        tot += (a1[i]-a2[i])**2
    return tot**0.5

def obj_mean(cluster):
    means = None
    for source in cluster:
        if source != None:
            if not means:
                 means = [0 for i in range(len(source))]
            for idx, coord in enumerate(source):
                means[idx] += coord
    return [mean/(len(cluster)-cluster.count(None)) for mean in means]

#Exact same function as in _scripts/collateutils, called at the end of running collate
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


#Collate function itself
def collate(hduls, name="CAT", tablename="OBJ", coords = ["x", "y", "flux"], algorithm="DBSCAN", minpts=4, eps=0.001, maxeps=np.inf, xi=0.85, collision_fix=False, show_collisions=False):
    """
    Clusters source data from HDULs to predict sky objects and then appends a TableHDU to each HDUL containing object/cluster data. 
    The nth index of each TableHDU for each HDUL refers to the same object with the value of the index referring to the
    (value)th source in the HDUL's source catalog.
    Collisions refer to when two sources in the same HDUL are clustered to the same object.  In the event of a collision,
    one source will overwrite the other in the cluster to ensure cluster sizes do not exceed the number of HDULs.
    Arguments:
        hduls -- list of fits HDULs
        name -- name of HDU for each HDUL containing catalog of sources
        tablename -- name of TableHDU to store cluster data
        coords -- list of dimensions to cluster sources over, multiple custom coordinates may be specified using: -c dim1 -c dim2 etc.
        algorithm -- choice of 'DBSCAN' or 'OPTICS' clustering algorithm
        minpts -- minimum cluster size
        eps -- point neighborhood radius parameter, used only in DBSCAN
        maxeps -- maximum neighborhood radius parameter, used only in OPTICS
        xi -- minimum steepness of cluster boundary, used only in OPTICS
        collision_fix -- allows for optional resolving of collisions which may result in increased runtime
        show_collisions -- Display colliding sources and their corresponding images and clusters
    """
    hduls = [hdul for hdul in hduls]
    collisions = []
    collision_count = 0
    labels = None
    irdf = [[] for i in range(len(coords))]
    for hdul in hduls:
        data = hdul[name].data
        for idx, coord in enumerate(coords):
            irdf[idx] += list(data[coord])
    normirdf = [_norm(array) for array in irdf]

    if algorithm == "DBSCAN":
        print('\n Collating using {} with {} coordinates, minpts={} and eps={}'.format(algorithm, coords, minpts, eps))
        cores, labels = dbscan(np.vstack(normirdf).T, eps=eps, min_samples=minpts)
    elif algorithm == "OPTICS":
        print('\n Collating using {} with {} coordinates, minpts={}, xi={}, and maxeps={}'.format(algorithm, coords, minpts, xi, maxeps))
        labels = OPTICS(min_samples=minpts, max_eps=maxeps, xi=xi).fit(np.vstack(normirdf).T).labels_
    else:
        print('Invalid clustering method, please use "DBSCAN" or "OPTICS"')
        return (hdul for hdul in hduls)
   
    obj_count = max(labels)
    if obj_count == -1:
        print('No clusters could be formed. Please try using more images, increasing eps, or changing algorithms.')
        return (hdul for hdul in hduls)

    for i, hdul in enumerate(hduls):
        source_count = len(hdul[name].data)
        cols = [-1 for j in range(obj_count+1)]
        for j in range(source_count):
            if labels[j] != -1:
                if cols[labels[j]] == -1:
                    cols[labels[j]] = j
                else:
                    collision_count += 1
                    #Store labels for collision image, cluster containing collision, and the two conflicting sources
                    collisions.append([i, labels[j], j, cols[labels[j]]])
                    if show_collisions:
                        print("collision at hdul {}, object '{}' at sources {} {}".format(i, labels[j], j, cols[labels[j]]))
        if len(labels) != source_count:
            labels = labels[source_count:]
        table = [fits.Column(name='objects', format="I", array=np.array(cols), ascii=True)]
        table = fits.TableHDU.from_columns(table)
        table.name = tablename
        hdul.append(table)

    print("{} collisions".format(collision_count))

    #Dealing with collisions, takes colliding sources and keeps the source closest to the mean coordinate of their corresponding cluster.
    if collision_fix == True:
        if len(collisions) > 0:
            print("Resolving collisions")
            all_sources = [hdul[name].data for hdul in hduls]
            count = len(hduls[0][tablename].data)
            clusters = [[] for i in range(count)]
            for hdul in hduls:
                sources = hdul[name].data
                for i in range(count):
                    index = hdul[tablename].data[i][0]
                    if index != -1:
                        clusters[i].append([sources[index][coord] for coord in coords])
                    else:
                        clusters[i].append(None)

            means = [obj_mean(cluster) for cluster in clusters]
            for collision in collisions:
                sources = all_sources[collision[0]]
                mean = means[collision[1]]
                #Gets coords of the two colliding sources before calculating their distance from the cluster mean to determine which source is kept.
                s1 = [sources[collision[2]][coord] for coord in coords]
                s2 = [sources[collision[3]][coord] for coord in coords]
                if dist_func(s1, mean) < dist_func(s2, mean):
                    hduls[collision[0]][tablename].data[collision[1]][0] = collision[2]

    print("Clustering complete, cluster data stored in '{}' TableHDU".format(tablename))
    cluster_ratio(hduls, name, tablename)
    return (hdul for hdul in hduls)

@cli.cli.command("collate")
@click.option("-n", "--name", default="CAT", help="name of HDU for each HDUL containing catalog of sources")
@click.option("-t", "--tablename", default="OBJ", help="name of TableHDU to store cluster data")
#@click.option("-c", "--coords", default="xy", type=click.Choice(["xy", "radec"], case_sensitive=False), help="choice of either 'xy' or 'radec' coordinate system used to calculate distance between sources for clustering")
@click.option("-c", "--coords", default=["x", "y", "flux"], multiple=True, help="list of dimensions to cluster sources over, multiple custom coordinates may be specified using: -c dim1 -c dim2 etc.")
@click.option("-a", "--algorithm", default="DBSCAN", type=click.Choice(["DBSCAN", "OPTICS"], case_sensitive=False), help="choice of 'DBSCAN' or 'OPTICS' clustering algorithm")
@click.option("-p", "--minpts", default=4, help="minimum cluster size")
@click.option("-e", "--eps", default=0.001, help="point neighborhood radius parameter, used only in DBSCAN")
@click.option("-m", "--maxeps", default=np.inf, help="maximum neighborhood radius parameter, used only in OPTICS")
@click.option("-x", "--xi", default=0.85, help="minimum steepness of cluster boundary, used only in OPTICS")
@click.option("-f", "--collision_fix", default=False, help="allows for optional resolving of collisions which may result in increased runtime")
@click.option("-s", "--show_collisions", default=False, help="Display colliding sources and their corresponding images and clusters")
@cli.operator

#collate function wrapper
def collate_cmd(hduls, name="CAT", tablename="OBJ", coords=["x", "y", "flux"], algorithm="DBSCAN", minpts=4, eps=0.001, maxeps=np.inf, xi=0.85, collision_fix=False, show_collisions=False):
    """
    Clusters source data from HDULs to predict sky objects and then appends a TableHDU to each HDUL containing object/cluster data. 
    The nth index of each TableHDU for each HDUL refers to the same object with the value of the index referring to the
    (value)th source in the HDUL's source catalog.
    Collisions refer to when two sources in the same HDUL are clustered to the same object.  In the event of a collision,
    one source will overwrite the other in the cluster to ensure cluster sizes do not exceed the number of HDULs.
    Arguments:
        hduls -- list of fits HDULs
        name -- name of HDU for each HDUL containing catalog of sources
        tablename -- name of TableHDU to store cluster data
        coords -- list of dimensions to cluster sources over, multiple custom coordinates may be specified using: -c dim1 -c dim2 etc.
        algorithm -- clustering algorithm of choice, either DBSCAN or OPTICS
        minpts -- minimum cluster size
        eps -- point neighborhood radius parameter, used only in DBSCAN
        maxeps -- maximum neighborhood radius parameter, used only in OPTICS
        xi -- minimum steepness of cluster boundary, used only in OPTICS
        collision_fix -- allows for optional resolving of collisions
        show_collisions -- Display colliding sources and their corresponding images and clusters
    """
    return  collate(hduls, name, tablename, coords, algorithm, minpts, eps, maxeps, xi, collision_fix, show_collisions)
