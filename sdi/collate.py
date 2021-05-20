import click
from . import _cli as cli
import numpy as np
from sklearn.cluster import dbscan, OPTICS
from astropy.io import fits


#Subfunctions used in collate function
def _norm(array):
    array -= min(array)
    array *= 1/max(array)
    return array

def dist_func(a1, a2):
    tot = 0
    for i in range(len(a1)):
        tot += (a1[i]-a2[i])**2
    return tot**0.5

def cluster_func(hduls, name="CAT", tablename="OBJ", coords="radec"):
    count = len(hduls[0][tablename].data[0])
    clusters = [[] for i in range(count)]
    for hdul in hduls:
        sources = hdul[name].data
        for i in range(count):
            index = hdul[tablename].data[0][i]
            if index != -1:
                if coords == "radec":
                    clusters[i].append([sources[index]["ra"], sources[index]["dec"], sources[index]["flux"]])
                elif coords == "xy":
                    clusters[i].append([sources[index]["x"], sources[index]["y"], sources[index]["flux"]])
            else:
                clusters[i].append(None)
    return clusters

def obj_mean(cluster):
    means = [0,0,0]
    for source in cluster:
        if source != None:
            for i in range(len(source)):
                means[i] += source[i]
    return [mean/(len(cluster)-cluster.count(None)) for mean in means]

#Collate function itself
def collate(hduls, name="CAT", tablename="OBJ", coords = "radec", algorithm="DBSCAN", minpts=4, eps=0.001, maxeps=np.inf, xi=0.85):
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
        coords -- coordinate system used to calculate distance between sources for clustering
        algorithm -- clustering algorithm of choice
        minpts -- minimum cluster size
        eps -- point neighborhood radius parameter, used only in DBSCAN
        maxeps -- maximum neighborhood radius parameter, used only in OPTICS
        xi -- minimum steepness of cluster boundary, used only in OPTICS
    """
    hduls = [hdul for hdul in hduls]
    collisions = []
    collision_count = 0
    irdf = [[] for i in range(3)]
    labels = None

    for hdul in hduls:
        data = hdul[name].data
        if coords == "radec":
            irdf[0] = irdf[0] + list(data["ra"])
            irdf[1] = irdf[1] + list(data["dec"])
        elif coords == 'xy':
            irdf[0] = irdf[0] + list(data["x"])
            irdf[1] = irdf[1] + list(data["y"])
        else:
            print('Invalid value for coords, please use "radec" or "xy"')
            return
        irdf[2] =  irdf[2] + list(data["flux"])
    normirdf = [_norm(array) for array in irdf]

    if algorithm == "DBSCAN":
        cores, labels = dbscan(np.vstack(normirdf).T, eps=eps, min_samples=minpts)
    elif algorithm == "OPTICS":
        labels = OPTICS(min_samples=minpts, max_eps=maxeps, xi=xi).fit(np.vstack(normirdf).T).labels_
    else:
        print('Invalid clustering method, please use "DBSCAN" or "OPTICS"')
        return
  
    obj_count = max(labels)
    for i in range(len(hduls)):
        source_count = len(hduls[i][name].data)
        cols = [-1 for j in range(obj_count+1)]
        for j in range(source_count):
            if labels[j] != -1:
                if cols[labels[j]] == -1:
                    cols[labels[j]] = j
                else:
                    collision_count += 1
                    collisions.append([i, labels[j], j, cols[labels[j]]])
                    print("collision at hdul {}, object '{}' at sources {} {}".format(i, labels[j], j, cols[labels[j]]))

        if len(labels) != source_count:
            labels = labels[source_count:]
        table = [fits.Column(name='{}'.format(j), format="I", array=np.array([cols[j]]), ascii=True) for j in range(obj_count+1)]
        table = fits.TableHDU.from_columns(table)
        table.name = tablename
        hduls[i].append(table)

    print("{} collisions".format(collision_count))

    #Dealing with collisions
    if len(collisions) > 0:
        print("Resolving collisions")
        all_sources = [hdul[name].data for hdul in hduls]
        clusters = cluster_func(hduls, name, tablename, coords)
        means = [obj_mean(cluster) for cluster in clusters]
        for collision in collisions:
            sources = all_sources[collision[0]]
            mean = means[collision[1]]
            if coords == 'radec':
                s1 = (sources[collision[2]]["ra"], sources[collision[2]]["dec"], sources[collision[2]]["flux"])
                s2 = (sources[collision[3]]["ra"], sources[collision[3]]["dec"], sources[collision[3]]["flux"])
            elif coords == 'xy':
                s1 = (sources[collision[2]]["x"], sources[collision[2]]["y"], sources[collision[2]]["flux"])
                s2 = (sources[collision[3]]["x"], sources[collision[3]]["y"], sources[collision[3]]["flux"])
            if dist_func(s1, mean) < dist_func(s2, mean):
                hduls[collision[0]][tablename].data[str(collision[1])] = np.array([collision[2]])
    print("Clustering complete, cluster data stored in '{}' TableHDU".format(tablename))
    return (hdul for hdul in hduls)

@cli.cli.command("collate")
@click.option("-n", "--name", default="CAT", help="name of HDU for each HDUL containing catalog of sources")
@click.option("-t", "--tablename", default="OBJ", help="name of TableHDU to store cluster data")
@click.option("-c", "--coords", default="radec", help="coordinate system used to calculate distance between sources for clustering")
@click.option("-a", "--algorithm", default="DBSCAN", help="clustering algorithm of choice")
@click.option("-p", "--minpts", default=4, help="minimum cluster size")
@click.option("-e", "--eps", default=0.001, help="point neighborhood radius parameter, used only in DBSCAN")
@click.option("-m", "--maxeps", default=np.inf, help="maximum neighborhood radius parameter, used only in OPTICS")
@click.option("-x", "--xi", default=0.85, help="minimum steepness of cluster boundary, used only in OPTICS")
@cli.operator

#collate function wrapper
def collate_cmd(hduls, name="CAT", tablename="OBJ", coords="radec", algorithm="DBSCAN", minpts=4, eps=0.001, maxeps=np.inf, xi=0.85):
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
        coords -- coordinate system used to calculate distance between sources for clustering
        algorithm -- clustering algorithm of choice
        minpts -- minimum cluster size
        eps -- point neighborhood radius parameter, used only in DBSCAN
        maxeps -- maximum neighborhood radius parameter, used only in OPTICS
        xi -- minimum steepness of cluster boundary, used only in OPTICS
    """
    return  collate(hduls, name, tablename, coords, algorithm, minpts, eps, maxeps, xi)

#Function to return clustered source as a list from collated hduls
#nth index of clusters refers to nth cluster/object
def hdultocluster(hduls, name="CAT", tablename="OBJ"):
    """
    Takes list of collated HDULs and returns a list of predicted objects.
    Each entry contains an ndarray of sources associated with a specific object.
    ndarrays contain source data taken from the HDU used to cluster sources.
    Arguments:
        hduls -- list of fits HDULs
        name -- name of HDU for each HDUL containing catalog of sources and used to collate sources
        tablename -- name of TableHDU appended to HDUL by collate function
    """
    data = [hdul[name].data for hdul in hduls]
    datatype = data[0].dtype
    clusters = [np.zeros(len(hduls), dtype=datatype) for i in range(len(hduls[0][tablename].data[0]))]
    for i in range(len(hduls)):
        indexes = hduls[i][tablename].data[0]
        for j in range(len(indexes)):
            if indexes[j] != -1:
                for k in range(len(datatype)):
                    clusters[j][i][k] = data[i][indexes[j]][k]
    return clusters


