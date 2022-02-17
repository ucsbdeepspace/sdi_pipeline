import sys
import numpy as np
import matplotlib.pyplot as plt
import time
from sklearn.cluster import dbscan, OPTICS
from astropy.io import fits
from glob import glob
from vis import image, circle_all

images = [fits.open(im) for im in sys.argv[1:]]

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

def cluster_func(hduls, name="CAT", tablename="OBJ", coords="radec"):
    count = len(hduls[0][tablename].data)
    clusters = [[] for i in range(count)]
    for hdul in hduls:
        sources = hdul[name].data
        for i in range(count):
            index = hdul[tablename].data[i][0]
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

#Clustering function itself
def collate(hduls, name="CAT", tablename="OBJ", algorithm="DBSCAN", minpts=4, eps=0.001, maxeps=np.inf, xi=0.85, coords="radec"):
    t1 = time.perf_counter()
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
        print('Collating using {} with {} coordinates, minpts={} and eps={}'.format(algorithm, coords, minpts, eps))
        cores, labels = dbscan(np.vstack(normirdf).T, eps=eps, min_samples=minpts)
    elif algorithm == "OPTICS":
        print('Collating using {} with {} coordinates, minpts={}, xi={}, and maxeps={}'.format(algorithm, coords, minpts, xi, maxeps))
        labels = OPTICS(min_samples=minpts, max_eps=maxeps, xi=xi).fit(np.vstack(normirdf).T).labels_
    else:
        print('Invalid clustering method, please use "DBSCAN" or "OPTICS"')
        return
    t2 = time.perf_counter()
    
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
        table = [fits.Column(name='objects'.format(j), format="I", array=np.array(cols), ascii=True)]
        table = fits.TableHDU.from_columns(table)
        table.name = tablename
        hduls[i].append(table)

    print("{} collisions".format(collision_count))
    t3 = time.perf_counter()

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
                hduls[collision[0]][tablename].data[collision[1]][0] = collision[2]
    t4 = time.perf_counter()
    print("Clustering complete, cluster data stored in '{}' TableHDU".format(tablename))
    print(t2-t1, t3-t1, t4-t1)


#cluster analysis tools
def hdultocluster(hduls, name="CAT", tablename="OBJ"):
    t1 = time.perf_counter()
    data = [hdul[name].data for hdul in hduls]
    datatype = data[0].dtype
    clusters = [np.zeros(len(hduls), dtype=datatype) for i in range(len(hduls[0][tablename].data))]
    for i in range(len(hduls)):
        imdata = data[i]
        indexes = hduls[i][tablename].data
        for j in range(len(indexes)):
            index = indexes[j][0]
            if index != -1:
                for k in range(len(datatype)):
                    clusters[j][i][k] = imdata[index][k]
    t2 = time.perf_counter()
    print(t2-t1)
    return clusters

def cluster_mean(cluster, dim='flux'):
    sources = len(cluster[dim]) - list(cluster[dim]).count(0)
    return sum(cluster[dim])/sources

def cluster_std(cluster, dim='flux'):
    values = [value for value in cluster[dim] if value != 0]
    return np.std(values)

def cluster_sm_ratio(cluster, dim='flux'):
    return cluster_std(cluster, dim)/cluster_mean(cluster, dim)

def ccluster_sm_ratio(cluster, dim='FLUX'):
    return cluster_std(cluster, dim)/cluster_mean(cluster, dim)

def cluster_search(clusters, x, y):
    means = [[cluster_mean(cluster, 'ra'), cluster_mean(cluster, 'dec')] for cluster in clusters]
    dist = 10000000
    smallest = None
    for i in range(len(means)):
        if dist_func(means[i], [x,y]) < dist:
            dist = dist_func(means[i], [x,y])
            smallest = i
    return clusters[smallest]

def cluster_search_radec(clusters, ra, dec):
    means = [[cluster_mean(cluster, 'ra'), cluster_mean(cluster, 'dec')] for cluster in clusters]
    dist = np.inf
    smallest = None
    for i in range(len(means)):
        if dist_func(means[i], [ra,dec]) < dist:
            dist = dist_func(means[i], [ra,dec])
            smallest = i
    return clusters[smallest]

def cluster_coords(cluster):
    coords = [[[cluster['ra'][i], cluster['dec'][i]]] for i in range(len(cluster))]
    return coords

def flux_func(cluster):
    return sum(list(cluster['flux']))/(len(cluster)-list(cluster['flux']).count(0))

def field_func(cluster, field='flux'):
    return sum(list(cluster[field]))/(len(cluster)-list(cluster[field]).count(0))

def source_search(hdul, coords, dims = ['x','y','flux'], name="CAT"):
    sources = np.array([hdul[name].data[dim] for dim in dims]).T
    dist = 10000000
    closest = None
    for source in sources:
        if dist_func(coords, source[:len(coords)]) < dist:
            dist = dist_func(coords, source)
            closest = source
    return closest

def cluster_plot(clusters, field='flux', label='record flux avg'):
    x = [field_func(cluster, field) for cluster in clusters]
    y = [len(cluster)-list(cluster[field]).count(0) for cluster in clusters]
    plt.plot(x,y,'go')
    plt.xscale('log')
    plt.xlabel(label)
    plt.ylabel('record sources count')
    plt.show()

def cluster_stat_plot(clusters, dim='flux'):
    x = [cluster_sm_ratio(cluster, dim) for cluster in clusters]
    y = [len(cluster)-list(cluster[dim]).count(0) for cluster in clusters]
    plt.plot(x,y,'go')
    plt.xscale('log')
    plt.xlabel('std to mean ratio')
    plt.ylabel('record sources count')
    plt.show()

def cluster_ratio(hduls, name="CAT", tablename="OBJ"):
    points = 0
    clustered_points = 0
    for hdul in hduls:
        points += len(hdul[name].data)
        noise = 0
        for point in hdul[tablename].data:
            if point[0] == -1:
                noise += 1
        clustered_points += len(hdul[tablename].data) - noise
    print(points, clustered_points, clustered_points/points)

#misc tools
def time_sort(image):
    return image["PRIMARY"].header["DATE-OBS"]

def all_coords(hduls):
    coords = [[] for hdul in hduls]
    for i in range(len(hduls)):
        data = hduls[i]["CAT"].data
        for coord in zip(data['x'], data['y']):
            coords[i].append(coord)
    return coords

def source_pixel_ratio(image):
    return len(image['CAT'].data)/(image[0].header['NAXIS1']*image[0].header['NAXIS2'])
        
