import pickle
import numpy as np
from scipy.spatial import distance

test1 = pickle.load(open('/home/williamwang/sources_Bramich.pkl', 'rb'))
# This line load the first pickle file
test2 = pickle.load(open('/home/williamwang/sources_hotpants.pkl', 'rb'))

x1 = test1[0][0]['x']  #read the x coordinate of the first list
y1 = test1[0][0]['y']  
x2 = test2[0][0]['x']  
y2 = test2[0][0]['y']  

coord1 = np.stack((x1,y1),axis=-1)  #stack the x,y coordinate of the first file to a list of coordniates in the form (x,y)
coord2 = np.stack((x2,y2),axis=-1)

dist = distance.cdist(coord1,coord2).min(axis=1)
# save the closest euclidean distance for every coordinate in coord1 from coord2

print(len(dist[dist<=1]))
# print the number of closest distance less than 1.
