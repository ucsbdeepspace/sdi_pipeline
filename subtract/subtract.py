import ois
import os
from astropy.io import fits

def subtract(sources, template):
    output = []

    if type(sources) == list:
        for element in sources:
            diff = ois.optimal_system(image=element.data, refimage=template, method='Bramich')
            output.append(diff)     
            
    return output
