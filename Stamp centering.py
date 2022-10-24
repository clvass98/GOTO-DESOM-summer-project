# This code is made for centering the random junk detections from the goto pipeline as a part of my summer project 2022 (GOTO). The centering was need fot the test set so the Deep Embedded Self Organizing Map (DESOM) can cluster the detections correctly.


from astropy.wcs import WCS
from astropy.io import fits
from random import randrange
from astropy.table import Table
from astropy.nddata import Cutout2D
import matplotlib.pyplot as plt
import pandas as pd
from pandas import concat
import numpy as np
from glob import glob
import math
from astropy.coordinates import SkyCoord
import astropy.units as u
import time
import warnings
import seaborn as sns
warnings.filterwarnings("ignore")


# scaling function is for normalization of the stamps

def scaling(stamps):
    stamps = np.nan_to_num(stamps, 1e-5)
    normalized_thumbnails = np.zeros(stamps.shape)
    for i, s in enumerate(stamps):
        th = s - np.median(s)
        pos_norm = np.abs(np.max(th[th>0]))
        neg_norm = np.abs(np.percentile(th[th<0], 0.05))
        th[th>0] /= pos_norm
        th[th<0] /= neg_norm
        th[th<-1] = -1
        normalized_thumbnails[i] = (th + 1)*0.5
    normalized_thumbnails = np.nan_to_num(normalized_thumbnails, 0.5)
    return normalized_thumbnails

lol = pd.read_csv('data_labels_copy.csv')
x_y_lista = [np.asarray(i[[1,-2,-1]]) for i in lol[lol['metalabel']=='randjunk'].values]

#closest_node function finds the closest coordinates from a list

def closest_node(node, nodes):
    nodes = np.asarray(nodes)
    dist_2 = np.sum((nodes - node)**2, axis=1)
    return np.argmin(dist_2)

new_x_coord = np.empty(0)
new_y_coord = np.empty(0)
x=0
y1 = []
np.savetxt('okay', np.array(y1), delimiter = ',')
stamps = np.empty([32,32],float)
stamps = stamps.reshape(1,32,32)
s_stamps = np.empty([32,32],float)
s_stamps = s_stamps.reshape(1,32,32)
error1 = []
for i in range(len(x_y_lista)):
    new_x_coord = np.empty(0)
    new_y_coord = np.empty(0)
    path = x_y_lista[i][0]
    x_y = [x_y_lista[i][1],x_y_lista[i][2]]
    print(np.shape(stamps))
    if path==x_y_lista[i-1][0]:
        new_y_coord = np.append(new_y_coord,sci[closest_node(x_y,sci)][1])
        new_x_coord = np.append(new_x_coord,sci[closest_node(x_y,sci)][0])
        stamps = np.concatenate((stamps, np.asarray([Cutout2D(diffimage, (x-1, y-1), (32,32), mode='partial').data for (x, y) in zip(new_x_coord, new_y_coord)])))
        if (abs(sci[closest_node(x_y,sci)][0] - x_y_lista[i][1]) > 15) or (abs(sci[closest_node(x_y,sci)][1] - x_y_lista[i][2]) > 15): 
            print('No science detection found in:')
            print(x)
            y1.append(x)
            continue
        try:
            s_stamps = np.concatenate((s_stamps,scaling(np.asarray([Cutout2D(diffimage, (x-1, y-1), (32,32), mode='partial').data for (x, y) in zip(new_x_coord, new_y_coord)]))))
        except ValueError:
            print('Error at')
            print(x)
            error1.append(x)
            pass
        except IndexError:
            print('Error at')
            print(x)
            error1.append(x)
            pass
    else:
        hdul = fits.open(path)
        diffphoto = Table(hdul['PHOTOMETRY_DIFF'].data).to_pandas()
        sciphoto = Table(hdul['PHOTOMETRY'].data).to_pandas()
        diffimage = hdul['DIFFERENCE'].data
        sci = sciphoto[['x','y']].to_numpy()
        new_y_coord = np.append(new_y_coord,sci[closest_node(x_y,sci)][1])
        new_x_coord = np.append(new_x_coord,sci[closest_node(x_y,sci)][0])
        stamps = np.concatenate((stamps, np.asarray([Cutout2D(diffimage, (x-1, y-1), (32,32), mode='partial').data for (x, y) in zip(new_x_coord, new_y_coord)])))
        if (abs(sci[closest_node(x_y,sci)][0] - x_y_lista[i][1]) > 15) or (abs(sci[closest_node(x_y,sci)][1] - x_y_lista[i][2]) > 15):
            print('No science detection found in:')
            print(x)
            y1.append(x)
            continue
        try:
            s_stamps = np.concatenate((s_stamps,scaling(np.asarray([Cutout2D(diffimage, (x-1, y-1), (32,32), mode='partial').data for (x, y) in zip(new_x_coord, new_y_coord)]))))
        except ValueError:
            print('Error at')
            print(x)
            error1.append(x)
            pass
        except IndexError:
            print('Error at')
            print(x)
            error1.append(x)
            pass
        hdul.close()
    x=x+1
#s_stamps = scaling(stamps[1:)
print(len(s_stamps))
tmp_df = pd.DataFrame(s_stamps[1:].reshape(-1, 1024))
tmp_df.to_csv("bogus_test_tom_2.csv", index=False)
np.savetxt('errors', np.array(error1), delimiter = ',')
np.savetxt('y',np.array(y1), delimiter = ',')




