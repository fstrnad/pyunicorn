import numpy as np

import matplotlib.pyplot as plt

import random

from pyunicorn.core import network, geo_network
from pyunicorn.core import spatial_network as sp
from pyunicorn.core.grid import Grid
import sys


#net=network.Network.ErdosRenyi(n_nodes=500, link_probability=0.05)
# Create Random Grids
rect_grid_num=10
grid=Grid.RegularGrid( time_seq=np.arange(2), lat_grid=np.random.randint(low=0, high=40, size=rect_grid_num), 
                  lon_grid=np.random.randint(low=0, high=40, size=rect_grid_num) )

erdos_renyi=sp.SpatialNetwork.ErdosRenyi(grid=grid, n_nodes=int(rect_grid_num**2),link_probability=0.1 )
#print(erdos_renyi)
geo_model=sp.SpatialNetwork(grid=erdos_renyi.grid, adjacency=erdos_renyi.adjacency )


fast_link_list,dic_fast=geo_model.GeoModel1(n_steps=int(5e3), tolerance=0.2, grid_type='euclidean', verbose=False)

print(fast_link_list)
sys.exit()

link_list1, dic1=geo_model.GeoModel1_slow(n_steps=int(5e3), tolerance=0.2, grid_type='euclidean', verbose=False)

erdos_renyi=sp.SpatialNetwork.ErdosRenyi(grid=grid, n_nodes=int(rect_grid_num**2),link_probability=0.1 )
geo_model=sp.SpatialNetwork(grid=erdos_renyi.grid, adjacency=erdos_renyi.adjacency)
link_list2, dic2=geo_model.GeoModel1_slow(n_steps=int(5e3), tolerance=0.05, grid_type='euclidean', verbose=False)

erdos_renyi=sp.SpatialNetwork.ErdosRenyi(grid=grid, n_nodes=int(rect_grid_num**2),link_probability=0.1 )
geo_model=sp.SpatialNetwork(grid=erdos_renyi.grid, adjacency=erdos_renyi.adjacency)
link_list3, dic3=geo_model.GeoModel1_slow(n_steps=int(5e3), tolerance=0.01, grid_type='euclidean', verbose=False)


#link_list2, dic2 = geo_model.GeoModel2(n_steps=100, grid_type='spherical', verbose=True)

print(dic1['x'])

# Test results
plt.figure(figsize=(6,4))

plt.xlabel('Steps')
plt.ylabel('Values')

plt.plot(dic1['x'], dic1['H'], label='H-Distance, GeoModel1')
plt.plot(dic2['x'], dic2['H'], label='H-Distance, GeoModel2')
plt.plot(dic3['x'], dic3['H'], label='H-Distance, GeoModel3')

#plt.plot(dic1['x'], dic1['L'], label='L-Distance, GeoModel1')
# plt.plot(dic1['x'], dic1['T'], label='Clustering coefficient, GeoModel1')
#plt.plot(dic2['x'], dic2['H'], label='H-Distance, GeoModel2')
plt.legend()
plt.tight_layout()
plt.show()