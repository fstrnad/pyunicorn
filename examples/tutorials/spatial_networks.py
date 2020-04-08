import numpy as np

import matplotlib.pyplot as plt

from pyunicorn.core import spatial_network as sp
from pyunicorn.core.grid import Grid


"""
This example code offers an overview of the spatial network code for GeoModel1 
and GeoModel2 can be called. Furthermore, with the python versions GeoModel1_py
and GeoModel2_py further analysis can be done, to visualize the working of the
network operations.
"""

# Create Random Grids
rect_grid_num=30
grid=Grid.RegularGrid( time_seq=np.arange(2), lat_grid=np.random.randint(low=0, high=40, size=rect_grid_num), 
                  lon_grid=np.random.randint(low=0, high=40, size=rect_grid_num) )

erdos_renyi=sp.SpatialNetwork.ErdosRenyi(grid=grid, n_nodes=int(rect_grid_num**2),link_probability=0.1 )
geo_model=sp.SpatialNetwork(grid=erdos_renyi.grid, adjacency=erdos_renyi.adjacency )

new_link_list=geo_model.GeoModel1(n_steps=int(5e4), tolerance=1, grid_type='euclidean', verbose=False)

print("New link list", new_link_list)

# Here the python version is used for visualizing the Hamming-distance
erdos_renyi=sp.SpatialNetwork.ErdosRenyi(grid=grid, n_nodes=int(rect_grid_num**2),link_probability=0.1 )
geo_model=sp.SpatialNetwork(grid=erdos_renyi.grid, adjacency=erdos_renyi.adjacency )

link_list1, dic1=geo_model.GeoModel1_py(n_steps=int(5e4), tolerance=0.1, grid_type='euclidean', verbose=False)
print("New link list" , link_list1)

erdos_renyi=sp.SpatialNetwork.ErdosRenyi(grid=grid, n_nodes=int(rect_grid_num**2),link_probability=0.1 )
geo_model=sp.SpatialNetwork(grid=erdos_renyi.grid, adjacency=erdos_renyi.adjacency)
link_list2, dic2=geo_model.GeoModel1_py(n_steps=int(5e4), tolerance=0.5, grid_type='euclidean', verbose=False)
 
erdos_renyi=sp.SpatialNetwork.ErdosRenyi(grid=grid, n_nodes=int(rect_grid_num**2),link_probability=0.1 )
geo_model=sp.SpatialNetwork(grid=erdos_renyi.grid, adjacency=erdos_renyi.adjacency)
link_list3, dic3=geo_model.GeoModel1_py(n_steps=int(5e4), tolerance=1, grid_type='euclidean', verbose=False)


#link_list2, dic2 = geo_model.GeoModel2(n_steps=100, grid_type='spherical', verbose=True)

#print(dic1['x'])

# Test results
plt.figure(figsize=(6,4))
 
plt.xlabel('Steps')
plt.ylabel('Values')
 
plt.plot(dic1['x'], dic1['H'], label='H-Distance, GeoModel1, tolerance=0.1')
plt.plot(dic2['x'], dic2['H'], label='H-Distance, GeoModel1, tolerance=0.5')
plt.plot(dic3['x'], dic3['H'], label='H-Distance, GeoModel1, tolerance=1')
 
plt.legend()
plt.tight_layout()
plt.show()