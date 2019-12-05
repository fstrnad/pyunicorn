#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# This file is part of pyunicorn.
# Copyright (C) 2008--2019 Jonathan F. Donges and pyunicorn authors
# URL: <http://www.pik-potsdam.de/members/donges/software>
# License: BSD (3-clause)
import igraph
import sys
from .geo_network import GeoNetwork
from .network import NetworkError
from pyunicorn.core._ext.numerics import _geo_model1
#from ._ext.numerics import _geo_model1



"""
Provides classes for analyzing spatial complex networks, that allow to generate 
random surrogates from a given spatially embedded network which can preserve 
global and local statistics associated with the nodes' embedding in a metric
space. 
This code is based on the work of Wiedermann et al. 2016 . 
"""

# array object and fast numerics
import numpy as np
from numpy.random import choice # random number generation

class SpatialNetwork(GeoNetwork):
    
    def __init__(self, grid, adjacency=None, edge_list=None, directed=False,
                 node_weight_type="surface", silence_level=0): 
        
        
        GeoNetwork.__init__(self, grid=grid, adjacency=adjacency, 
                            edge_list=edge_list, directed=directed, node_weight_type=node_weight_type, 
                            silence_level=silence_level)
           
    
    """
    Gives distances in spherical/euclidian networks from i to all surrounding nodes
    """
    def distance(self,  D, node,  list_of_neighbors=None):
        
        d_i_all=[]
        if list_of_neighbors is None:
            list_of_neighbors=range(len(D[node]))
        for i in list_of_neighbors:
            d_i_all.append(D[node][i])

        return np.array(d_i_all)
    
    def link_id(self, node1, node2, N):
        
        id = -0.5 * node1 * (node1 - 2 * N + 1)
        id += node2 - node1 - 1
        
        return id, node1, node2 
        
    
    def check_requirements(self, link_list, n_steps, ):
    
        return True
    
    def check_GeoModel_requirements(self, grid, tolerance, grid_type):
        return True
    
    def _Hamming_Distance(self, A, B):
        #  Check if the graphs have the same number of vertices
        if len(A) == len(B):
            #  Calculate the hamming distance
            hamming = (A != B).sum()

            #  Return the normalized hamming distance
            return hamming / float(len(A) * (len(A) - 1))
        else:
            raise NetworkError("Only defined for networks with same number of \
                               nodes.")  
    
    """
    GeoModel1 preserves the global link-length distribution P(l). Hence,
    the potentially newly established links are of the same
    length as those that are removed from the network. This means
    that the four randomly drawn nodes i, j, k, and l form
    a kite with exactly one link present at each of the two sides
    of the same length 
    """ 
    def GeoModel1(self, n_steps, tolerance, grid_type="spherical", verbose=False):
        
        if grid_type == "spherical":
            D = self.grid.angular_distance()
        elif grid_type == "euclidean":
            D = self.grid.euclidean_distance()
        else:
            print("ERROR, This grid type is not recognized: ",  grid_type)
            sys.exit(1)
       
        #  Get number of nodes
        N = self.N
        #  Get number of links
        E = self.n_links 
        #check_requirements(link_list, n_steps, n_steps)
        #check_GeoModel_requirements(grid, tolerance, grid_type)
        
        #grid=self.grid.coord_sequence_from_rect_grid(lat_grid, lon_grid)
        # Needs to be done to prevent error in computation!
        link_list=np.array(self.graph.get_edgelist(), np.int).copy(order='c')        
#         link_list = link_list.copy()
        
        A= self.adjacency.copy(order='c')
        
        n_sampling_points=self.n_links
        
       
    
        original_link_ids = -0.5 * link_list[:, 0] * (link_list[:, 0] - 2 * N + 1)
        original_link_ids += link_list[:, 1] - link_list[:, 0] - 1
        
        #print("Original_link_ids", link_list[:,0], link_list[:,1], original_link_ids)
        sur_link_ids = original_link_ids.copy()
    
        compute_at = int(n_steps / n_sampling_points)
        print(compute_at)
        g = igraph.Graph(link_list.tolist())
    
        T = [g.transitivity_avglocal_undirected()]
        x = [0]
        H = [0]
        L = [g.average_path_length()]
        ass = [g.assortativity_degree()]
    
        c = 0
        
        #print(link_list.shape)
        np.savetxt('link_list.txt', link_list, delimiter=' ', fmt='%d')
        
        print(A.shape)
        np.savetxt('A.txt', A, delimiter=' ', fmt='%d')
        np.savetxt('D.txt', D, delimiter=' ', fmt='%1f')

        #link_list = np.reshape(link_list, E*2)
        #print(link_list.shape)
      
        _geo_model1(n_steps, tolerance, A, D, link_list, N, E)
        
        
# #             print("Second Link List", second_link_index, k,l)
#         
#         # gives id for link i<->k resp. j<->l in original_link_ids
#         id1, i, k = self.link_id(i, k, N)
#         id2, j, l = self.link_id(j, l, N)
#         
#         link_list[first_link_index] = [i, k]
#         link_list[second_link_index] = [j, l]
#         sur_link_ids[first_link_index] = id1
#         sur_link_ids[second_link_index] = id2
#         
#         c += 1
#         if c == compute_at:
#             g = igraph.Graph(link_list.tolist())
#             x.append(u)
#             T.append(g.transitivity_avglocal_undirected())
#             H.append(self._Hamming_Distance(original_link_ids, sur_link_ids))
#             L.append(g.average_path_length())
#             ass.append(g.assortativity_degree())
#             c = 0
#             print(c,compute_at)
    
    
        dic = {"x": x, "T": T, "L": L, "H": H, "links": link_list,
               "assortativity": ass}
        
        # Update Adjcancy matrix
        self.adjacency=A
            
        
        return link_list, dic

    
    
    
    def GeoModel1_slow(self, n_steps, tolerance, grid_type="spherical", verbose=False):
        
        if grid_type == "spherical":
            D = self.grid.angular_distance()
        elif grid_type == "euclidean":
            D = self.grid.euclidean_distance()
        else:
            print("ERROR, This grid type is not recognized: ",  grid_type)
            sys.exit(1)
        
        #check_requirements(link_list, n_steps, n_steps)
        #check_GeoModel_requirements(grid, tolerance, grid_type)
        
        # Get an array of all links between all nodes in terms of [...,[i,j],...] for all Nodes
        link_list=np.array(self.graph.get_edgelist(), np.int32).copy(order='c')        
        
        #print(link_list)
        
        A= self.adjacency.copy(order='c')
        
        n_sampling_points=self.n_links
        
        #  Get number of nodes
        N = self.N
        #  Get number of links
        E = self.n_links
        print("len link list:" , len(link_list), E)

    
        original_link_ids = -0.5 * link_list[:, 0] * (link_list[:, 0] - 2 * N + 1)
        original_link_ids += link_list[:, 1] - link_list[:, 0] - 1
        
        #print("Original_link_ids", link_list[:,0], link_list[:,1], original_link_ids)
        sur_link_ids = original_link_ids.copy()
    
        g = igraph.Graph(link_list.tolist())
    
        T = [g.transitivity_avglocal_undirected()]
        x = [0]
        H = [0]
        L = [g.average_path_length()]
        ass = [g.assortativity_degree()]
    
        c = 0
        q = 0
        for u in range(n_steps):
            if verbose:
                print ("Step {0}/{1}".format(u+1, n_steps))
    
            cond = True
            while cond:
                q += 1
    
                first_link_index = np.random.randint(E)
                active_link = link_list[first_link_index]
                active_link = np.random.permutation(active_link)
                i, j = active_link
#                 print('First Link index', first_link_index, active_link, i, j)

                # If second argument is None, distance to any neighbor node is calculated
                d_i_all = self.distance( D, i, None  )
                D_tot = d_i_all - d_i_all[j]
    
                mask = np.abs(D_tot) < tolerance * d_i_all[j]
                mask[i] = False
                mask[j] = False
    
                possible_nbs = np.arange(N)[mask]
                possible_nbs = np.random.permutation(possible_nbs)
    
                l = None
    
                for k in possible_nbs:
                    #print(link_list == k)
                    nbs_of_k = np.fliplr(link_list == k)
                    nbs_of_k = link_list[nbs_of_k]
                    nbs_of_k = np.array(list(set(nbs_of_k) - set([k])))
                    
                    #print(nbs_of_k)
                    if i in nbs_of_k or len(nbs_of_k) == 0:
                        continue
    
                    d_k_all = self.distance(D, k, nbs_of_k)
                    d_j_all = self.distance(D, j, nbs_of_k)
                    
                    #print('Lengths', len(d_k_all), len(d_j_all))
                    D_tot = d_k_all - d_j_all
                    mask = np.abs(D_tot) < tolerance * d_k_all
    
                    if mask.any():
                        l_candidate = choice(nbs_of_k[mask])
                        nbs_of_l = np.fliplr(link_list == l_candidate)
                        nbs_of_l = link_list[nbs_of_l]
                        if j not in nbs_of_l:
                            l = l_candidate
                            break
    
                if l is None:
                    continue
                
                cond = False
    
            second_link_index = ((link_list == k) | (link_list == l))
            second_link_index = second_link_index.sum(axis=1) == 2
            second_link_index = np.arange(E)[second_link_index]
            
            A[i,j] =  A[j,i] = 0  # Delete link i<->j
            A[k,l] =  A[l,k] = 0  # Delete link k<->l
            A[i,k] =  A[k,i] = 1  # Add link i<->k
            A[j,l] =  A[l,j] = 1  # Add link j<->l
            
            
#             print("Second Link List", second_link_index, k,l)
            
            # gives id for link i<->k resp. j<->l in original_link_ids
            id1, i, k = self.link_id(i, k, N)
            id2, j, l = self.link_id(j, l, N)
            
            link_list[first_link_index] = [i, k]
            link_list[second_link_index] = [j, l]
            sur_link_ids[first_link_index] = id1
            sur_link_ids[second_link_index] = id2
            
            c += 1

            # compute_at = int(n_steps / n_sampling_points)
            if c % n_sampling_points == 0:
                g = igraph.Graph(link_list.tolist())
                x.append(u)
                T.append(g.transitivity_avglocal_undirected())
                H.append(self._Hamming_Distance(original_link_ids, sur_link_ids))
                L.append(g.average_path_length())
                ass.append(g.assortativity_degree())
                print(c)
    
        print ("# Total steps:", q)
    
        dic = {"x": x, "T": T, "L": L, "H": H, "links": link_list,
               "assortativity": ass}
    
        return link_list, dic
    
        print ("# Total steps:", q)
    
        dic = {"x": x, "T": T, "L": L, "H": H, "links": link_list,
               "assortativity": ass}
        
        # Update Adjcancy matrix
        self.adjacency=A
            
        
        return link_list, dic
    
    
        
    def GeoModel2(self, n_steps, grid_type="spherical", verbose=False):
    
#         check_requirements(link_list, n_steps, n_steps)
#         check_GeoModel_requirements(grid, tolerance, grid_type)
        if grid_type == "spherical":
            D = self.grid.angular_distance()
        elif grid_type == "euclidean":
            D = self.grid.euclidean_distance()
        else:
            print("ERROR, This grid type is not recognized: ",  grid_type)
            sys.exit(1)
     
     
        tolerance=self.tolerance
        # Needs to be done to prevent error in computation!
        link_list=np.array(self.graph.get_edgelist()).copy(order='c')        
#         link_list = link_list.copy()
    
        n_sampling_points=self.n_links    
        M = len(link_list)
        N = link_list.max() + 1
    
        original_link_ids = -0.5 * link_list[:, 0] * (link_list[:, 0] - 2 * N + 1)
        original_link_ids += link_list[:, 1] - link_list[:, 0] - 1
        sur_link_ids = original_link_ids.copy()
    
        compute_at = int(n_steps / n_sampling_points)
        g = igraph.Graph(link_list.tolist())
    
        T = [g.transitivity_avglocal_undirected()]
        x = [0]
        H = [0]
        L = [g.average_path_length()]
        ass = [g.assortativity_degree()]
    
        c = 0
        q = 0
        for u in range(n_steps):
            if verbose:
                print ("Step {0}/{1}".format(u+1, n_steps))
    
            cond = True
            while cond:
                q += 1
    
                first_link_index = np.random.randint(M)
                active_link = link_list[first_link_index]
                active_link = np.random.permutation(active_link)
                i, j = active_link
    
                d_i_all = self.distance(D, i, None )
                D_tot = d_i_all - d_i_all[j]
    
                mask = np.abs(D_tot) < tolerance * d_i_all[j]
                mask[i] = False
                mask[j] = False
    
                possible_nbs = np.arange(N)[mask]
                possible_nbs = np.random.permutation(possible_nbs)
    
                l = None
    
                for k in possible_nbs:
                    nbs_of_k = np.fliplr(link_list == k)
                    nbs_of_k = link_list[nbs_of_k]
                    nbs_of_k = np.array(list(set(nbs_of_k) - set([k])))
    
                    if i in nbs_of_k or len(nbs_of_k) == 0:
                        continue
    
                    d_k_all = self.distance(D, k, nbs_of_k)
                    d_j_all = self.distance(D, j, nbs_of_k)
                    D_tot = d_k_all - d_j_all
    
                    mask1 = np.abs(D_tot) < tolerance * d_k_all
                    mask2 = np.abs(d_k_all - d_i_all[j]) < tolerance * d_i_all[j]
                    mask = mask1 & mask2
    
                    if mask.any():
                        l_candidate = choice(nbs_of_k[mask])
                        nbs_of_l = np.fliplr(link_list == l_candidate)
                        nbs_of_l = link_list[nbs_of_l]
                        if j not in nbs_of_l:
                            l = l_candidate
                            break
    
                if l is None:
                    continue
    
                cond = False
    
            second_link_index = ((link_list == k) | (link_list == l))
            second_link_index = second_link_index.sum(axis=1) == 2
            second_link_index = np.arange(M)[second_link_index]
    
            id1, i, k = self.link_id(i, k, N)
            id2, j, l = self.link_id(j, l, N)
    
            link_list[first_link_index] = [i, k]
            link_list[second_link_index] = [j, l]
            sur_link_ids[first_link_index] = id1
            sur_link_ids[second_link_index] = id2
            c += 1
            if c == compute_at:
                g = igraph.Graph(link_list.tolist())
                x.append(u)
                T.append(g.transitivity_avglocal_undirected())
                H.append(self._Hamming_Distance(original_link_ids, sur_link_ids))
                L.append(g.average_path_length())
                ass.append(g.assortativity_degree())
                c = 0
    
        print ("# Total steps:", q)

        dic = {"x": x, "T": T, "L": L, "H": H, "links": link_list,
               "assortativity": ass}
    
        return link_list, dic
    
    