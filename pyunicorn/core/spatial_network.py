#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# This file is part of pyunicorn.
# Copyright (C) 2008--2019 Jonathan F. Donges and pyunicorn authors
# URL: <http://www.pik-potsdam.de/members/donges/software>
# License: BSD (3-clause)
from pyunicorn.core import geo_network
import igraph.Graph as Graph
"""
Provides classes for analyzing spatial complex networks, that allow to generate 
random surrogates from a given spatially embedded network which can preserve 
global and local statistics associated with the nodes' embedding in a metric
space. 
This code is based on the work of Wiedermann et al. 2016 . 
"""

# array object and fast numerics
import numpy as np
# random number generation
from numpy import random

class spatialNetwork(geo_network):
    
    """
    Gives distances in spherical/euclidian networks from i to all surrounding nodes
    """
    def distance(self, node, grid, grid_type):
        if grid_type == "spherical":
            D = grid.angular_distance()
        elif grid_type == "euclidean":
            D = grid.euclidean_distance()
        else:
            print("ERROR, This grid type is not recognized: ",  grid_type)
        
        d_i_all=[]
        for i in range(len(D)):
            d_i_all.append(D[node][i])
        
        return d_i_all
                
    def GeoModel1(self, link_list, grid, n_steps, n_sampling_points, tolerance,
                  grid_type, verbose=False):
    
        check_requirements(link_list, n_steps, n_steps)
        check_GeoModel_requirements(grid, tolerance, grid_type)
    
        # Needs to be done to prevent error in computation!
        link_list = link_list.copy()
    
        M = len(link_list)
        N = link_list.max() + 1
    
        original_link_ids = -0.5 * link_list[:, 0] * (link_list[:, 0] - 2 * N + 1)
        original_link_ids += link_list[:, 1] - link_list[:, 0] - 1
        sur_link_ids = original_link_ids.copy()
    
        compute_at = int(n_steps / n_sampling_points)
        g = Graph(link_list.tolist())
    
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
    
                d_i_all = self.distance(grid[i], grid, grid_type=grid_type)
                D = d_i_all - d_i_all[j]
    
                mask = np.abs(D) < tolerance * d_i_all[j]
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

                d_k_all = self.distance(grid[k], grid[nbs_of_k], grid_type)
                d_j_all = self.distance(grid[j], grid[nbs_of_k], grid_type)
                D = d_k_all - d_j_all
                mask = np.abs(D) < tolerance * d_k_all

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

        id1, i, k = link_id(i, k, N)
        id2, j, l = link_id(j, l, N)

        link_list[first_link_index] = [i, k]
        link_list[second_link_index] = [j, l]
        sur_link_ids[first_link_index] = id1
        sur_link_ids[second_link_index] = id2
        c += 1
        if c == compute_at:
            g = Graph(link_list.tolist())
            x.append(u)
            T.append(g.transitivity_avglocal_undirected())
            H.append(_Hamming_Distance(original_link_ids, sur_link_ids))
            L.append(g.average_path_length())
            ass.append(g.assortativity_degree())
            c = 0

    print ("# Total steps:", q)

    dic = {"x": x, "T": T, "L": L, "H": H, "links": link_list,
           "assortativity": ass}

    return link_list, dic


def GeoModel2alt(link_list, grid, n_steps, n_sampling_points, tolerance,
                 grid_type, verbose=False):

    check_requirements(link_list, n_steps, n_steps)
    check_GeoModel_requirements(grid, tolerance, grid_type)

    # Needs to be done to prevent error in computation!
    link_list = link_list.copy()

    M = len(link_list)
    N = link_list.max() + 1

    original_link_ids = -0.5 * link_list[:, 0] * (link_list[:, 0] - 2 * N + 1)
    original_link_ids += link_list[:, 1] - link_list[:, 0] - 1
    sur_link_ids = original_link_ids.copy()

    compute_at = int(n_steps / n_sampling_points)
    g = Graph(link_list.tolist())

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

            d_i_all = distance(grid[i], grid, grid_type=grid_type)
            D = d_i_all - d_i_all[j]

            mask = np.abs(D) < tolerance * d_i_all[j]
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

                d_k_all = distance(grid[k], grid[nbs_of_k], grid_type)
                d_j_all = distance(grid[j], grid[nbs_of_k], grid_type)
                D = d_k_all - d_j_all

                mask1 = np.abs(D) < tolerance * d_k_all
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

        id1, i, k = link_id(i, k, N)
        id2, j, l = link_id(j, l, N)

        link_list[first_link_index] = [i, k]
        link_list[second_link_index] = [j, l]
        sur_link_ids[first_link_index] = id1
        sur_link_ids[second_link_index] = id2
        c += 1
        if c == compute_at:
            g = Graph(link_list.tolist())
            x.append(u)
            T.append(g.transitivity_avglocal_undirected())
            H.append(_Hamming_Distance(original_link_ids, sur_link_ids))
            L.append(g.average_path_length())
            ass.append(g.assortativity_degree())
            c = 0

    print ("# Total steps:", q)

    dic = {"x": x, "T": T, "L": L, "H": H, "links": link_list,
           "assortativity": ass}

    return link_list, dic