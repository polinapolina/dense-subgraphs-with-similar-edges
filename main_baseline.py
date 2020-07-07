from networkx import *
from networkx.generators.random_graphs import *
import random
import numpy as np
import matplotlib.pyplot as plt
from ctypes import cdll, c_double, c_int, byref
import sys
import os
import pickle
import time    

import argparse
    
INFINITY_NUMBER = 1.79769e+308
ALLOWED_ERROR = 1e-10
ALLOWED_DIFFERENCE = 1e-5

def print_edgeset(reverse_edge_map, srcSet):
    num_edges = len(reverse_edge_map)
    count = 0
    print('Edgelist: ', end='')
    for i in range(len(srcSet)):
        if srcSet[i] == 1:
            print(reverse_edge_map[i+1], end=' ')
            count+=1
    print('\nSize:', count)

def get_edgeset(reverse_edge_map, srcSet):
    num_edges = len(reverse_edge_map)
    count = 0
    out_edges = []
    #print('Edgelist: ', end='')
    for i in range(len(srcSet)):
        if srcSet[i] == 1:
            #if kind == 'BLDen':
            out_edges.append(reverse_edge_map[i+1])
            #elif kind == 'BLSim':
            #    out_edges.append(reverse_edge_map[i+1])
            #print(reverse_edge_map[i+1], end=' ')
            count+=1
    #print('\nSize:', count)
    return out_edges

def solve_MinCut_BL(num_elems, total_sim, precision = ALLOWED_ERROR, max_iters = 1000):    
    c = 0
    counter = 0
    
    srcSet_odd = (c_int*num_elems)()
    srcSet_even = (c_int*num_elems)()
    total_time = 0
    while counter < max_iters:
        ts = time.time()
        lib.c_pseudoflowPhase1()
        mincut_c = lib.c_getMinCutValue()
        #print('MC', mincut_c)
                
        F_elems = c_int()
        tmp = c_int()
        lib.c_getSizeOfMinCutSet(c_int(num_elems), byref(F_elems), byref(tmp))
        #print(F_elems, num_elems)
                
        if counter%2 == 0:           
            lib.c_getMinCutEdgeSet(c_int(num_elems), byref(srcSet_even))
        else:           
            lib.c_getMinCutEdgeSet(c_int(num_elems), byref(srcSet_odd))
        
        Q = -mincut_c + .5*total_sim

        #print(Q, mincut_c, .5*total_sim)
        if (Q < precision) or (F_elems.value == 0) or counter == max_iters-1:
            if  counter%2 == 0:
                total_time += time.time() - ts
                #print('Q-iterations:', counter, total_time)
                return srcSet_odd
            else:
                total_time += time.time() - ts
                #print('Q-iterations:', counter, total_time)
                return srcSet_even
            break
        
        newc_c = Q/F_elems.value + c
        #print('C:', newc_c)
        
        lib.c_updateSrcCapacities(c_double(newc_c), c_int(num_elems))
        
        c = newc_c 
        counter += 1
        total_time += time.time() - ts


def get_values_BLN(nodes, edge_map, edge_sim):
    edges = []
    for i1 in range(len(nodes)):
        for i2 in range(i1+1, len(nodes)):
            n1, n2 = nodes[i1], nodes[i2]
            pair = (n1, n2) if n1 < n2 else (n2, n1)
            if pair in edge_map:
                edges.append(pair)
    sim_res = 0
    for i1 in range(len(edges)):
        for i2 in range(i1+1, len(edges)):
            n1, n2 = edges[i1], edges[i2]
            pair = (n1, n2) if n1 < n2 else (n2, n1)
            if pair in edge_sim:
                sim_res += edge_sim[pair]
    return sim_res/len(edges), len(edges)/len(nodes)
            
def get_values_BLE(edges, edge_sim):
    #print(edges)
    nodes = set()
    for e in edges:        
        nodes.add(e[0])
        nodes.add(e[1])
    sim_res = 0
    for i1 in range(len(edges)):
        for i2 in range(i1+1, len(edges)):
            n1, n2 = edges[i1], edges[i2]
            pair = (n1, n2) if n1 < n2 else (n2, n1)
            if pair in edge_sim:
                sim_res += edge_sim[pair]
    return sim_res/len(edges), len(edges)/len(nodes)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='Run baseline for 11 gamma values in range [0,1]')
    parser.add_argument('--dataset', '-d', type=str, default='CS-Aarhus_multiplex', help='dataset filename (options: CS-Aarhus_multiplex)')
    parser.add_argument('--indir', type=str, default='Data', help='input data folder')
    parser.add_argument('--baseline', '-b', type=str, default='BLSim', help='baseline name: BLDen or BLSim')

    args = parser.parse_args()
    filename = args.dataset
    kind = args.baseline

    pickled_file = kind + '_' + filename + '.p'
    meta_file = 'metagraph_' + kind + '_' + filename + '.txt'
    if kind == 'BLDen':
        node_sim, node_sim_deg, node_link, node_link_deg, node_map, reverse_node_map = pickle.load(open(os.path.join('.', args.indir, pickled_file), "rb" ))
        reverse_map = reverse_node_map
        total_sim1 = sum(node_link_deg.values())
        total_sim2 = sum(node_sim_deg.values())
        num_elems = len(node_map)  
        max_mu = len(node_map)
        
    elif kind == 'BLSim':        
        edge_sim, edge_sim_deg, edge_link, edge_link_deg, edge_map, reverse_edge_map = pickle.load(open(os.path.join('.', args.indir, pickled_file), "rb" ))
        reverse_map = reverse_edge_map
        total_sim1 = sum(edge_sim_deg.values())
        total_sim2 = sum(edge_link_deg.values())
        num_elems = len(edge_map)  
        max_mu = len(edge_map)        
    
    lib = cdll.LoadLibrary('./bin/lib_pseudopar_baseline.so')
    lib.c_readDimacsFileCreateList(bytes(os.path.join('.', args.indir, meta_file), "utf8"))
    print('read done')
    lib.c_simpleInitialization()    
    mincut_func = lib.c_getMinCutValue
    mincut_func.restype = c_double
    srcSetSize_func = lib.c_getSizeOfMinCutSet
    srcSetSize_func.restype = c_int
    
    
    c = 0
    max_mu = 10
    for mu in np.linspace(0, max_mu, 101):
        lib.c_reCreateGraph(c_int(num_elems), c_double(mu), c_double(c))
        srcSet = solve_MinCut_BL(num_elems, total_sim1 + mu*total_sim2, precision = ALLOWED_ERROR, max_iters = 1000)
        #print_edgeset(reverse_map, srcSet)
        out_elems = get_edgeset(reverse_map, srcSet)

        sim, _ , edge_map, _, _ = pickle.load(open(os.path.join('.', args.indir, 'metagraph_' + filename + '.p'), "rb" ))
        if kind == 'BLDen':
            s, d = get_values_BLN(out_elems, edge_map, sim)
        elif kind == 'BLSim':
            s, d = get_values_BLE(out_elems, sim)
        print('gamma: {:10}, similarity: {:10}, density: {:10}'.format(mu, s, d))
    
    lib.c_finalfreeMemory()