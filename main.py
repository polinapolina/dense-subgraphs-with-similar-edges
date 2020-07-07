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

def lambda_search(l_max, l_min, l_delta, max_iters, reverse_edge_map, printedeges=False):
    P = []
    l_l = l_min
    l_u = l_max
    queue = []
    total_time_search = 0
    t1 = time.time()
    sim_l, den_l, srcSet = solve_MinCut(l_l) 
    total_time_search += time.time() - t1
    print('New solution found:')
    print('lambda: {:10}, similarity: {:10}, density: {:10}'.format(l_l, sim_l, den_l))
    if printedeges:
        print_edgeset(reverse_edge_map, srcSet)
    t1 = time.time()    
    sim_u, den_u, srcSet = solve_MinCut(l_u)  
    total_time_search += time.time() - t1
    
    iterations = 2
    if abs(sim_l - sim_u) > ALLOWED_DIFFERENCE or abs(den_l - den_u) > ALLOWED_DIFFERENCE:
        print('New solution found:')
        print('lambda: {:10}, similarity: {:10}, density: {:10}'.format(l_u, sim_u, den_u))
        if printedeges:
            print_edgeset(reverse_edge_map, srcSet)
        queue.append((l_l, l_u, sim_l, sim_u, den_l, den_u))
    while queue and iterations < max_iters:
        l_l, l_u, sim_l, sim_u, den_l, den_u = queue.pop(0)
        l_m = (l_l + l_u)/2
        t1 = time.time()
        sim_m, den_m, srcSet = solve_MinCut(l_m)
        total_time_search += time.time() - t1
        iterations += 1
        
        distinct_l = abs(sim_m - sim_l) > ALLOWED_DIFFERENCE or abs(den_m - den_l) > ALLOWED_DIFFERENCE        
        if distinct_l and l_m - l_l > l_delta:
            queue.append((l_l, l_m, sim_l, sim_m, den_l, den_m))
            
        distinct_u = abs(sim_m - sim_u) > ALLOWED_DIFFERENCE or abs(den_m - den_u) > ALLOWED_DIFFERENCE           
        if distinct_u and l_u - l_m > l_delta:
            queue.append((l_m, l_u, sim_m, sim_u, den_m, den_u))
            
        if distinct_l and distinct_u:
            print('New solution found:')
            print('lambda: {:10}, similarity: {:10}, density: {:10}'.format(l_m, sim_m, den_m))
            if printedeges:
                print_edgeset(reverse_edge_map, srcSet)
            
        sys.stdout.flush()
    print('Lambda-search stats:', iterations, 'iterations', total_time_search)
    print('iterations:', iterations, '; total time:', total_time_search)

    
def solve_MinCut(lmbda, precision = ALLOWED_ERROR, max_iters = 1000):   
    #print('Lambda:', lmbda)
    c = -lmbda*num_nodes
    lib.c_reCreateGraph(c_int(num_edges), c_double(lmbda), c_double(c))
    
    counter = 0
    srcSet_odd = (c_int*num_edges)()
    srcSet_even = (c_int*num_edges)()
    cur_sim, cur_den = -1, -1
    old_sim, old_den = -10, -10
    total_time = 0
    while counter < max_iters:
        ts = time.time()
        lib.c_pseudoflowPhase1()
        mincut_c = lib.c_getMinCutValue()       
                
        F_edges = c_int()
        F_nodes = c_int()
        lib.c_getSizeOfMinCutSet(c_int(num_edges), byref(F_edges), byref(F_nodes))
                
        if counter%2 == 0:           
            lib.c_getMinCutEdgeSet(c_int(num_edges), byref(srcSet_even))
        else:           
            lib.c_getMinCutEdgeSet(c_int(num_edges), byref(srcSet_odd))
        
        Q = -mincut_c + .5*total_sim

        if (Q < precision) or (F_edges.value == 0) or counter == max_iters-1:
            if  counter%2 == 0:
                total_time += time.time() - ts
                #print('Q-iterations:', counter, total_time)
                return cur_sim, cur_den, srcSet_odd
            else:
                total_time += time.time() - ts
                #print('Q-iterations:', counter, total_time)
                return cur_sim, cur_den, srcSet_even
            break         
     
        newc_c = Q/F_edges.value + c
        cur_sim, cur_den = c + (Q + lmbda*F_nodes.value)/F_edges.value, F_edges.value/F_nodes.value
        
        lib.c_updateSrcCapacities(c_double(newc_c), c_int(num_edges))
        
        c = newc_c 
        counter += 1
        total_time += time.time() - ts

    
if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='Run DenSim')
    parser.add_argument('--dataset', '-d', type=str, default='CS-Aarhus_multiplex', help='dataset filename (options: CS-Aarhus_multiplex)')
    parser.add_argument('--indir', '-i', type=str, default='Data', help='data folder of a metagraph')
    parser.add_argument('--printedeges', '-p', action='store_true', help='print solutions edgesets')
    
    args = parser.parse_args()
    filename = args.dataset
    
    sim, simdegree, edge_map, node_map, reverse_edge_map = pickle.load(open(os.path.join('.', args.indir, 'metagraph_' + filename + '.p'), "rb" ))
    
    max_num_lambdas = INFINITY_NUMBER
   
    total_sim = sum(simdegree.values())
    num_edges = len(edge_map)
    num_nodes = len(node_map)
    
    lib = cdll.LoadLibrary('./bin/lib_pseudopar.so')
    filename = 'metagraph_' + filename + '.txt'
    lib.c_readDimacsFileCreateList(bytes(os.path.join('.', args.indir, filename), "utf8"))
    lib.c_simpleInitialization()
    mincut_func = lib.c_getMinCutValue
    mincut_func.restype = c_double
    srcSetSize_func = lib.c_getSizeOfMinCutSet
    srcSetSize_func.restype = c_int    
    
    ds = 0.001
    l_delta = ds/(num_edges*num_edges)
    lambda_search(l_max = 1000000, l_min = 0, l_delta = l_delta, max_iters = max_num_lambdas, reverse_edge_map = reverse_edge_map)
    
    
    lib.c_finalfreeMemory()