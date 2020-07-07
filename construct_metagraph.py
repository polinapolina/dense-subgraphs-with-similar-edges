from networkx import *
from networkx.generators.random_graphs import *
import random
import numpy as np
import matplotlib.pyplot as plt
from ctypes import cdll, c_double, c_int, byref
import sys
import os
import pickle

import argparse
    
infinitynumber = 1.79769e+308

def construct_graph(sim, simdegree,  edge_map, node_map, lmbda = 1.0, c = 1.0):
    G = nx.DiGraph()
    maxcount = len(edge_map) + len(node_map)
    src = maxcount + 1
    sink = maxcount + 2
    maxcap = infinitynumber
    for (e1, e2), s in sim.items():
        G.add_edge(edge_map[e1], edge_map[e2], cap = s/2)
        G.add_edge(edge_map[e2], edge_map[e1], cap = s/2)    
    for e, ie in edge_map.items():
        G.add_edge(src, ie, cap = c)
        G.add_edge(ie, sink, cap = simdegree.get(e, 0)/2)
        G.add_edge(node_map[e[0]], edge_map[e], cap = maxcap)
        G.add_edge(node_map[e[1]], edge_map[e], cap = maxcap)
    for u, iu in node_map.items():
        G.add_edge(src, iu, cap = lmbda)
    return G, src, sink
    
def write_graph_file(G, s, t, filename = 'test_graph.txt'):
    file = open(filename, "w+")
    file.write('p par-max '+ str(G.number_of_nodes()) + ' ' + str(G.number_of_edges()) + '\n')
    file.write('n ' + str(s) + ' s\n')
    file.write('n ' + str(t) + ' t\n')
    for u,v in G.edges():
        file.write('a ' + str(u) + ' ' + str(v) + ' ' + str(G[u][v]['cap']) + '\n')
    file.close() 
    
    
def read_multi_layer(filename):
    edge_map = {}
    node_map = {}
    reverse_edge_map = {}
    reverse_node_map = {}
    edge_labels = {}
    
    sum = {}
    stddev = {}
    exp = {}
    count = 1
    sumxy = {} 
    sim = {}
    simdegree = {}
    nodes = set()
    with open(filename) as file:
        for line in file:
            chars = line.strip().split(' ')
            layer = int(chars[0])
            n1, n2 = int(chars[1]), int(chars[2])            
            e1 = (n1, n2) if n1 < n2 else (n2, n1)            
            if e1 not in edge_map:
                edge_labels[e1] = set()
                edge_map[e1] = count
                reverse_edge_map[count] = e1
                count += 1
                node_map[e1[0]] = -1
                node_map[e1[1]] = -1
            edge_labels[e1].add(layer)
                
    for u in node_map:
        node_map[u] = count
        count += 1
                
    edges = list(edge_map.keys())
    for ie1 in range(len(edges)):
        for ie2 in range(ie1+1, len(edges)):
            e1, e2 = edges[ie1], edges[ie2]
            js = len(edge_labels[e1] & edge_labels[e2])/len(edge_labels[e1] | edge_labels[e2])
            if js > 0:
                pair = (e1, e2) if e1 < e2 else (e2, e1)
                sim[pair] = js
                simdegree[e1] = simdegree.get(e1, 0) + js
                simdegree[e2] = simdegree.get(e2, 0) + js
                
    return sim, simdegree, edge_map, node_map, reverse_edge_map
    
    
def write_metagraph(filename, indir='Data', outdir='Data'):
    print('read input file')
    sim, simdegree, edge_map, node_map, reverse_edge_map = read_multi_layer(os.path.join('.', indir, filename + '.edges'))
        
    lmbda = 0.0
    print('construct metagraph')
    H, src, sink = construct_graph(sim, simdegree,  edge_map, node_map, lmbda = lmbda, c = 0.0)
    filename_meta = 'metagraph_' + filename + '.txt'
    print('write metagraph')
    write_graph_file(H, src, sink, os.path.join('.', 'Data', filename_meta))
    filename_pickle = 'metagraph_' + filename + '.p'
    pickle.dump((sim, simdegree, edge_map, node_map, reverse_edge_map), open( os.path.join('.', outdir, filename_pickle), "wb" ))
    
    
if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='Construct a metagraph from an input graph in format: "layerID nodeID nodeID weight"')
    parser.add_argument('--dataset', '-d', type=str, default='CS-Aarhus_multiplex', help='dataset filename (options: CS-Aarhus_multiplex)')
    parser.add_argument('--indir', '-i', type=str, default='Data', help='input data folder')
    parser.add_argument('--outdir', '-o', type=str, default='Data', help='output data folder')
    
    args = parser.parse_args()
    write_metagraph(args.dataset, args.indir, args.outdir)