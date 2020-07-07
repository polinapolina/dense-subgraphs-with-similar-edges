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
    
def construct_write_graph_BL(node_link, node_link_deg, node_sim, node_sim_deg, node_map, filename = 'test_graph.txt', c = 0.0):
    maxcount = len(node_map)
    src = maxcount + 1
    sink = maxcount + 2
    
    pairs = set(node_sim.keys()) | set(node_link.keys())
    ne = len(pairs)*2 + len(node_map)*2
    
    file = open(filename, "w+")
    file.write('p par-max '+ str(len(node_map)+2) + ' ' + str(ne) + '\n')
    file.write('n ' + str(src) + ' s\n')
    file.write('n ' + str(sink) + ' t\n')  
    
    for (n1, n2) in pairs:
        s1 = node_link.get((n1, n2), 0)
        s2 = node_sim.get((n1, n2), 0)
        file.write('a ' + str(node_map[n1]) + ' ' + str(node_map[n2]) + ' ' + str(s1/2) + ' ' + str(s2/2) + '\n')
        file.write('a ' + str(node_map[n2]) + ' ' + str(node_map[n1]) + ' ' + str(s1/2) + ' ' + str(s2/2) + '\n')
    for n, i in node_map.items():
        file.write('a ' + str(src) + ' ' + str(i) + ' ' + str(c) + ' ' + str(c) + '\n')
        file.write('a ' + str(i) + ' ' + str(sink) + ' ' + str(node_link_deg.get(n, 0)/2) + ' ' + str(node_sim_deg.get(n, 0)/2) + '\n')
    total_sim1 = sum(node_link_deg.values())
    total_sim2 = sum(node_sim_deg.values())
    file.close()
    return total_sim1, total_sim2

    
def read_multi_layer_EOG(filename):
    edge_map = {}
    reverse_edge_map = {}
    edge_labels = {}
    
    edge_sim = {}
    edge_link = {}
    edge_sim_deg = {}
    edge_link_deg = {} 
    
    count = 1
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
            edge_labels[e1].add(layer)                

    edges = list(edge_map.keys())
    for in1 in range(len(edges)):
        for ie2 in range(in1+1, len(edges)):
            e1, e2 = edges[in1], edges[ie2]
            js = len(edge_labels[e1] & edge_labels[e2])/len(edge_labels[e1] | edge_labels[e2])
            if js > 0:
                pair = (e1, e2) if e1 < e2 else (e2, e1)
                edge_sim[pair] = js
                edge_sim_deg[e1] = edge_sim_deg.get(e1, 0) + js
                edge_sim_deg[e2] = edge_sim_deg.get(e2, 0) + js
            if set(e1) & set(e2):
                edge_link[pair] = 1
                edge_link_deg[e1] = edge_link_deg.get(e1, 0) + 1
                edge_link_deg[e2] = edge_link_deg.get(e2, 0) + 1
                
    return edge_sim, edge_sim_deg, edge_link, edge_link_deg, edge_map, reverse_edge_map     
    
 
    
def read_multi_layer_NOG(filename):
    node_map = {}
    reverse_node_map = {}
    
    node_sim = {}
    node_link = {}
    node_sim_deg = {}
    node_link_deg = {}    
    
    node_labels = {}  
   
    count = 1 
    with open(filename) as file:
        for line in file:
            chars = line.strip().split(' ')
            layer = int(chars[0])
            n1, n2 = int(chars[1]), int(chars[2])            
            
            if n1 not in node_map:
                node_labels[n1] = set()
                node_map[n1] = count
                reverse_node_map[count] = n1
                count += 1
            if n2 not in node_map:
                node_labels[n2] = set()
                node_map[n2] = count
                reverse_node_map[count] = n2
                count += 1
            pair = (n1, n2) if n1 < n2 else (n2, n1)
            node_link[pair] = 1
            node_labels[n1].add(layer)
            node_labels[n2].add(layer)
                
    nodes = list(node_map.keys())
    for in1 in range(len(nodes)):
        for in2 in range(in1+1, len(nodes)):
            n1, n2 = nodes[in1], nodes[in2]
            pair = (n1, n2) if n1 < n2 else (n2, n1)
            js = len(node_labels[n1] & node_labels[n2])/len(node_labels[n1] | node_labels[n2])
            if js > 0:                
                node_sim[pair] = js
                node_sim_deg[n1] = node_sim_deg.get(n1, 0) + js
                node_sim_deg[n2] = node_sim_deg.get(n2, 0) + js
            if pair in node_link:
                node_link_deg[n1] = node_link_deg.get(n1, 0) + 1
                node_link_deg[n2] = node_link_deg.get(n2, 0) + 1
                
    return node_sim, node_sim_deg, node_link, node_link_deg, node_map, reverse_node_map


if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='Construct a metagraph for the baselines from an input graph in format: "layerID nodeID nodeID weight"')
    parser.add_argument('--dataset', '-d', type=str, default='CS-Aarhus_multiplex', help='dataset filename (options: CS-Aarhus_multiplex)')
    parser.add_argument('--indir', '-i', type=str, default='Data', help='input data folder')
    parser.add_argument('--outdir', '-o', type=str, default='Data', help='output data folder')
    parser.add_argument('--baseline', '-b', type=str, default='BLDen', help='baseline name: BLDen or BLSim')

    args = parser.parse_args()
    filename = args.dataset
    kind = args.baseline

    pickled_file = kind + '_' + filename + '.p'
    meta_file = 'metagraph_' + kind + '_' + filename + '.txt'
    if kind == 'BLDen':
        print('read input file')
        node_sim, node_sim_deg, node_link, node_link_deg, node_map, reverse_node_map = read_multi_layer_NOG(os.path.join('.', args.indir, filename + '.edges'))
        pickle.dump((node_sim, node_sim_deg, node_link, node_link_deg, node_map, reverse_node_map), open(os.path.join('.', args.outdir, pickled_file), "wb" ))        
        print('construct and write metagraph')
        construct_write_graph_BL(node_link, node_link_deg, node_sim, node_sim_deg, node_map, filename = os.path.join('.', args.outdir, meta_file))
        
    elif kind == 'BLSim':
        print('read input file')
        edge_sim, edge_sim_deg, edge_link, edge_link_deg, edge_map, reverse_edge_map = read_multi_layer_EOG(os.path.join('.', args.indir, filename + '.edges'))        
        pickle.dump((edge_sim, edge_sim_deg, edge_link, edge_link_deg, edge_map, reverse_edge_map), open(os.path.join('.', args.outdir, pickled_file), "wb" ))  
        print('construct and write metagraph')
        construct_write_graph_BL(edge_sim, edge_sim_deg, edge_link, edge_link_deg, edge_map, filename = os.path.join('.', args.outdir, meta_file))
        
      
