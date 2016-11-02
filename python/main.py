import MultiTensor as mt
import numpy as np
import networkx as nx
from argparse import ArgumentParser
import sys
import tools as tl
import time

def main():
	inf=10000000000000
	err_max=0.0000001
	p = ArgumentParser()
	p.add_argument('-f', '--folder', type=str, default='')
	p.add_argument('-a', '--adj', type=str, default='adjacency.dat')
	p.add_argument('-E', '--end_file', type=str, default='.dat')
	p.add_argument('-w', '--w_file', type=str, default='w.dat')
	p.add_argument('-l', '--L', type=int, default=4)
	p.add_argument('-i', '--initialization', type=int, default=0)
	p.add_argument('-k', '--K', type=int, default=5)
	p.add_argument('-r', '--N_real', type=int, default=1)
	p.add_argument('-t', '--maxit', type=int, default=500)
	p.add_argument('-e', '--tolerance', type=float,default=0.1)
	p.add_argument('-g', '--err', type=float,default=0.1)
	p.add_argument('-o','--out_adjacency',type=int,default=0)
	p.add_argument('-A','--assortative',type=int,default=0)
	p.add_argument('-u','--undirected',type=int,default=0)
	p.add_argument('-z', '--rseed', type=int, default=0)
	p.add_argument('-y','--decision',type=int,default=2)
	args = p.parse_args()
	
	folder="../data/"+args.folder
	if(args.undirected==True):A=[ nx.MultiGraph() for l in range(args.L) ]   # list of graphs
	else:A=[ nx.MultiDiGraph() for l in range(args.L) ]   # list of graphs

	tl.read_graph(folder,args.adj,A)
	print "Undirected=",bool(args.undirected)
	print "Assortative=",bool(args.assortative)
	tl.print_graph_stat(A)

	if(args.out_adjacency):tl.out_graph(folder,A)

	if(args.undirected==True): 
		u_list=v_list=tl.remove_zero_entries_undirected(A) 
	else:
		u_list=tl.remove_zero_entries_u(A)   # list of INT INDECES of nodes with zero out degree
		v_list=tl.remove_zero_entries_v(A)   # list of INT INDECES of nodes with zero in degree

	MT=mt.MultiTensor(  N=A[0].number_of_nodes(),
			L=args.L,K=args.K, 
			N_real=args.N_real,
			tolerance=args.tolerance,
			decision=args.decision,
			maxit=args.maxit,
			rseed=args.rseed,
			out_adjacency=bool(args.out_adjacency),
			inf=inf,
			err_max=err_max,
			err=args.err,
			initialization=args.initialization,
			undirected=bool(args.undirected),
			folder=folder,
			end_file=args.end_file,
			adj=args.adj,
			w_file=args.w_file,
			assortative=bool(args.assortative)
			)

	tic = time.clock()
	N=A[0].number_of_nodes()
	B=np.empty(shape=[args.L,N,N])

	for l in range(args.L):B[l,:,:]=nx.to_numpy_matrix(A[l],weight='weight')
	
	MT.cycle_over_realizations(A,B,u_list,v_list)		

	#tl.print_graph_stat(A)	

	toc = time.clock()	  
	print "It took ",toc-tic," seconds."; 	


if __name__ == '__main__':
	main()
