
# -----------------------------------------------------------------
#	Functions needed to perform different tasks inside the update routine
# -----------------------------------------------------------------
import networkx as nx                                            
import os
import numpy as np

def remove_zero_entries_v(A):
	"INPUT:  Multilayer graph A"
	"OUTPUT: list with INT INDECES of nodes having zero total in_degree over the all layers "
	L=len(A)  # number of layers
	zero_in=[]
	nodes=list(A[0].nodes())
	for i in nodes :  # cycle over nodes
		k=0 
		for l in range(L):
			if(type(A[l].in_degree(i))!=dict):k+=A[l].in_degree(i)	
		if(k>0):zero_in.append(nodes.index(i))	
	return zero_in	

def remove_zero_entries_u(A):
	"INPUT:  Multilayer graph A"
	"OUTPUT: list with INT INDECES of nodes having zero total out_degree over the all layers "
	L=len(A)  # number of layers
	zero_out=[]
	nodes=list(A[0].nodes())
	for i in nodes:  # cycle over nodes
		k=0 
		for l in range(L):
			if(type(A[l].out_degree(i))!=dict):k+=A[l].out_degree(i)
		if(k>0):zero_out.append(nodes.index(i))	
	return zero_out

def remove_zero_entries_undirected(A):
	"INPUT:  Multilayer UNDIRECTED graph A"
	"OUTPUT: list with INT INDECES of nodes having zero total degree over the all layers "
	L=len(A)  # number of layers
	zero_in=[]
	nodes=list(A[0].nodes())
	for i in nodes:  # cycle over nodes
		k=0 
		for l in range(L):
			if(type(A[l].degree(i))!=dict):k+=A[l].degree(i)
			
		if(k>0):zero_in.append(nodes.index(i))	
	return zero_in	

def idx(i,A):
	" Adds node i to all layers"
	" returns node index "
	L=len(A)
	if(i not in list(A[0].nodes())):
		for l in range(L):A[l].add_node(i)
	
	#return A[0].nodes().index(i)

def read_graph(folder,adjacency_file,A):
	assert( os.path.isfile(folder+adjacency_file) and os.access(folder+adjacency_file, os.R_OK))
   	print "Adjacency file =",folder+adjacency_file
   	infile=open(folder+adjacency_file,'r');
   	nr=0
   	L=len(A)
   	for line in infile:
   		a=line.strip('\n').split()
   		if(a[0]=="E"):  # Flag to check the entry is an edge
	   		if(nr==0): # check format file is ok
	   			l=len(a)-3
	   			assert(l==L)
	   		v1=a[1]
	   		v2=a[2]
	   		idx(v1,A)
	   		idx(v2,A)
	   		for l in range(L):
	   			is_edge=int(a[l+3])
	   			if(is_edge>0):A[l].add_edge(v1, v2, weight=is_edge)
	infile.close()   			

def print_graph_stat(A, undirected=False):
	L=len(A);N=A[0].number_of_nodes()
	print "N=",N
	for l in range(L):
		B=nx.to_numpy_matrix(A[l],weight='weight')
		if undirected==False: E=np.sum(B)
		else:E=0.5*np.sum(B)
		print 'E[',l,']=',E," density=",100*float(E)/float(N*(N-1))

def out_graph(folder,A):
	L=len(A)
	for a in range(L):
		outfile=folder+"out_adjacency_"+str(a)+".dat";
		outf=open(outfile,'w')
		print "Adjacency of layer ",a," output in: ",outfile
		for e in A[a].edges():
			i=e[0]
			j=e[1]
			print >> outf,i,j
		outf.close()
		
def can_cast(string):
    try:
        int(string)
        return True
    except ValueError:
        return False
