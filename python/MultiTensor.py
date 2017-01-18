"""
Poisson Tensor factorization for Multi-layer networks.
"""
import time
import sys
import numpy as np
from numpy.random import RandomState
import tools as tl

class MultiTensor :
	def __init__(self,N=100,L=1,K=2, N_real=1,tolerance=0.1,decision=10,maxit=500,rseed=0,out_adjacency=False,inf=1e10,err_max=0.00001,err=0.1,initialization=0,undirected=False,folder="data/",end_file="",adj="adjacency.dat",w_file="w.dat",assortative=False):	
		self.N=N
		self.L=L
		self.K=K
		self.N_real=N_real
		self.tolerance=tolerance
		self.decision=decision
		self.maxit=maxit
		self.rseed=rseed
		self.out_adjacency=out_adjacency
		self.inf=inf
		self.err_max=err_max
		self.err=err
		self.initialization=initialization
		self.undirected=undirected
		self.folder=folder
		self.end_file=end_file
		self.adj=adj
		self.w_file=w_file
		self.assortative=assortative

		# Values used inside the update
		self.u=np.zeros((self.N,self.K),dtype=float)  # Out-going membership
		self.v=np.zeros((self.N,self.K),dtype=float)  # In-going membership

		# Old values
		self.u_old=np.zeros((self.N,self.K),dtype=float)  # Out-going membership
		self.v_old=np.zeros((self.N,self.K),dtype=float)  # In-going membership
		# Final values after convergence --> the ones that maximize Likelihood
		self.u_f=np.zeros((self.N,self.K),dtype=float)  # Out-going membership
		self.v_f=np.zeros((self.N,self.K),dtype=float)  # In-going membership

		if(self.assortative==True): # Purely diagonal matrix
			self.w=np.zeros((self.K,self.L),dtype=float)  # Affinity Matrix
			self.w_old=np.zeros((self.K,self.L),dtype=float)  # Affinity Matrix
			self.w_f=np.zeros((self.K,self.L),dtype=float)  # Affinity Matrix
		else:
			self.w=np.zeros((self.K,self.K,self.L),dtype=float)  # Affinity Matrix
			self.w_old=np.zeros((self.K,self.K,self.L),dtype=float)  # Affinity Matrix
			self.w_f=np.zeros((self.K,self.K,self.L),dtype=float)  # Affinity Matrix


	def _randomize_w(self,rng):
		" Assign a random number in (0,1.) to each entry"
		for i in range(self.L):
			for k in range(self.K):
				if(self.assortative==True):self.w[k,i]=rng.random_sample(1)
				else:
					for q in range(k,self.K):
						if(q==k):self.w[k,q,i]=rng.random_sample(1)
						else: self.w[k,q,i]=self.w[q,k,i]=self.err*rng.random_sample(1)


	def _randomize_u_v(self,rng,u_list,v_list):
		" Randomize the memberships' entries different from zero"				
		rng=np.random.RandomState(self.rseed)   # Mersenne-Twister random number generator
		for k in range(self.K):
			for i in range(len(u_list)):
				j=u_list[i]
				self.u[j][k]=rng.random_sample(1)
				if(self.undirected==True):self.v[j][k]=self.u[j][k]
			if(self.undirected==False):
				for i in range(len(v_list)):
					j=v_list[i]
					self.v[j][k]=rng.random_sample(1)
			
				

	def _initialize_w(self,rng,infile_name):
		" Initialize affinity matix from diagonal one extracted from file"	
		infile=open(infile_name,'r')
		nr=0
		for line in infile:
			if(nr>0):
				a=line.strip('\n').split()
				l=a[0]  # layer index
				assert(len(a)==self.K+1)
				for k in range(self.K):
					if(assortative==False):self.w[k][k][l]=float(a[k+1])
					else:self.w[k][l]=float(a[k+1])
		infile.close()			
		for l in range(self.L):
			for k in range(self.K):
				if(self.assortative==True):self.w[k][l]+=self.err*rng.random_sample(1)	
				else:
					for q in range(self.K):
						self.w[k][q][l]+=self.err*rng.random_sample(1)			



	def _initialize_u(self,rng,infile_name,nodes):
		" Initialize membership from file and nodes' list"	
		" INPUT 'nodes' is graph node list G.nodes() containing the labels"		
		infile=open(infile_name,'r')
		nr=0;
		max_entry=0.
		assert(len(nodes)==self.N)
		for line in infile:
			a= line.strip('\n').split(); # removes \n and split words in a list of lenght 1+K
			if(nr>0 and len(a)>0):
				assert(self.K==len(a)-1)
				if(a[0] in nodes):
					i=nodes.index(a[0])
					for k in range(self.K):
						z=float(a[k+1]);   # Value of the memebership for node a[0] and group k
						self.u[i][k]=z;
						max_entry=max(max_entry,z)
			nr+=1;  
		for n in range(self.N):
			for k in range(self.K):
				self.u[n][k]+=max_entry*self.err*rng.random_sample(1)    
		infile.close()    

	def _initialize_v(self,rng,infile_name,nodes):
		" Initialize membership from file and nodes' list"	
		" INPUT 'nodes' is graph node list G.nodes() containing the labels"		
		if(self.undirected==True):self.v=self.u;
		else:
			infile=open(infile_name,'r')
			nr=0;max_entry=0.;
			assert(len(nodes)==self.N)
			for line in infile:
				a= line.strip('\n').split(); # removes \n and split words in a list of lenght 1+K
				if(nr>0 and len(a)>0):
					assert(self.K==len(a)-1)
					if(a[0] in nodes):
						i=nodes.index(a[0])
						for k in range(self.K):
							z=float(a[k+1]);   # Value of the memebership for node a[0] and group k
							self.v[i][k]=z;
							max_entry=max(max_entry,z)
				nr+=1;  
			for n in range(self.N):
				for k in range(self.K):
					self.v[n][k]+=max_entry*self.err*rng.random_sample(1)      
			infile.close()   

	def _initialize(self,u_list,v_list,nodes):         
		
		rng=np.random.RandomState(self.rseed)   # Mersenne-Twister random number generator

		infile1=self.folder+'u_K'+str(self.K)+self.w_file
		infile2=self.folder+'v_K'+str(self.K)+self.w_file
		w_infile=self.folder+'w_K'+str(self.K)+self.w_file

		if(self.initialization==0):
			print " Random initializations"
			self._randomize_w(rng)	
			self._randomize_u_v(rng,u_list,v_list)			

		elif(self.initialization==1):
			print " W, U and V are initialized using: ";
			print infile1;
			print infile2;
			print w_infile;
			self._initialize_u(rng,infile1,nodes)
			self._initialize_v(rng,infile2,nodes)
			self._initialize_w(w_infile)

		elif(self.initialization==2):
			print " W initialized using: ";
			print w_infile;
			self._initialize_w(rng,w_infile)	
			self._randomize_u_v(rng,u_list,v_list)	
		
		elif(initialization==3):
			print " U and V are initialized using: ";
			print infile1;
			print infile2;
			self._randomize_w(rng)	
			self._initialize_u(rng,infile1,nodes)
			self._initialize_v(rng,infile2,nodes)

	def output_membership(self,nodes):		
		" INPUT 'nodes' is graph node list G.nodes() containing the labels"		
		print " u : ";
		for i in range(self.N):
			print nodes[i],
			for k in range(self.K):
				print self.u[i][k],
			print;	
		print;
		if(self.undirected==False):
			print " v : ";
			for i in range(self.N):
				print nodes[i],
				for k in range(self.K):
					print self.v[i][k],
				print;	
			
	def _output_affinity_matrix(self):
		print " W:";
		for l in range(self.L):
			if(self.assortative==False):
				print "a=",l;
				for k in range(self.K):
					for q in range(self.K):
						print self.w[k][q][l],
					print;
				print;	
			else:
				print l,	
				for k in range(self.K):	print self.w[k][l],	
				print;
		print;	

	def _update_old_variables(self,u_list,v_list):
		for i in range(len(u_list)):
			for k in range(self.K):
				self.u_old[u_list[i]][k]=self.u[u_list[i]][k]	
		for i in range(len(v_list)):
			for k in range(self.K):self.v_old[v_list[i]][k]=self.v[v_list[i]][k]	
		for l in range(self.L):
			for k in range(self.K):	
				if(self.assortative==True):self.w_old[k][l]=self.w[k][l]	
				else:
					for q in range(self.K):
						self.w_old[k][q][l]=self.w[k][q][l]	


	def _update_optimal_parameters(self):
		self.u_f=self.u;				
		self.v_f=self.v;				
		self.w_f=self.w;	

	def output_results(self,maxL,nodes):
		" Output results after convergence "

		# SORT node list if possible
		sorting=tl.can_cast(nodes[0])
		if(sorting==True):
			node_list=np.sort( [int(i) for i in nodes] )
			print "Sorting the membership vectors..."
		infile1=self.folder+"u_K"+str(self.K)+self.end_file
		infile3=self.folder+"w_K"+str(self.K)+self.end_file
		in1=open(infile1,'w')				
		in3=open(infile3,'w')
		print >>in1,"# Max Likelihood= ",maxL," N_real=",self.N_real				
		print >>in3,"# Max Likelihood= ",maxL," N_real=",self.N_real				
		if(self.undirected==False):
			infile2=self.folder+"v_K"+str(self.K)+self.end_file
			in2=open(infile2,'w')				
			print >>in2,"# Max Likelihood= ",maxL," N_real=",self.N_real				


		# Output membership

		if(sorting==True):
			for u in node_list:	
				i=nodes.index(str(u))
				print >>in1,u,
				if(self.undirected==False):print >> in2, u,
				for k in range(self.K):
					print >> in1,self.u_f[i][k],	
					if(self.undirected==False):print >> in2,self.v_f[i][k],	
				print >> in1;
				if(self.undirected==False):print >> in2;	
		else:	
			for i in range(self.N):
				print >> in1,nodes[i],		
				if(self.undirected==False):print >> in2, nodes[i],
				for k in range(self.K):
					print >> in1,self.u_f[i][k],	
					if(self.undirected==False):print >> in2,self.v_f[i][k],	
				print >> in1;
				if(self.undirected==False):print >> in2;	

		in1.close();

		if(self.undirected==False):in2.close();
		
		# Output affinity matrix
		for l in range(self.L):
			if(self.assortative==False):
				print >> in3, "a=",l;
				for k in range(self.K):
					for q in range(self.K):
						print >> in3, self.w_f[k][q][l],
					print >> in3;
				print >> in3;
			else:
				print >> in3,l,
				for k in range(self.K):	print >> in3, self.w_f[k][l],
				print >> in3;
		in3.close();
		
		self._output_affinity_matrix()	# output on screen				 

		print "Data saved in:";
		print infile1;print infile3;
		if(self.undirected==False):print infile2;


	# ----------	----------	----------	----------	----------	
	# ----------  Functions needed in the update_EM routine ----------	
	# ----------	----------	----------	----------	----------	

	def _update_U(self,A):

		Du=np.einsum('iq->q',self.v_old)
		if(self.assortative==False):
			w_k=np.einsum('kqa->kq',self.w_old)
			Z_uk=np.einsum('q,kq->k',Du,w_k)
			rho_ijka=np.einsum('jq,kqa->jka',self.v_old,self.w_old)
		else:
			w_k=np.einsum('ka->k',self.w_old)
			Z_uk=np.einsum('k,k->k',Du,w_k)
			rho_ijka=np.einsum('jk,ka->jka',self.v_old,self.w_old)
		
		rho_ijka=np.einsum('ik,jka->ijka',self.u,rho_ijka)

		Z_ija=np.einsum('ijka->ija',rho_ijka)
		Z_ijka=np.einsum('k,ija->ijka',Z_uk,Z_ija)

		non_zeros=Z_ijka>0.
		
		rho_ijka[non_zeros]/=Z_ijka[non_zeros]

		self.u=np.einsum('aij,ijka->ik',A,rho_ijka)
		low_values_indices = self.u < self.err_max  # Where values are low
		self.u[low_values_indices] = 0.  # All low values set to 0
		dist_u=np.amax(abs(self.u-self.u_old))	
		self.u_old=self.u

		return dist_u	

	def _update_V(self,A):

		Dv=np.einsum('iq->q',self.u_old)
		if(self.assortative==False):
			w_k=np.einsum('qka->qk',self.w_old)
			Z_vk=np.einsum('q,qk->k',Dv,w_k)
			rho_jika=np.einsum('jq,qka->jka',self.u_old,self.w_old)

		else:	
			w_k=np.einsum('ka->k',self.w_old)
			Z_vk=np.einsum('k,k->k',Dv,w_k)
			rho_jika=np.einsum('jk,ka->jka',self.u_old,self.w_old)

		rho_jika=np.einsum('ik,jka->jika',self.v,rho_jika)
		
		Z_jia=np.einsum('jika->jia',rho_jika)
		Z_jika=np.einsum('k,jia->jika',Z_vk,Z_jia)
		non_zeros=Z_jika>0.

		rho_jika[non_zeros]/=Z_jika[non_zeros]

		self.v=np.einsum('aji,jika->ik',A,rho_jika)

		low_values_indices = self.v < self.err_max  # Where values are low
		self.v[low_values_indices] = 0.  # All low values set to 0

		dist_v=np.amax(abs(self.v-self.v_old))	
		self.v_old=self.v

		return dist_v		

	def _update_W(self,A):

		if(self.assortative==False):
			uk=np.einsum('ik->k',self.u)
			vk=np.einsum('ik->k',self.v)
			Z_kq=np.einsum('k,q->kq',uk,vk)
			#Z_kq=np.einsum('ik,jq->kq',self.u,self.v)
			Z_ija=np.einsum('jq,kqa->jka',self.v,self.w_old)
		else:
			uk=np.einsum('ik->k',self.u)
			vk=np.einsum('ik->k',self.v)
			Z_k=np.einsum('k,k->k',uk,vk)
			#Z_k=np.einsum('ik,jk->k',self.u,self.v)
			Z_ija=np.einsum('jk,ka->jka',self.v,self.w_old)
		
		Z_ija=np.einsum('ik,jka->ija',self.u,Z_ija)

		B=np.einsum('aij->ija',A)
		non_zeros=Z_ija>0.
		Z_ija[non_zeros]=B[non_zeros]/Z_ija[non_zeros]

		rho_ijkqa=np.einsum('ija,ik->jka',Z_ija,self.u)
		
		if(self.assortative==False):
			rho_ijkqa=np.einsum('jka,jq->kqa',rho_ijkqa,self.v)
			rho_ijkqa=np.einsum('kqa,kqa->kqa',rho_ijkqa,self.w_old)
			self.w=np.einsum('kqa,kq->kqa',rho_ijkqa,1./Z_kq)
		else: 
			rho_ijkqa=np.einsum('jka,jk->ka',rho_ijkqa,self.v)
			rho_ijkqa=np.einsum('ka,ka->ka',rho_ijkqa,self.w_old)
			self.w=np.einsum('ka,k->ka',rho_ijkqa,1./Z_k)
		
		low_values_indices = self.w < self.err_max  # Where values are low
		self.w[low_values_indices] = 0.  # All low values set to 0

		dist_w=np.amax(abs(self.w-self.w_old))	
		self.w_old=self.w

		return dist_w		

	
	def _update_em(self,B):

		d_u=self._update_U(B)
		d_v=self._update_V(B)
		d_w=self._update_W(B)

		return d_u,d_v,d_w


	# --------------------------------------------------
	# Function needed to iterate
	# --------------------------------------------------				

	def _Likelihood(self,A):
		if(self.assortative==False):
			mu_ija=np.einsum('kql,jq->klj',self.w,self.v);
		else:	
			mu_ija=np.einsum('kl,jk->klj',self.w,self.v);
		mu_ija=np.einsum('ik,klj->lij',self.u,mu_ija);   
		l=-mu_ija.sum()
		non_zeros=A>0
		logM=np.log(mu_ija[non_zeros])
		Alog=A[non_zeros]*logM
		l+=Alog.sum()
		
		if(np.isnan(l)):
			print "Likelihood is NaN!!!!"
			sys.exit(1)
		else:return l			


	def _check_for_convergence(self,B,it,l2,coincide,convergence):
		if(it % 10 ==0):
			old_L=l2
			l2=self._Likelihood(B)	
			if(abs(l2-old_L)<self.tolerance): coincide+=1
			else: coincide=0
		if(coincide>self.decision):convergence=True	
		it+=1

		return it,l2,coincide,convergence	

	def cycle_over_realizations(self,A,B,u_list,v_list):
		maxL=-self.inf
		nodes=A[0].nodes()

		for r in range(self.N_real):
				
			self._initialize(u_list,v_list,nodes)
			
			self._update_old_variables(u_list,v_list)

			# Convergence local variables
			coincide=0
			convergence=False
			it=0
			l2=self.inf	
			delta_u=delta_v=delta_w=self.inf

			print "Updating r=",r," ..."
			tic=time.clock()
			# ------------------- Single step iteration update ------------------*/
			while(convergence==False and it<self.maxit):
				# Main EM update: updates membership and calculates max difference new vs old
				delta_u,delta_v,delta_w=self._update_em(B)

				it,l2,coincide,convergence=self._check_for_convergence(B,it,l2,coincide,convergence)
			print "r=",r," Likelihood=",l2," iterations=",it,' time=',time.clock()-tic,'s';
			if(l2>maxL): 
				self._update_optimal_parameters()
				maxL=l2
			self.rseed+=1	
		# end cycle over realizations
	
		print "Final Likelihood=",maxL

		self.output_results(maxL,A[0].nodes())	
				
	
