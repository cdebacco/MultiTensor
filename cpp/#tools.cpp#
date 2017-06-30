
#include "mlg.hpp"

#include <boost/random.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <math.h> 
#include <stdio.h> 
#include <fstream>

using namespace boost;
using namespace std;


boost::mt19937 gen ;   // random real and integer number generators
uniform_real<double> const uni_dist(0,1);    // uniform distribution for real numbers
variate_generator<boost::mt19937 &, uniform_real<> > real01(gen, uni_dist);   // picks a random real number in (0,1)
extern double err;  // defined in main.hpp

// Assign a random real number in (0,1) to each entry

void randomize_w(vector<vector< vector<double> > > & w) {
	int L=(int)w[0][0].size();
	int K=(int)w.size();	
  for(int i=0;i<L;i++)for(int k=0; k<K; k++ )for(int q=k; q<K; q++ ){
    if(q==k)w[k][q][i]=real01();  // for better convergence
    else w[k][q][i]=w[q][k][i]=real01();
  }
}


void randomize_v_in_i(vector< vector<double> > & v, vector<int> v_list) {
	int K=(int)v[0].size();	
 	int list_length=v_list.size();
    for (int k=0; k<K; k++ )for(int i=0;i<list_length;i++){ 
    	int j=v_list[i];
    	v[j][k] =real01();
    }
}


// Calculate in-degree, summing over all layers
int in_neigh(Vertex i, vector< Graph> A){
	int L=A.size();  // number of layers
	int k=0;
	in_edge_iterator eit, eend;   // Cycle over in-neighbors of i in layer a
	for(int a=0;a<L;a++)for(tie(eit, eend) = in_edges(i,A[a]); eit != eend; ++eit) k++;
	return k;
}

/*
int in_neigh(Vertex i, vector< Graph_undirected> A){
  int L=A.size();  // number of layers
  int k=0;
  in_edge_iterator_un eit, eend;   // Cycle over in-neighbors of i in layer a
  for(int a=0;a<L;a++)for(tie(eit, eend) = in_edges(i,A[a]); eit != eend; ++eit) k++;
  return k;
}*/

// Calculate out-degree, summing over all layers
int out_neigh(Vertex i, vector< Graph> A){
	int L=A.size();
	int k=0;
	edge_iterator eit, eend;   // Cycle over out-neighbors of i in layer a
	for(int a=0;a<L;a++)for(tie(eit, eend) = out_edges(i,A[a]); eit != eend; ++eit) k++;
	return k;
}

int out_neigh(Vertex i, vector< Graph_undirected> A){
  int L=A.size();
  int k=0;
  g_edge_iterator_und eit, eend;   // Cycle over out-neighbors of i in layer a
  for(int a=0;a<L;a++)for(tie(eit, eend) = out_edges(i,A[a]); eit != eend; ++eit) k++;
  return k;
}

// If out/in -degree is zero --> membership has zero entries
void set_zero_entries(vector< vector<double> > & u, vector< vector<double> > & v,  vector<Graph>  A){
	int N=(int)num_vertices(A[0]);
	int K=(int)v[0].size();
	for(int i=0;i<N;i++){
	  if(out_neigh(i,A)==0)for(int k=0;k<K;k++)u[i][k]=0.;
	  if(in_neigh(i,A)==0)for(int k=0;k<K;k++)v[i][k]=0.;
	}
}
// Outputs a vector of nodes that have at least one in-coming neighbor
vector<int> remove_zero_entries_v(  vector<Graph>  A){
  vector<int> v;
  int N=(int)num_vertices(A[0]);
  for(int i=0;i<N;i++)if(in_neigh(i,A)>0)v.push_back(i);
  return v;
}
// Outputs a vector of nodes that do NOT have any in-coming neighbor
vector<int> collect_zero_entries_v(  vector<Graph>  A){
  vector<int> v;
  int N=(int)num_vertices(A[0]);
  for(int i=0;i<N;i++)if(in_neigh(i,A)==0)v.push_back(i);
  return v;
}

// Outputs a vector of nodes that have at least one out-going neighbor
vector<int> remove_zero_entries_u(  vector<Graph>  A){
  vector<int> u;
  int N=(int)num_vertices(A[0]);
  for(int i=0;i<N;i++)if(out_neigh(i,A)>0)u.push_back(i);
  return u;	
}
vector<int> remove_zero_entries_u(  vector<Graph_undirected>  A){
  vector<int> u;
  int N=(int)num_vertices(A[0]);
  for(int i=0;i<N;i++)if(out_neigh(i,A)>0)u.push_back(i);
  return u; 
}
// Outputs a vector of nodes that do NOT have any out-going neighbor
vector<int> collect_zero_entries_u(  vector<Graph>  A){

  vector<int> u;
  int N=(int)num_vertices(A[0]);
  for(int i=0;i<N;i++)if(out_neigh(i,A)==0)u.push_back(i);
  return u;
}

vector<int> collect_zero_entries_u(  vector<Graph_undirected>  A){

  vector<int> u;
  int N=(int)num_vertices(A[0]);
  for(int i=0;i<N;i++)if(out_neigh(i,A)==0)u.push_back(i);
  return u;
}
/* Initialize the u,v, or w from input file ------------------------------------------------ */

// Input file is from a symmetric (assortative) matrix W
void initialize_w(vector<vector< vector<double> > > & w, istream & file) {
    int a=-1;
    int L=(int)w[0][0].size();
    int K=(int)w.size();
    for(string line; getline( file, line ); ){
      std::istringstream is( line );
      if(a>=0){
      double l;
      is>>l; // first number is the layer index
      for(int k=0;k<K;k++)is>>w[k][k][a];
      }
      a++;
    }
    for(int a=0;a<L;a++)for(int k=0;k<K;k++)for(int q=0;q<K;q++)w[k][q][a]+=err*real01();  // add some noise so to have non-zero entries out of the diagonal
}

//  Input file is the membership vector v or u
void initialize_v(vector< vector<double> > & v, istream & file) {
    int a=-1;double max_entry=0.;
    int K=(int)v[0].size();
    int N=v.size();
    for(string line; getline( file, line ); ){
      std::istringstream is( line );
      if(a>=0){
      string i;
      is>>i; // node label, is different from the node index a!
      for(int k=0;k<K;k++){is>>v[a][k];max_entry=std::max(max_entry,v[a][k]);}  // need to know magnitude of max_entry to add noise consequently
      }
      a++; // node index from 0 to N-1
    }
    // Add noise (perturbation)
    for(int i=0;i<N;i++)for(int k=0;k<K;k++)v[i][k]+=max_entry*err*real01();
}

void initialize(int initialization,  vector<int> u_list,  vector<int> v_list,vector< vector<double> > & u,vector< vector<double> > & v, vector< vector< vector<double> > > & w, string file, string w_file,vector<Graph> A){
	
	int K=(int)w.size();  // number of communities

	if(initialization==0){  // everything is random
	      cout<<" Random initializations"<<endl;
	      randomize_w(w);
	      randomize_v_in_i(v,v_list);  
	      randomize_v_in_i(u,u_list);  
    }
    else if(initialization==1){ // everything is initialized
      cout<< " W, U and V are initialized using: "<< file+"v_K"+to_string(K)+w_file<<" ";
      cout<<file+"u_K"+to_string(K)+w_file<<" "; 
      cout<<file+"w_K"+to_string(K)+w_file<<endl; 
      std::ifstream in_w;
      in_w.open((file+"w_K"+to_string(K)+w_file).c_str());
      initialize_w(w,in_w);
      in_w.close();
      in_w.open((file+"u_K"+to_string(K)+w_file).c_str());
      initialize_v(u,in_w);
      in_w.close();
      in_w.open((file+"v_K"+to_string(K)+w_file).c_str());
      initialize_v(v,in_w);
      in_w.close();
      
      vector<int> u_to_be_removed=collect_zero_entries_u(A); // Creates a vector of nodes that do NOT have any out-going neighbor
      vector<int> v_to_be_removed=collect_zero_entries_v(A); // Creates a vector of nodes that do NOT have any in-coming neighbor
      for(int i=0;i<u_to_be_removed.size();i++)for(int k=0;k<K;k++)u[u_to_be_removed[i]][k]=0.;
      for(int i=0;i<v_to_be_removed.size();i++)for(int k=0;k<K;k++)v[v_to_be_removed[i]][k]=0.;
    }
    else if(initialization==2){  // only w is initialized
      std::ifstream in_w;
      cout<< "W is initialized using: "<< file+"w_K"+to_string(K)+w_file<<endl; 
      in_w.open((file+"w_K"+to_string(K)+w_file).c_str());
      initialize_w(w,in_w);
      in_w.close();
      randomize_v_in_i(v,v_list);  
      randomize_v_in_i(u,u_list);
    }
    if(initialization==3){  // only u and v are initialized
      std::ifstream in_w;
      cout<< " U and V are initialized using: "<< file+"v_K"+to_string(K)+w_file<<" ";
      cout<<file+"u_K"+to_string(K)+w_file<<endl;
      in_w.open((file+"v_K"+to_string(K)+w_file).c_str());
      initialize_v(v,in_w);
      in_w.close();
      in_w.open((file+"u_K"+to_string(K)+w_file).c_str());
      initialize_v(u,in_w);
      in_w.close();
      vector<int> u_to_be_removed=collect_zero_entries_u(A); // Creates a vector of nodes that do NOT have any out-going neighbor
      vector<int> v_to_be_removed=collect_zero_entries_v(A); // Creates a vector of nodes that do NOT have any in-coming neighbor
      for(int i=0;i<u_to_be_removed.size();i++)for(int k=0;k<K;k++)u[u_to_be_removed[i]][k]=0.;
      for(int i=0;i<v_to_be_removed.size();i++)for(int k=0;k<K;k++)v[v_to_be_removed[i]][k]=0.;
      randomize_w(w);
    }

}

void initialize_undirected(int initialization,  vector<int> u_list, vector< vector<double> > & u,vector< vector< vector<double> > > & w, string file, string w_file,vector<Graph_undirected> A){
  
  int K=(int)w.size();  // number of communities

  if(initialization==0){  // everything is random
        cout<<" Random initializations"<<endl;
        randomize_w(w); 
        randomize_v_in_i(u,u_list);  
    }
    else if(initialization==1){ // everything is initialized
      cout<< " W and U  are initialized using: "<< file+"u_K"+to_string(K)+w_file<<" ";
      cout<<file+"w_K"+to_string(K)+w_file<<endl; 
      std::ifstream in_w;
      in_w.open((file+"w_K"+to_string(K)+w_file).c_str());
      initialize_w(w,in_w);
      in_w.close();
      in_w.open((file+"u_K"+to_string(K)+w_file).c_str());
      initialize_v(u,in_w);
      in_w.close();
      
      vector<int> u_to_be_removed=collect_zero_entries_u(A); // Creates a vector of nodes that do NOT have any out-going neighbor
      for(int i=0;i<u_to_be_removed.size();i++)for(int k=0;k<K;k++)u[u_to_be_removed[i]][k]=0.;
    }
    else if(initialization==2){  // only w is initialized
      std::ifstream in_w;
      cout<< "W is initialized using: "<< file+"w_K"+to_string(K)+w_file<<endl; 
      in_w.open((file+"w_K"+to_string(K)+w_file).c_str());
      initialize_w(w,in_w);
      in_w.close();
      randomize_v_in_i(u,u_list);
    }
    if(initialization==3){  // only u and v are initialized
      std::ifstream in_w;
      cout<< " U is initialized using: "<< file+"u_K"+to_string(K)+w_file<<endl;
      in_w.open((file+"u_K"+to_string(K)+w_file).c_str());
      initialize_v(u,in_w);
      in_w.close();
      vector<int> u_to_be_removed=collect_zero_entries_u(A); // Creates a vector of nodes that do NOT have any out-going neighbor
      for(int i=0;i<u_to_be_removed.size();i++)for(int k=0;k<K;k++)u[u_to_be_removed[i]][k]=0.;
      randomize_w(w);
    }

}

/* ---------- Graph related functions ----------------------------------------*/

// Add vertex with label 'id'; only if the vertex was not present in the graph yet
// returns the int index corresponding to the string input vertex label
int idx(string const & id, vector<Graph> & g)
{
  static map <string, int> idx_map;
  map<string, int>::iterator mit = idx_map.find(id);
  if (mit == idx_map.end()){  // Case where the node has not been already added into the graph
    for(int a=0;a<g.size();a++){   // Run through layers
      if(a==0)idx_map[id] = add_vertex(VertexProperties(id), g[a]);   // Update the map just once 
      else add_vertex(VertexProperties(id), g[a]);
    } 
    return idx_map[id] ;
  }
  return mit->second;
}

int idx(string const & id, vector<Graph_undirected> & g)
{
  static map <string, int> idx_map;
  map<string, int>::iterator mit = idx_map.find(id);
  if (mit == idx_map.end()){  // Case where the node has not been already added into the graph
    for(int a=0;a<g.size();a++){   // Run through layers
      if(a==0)idx_map[id] = add_vertex(VertexProperties(id), g[a]);   // Update the map just once 
      else add_vertex(VertexProperties(id), g[a]);
    } 
    return idx_map[id] ;
  }
  return mit->second;
}
//-------------------------------------------------------------------------------

/* INPUT multilayer adjacency file format:
 E 102 231 1 0 0 1     
 "102" and "231" are the source and target node (string) labels; 
" 1 0 0 1"  stand for: the direct edge  102 -> 231  exists in layer 1 and 4, not in layer 2 and 3; 
*/

int read_graph(string folder,string adj, vector<Graph> & A)
{  
	string file="../data/"+folder+adj;
	std::ifstream in1;
	in1.open((file).c_str()); 
	if(!in1.is_open()){
		cout<< " Error opening adjacency matrix file:"<<file<<endl;
		return -1;
	}
	cout<<"Infile = "<<file<<endl;

	int L=(int)A.size();
	string tok;
  	while (in1 >> tok){
	    if (tok == "E") {  // control variable "E"; to check that the input file containes the edge information
	      string i, j;  // source and target node labels
	      in1 >> i >> j ;
	      int v1=idx(i,A);  // Adds node i to all layers and returns corresponding int index
	      int v2=idx(j,A);  // Adds node j to all layers and returns corresponding int index
	      int is_edge;      // Flag saying if the edge exists
	      for(int alpha=0;alpha<L;alpha++){
	       in1>>is_edge;  // if this is 1 then creates an edge in the correspending layer
	        if(is_edge>0)for(int weigth=0;weigth<is_edge;weigth++)add_edge(vertex(v1, A[alpha]), vertex(v2, A[alpha]), A[alpha]);
	      }  // end cycle for over layers
	    }  // end if(tok=="E")
	  }  // end while
	in1.close();

	int N=(int)num_vertices(A[0]);
	return N;
 } 


int read_graph(string folder,string adj, vector<Graph_undirected> & A)
{  
  string file="../data/"+folder+adj;
  std::ifstream in1;
  in1.open((file).c_str()); 
  if(!in1.is_open()){
    cout<< " Error opening adjacency matrix file:"<<file<<endl;
    return -1;
  }
  cout<<"Infile = "<<file<<endl;

  int L=(int)A.size();
  string tok;
    while (in1 >> tok){
      if (tok == "E") {  // control variable "E"; to check that the input file containes the edge information
        string i, j;  // source and target node labels
        in1 >> i >> j ;
        int v1=idx(i,A);  // Adds node i to all layers and returns corresponding int index
        int v2=idx(j,A);  // Adds node j to all layers and returns corresponding int index
        int is_edge;      // Flag saying if the edge exists
        for(int alpha=0;alpha<L;alpha++){
         in1>>is_edge;  // if this is 1 then creates an edge in the correspending layer
          if(is_edge>0)for(int weigth=0;weigth<is_edge;weigth++)add_edge(vertex(v1, A[alpha]), vertex(v2, A[alpha]), A[alpha]);
        }  // end cycle for over layers
      }  // end if(tok=="E")
    }  // end while
  in1.close();

  int N=(int)num_vertices(A[0]);
  return N;
 } 

// Outputs single layers networks in the format of a directed edge list
void out_graph(string folder, vector<Graph> A)
{
	int L=(int)A.size();
	for(int a=0;a<L;a++){
	std::ofstream out;
	string outfile="../data/"+folder+"out_adjacency_"+to_string(a)+".dat";
	cout<< "Adjacency of layer "<<a<<" output in "<<outfile<<endl;
	out.open((outfile).c_str());  // can change it if you want a different output file
	Graph gr=A[a];
	graph_edge_iterator eit, eend;
	for (tie(eit, eend) = edges(gr); eit != eend; ++eit) {
	  Edge e=*eit;
	  Vertex j = target(e, gr); 
	  Vertex i = source(e, gr);  
	  out << gr[i].name << " " << gr[j].name << endl;
	} 
	out.close();
  }
}

void out_graph(string folder, vector<Graph_undirected> A)
{
  int L=(int)A.size();
  for(int a=0;a<L;a++){
  std::ofstream out;
  string outfile="../data/"+folder+"out_adjacency_"+to_string(a)+".dat";
  cout<< "Adjacency of layer "<<a<<" output in "<<outfile<<endl;
  out.open((outfile).c_str());  // can change it if you want a different output file
  Graph_undirected gr=A[a];
//  graph_edge_iterator eit, eend;
  edge_iterator_und eit, eend;
  for (tie(eit, eend) = edges(gr); eit != eend; ++eit) {
    Edge_und e=*eit;
    Vertex j = target(e, gr); 
    Vertex i = source(e, gr);  
    out << gr[i].name << " " << gr[j].name << endl;
  } 
  out.close();
  }
}
void output_membership(vector< vector<double> > v,vector<Graph> A)
{
	int N=(int)v.size();
	int K=(int)v[0].size();
	for(int i=0;i<N;i++){
		cout<<A[0][i].name<<" ";
		for(int k=0;k<K;k++)cout<<v[i][k]<<" ";
		cout<<endl;	
	}

}

void output_membership(vector< vector<double> > v,vector<Graph_undirected> A)
{
  int N=(int)v.size();
  int K=(int)v[0].size();
  for(int i=0;i<N;i++){
    cout<<A[0][i].name<<" ";
    for(int k=0;k<K;k++)cout<<v[i][k]<<" ";
    cout<<endl; 
  }

}

void output_affinity_matrix( vector< vector< vector<double> > > w)
{
	int L=(int)w[0][0].size();
	int K=(int)w.size();
	for(int l=0;l<L;l++){
		cout<<"a="<<l<<endl;
		for(int k=0;k<K;k++){
			for(int q=0;q<K;q++)cout<<w[k][q][l]<<" ";
			cout<<endl;	
		}
		cout<<endl;
	}
}
//  -------   Calculate Likelihood  ------------------------------------------------------------------------------------
double Likelihood(vector< vector<double> > & u, vector< vector<double> > & v,vector<vector< vector<double> > > & w, vector<Graph> & A){

  double l=0.;
  int L=(int)A.size();
  int K=(int)w.size();
  int N=(int)num_vertices(A[0]);
  // Calculate argument inside log
  for(int a=0;a<L;a++){
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
		    double log_arg=0.;
        for(int k=0;k<K;k++)for(int q=0;q<K;q++){
          double uvw=u[i][k]*v[j][q]*w[k][q][a];
          l-=uvw;   // Add this term regardeles of the value of A_ijk
    	    if(edge(i, j, A[a]).second)log_arg+=uvw;  // if edge exists, consider this term inside the log argument
  	    }// end cycle over k and q
  	    if(log_arg>0.){
          int c=0;  // count parallel edges
          edge_iterator eit, eend;   // Cycle over out-neighbors of i in layer a,  --> count parallel edges
          for (tie(eit, eend) = out_edges(i,A[a]); eit != eend; ++eit)if(target(*eit,A[a])==j)c++;
          l+=c*log(log_arg);
        } // end cycle to count parallel edges
      }// end cycle over j      
    }// end cycle over i
  }// end cycle over a
  if(isnan(l)){
  	cout<<"Likelihood is Nan!!!"<<endl;
  	exit(-1);
  }
  else return l;
}

double Likelihood_undirected(vector< vector<double> > & u,vector<vector< vector<double> > > & w, vector<Graph_undirected> & A){

  double l=0.;
  int L=(int)A.size();
  int K=(int)w.size();
  int N=(int)num_vertices(A[0]);
  // Calculate argument inside log
  for(int a=0;a<L;a++){
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
        double log_arg=0.;
        for(int k=0;k<K;k++)for(int q=0;q<K;q++){
          double uvw=u[i][k]*u[j][q]*w[k][q][a];
          l-=uvw;   // Add this term regardeles of the value of A_ijk
          if(edge(i, j, A[a]).second)log_arg+=uvw;  // if edge exists, consider this term inside the log argument
       }// end cycle over k and q
       if(log_arg>0.){
          int c=0;  // count parallel edges
          g_edge_iterator_und eit, eend;   // Cycle over out-neighbors of i in layer a,  --> count parallel edges
          for (tie(eit, eend) = out_edges(i,A[a]); eit != eend; ++eit)if(target(*eit,A[a])==j)c++;
          l+=c*log(log_arg);
        } // end cycle to count parallel edges
      }// end cycle over j      
    }// end cycle over i
  }// end cycle over a
  if(isnan(l)){
    cout<<"Likelihood is Nan!!!"<<endl;
    exit(-1);
  }
  else return l;
}
void update_old_variables(  vector< vector<double> > & u_old,  vector< vector<double> > & v_old,  vector< vector<double> > u,  vector< vector<double> > v,  vector< vector< vector<double> > > & w_old,  vector< vector< vector<double> > > w, vector<int> u_list,vector<int> v_list){
    int K=(int) w.size();
    int L=(int) w[0][0].size();
    int v_length=(int)v_list.size();
    int u_length=(int)u_list.size();
    for(int i=0;i<v_length;i++)for(int k=0;k<K;k++) v_old[v_list[i]][k]=v[v_list[i]][k]; 
    for(int i=0;i<u_length;i++)for(int k=0;k<K;k++) u_old[u_list[i]][k]=u[u_list[i]][k]; 
    for(int k=0;k<K;k++)for(int q=0;q<K;q++)for(int a=0;a<L;a++)w_old[k][q][a]=w[k][q][a];  
}

void update_old_variables_undirected(  vector< vector<double> > & u_old,  vector< vector<double> > u,  vector< vector< vector<double> > > & w_old,  vector< vector< vector<double> > > w, vector<int> u_list){
    int K=(int) w.size();
    int L=(int) w[0][0].size();
    int u_length=(int)u_list.size();
    for(int i=0;i<u_length;i++)for(int k=0;k<K;k++) u_old[u_list[i]][k]=u[u_list[i]][k]; 
    for(int k=0;k<K;k++)for(int q=0;q<K;q++)for(int a=0;a<L;a++)w_old[k][q][a]=w[k][q][a];  
}

void check_for_convergence(vector<Graph> A,int & it,double & l2,double tolerance,int & coincide, int decision, bool & convergence,vector< vector<double> > u,vector< vector<double> > v, vector< vector< vector<double> > > w){
    // Check new vs old every 10 iterations
    if(it % 10 ==0){    
    	double old_L=l2; 
    	l2=Likelihood(u,v,w,A);   // Calculate likelihood
    	if(abs(l2-old_L)<tolerance)coincide++;  // can use maxD or abs(l2-old_L) to check convergence
    	else coincide=0;
      }   // end if it multiple of 10
    if(coincide==decision)convergence=true;
    it++;
}

void check_for_convergence_undirected(vector<Graph_undirected> A,int & it,double & l2,double tolerance,int & coincide, int decision, bool & convergence,vector< vector<double> > u, vector< vector< vector<double> > > w){
    // Check new vs old every 10 iterations
    if(it % 10 ==0){    
      double old_L=l2; 
      l2=Likelihood_undirected(u,w,A);   // Calculate likelihood
      if(abs(l2-old_L)<tolerance)coincide++;  // can use maxD or abs(l2-old_L) to check convergence
      else coincide=0;
      }   // end if it multiple of 10
    if(coincide==decision)convergence=true;
    it++;
}
void update_optimal_parameters(double & maxL, double l2,vector< vector<double> > & u_f,vector< vector<double> > & v_f,vector< vector<double> > u,vector< vector<double> > v, vector< vector< vector<double> > > & w_f, vector< vector< vector<double> > > w){
  	
  	if(maxL<l2){
        maxL=l2;
        u_f=u;v_f=v;w_f=w;
  	}
}

void update_optimal_parameters_undirected(double & maxL, double l2,vector< vector<double> > & u_f,vector< vector<double> > u, vector< vector< vector<double> > > & w_f, vector< vector< vector<double> > > w){
    
    if(maxL<l2){
        maxL=l2;
        u_f=u;w_f=w;
    }
}

void output_results(std::vector<Graph> A,string file,string end_file,int N_real,double maxL,vector< vector<double> > u_f,vector< vector<double> > v_f,vector< vector< vector<double> > > w_f){
	
	int L=(int)A.size();
	int K=(int) w_f.size();
	int N=(int)num_vertices(A[0]);

	std::ofstream out1,out2,out3;
	out1.open((file+"u_K"+to_string(K)+end_file).c_str());
	out2.open((file+"v_K"+to_string(K)+end_file).c_str());
	out3.open((file+"w_K"+to_string(K)+end_file).c_str());
	out1<<"# Max likelihood= "<<maxL<<" N_real="<<N_real<<endl;
	out2<<"# Max likelihood= "<<maxL<<" N_real="<<N_real<<endl;
	out3<<"# Max likelihood= "<<maxL<<" N_real="<<N_real<<endl;
	
	// Output membership vectors u and v ----------------
	
	for(int i=0;i<N;i++){   // output node following their int index order
	out1<<A[0][i].name<<" ";  // output node label, a string
	out2<<A[0][i].name<<" ";
	for(int k=0;k<K;k++){
	  out1<<u_f[i][k]<<" ";
	  out2<<v_f[i][k]<<" ";
	} 
	out1<<endl;out2<<endl;
	} 
	// Output affinity matrix W ----------------------------
	cout<<"Final affinity matrix:"<<endl;
	for(int a=0;a<L;a++){
	out3<<"a= "<<a<<endl;
	for(int k=0;k<K;k++){
	  for(int q=0;q<K;q++){out3<<w_f[k][q][a]<<" ";cout<<w_f[k][q][a]<<" ";}
	  out3<<endl; cout<<endl; 
	}
	out3<<endl; cout<<endl;
	}

	cout<< "Data saved in "<<file+"u_K"+to_string(K)+end_file;
	cout<< "  "<<file+"v_K"+to_string(K)+end_file;
	cout<< "  "<<file+"w_K"+to_string(K)+end_file<<endl;

	out1.close();out2.close();out3.close();

}

void output_results_undirected(std::vector<Graph_undirected> A,string file,string end_file,int N_real,double maxL,vector< vector<double> > u_f,vector< vector< vector<double> > > w_f){
  
  int L=(int)A.size();
  int K=(int) w_f.size();
  int N=(int)num_vertices(A[0]);

  std::ofstream out1,out2,out3;
  out1.open((file+"u_K"+to_string(K)+end_file).c_str());
  out3.open((file+"w_K"+to_string(K)+end_file).c_str());
  out1<<"# Max likelihood= "<<maxL<<" N_real="<<N_real<<endl;
  out3<<"# Max likelihood= "<<maxL<<" N_real="<<N_real<<endl;
  
  // Output membership vectors u and v ----------------
  
  for(int i=0;i<N;i++){   // output node following their int index order
  out1<<A[0][i].name<<" ";  // output node label, a string
  for(int k=0;k<K;k++){
    out1<<u_f[i][k]<<" ";
  } 
  out1<<endl;
  } 
  // Output affinity matrix W ----------------------------
  cout<<"Final affinity matrix:"<<endl;
  for(int a=0;a<L;a++){
  out3<<"a= "<<a<<endl;
  for(int k=0;k<K;k++){
    for(int q=0;q<K;q++){out3<<w_f[k][q][a]<<" ";cout<<w_f[k][q][a]<<" ";}
    out3<<endl; cout<<endl; 
  }
  out3<<endl; cout<<endl;
  }

  cout<< "Data saved in "<<file+"u_K"+to_string(K)+end_file;
  cout<< "  "<<file+"w_K"+to_string(K)+end_file<<endl;

  out1.close();out3.close();

}

void print_graph_stat(std::vector<Graph> A){

	int L=(int) A.size();  // number of layers
	vector<int> E(L,0);  // number of edges in each layer
	int N=(int)num_vertices(A[0]);
	if(N<2){
  		cout<<" N="<<N<<" Too few nodes to build a network!"<<endl;
  		exit(-2);
  	}
	
	cout<<"N= "<<N<<endl;
	for(int a=0;a<L;a++){
	    E[a]=(int)num_edges(A[a]);
	    cout<<"E["<<a<<"] = "<<E[a]<<"  density= "<<100.*((double)E[a]/((double)(N*(N-1.))))<< endl;
    	if(E[a]<2){
  		cout<<"E["<<a<<"] = "<<E[a]<<" Too few edges to build a network!"<<endl;
  		exit(-2);
	  	}
  	}
  
}

void print_graph_stat(std::vector<Graph_undirected> A){

  int L=(int) A.size();  // number of layers
  vector<int> E(L,0);  // number of edges in each layer
  int N=(int)num_vertices(A[0]);
  if(N<2){
      cout<<" N="<<N<<" Too few nodes to build a network!"<<endl;
      exit(-2);
    }
  
  cout<<"N= "<<N<<endl;
  for(int a=0;a<L;a++){
      E[a]=(int)num_edges(A[a]);
      cout<<"E["<<a<<"] = "<<E[a]<<"  density= "<<100.*((double)E[a]/((double)(N*(N-1.))))<< endl;
      if(E[a]<2){
      cout<<"E["<<a<<"] = "<<E[a]<<" Too few edges to build a network!"<<endl;
      exit(-2);
      }
    }
  
}
