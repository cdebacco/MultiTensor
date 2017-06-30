
#include "mlg.hpp"

#include <boost/random.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <math.h> 
#include <stdio.h> 
#include <fstream>

using namespace boost;
using namespace std;

extern boost::mt19937 gen ;   // random real and integer number generators
uniform_real<double> const uni_dist_as(0,1);    // uniform distribution for real numbers
variate_generator<boost::mt19937 &, uniform_real<> > real01_as(gen,uni_dist_as);   // picks a random real number in (0,1)

extern double err;  // defined in main.hpp

// Functions contained in tools.cpp
void randomize_v_in_i(vector< vector<double> > & v, vector<int> v_list);
vector<int> collect_zero_entries_v(  vector<Graph>  A);
vector<int> collect_zero_entries_u(  vector<Graph>  A);
vector<int> collect_zero_entries_u(  vector<Graph_undirected>  A);
void initialize_v(vector< vector<double> > & v, istream & file);

// Assign a random real number in (0,1) to each entry
void randomize_w(vector< vector<double> > & w) {

	int L=(int)w[0].size();
	int K=(int)w.size();	
	for(int i=0;i<L;i++)for(int k=0; k<K; k++ )w[k][i]=real01_as();
}

/* Initialize the u,v, or w from input file ------------------------------------------------ */

// Input file is from a symmetric (assortative) matrix W
void initialize_w(vector< vector<double> > & w, istream & file) {
	
	/*  INPUT Format:
		l_index w_l1 w_l2 ... w_lK

		First line has the Likelihood, donc skip it.
	 */

    int a=-1;
    int L=(int)w[0].size();
    int K=(int)w.size();
    for(string line; getline( file, line ); ){
      std::istringstream is( line );
      if(a>=0){
      double l;
      is>>l; // first number is the layer index
      for(int k=0;k<K;k++)is>>w[k][a];
      }
      a++;
    }
    for(int a=0;a<L;a++)for(int k=0;k<K;k++)w[k][a]+=w[k][a]*err*real01_as();  // add some noise so to have non-zero entries out of the diagonal
}

void initialize(int initialization,  vector<int> u_list,  vector<int> v_list,vector< vector<double> > & u,vector< vector<double> > & v, vector< vector<double> > & w, string file, string w_file,vector<Graph> A){
	
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
      cout<< "W is initializaed using: "<< file+"w_K"+to_string(K)+w_file<<endl; 
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

void initialize_undirected(int initialization,  vector<int> u_list,vector< vector<double> > & u, vector< vector<double> > & w, string file, string w_file,vector<Graph_undirected> A){
  
  int K=(int)w.size();  // number of communities

  if(initialization==0){  // everything is random
        cout<<" Random initializations"<<endl;
        randomize_w(w);
        randomize_v_in_i(u,u_list);  
    }
    else if(initialization==1){ // everything is initialized
      cout<< " W and U are initialized using: "<< file+"u_K"+to_string(K)+w_file<<endl;
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
      cout<< "W is initializaed using: "<< file+"w_K"+to_string(K)+w_file<<endl; 
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

void output_affinity_matrix( vector< vector<double> > w)
{	/* OUTPUT Format:
		l w_l1 w_l2 ... w_lK
		a w_a1 w_a2 ... w_aK
	 */
	int L=(int)w[0].size();
	int K=(int)w.size();
	for(int l=0;l<L;l++){
		cout<<l<<" ";
		for(int k=0;k<K;k++)cout<<w[k][l]<<" ";
		cout<<endl;	
	}
}
//  -------   Calculate Likelihood  ------------------------------------------------------------------------------------
double Likelihood(vector< vector<double> > & u, vector< vector<double> > & v,vector< vector<double> >  & w, vector<Graph> & A){

  double l=0.;
  int L=(int)A.size();
  int K=(int)w.size();
  int N=(int)num_vertices(A[0]);
  // Calculate argument inside log
  for(int a=0;a<L;a++){
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
		    double log_arg=0.;
        for(int k=0;k<K;k++){
	       double uvw=u[i][k]*v[j][k]*w[k][a];
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

double Likelihood_undirected(vector< vector<double> > & u,vector< vector<double> >  & w, vector<Graph_undirected> & A){

  double l=0.;
  int L=(int)A.size();
  int K=(int)w.size();
  int N=(int)num_vertices(A[0]);
  // Calculate argument inside log
  for(int a=0;a<L;a++){
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
    double log_arg=0.;
        for(int k=0;k<K;k++){
      double uvw=u[i][k]*u[j][k]*w[k][a];
      l-=uvw;   // Add this term regardeles of the value of A_ijk
      if(edge(i, j, A[a]).second)log_arg+=uvw;  // if edge exists, consider this term inside the log argument
    }// end cycle over k and q
    if(log_arg>0.)l+=log(log_arg);
      }// end cycle over j      
    }// end cycle over i
  }// end cycle over a
  if(isnan(l)){
    cout<<"Likelihood is Nan!!!"<<endl;
    exit(-1);
  }
  else return l;
}

void update_old_variables(  vector< vector<double> > & u_old,  vector< vector<double> > & v_old,  vector< vector<double> > u,  vector< vector<double> > v, vector< vector<double> > & w_old, vector< vector<double> > w, vector<int> u_list,vector<int> v_list){
    int K=(int) w.size();
    int L=(int) w[0].size();
    int v_length=(int)v_list.size();
    int u_length=(int)u_list.size();
    for(int i=0;i<v_length;i++)for(int k=0;k<K;k++) v_old[v_list[i]][k]=v[v_list[i]][k]; 
    for(int i=0;i<u_length;i++)for(int k=0;k<K;k++) u_old[u_list[i]][k]=u[u_list[i]][k]; 
    for(int k=0;k<K;k++)for(int a=0;a<L;a++)w_old[k][a]=w[k][a];  
}

void update_old_variables_undirected(  vector< vector<double> > & u_old,  vector< vector<double> > u, vector< vector<double> > & w_old, vector< vector<double> > w, vector<int> u_list){
    int K=(int) w.size();
    int L=(int) w[0].size();
    int u_length=(int)u_list.size();
    for(int i=0;i<u_length;i++)for(int k=0;k<K;k++) u_old[u_list[i]][k]=u[u_list[i]][k]; 
    for(int k=0;k<K;k++)for(int a=0;a<L;a++)w_old[k][a]=w[k][a];  
}

void check_for_convergence(vector<Graph> A,int & it,double & l2,double tolerance,int & coincide, int decision, bool & convergence,vector< vector<double> > u,vector< vector<double> > v, vector< vector<double> > w){
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

void check_for_convergence_undirected(vector<Graph_undirected> A,int & it,double & l2,double tolerance,int & coincide, int decision, bool & convergence,vector< vector<double> > u, vector< vector<double> > w){
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

void update_optimal_parameters(double & maxL, double l2,vector< vector<double> > & u_f,vector< vector<double> > & v_f,vector< vector<double> > u,vector< vector<double> > v, vector< vector<double> >  & w_f, vector< vector<double> >  w){
  	
  	if(maxL<l2){
        maxL=l2;
        u_f=u;v_f=v;w_f=w;
  	}
}

void update_optimal_parameters_undirected(double & maxL, double l2,vector< vector<double> > & u_f,vector< vector<double> > u, vector< vector<double> >  & w_f, vector< vector<double> >  w){
    
    if(maxL<l2){
        maxL=l2;
        u_f=u;w_f=w;
    }
}

void output_results(std::vector<Graph> A,string file,string end_file,int N_real,double maxL,vector< vector<double> > u_f,vector< vector<double> > v_f,vector< vector<double> >  w_f){
	
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
		out3<<a<<" ";cout<<a<<" ";
		for(int k=0;k<K;k++){out3<<w_f[k][a]<<" ";cout<<w_f[k][a]<<" ";}
		out3<<endl; cout<<endl; 
	} // end cycle over layers a
	cout<< "Data saved in "<<file+"u_K"+to_string(K)+end_file;
	cout<< "  "<<file+"v_K"+to_string(K)+end_file;
	cout<< "  "<<file+"w_K"+to_string(K)+end_file<<endl;

	out1.close();out2.close();out3.close();

}

void output_results_undirected(std::vector<Graph_undirected> A,string file,string end_file,int N_real,double maxL,vector< vector<double> > u_f,vector< vector<double> >  w_f){
  
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
    out3<<a<<" ";cout<<a<<" ";
    for(int k=0;k<K;k++){out3<<w_f[k][a]<<" ";cout<<w_f[k][a]<<" ";}
    out3<<endl; cout<<endl; 
  } // end cycle over layers a
  cout<< "Data saved in "<<file+"u_K"+to_string(K)+end_file;
  cout<< "  "<<file+"v_K"+to_string(K)+end_file;
  cout<< "  "<<file+"w_K"+to_string(K)+end_file<<endl;

  out1.close();out3.close();

}


