
#include <boost/config.hpp>


#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/breadth_first_search.hpp>
//#include <boost/graph/visitors.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/random.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/property_map/property_map.hpp>
//#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
//#include <boost/graph/connected_components.hpp>

#include <fstream>
#include <vector>
#include <list>
#include <utility>
#include <string>
#include <math.h>
#include <iomanip>
#include <boost/limits.hpp>
#include <queue>
#include <time.h>
#include <ctime> 


using namespace boost;
using namespace std;

/* ----------------GRAPH VARIABLES ----------------*/
int N=500;   // Number of nodes
int L=12; 	// Number of types of edges/layers
int K=4;	// Number of groups/communities
int N_real=1;

/* ----------------CONVERGENCE VARIABLES ----------------*/
double tolerance =0.00001;
int decision = 15; 
int maxit=10;

unsigned int seed=0;
unsigned int rseed=0;

double inf=1e10; 
double err_max=0.00000000001; // max err in the converge() routine


/* ----------------I/O FILES ----------------*/
std::ofstream out1,out2,out3,out_d;
std::ifstream in1, in2;
std::string folder="Elly";
std::string end_file="";
std::string adj="adjacency1.dat";

//----------------RANDOM NUMBER GENERATORS----------------------------------------------------

boost::mt19937 gen,gen_int ;

uniform_real<double> const uni_dist(0,1);
variate_generator<boost::mt19937 &, uniform_real<> > real01(gen, uni_dist);

// Pick an int random number btw [min,max] 
int roll_die(int min, int max) {
  boost::uniform_int<> dist(min, max );
  boost::variate_generator<boost::mt19937&, boost::uniform_int<> > die(gen_int, dist);
  return die();
  //return (int) rand()%(max-min)+min;
}

//------------------------------------------------------------


/* ---------------- GRAPH TOOLS -----------------------------*/



void randomize_v(vector<double> & v) {
  double Z=0.;
  // Assign a random weight to each entry
  for (int i=0; i<v.size() ; i++ ) {
    //m.E[d] =roll_die(0,1);
    v[i] =real01();

    Z+=v[i];
  }
  // Normalize such that the sum over the entries is 1
  for (int i=0; i<v.size() ; i++) v[i]/=Z;
}

void randomize_w(vector< vector<double> > & w) {

   for(int i=0;i<L;i++)for (int k=0; k<K; k++ ){
      w[k][i] =real01();
    }
}

void randomize_v_in_i(vector< vector<double> > & v, vector<int> v_list) {
  // Assign a random weight to each entry
  int list_length=v_list.size()	;
    for (int k=0; k<K; k++ ){
    double Z=0.; 
//    double r=real01();
    for(int i=0;i<list_length;i++){ 
    	int j=v_list[i];
      v[j][k] =real01();
      //v[i][k] =(double)roll_die(1,4);
      Z+=v[j][k];
    }
     for(int i=0;i<list_length;i++)v[v_list[i]][k]/=Z; // Normalize such that the sum over the entries is 1
  }
}



//-------------------------------------------------------------------------------

/*
struct EdgeProperties {
  EdgeProperties() :  {int i=1;}
  };
*/
//--------------------------------------------------------------------

struct VertexProperties  {
  VertexProperties() : name("0"){}
  VertexProperties(string const & name) : name(name){}
  string name;
};
//--------------------------------------------------------------------
typedef adjacency_list<setS, vecS,/*directedS,undirectedS bidirectionalS*/bidirectionalS,VertexProperties
		       /*VertexProperties, EdgeProperties*/> GraphBase;

typedef graph_traits<GraphBase>::vertex_iterator vertex_iterator;
typedef graph_traits<GraphBase>::out_edge_iterator edge_iterator;
typedef graph_traits<GraphBase>::in_edge_iterator in_edge_iterator;
typedef graph_traits<GraphBase>::edge_iterator graph_edge_iterator;
typedef graph_traits<GraphBase>::edge_descriptor Edge;
typedef graph_traits<GraphBase>::vertex_descriptor Vertex;

struct Graph : public GraphBase {
  Vertex rootid;
};

vector<Graph> A(L);  // a vector of graph, one for each of the L layers
//Graph g,g_f;


// Graph related functions
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

//-------------------------------------------------------------------------------



void read_graph(istream & file, vector<Graph> & A)
{  string tok;
  while (file >> tok){
    if (tok == "E") {
      string i, j;
      file >> i >> j ;
      int v1=idx(i,A);  // Adds node i to all layers
      int v2=idx(j,A);  // Adds node j to all layers
      int is_edge;
      for(int alpha=0;alpha<L;alpha++){
        file>>is_edge;
        if(is_edge==1)Edge e = add_edge(vertex(v1, A[alpha]), vertex(v2, A[alpha]), A[alpha]).first;
      }
    }
  }
 } 

//-------------------------------------------------------------------------------
void out_graph(string folder, vector<Graph> A)
{
  for(int a=0;a<L;a++){
    std::ofstream out;
    out.open((folder+"adjacency1_"+to_string(a)+".dat").c_str());
    Graph gr=A[a];
    graph_edge_iterator eit, eend;
    for (tie(eit, eend) = edges(gr); eit != eend; ++eit) {
      Edge e=*eit;
      Vertex j = target(e, gr); 
      Vertex i = source(e, gr);  
      out << gr[i].name << " " << gr[j].name << endl;
     // out<<"{ ";  
     // out << gr[i].name << " -> " << gr[j].name <<" , "<<a<<" } "<<","<< endl;
    } 
    out.close();
  }
}



int in_neigh(Vertex i, vector< Graph> A){
	int L=A.size();
	int k=0;
	in_edge_iterator eit, eend;   // Cycle over out-neighbors of i in layer a
	for(int a=0;a<L;a++){
       for(tie(eit, eend) = in_edges(i,A[a]); eit != eend; ++eit) k++;
	}
return k;
}

int out_neigh(Vertex i, vector< Graph> A){
	int L=A.size();
	int k=0;
	edge_iterator eit, eend;   // Cycle over out-neighbors of i in layer a
	for(int a=0;a<L;a++){
       for(tie(eit, eend) = out_edges(i,A[a]); eit != eend; ++eit) k++;
	}
return k;
}

void set_zero_entries(vector< vector<double> > & u, vector< vector<double> > & v,  vector<Graph>  A){

	for(int i=0;i<N;i++){
		if(out_neigh(i,A)==0)for(int k=0;k<K;k++)u[i][k]=0.;
		if(in_neigh(i,A)==0)for(int k=0;k<K;k++)v[i][k]=0.;
	}
}


// Normalize k-column
void normalize_w_over_a(vector<vector<double> > & w,int k){
  double Z=0.;
  for(int a=0;a<L;a++)Z+=w[k][a];
  for(int a=0;a<L;a++){
    if(Z>0)w[k][a]/=Z;
    if(w[k][a]<err_max)w[k][a]=0.;
  }
}
void normalize_w_over_k(vector<vector<double> > & w){
  for(int a=0;a<L;a++){
    double Z=0.;
    for(int k=0;k<K;k++)Z+=w[k][a];
    for(int k=0;k<K;k++){
      if(Z>0)w[k][a]/=Z;
      if(w[k][a]<err_max)w[k][a]=0.;
    }
  }  
}
// Normalize k-column
void normalize_over_i(vector<vector<double> > & u,int k){
  double Z=0.;
  for(int i=0;i<N;i++)Z+=u[i][k];
  for(int i=0;i<N;i++){
    if(Z>0)u[i][k]/=Z;
    if(u[i][k]<err_max)u[i][k]=0.;
  }
}
//Normalize in k
void normalize_over_k(vector<vector<double> > & u){
  for(int i=0;i<u.size();i++){
    double Zu=0.;
    for(int k=0;k<u[i].size();k++)Zu+=u[i][k];
    for(int k=0;k<u[i].size();k++){if(Zu>0)u[i][k]/=Zu;if(u[i][k]<err_max)u[i][k]=0.;}
  }
}

vector<int> remove_zero_entries_v(  vector<Graph>  A){

  vector<int> v;
  for(int i=0;i<N;i++)if(in_neigh(i,A)>0){v.push_back(i);}
  return v;
}


vector<int> remove_zero_entries_u(  vector<Graph>  A){

  vector<int> u;
  for(int i=0;i<N;i++)if(out_neigh(i,A)>0)u.push_back(i);else cout<<A[0][i].name<<endl;
  return u; 
}

// ------  Update in sequence  -------------
double update_sequential(vector< vector<double> > & u,vector< vector<double> > & u_old,vector< vector<double> > & v,vector< vector<double> > & v_old,vector< vector<double> > & w,vector< vector<double> > & w_old, vector<Graph> A, vector<int> u_list,vector<int> v_list,double & d_u,double & d_v,double & d_w){
 
 double deltaL=0.;  // Delta Likelihood
 double dist_u,dist_v,dist_w;
  /* ---------------------------------------------------------------------------
  --------------- 1.  Cycle over u_ik -----------------------------------------
  ---------------------------------------------------------------------------*/
  for(int k=0;k<K;k++){

    // Calculate Du and w_k, they are the same for all u_i
    double w_k=0.,Du=0.;
    for(int a=0;a<L;a++)w_k+=w_old[k][a];
    //for(int i=0;i<N;i++)Du+=v_old[i][k];
    for(int i=0;i<v_list.size();i++)Du+=v_old[v_list[i]][k];
     // cout<<"Du "<<Du<<" "<<" w_k "<<w_k<<" "<<k<<endl;
    
    // Cycle over i ----
    for(int z=0;z<u_list.size();z++)if(u_old[u_list[z]][k]>0.){     // Update only if u_ik >0 	
   // for(int i=0;i<N;i++)if(u_old[i][k]>0.){     // Update only if u_ik >0 
      int i=u_list[z];	
      u[i][k]=0.;
      double deltaU=0.;

      for(int a=0;a<L;a++){ // Cycle over layers a -----------
        edge_iterator eit, eend;   // Cycle over out-neighbors of i in layer a
        for (tie(eit, eend) = out_edges(i,A[a]); eit != eend; ++eit) {
          Vertex j=target(*eit,A[a]);   // OUT-EDGE

          double Zij_a=0.;
          for(int q=0;q<K;q++)Zij_a+=u_old[i][q]*v_old[j][q]*w_old[q][a];            
          if(Zij_a==0)out_d<<"Zij_a "<<A[a][i].name<<"--> "<<A[a][j].name<<" "<<Zij_a<<" "<<a<<endl;

          double rho_ijka=v_old[j][k]*w_old[k][a];   //  u_old[i][k] will be multiplied only at the end!!!!
          if(Zij_a>0)rho_ijka/=Zij_a;
          else rho_ijka=0.;
          u[i][k]+=rho_ijka;

        } //end cycle over out-j  in layer a
      } // end cycle over layer a - --- --- --- - - - - - - -- - - - - - - - - 

      deltaU=Du*w_k - u[i][k];///u_old[i][k];
      if((Du!=0.) and (w_k!=0.))u[i][k]*= u_old[i][k]/(Du*w_k);   // final update for u_ik
      else u[i][k]=0.;

     /* deltaU=1. - u[i][k];
      u[i][k]*= u_old[i][k];*/
      deltaU*=(u[i][k]-u_old[i][k]);
      deltaL+=deltaU;
      //cout<<" delta u "<<A[0][i].name<<" "<<deltaU<<" "<<u[i][k]<<" "<<u_old[i][k]<<" "<<deltaU/(u[i][k]-u_old[i][k])<<endl;
      if(u[i][k]<err_max){u[i][k]=0.;  /*cout<<" u=0 "<<A[0][i].name<<endl;*/}
      dist_u=std::max(abs(u[i][k]-u_old[i][k]),dist_u);
      u_old[i][k]= u[i][k];
 // < u_old[i][k]=u[i][k];  // update so that it will appear correctly in the normalization Zij_a

    } // cycle over i --------------------------------------------------
    //  normalize_over_k(u);
 // normalize_over_i(u,k);   // Normalize the k-th column u_k

  }  // end cycle over k --------------------------------------------------
  d_u=dist_u;
   // normalize_over_k(u);
   //normalize_over_k(v);
  //for(int k=0;k<K;k++)normalize_over_i(u,k); 

    /*--------------- 1.  end cycle  over u_ik --------------------------*/
 /* ---------------------------------------------------------------------------
  --------------- 2.  Cycle  over v_ik -----------------------------------------
  ---------------------------------------------------------------------------*/
  for(int k=0;k<K;k++){

    // Calculate Dv and w_k, they are the same for all v_i
 
    double w_k=0.,Dv=0.;
    for(int a=0;a<L;a++)w_k+=w_old[k][a];
    //for(int i=0;i<N;i++)Dv+=u[i][k];   // Use the new u_ik !!!!!!!
  for(int i=0;i<u_list.size();i++)Dv+=u[u_list[i]][k];   // Use the new u_ik !!!!!!!
  
   //       cout<<"Dv "<<Dv<<" "<<" w_k "<<w_k<<" "<<k<<endl;
    // Cycle over i ----
    //for(int i=0;i<N;i++)if(v_old[i][k]>0.){     // Update only if v_ik >0 
    for(int z=0;z<v_list.size();z++)if(v_old[v_list[z]][k]>0.){     // Update only if v_ik >0 
    
      int i=v_list[z];  
      v[i][k]=0.;
      double deltaV=0.;

      for(int a=0;a<L;a++){ // Cycle over layers a -----------
        in_edge_iterator eit, eend;   // Cycle over out-neighbors of i in layer a
        //edge_iterator eit, eend; 
        for (tie(eit, eend) = in_edges(i,A[a]); eit != eend; ++eit) {

          Vertex j=source(*eit,A[a]);  // IN-EDGE

          double Zij_a=0.;
          for(int q=0;q<K;q++)Zij_a+=u_old[j][q]*v_old[i][q]*w_old[q][a];            
          if(Zij_a==0)out_d<<"Zij_a "<<A[a][j].name<<"--> "<<A[a][i].name<<" "<<Zij_a<<" "<<a<<endl;

          double rho_jika=u_old[j][k]*w_old[k][a];   //  v_old[i][k] will be multiplied only at the end!!!!
          if(Zij_a>0)rho_jika/=Zij_a;
          else rho_jika=0.;
          v[i][k]+=rho_jika;
    
        } //end cycle over in-j  in layer a
      } // end cycle over layer a - --- --- --- - - - - - - -- - - - - - - - - 

      deltaV=Dv*w_k - v[i][k];///v_old[i][k];
      if((Dv!=0.) and (w_k!=0.))v[i][k]*= v_old[i][k]/(Dv*w_k);   // final update for u_ik
      else v[i][k]=0.;
     // cout<<A[0][i].name<<" "<<v[i][k]<<endl;
      
     /* deltaV=1. - v[i][k];
      v[i][k]*= v_old[i][k];*/

     deltaV*=(v[i][k]-v_old[i][k]);
      deltaL+=deltaV;

//cout<<" delta v "<<A[0][i].name<<" "<<deltaV<<" "<<v[i][k]<<" "<<v_old[i][k]<<" "<<deltaV/(v[i][k]-v_old[i][k])<<endl;
     if(v[i][k]<err_max){v[i][k]=0.; }//cout<<" v=0 "<<A[0][i].name<<endl;
    dist_v=std::max(abs(v[i][k]-v_old[i][k]),dist_v);
      v_old[i][k]= v[i][k];
    // v_old[i][k]=v[i][k];  // update so that it will appear correctly in the normalization Zij_a
//v[i][k]=u[i][k];
    } // cycle over i --------------------------------------------------
   //  normalize_over_k(v);

 //  normalize_over_i(v,k);   // Normalize the k-th column u_k
  

  }  // end cycle over k -------------------------------------------------- 
  d_v=dist_v;
 //normalize_over_k(v);
  //for(int k=0;k<K;k++)normalize_over_i(v,k); 

  /*--------------- 2.  end cycle  over v_ik --------------------------*/
   /* ---------------------------------------------------------------------------
  --------------- 3.  Cycle  over w_ka -----------------------------------------
  ---------------------------------------------------------------------------*/
  for(int k=0;k<K;k++){
    // Calculate Z_k, they are the same for all w_a
   double Z_k=0.;
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++) {Z_k+=u[i][k]*v[j][k];}// cout<<k<<" "<<Z_k<<" "<<u[i][k]<<" "<<v[j][k]<<endl;}
    }


    for(int a=0;a<L;a++)if(w_old[k][a]>0.){  // Update only if w_ka >0 
      w[k][a]=0.;
      double deltaW=0.;

      for(int i=0;i<N;i++){ // Cycle over i ----
        edge_iterator eit, eend;   // Cycle over out-neighbors of i in layer a
        double rho_ijka=0.;
        for (tie(eit, eend) = out_edges(i,A[a]); eit != eend; ++eit) {
          Vertex j=target(*eit,A[a]);
          double Zij_a=0.;
          for(int q=0;q<K;q++)Zij_a+=u_old[i][q]*v_old[j][q]*w_old[q][a];            
          if(Zij_a==0){out_d<<"Zij_a (in updating w_ka) "<<A[a][i].name<<"--> "<<A[a][j].name<<" "<<Zij_a<<" "<<a<<endl;}
          if(Zij_a>0)rho_ijka+=v_old[j][k]/Zij_a;   //  w_old[k][a] will be multiplied only at the end!!!!
          else rho_ijka=0.;
          //if(Zij_a>0)rho_ijka/=Zij_a;
        } //end cycle over out-j  of i
        w[k][a]+=u_old[i][k]*rho_ijka;
       // cout<<k<<" "<<a<<" "<<w[k][a]<<" "<<i<<" "<<u_old[i][k]<<" "<<rho_ijka<<endl;
      } // end cycle over i - --- --- --- - - - - - - -- - - - - - - - - 
      
      deltaW=Z_k - w[k][a];///w_old[k][a];
      if(Z_k!=0.)w[k][a]*=w_old[k][a]/Z_k;   // final update for w_ka
      else w[k][a]=0.;//*=w_old[k][a];

     /* deltaW=1. - w[k][a];
      w[k][a]*=w_old[k][a];*/
      deltaW*=(w[k][a]-w_old[k][a]);
      deltaL+=deltaW;
//cout<<" delta w "<<deltaW<<endl;
      //cout<<w[k][a]<<endl;
     if(w[k][a]<err_max)w[k][a]=0.;
      dist_w=std::max(abs(w[k][a]-w_old[k][a]),dist_w);
      w_old[k][a]= w[k][a];

    //  w_old[k][a]=w[k][a];  // update so that it will appear correctly in the normalization Zij_a
    } // cycle over a --------------------------------------------------
   // normalize_w_over_k(w);
  // normalize_w_over_a(w,k);   // Normalize the k-th column w_k
  }  // end cycle over k --------------------------------------------------

d_w=dist_w;
 // for(int k=0;k<K;k++)normalize_w_over_a(w,k); 
   //for(int a=0;a<L;a++)
 // normalize_w_over_k(w);
  return -deltaL;
} // end update
// ---------------------- end update sequentially

// ------  Update in parallel  -------------
void update_parallel(vector< vector<double> > & u,vector< vector<double> >  u_old,vector< vector<double> > & v,vector< vector<double> >  v_old,vector< vector<double> > & w,vector< vector<double> >  w_old, vector<Graph> A){

  //static std::vector< vector<int> > non_zero_u(K),non_zero_v(K);  // allows to select for update only the non zero components
  // Cycle over groups
for(int k=0;k<K;k++){
  
 // if(non_zero_u[k].size()==0)for(int i=0;i<N;i++)if(u_old[i][k]>err_max)non_zero_u[k].push_back(i);
  //if(non_zero_v[k].size()==0)for(int i=0;i<N;i++)if(v_old[i][k]>err_max)non_zero_v[k].push_back(i);

 
  // Initialize variables
  double w_k=0.,Z_k=0.,Du=0.,Dv=0.;
  for(int a=0;a<L;a++){
    w[k][a]=0.;
    w_k+=w_old[k][a];
  }

 for(int i=0;i<N;i++){
  Du+=v_old[i][k];
  Dv+=u_old[i][k];
  for(int j=0;j<N;j++)Z_k+=u_old[i][k]*v_old[j][k];
  }
  for(int i=0;i<N;i++){
    Du+=v_old[i][k];
    Dv+=u_old[i][k];
    for(int j=0;j<N;j++)Z_k+=u_old[i][k]*v_old[j][k];
  }
  
  //Cycle over nodes i
  for(int i=0;i<N;i++){
    u[i][k]=v[i][k]=0.;
    double Bijk=0.,Bjik=0.;

    //Cycle over layers
    for(int a=0;a<L;a++){
      //Cycle over out-neighbors of i in layer a
      
      // Update only if u_ik >0 
      if(u_old[i][k]>0.){
        edge_iterator eit, eend;
        for (tie(eit, eend) = out_edges(i,A[a]); eit != eend; ++eit) {
          Vertex j=target(*eit,A[a]);
          double Zij_a=0.;
          for(int q=0;q<K;q++){ 
            Zij_a+=u_old[i][q]*v_old[j][q]*w_old[q][a];            
          }
          if(Zij_a==0)out_d<<"Zij_a "<<A[a][i].name<<"--> "<<A[a][j].name<<" "<<Zij_a<<" "<<a<<endl;
          double rho_ijka=u_old[i][k]*v_old[j][k]*w_old[k][a];
          if(Zij_a>0)rho_ijka/=Zij_a;
          w[k][a]+=rho_ijka;
          Bijk+=rho_ijka;
        } //end cycle over out-j
        u[i][k]+=Bijk;
       } // end only if u_ik > 0 

       //Cycle over in-neighbors of i in layer a
      // Update only if v_jk >0 
      if(v_old[i][k]>0.){
        in_edge_iterator in_eit, in_eend;
        for (tie(in_eit, in_eend) = in_edges(i,A[a]); in_eit != in_eend; ++in_eit) {
          Vertex j=source(*in_eit,A[a]);
          double Zji_a=0.;
          for(int q=0;q<K;q++){      
            Zji_a+=u_old[j][q]*v_old[i][q]*w_old[q][a];
          }
          if(Zji_a==0)out_d<<"Zji_a "<<A[a][i].name<<"--> "<<A[a][j].name<<" "<<Zji_a<<" "<<a<<endl;
          double rho_jika=u_old[j][k]*v_old[i][k]*w_old[k][a];
          if(Zji_a>0)rho_jika/=Zji_a;
          w[k][a]+=rho_jika;
          Bjik+=rho_jika;
        } //end cycle over in-j
        v[i][k]+=Bjik;
      } // end only if v_jk > 0  

    } // end cycle over layers a
     if((Du!=0.) and (Dv!=0.) and (w_k!=0.)){
     u[i][k]/=Du*w_k; //Normalize
     v[i][k]/=Dv*w_k; //Normalize 
     }

  } // end cycle over i

  for(int a=0;a<L;a++)w[k][a]/=Z_k; //Normalize
  
} //end cycle over k groups
  
  // Normalize vectors

  
  for(int a=0;a<L;a++)
  {
    double Za=0.;
    for(int k=0;k<K;k++)Za+=w[k][a];
    if(Za>0)for(int k=0;k<K;k++)w[k][a]/=Za;
    for(int k=0;k<K;k++)if(w[k][a]<err_max)w[k][a]=0.;
  }

  for(int k=0;k<K;k++){
    double Zu=0.,Zv=0.,Zw=0.;
  
    for(int i=0;i<N;i++){
      Zu+=u[i][k];Zv+=v[i][k];
    }
    for(int a=0;a<L;a++)Zw+=w[k][a];
    if(Zw>0)for(int a=0;a<L;a++)w[k][a]/=Zw;
    for(int a=0;a<L;a++)if(w[k][a]<err_max)w[k][a]=0.;

    for(int i=0;i<N;i++){
      if(Zu>0)u[i][k]/=Zu;if(u[i][k]<err_max)u[i][k]=0.;
      if(Zv>0)v[i][k]/=Zv;if(v[i][k]<err_max)v[i][k]=0.;
    }
  }
  
   // normalize_over_k(v);normalize_over_k(u);
} // end update
// --------------------------------------------------

 //-------Calculate Likelihood ---------------------
double Likelihood(vector< vector<double> > u, vector< vector<double> > v,vector< vector<double> > w, vector<Graph>  A){

double l=0.;
  // Calculate argument inside log

  for(int a=0;a<L;a++){
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
         double log_arg=0.;double x=0.;
        for(int k=0;k<K;k++){
          double uvw=u[i][k]*v[j][k]*w[k][a];
          x+=uvw;
          l-=uvw;
          if(edge(i, j, A[a]).second)log_arg+=uvw;
        }// end cycle over k
        if(log_arg>0.)l+=log(log_arg);
        if((x==log_arg==0.) ){
        }
      }// end cycle over j
       //log_arg+=u[i];
    }// end cycle over i
      
   

  }// end cycle over a

  return l;
}
//--------------------------------------------------
// 2-norm of a K-dimensional vector
double distance(vector< vector<double> > u,vector< vector<double> > u_old){
  double maxD=0.;
  for(int k=0;k<u[0].size();k++){
    double d=0.;
    for(int i=0;i<u.size();i++)d+=(u[i][k]-u_old[i][k])*(u[i][k]-u_old[i][k]);
    maxD=std::max(maxD,d);  
  }
  return std::sqrt(maxD);
}
//--------------------------------------------------

namespace po = boost::program_options;

po::variables_map parse_command_line(int ac, char ** av)
{
  po::options_description desc("Usage: " + string(av[0]) + " <option> ... \n\twhere <option> is one or more of");
  desc.add_options()
    ("help", "produce help message")
    ("folder,f", po::value(&folder)->default_value("Elly/"), "set folder name")
    ("adj,a", po::value(&adj)->default_value("adjacency1.dat"), "set adjacency file name") 
    ("end_file,E", po::value(&end_file)->default_value(".dat"), "set end of file name")
    ("L,l", po::value(&L)->default_value(12), "set number of layers")
    ("K,k", po::value(&K)->default_value(4), "set number of groups")
    ("N_real,r", po::value(&N_real)->default_value(1), "set number of realizations")
    ("N,N", po::value(&N)->default_value(100), "set number of nodes")
    ("maxit,t", po::value(&maxit)->default_value(10), "set maximum number of iterations")
    ("tolerance,e", po::value(&tolerance)->default_value(1.), "set convergence tolerance")
    ("seed,s", po::value<unsigned>(), "sets instance seed")
    ("rseed,z", po::value<unsigned>(), "sets biases seed")
    ("decision,y", po::value<int>(&decision)->default_value(1), "program converges after this # repeats of the decision variables");

  po::variables_map vm;
  po::store(po::parse_command_line(ac, av, desc), vm);
  po::notify(vm);

  if (vm.count("seed")) {
    unsigned s = vm["seed"].as<unsigned>();
    gen_int.seed(s);
  }
  if (vm.count("rseed")) {
    unsigned s = vm["rseed"].as<unsigned>();
    gen.seed(s);
  }
  
  if (vm.count("help")) {
    cerr << desc << "\n";
    exit(1);
  }

  return vm;
}

/* ---------------- MAIN ---------------------------------------------*/
//--------------------------------------------------------------------
int main(int ac, char** av)
{
  
  cout.setf(ios_base::fixed, ios_base::floatfield);
  po::variables_map vm = parse_command_line(ac, av);
  
  // FILE
  string file="data/"+folder;

  //Variables
  //vector<Graph> A(L);
  
  clock_t t=0.;//time to process one m (N_realizations)
  
  /*-------------  FILE I/O ------------------------------*/

  //in1.open(("/Users/Caterina/Dropbox/Lavoro/Noverlap/Data/"+folder+"/k3/Single/adjacency.dat").c_str());
  in1.open(("data/"+folder+adj).c_str());
  cout<<"Infile = "<<"data/"+folder+adj<<endl;
  out_d.open((file+"debug_K"+to_string(K)+end_file).c_str()/*, ios::app*/);

  /*  gen.seed(rseed);
  mes_gen.seed(rseed);
  gen_int.seed(seed);*/

  

  read_graph(in1,A);
  in1.close();

  N=(int)num_vertices(A[0]);
  cout<<"N= "<<N<<endl;
  vector<int> E(L,0);
  for(int a=0;a<L;a++){
    E[a]=(int)num_edges(A[a]);
    cout<<"E["<<a<<"] = "<<E[a]<<"  density= "<<100.*0.5*((double)E[a]/((double)(N*(N-1.))))<< endl;
  }
  out_graph("data/"+folder+"Graph/",A);

   vector<int> u_list=remove_zero_entries_u(A);  // Contains vertex indices of vertices having  out-edges
  vector<int> v_list=remove_zero_entries_v(A); // Contains vertex indices of vertices having in-edges

  /*--------------------------------------------------*/

/* ----------------RESULTS VARIABLES ----------------*/
vector< vector<double> > u_f(N, vector<double>(K,0.)),v_f(N, vector<double>(K,0.)),w_f(K, vector<double>(L));   // Output results
double maxL=-inf;

// Cycle over different realizations ---> try to avoid local minima
out1.open((file+"Likelihood_K"+to_string(K)+end_file).c_str());

 vector< vector<double> > u(N, vector<double>(K,0.)),v(N, vector<double>(K,0.)),w(K, vector<double>(L,0.));   // Output results

 //vector< vector<double> > u(N, vector<double>(K)),v(N, vector<double>(K)),w(K, vector<double>(L));   // Output results
 vector< vector<double> > u_old(N, vector<double>(K,0.)),v_old(N, vector<double>(K,0.)),w_old(K, vector<double>(L));   // Previous time step values 


 int v_length=v_list.size();
 int u_length=u_list.size();

for(int r=0;r<N_real;r++){

 
  // Randomize initial vectors
  // for(int i=0;i<N;i++){
  //   randomize_v(u[i]);
  //   randomize_v(v[i]);
  //   if(i<K)randomize_v(w[i]);
  // }

   randomize_v_in_i(v,v_list);  
   randomize_v_in_i(u,u_list); 
   randomize_w(w);

  for(int i=0;i<v_length;i++)for(int k=0;k<K;k++) v_old[v_list[i]][k]=v[v_list[i]][k]; 
  for(int i=0;i<u_length;i++)for(int k=0;k<K;k++) u_old[u_list[i]][k]=u[u_list[i]][k]; 

  for(int k=0;k<K;k++)for(int a=0;a<L;a++)w_old[k][a]=w[k][a];  

  	
  out1<<-1<< " "<<Likelihood(u,v,w,A)<<endl;
   /* cout<< " Resetting ..."<<endl;
  set_zero_entries(u,v,A);
  cout<< " Finished!"<<endl;*/
//double l;
 // double l=Likelihood(u,v,w,A);
  //cout<<"r= "<<r<<" initial Likelihood= "<<l<<endl;
  //out_d<<"r= "<<r<<" initial Likelihood= "<<l<<endl;

  // Convergence variables
  int coincide=0;  // Number of consecutive times error is lower than the threshold
  bool convergence=false; 
  int it=0; // Iteration time

  double l2=-inf;
  double delta_u,delta_v,delta_w;
   //for(int i=0;i<N;i++) u[i]=v[i];
/* ------------------- Single step iteration update ------------------*/
  while(!convergence and it<maxit){

  
  // Re-set variables

	/*
  for(int i=0;i<N;i++){
  //u_old[i]=u[i];   
  //v_old[i]=v[i];  
  }  */

  
  //	cout<<" Updating ...."<<endl;
  double delta_l=	update_sequential(u,u_old,v,v_old,w,w_old,A,u_list,v_list,delta_u,delta_v,delta_w);
 // cout<<"End update!"<<endl;
  // update_parallel(u,u_old,v,v_old,w,w_old,A);
  
  if(it % 10 ==0){
  double old_L=l2; 
  double maxD=0.;
  /*
  double delta_u=distance(u,u_old);//cout<<" delta_u= "<<delta_u<<endl;
  double delta_v=distance(v,v_old);//cout<<" delta_v= "<<delta_v<<endl;
  double delta_w=distance_w(w,w_old);//cout<<" delta_w= "<<delta_w<<endl;
  */
  maxD=std::max(maxD,delta_u);maxD=std::max(maxD,delta_v);std::max(maxD,delta_w);

  //l=Likelihood(u,v,w,A);
 // if(it==0)l2=Likelihood(u,v,w,A);
 // else l2+=delta_l; 

  l2=Likelihood(u,v,w,A);
  

  out_d<<it<<" "<<maxD<<" "<<l2<<endl;
  
  
  if(/*maxD*/abs(l2-old_L)<tolerance)coincide++;
  else coincide=0;
  }
  out1<<it<< " "<<l2<<" "<<delta_u<<" "<<delta_v<<" "<<delta_w<< endl;


  if(coincide==decision)convergence=true;

  it++;
  }// end while convergence

  l2=Likelihood(u,v,w,A);
  cout<<"r= "<<r<<" Final Likelihood= "<<l2<<endl;
  out_d<<"r= "<<r<<" Final Likelihood= "<<l2<<endl;
  out_d<<"----------------------"<<endl;

  // Keep track of the realization that achieves the max likelihood
  

  if(maxL<l2){
    maxL=l2;
    u_f=u;v_f=v;w_f=w;
  }

 } // end cycle over realizations

 // Normalize vectors over k
 /*
 normalize_over_k(u_f);
 normalize_over_k(v_f);
 normalize_w_over_k(w_f);
 */
 maxL=Likelihood(u_f,v_f,w_f,A); 
 cout<<" Final Likelihood= "<<maxL<<endl;
 out1.close();
 // Output results
 out1.open((file+"u_K"+to_string(K)+end_file).c_str());
 out2.open((file+"v_K"+to_string(K)+end_file).c_str());
 out3.open((file+"w_K"+to_string(K)+end_file).c_str());
 out1<<"# Max likelihood= "<<maxL<<endl;
 out2<<"# Max likelihood= "<<maxL<<endl;
 out3<<"# Max likelihood= "<<maxL<<endl;
 cout<< " Data saved in "<<file+"u_K"+to_string(K)+end_file;
 cout<< "  "<<file+"v_K"+to_string(K)+end_file;
 cout<< "  "<<file+"w_K"+to_string(K)+end_file<<endl;
 for(int i=0;i<N;i++){
  out1<<A[0][i].name<<" ";
  out2<<A[0][i].name<<" ";
  for(int k=0;k<K;k++){
  out1<<u_f[i][k]<<" ";
  out2<<v_f[i][k]<<" ";
  } 
  out1<<endl;out2<<endl;
 } 

 for(int a=0;a<L;a++){
  out3<<a<<" ";
  for(int k=0;k<K;k++){
   out3<<w_f[k][a]<<" ";
  }
  out3<<endl;
 }

 out1.close();out2.close();out3.close();
 out_d.close();

  cout<<"N= "<<N<<endl;

  for(int a=0;a<L;a++){
    E[a]=(int)num_edges(A[a]);
    cout<<"E["<<a<<"] = "<<E[a]<<"  density= "<<100.*0.5*((double)E[a]/((double)(N*(N-1.))))<< endl;
  }

return 1;
}


