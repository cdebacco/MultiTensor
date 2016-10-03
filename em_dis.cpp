
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
#include <sstream>
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
int N_real=1;   // Number of different initial realizations --> to tackle issue of getting stuck in local minima

/* ----------------CONVERGENCE VARIABLES ----------------*/
double tolerance =0.00001;    // difference btw old and new values of the convergence variables
int decision = 15;            // number of times the convergence criterium is satisfied before stopping the simulation
int maxit=10;                
int teq=500;                   // equilibration time: if Likelihood is lower than the previous realizations' ones within teq number of realization => convergence is satisfied

unsigned int seed=0;   // seed for gen_int --> used only if you pick a random int using roll_die()
unsigned int rseed=0;  // seed for gen ---> used to pick a random real number using real01(). E.g. to initialize the variables u,v and w

double inf=1e10; 
double err_max=0.000001; // max err in the converge() routine
double err=0.01;              // magnitude of the noise added for the initialization form file input
int initialization=0;    // 0-> random initialization; 1-> u,v,w are initializes; 2 --> only w is initialized; 3 --> only u and v are initialized 
/* ----------------I/O FILES ----------------*/
std::ofstream out1,out2,out3,out_d;
std::ifstream in1, in2, in_w;
std::string folder="Elly";   // input and output folder
std::string end_file="";    // output end of file
std::string adj="adjacency1.dat";   // input adjacency matrix file
std::string w_file="w.dat";   // input file for the W matrix
//----------------RANDOM NUMBER GENERATORS----------------------------------------------------

boost::mt19937 gen,gen_int ;   // random real and integer number generators

uniform_real<double> const uni_dist(0,1);    // uniform distribution for real numbers
variate_generator<boost::mt19937 &, uniform_real<> > real01(gen, uni_dist);   // picks a random real number in (0,1)

// Pick an int random number btw [min,max] 
int roll_die(int min, int max) {
  boost::uniform_int<> dist(min, max );
  boost::variate_generator<boost::mt19937&, boost::uniform_int<> > die(gen_int, dist);
  return die();
}

//------------------------------------------------------------


/* ---------------- GRAPH TOOLS -----------------------------*/

// Assign a random real number in (0,1) to each entry

void randomize_w(vector<vector< vector<double> > > & w) {
  for(int i=0;i<L;i++)for(int k=0; k<K; k++ )for(int q=k; q<K; q++ ){
    if(q==k)w[k][q][i]=real01();  // for better convergence
    else w[k][q][i]=w[q][k][i]=err*real01();
  }
}


void randomize_v_in_i(vector< vector<double> > & v, vector<int> v_list) {
  int list_length=v_list.size()	;
    for (int k=0; k<K; k++ )for(int i=0;i<list_length;i++){ 
    	int j=v_list[i];
	v[j][k] =real01();
    }
}

/* Initialize the u,v, or w from input file ------------------------------------------------ */

// Input file is from a symmetric (assortative) matrix W
void initialize_w(vector<vector< vector<double> > > & w, istream & file) {
    int a=-1;
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

/* ------------------------------------------------------------------------------------------*/

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

vector<Graph> A(L);  // a vector of graph, one for each of the L layers; Global variable

/* ---------- Graph related functions ----------------------------------------*/

// add vertex with label 'id'; only if the vertex was not present in the graph yet
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

//-------------------------------------------------------------------------------

/* INPUT multilayer adjacency file format:
 E 102 231 1 0 0 1     
 "102" and "231" are the source and target node string labels; 
" 1 0 0 1"  stand for the direct edge  102 -> 231  exists in layer 1 and 4, not in layer 2 and 3; 
*/

void read_graph(istream & file, vector<Graph> & A)
{  string tok;
  while (file >> tok){
    if (tok == "E") {  // control variable "E"; to check that the input file containes the edge information
      string i, j;  // source and target node labels
      file >> i >> j ;
      int v1=idx(i,A);  // Adds node i to all layers and returns corresponding int index
      int v2=idx(j,A);  // Adds node j to all layers and returns corresponding int index
      int is_edge;
      for(int alpha=0;alpha<L;alpha++){
        file>>is_edge;  // if this is 1 then created an edge in the correspending layer
        if(is_edge==1)Edge e = add_edge(vertex(v1, A[alpha]), vertex(v2, A[alpha]), A[alpha]).first;
      }  // end cycle for over layers
    }  // end if(tok=="E")
  }  // end while
 } 

// Outputs single layers networks in the format of a directed edge list
void out_graph(string folder, vector<Graph> A)
{
  for(int a=0;a<L;a++){
    std::ofstream out;
    out.open((folder+"adjacency_"+to_string(a)+".dat").c_str());  // can change it if you want a different output file
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

// Calculate in-degree, summing over all layers
int in_neigh(Vertex i, vector< Graph> A){
	int L=A.size();  // number of layers
	int k=0;
	in_edge_iterator eit, eend;   // Cycle over in-neighbors of i in layer a
	for(int a=0;a<L;a++)for(tie(eit, eend) = in_edges(i,A[a]); eit != eend; ++eit) k++;
	return k;
}

// Calculate out-degree, summing over all layers
int out_neigh(Vertex i, vector< Graph> A){
	int L=A.size();
	int k=0;
	edge_iterator eit, eend;   // Cycle over out-neighbors of i in layer a
	for(int a=0;a<L;a++)for(tie(eit, eend) = out_edges(i,A[a]); eit != eend; ++eit) k++;
	return k;
}

// If out/in -degree is zero --> membership has zero entries
void set_zero_entries(vector< vector<double> > & u, vector< vector<double> > & v,  vector<Graph>  A){
	for(int i=0;i<N;i++){
	  if(out_neigh(i,A)==0)for(int k=0;k<K;k++)u[i][k]=0.;
	  if(in_neigh(i,A)==0)for(int k=0;k<K;k++)v[i][k]=0.;
	}
}
// Outputs a vector of nodes that have at least one in-coming neighbor
vector<int> remove_zero_entries_v(  vector<Graph>  A){
  vector<int> v;
  for(int i=0;i<N;i++)if(in_neigh(i,A)>0)v.push_back(i);
  return v;
}
// Outputs a vector of nodes that do NOT have any in-coming neighbor
vector<int> collect_zero_entries_v(  vector<Graph>  A){
  vector<int> v;
  for(int i=0;i<N;i++)if(in_neigh(i,A)==0)v.push_back(i);
  return v;
}

// Outputs a vector of nodes that have at least one out-going neighbor
vector<int> remove_zero_entries_u(  vector<Graph>  A){
  vector<int> u;
  for(int i=0;i<N;i++)if(out_neigh(i,A)>0)u.push_back(i);
  return u;	
}
// Outputs a vector of nodes that do NOT have any out-going neighbor
vector<int> collect_zero_entries_u(  vector<Graph>  A){

  vector<int> u;
  for(int i=0;i<N;i++)if(out_neigh(i,A)==0)u.push_back(i);
  return u;
}

/* -----------------------------------------------------------------------------------------
   -------------------- Main routine for the EM algorithm ----------------------------------
   ----------------------------------------------------------------------------------------- */

// ------  Update in sequence  -------------
double update_sequential(vector< vector<double> > & u,vector< vector<double> > & u_old,vector< vector<double> > & v,vector< vector<double> > & v_old,vector<vector< vector<double> > > & w,vector< vector< vector<double> > > & w_old, vector<Graph> A, vector<int> u_list,vector<int> v_list,double & d_u,double & d_v,double & d_w){
 
 double deltaL=0.;  // Delta Likelihood
 double dist_u=0.,dist_v=0.,dist_w=0.;  // max difference btw old and nex membership entries --> for convergence

  /* ---------------------------------------------------------------------------
  --------------- 1.  Cycle over u_ik -----------------------------------------
  ---------------------------------------------------------------------------*/
  for(int k=0;k<K;k++){
    // Calculate Z_u, is the same for all u_ik
    double Z_u=0.;
    for(int q=0;q<K;q++){ 
        double w_k=0.,Du=0.;
        for(int a=0;a<L;a++)w_k+=w_old[k][q][a];
        for(int i=0;i<v_list.size();i++)Du+=v_old[v_list[i]][q];
        Z_u+=w_k*Du;  
    }
    // Cycle over i ----
    for(int z=0;z<u_list.size();z++)if(u_old[u_list[z]][k]>0.){     // Update only if u_ik >0 	
	int i=u_list[z];	
	u[i][k]=0.;
	double deltaU=0.;
	for(int a=0;a<L;a++){ // Cycle over layers a -----------
	  edge_iterator eit, eend;   // Cycle over out-neighbors of i in layer a
	  for (tie(eit, eend) = out_edges(i,A[a]); eit != eend; ++eit) {
	    Vertex j=target(*eit,A[a]);   // OUT-EDGE
	    double rho_ijka=0.;
	    double Zij_a=0.;
	    for(int m=0;m<K;m++)for(int l=0;l<K;l++)Zij_a+=u_old[i][m]*v_old[j][l]*w_old[m][l][a];
	    if(Zij_a==0)out_d<<"Zij_a in u "<<A[a][i].name<<"--> "<<A[a][j].name<<" "<<Zij_a<<" "<<a<<endl;
	    if(Zij_a!=0.){
	      for(int q=0;q<K;q++)rho_ijka+=v_old[j][q]*w_old[k][q][a];
	      rho_ijka/=Zij_a;   //  u_old[i][k] will be multiplied only at the end!!!!
	    }  // end if Zij_a
             
	    u[i][k]+=rho_ijka;

	  } //end cycle over out-j  in layer a
	} // end cycle over layer a - --- --- --- - - - - - - -- - - - - - - - - 

	deltaU=Z_u - u[i][k];///u_old[i][k];
	if(Z_u!=0.)u[i][k]*= u_old[i][k]/Z_u;   // final update for u_ik
	else{
	  cout<<" Z_u==0 "<<A[0][i].name<<" "<<u[i][k]<<" "<<u_old[i][k]<<endl;
	  u[i][k]=0.;
	} 
    
	deltaU*=(u[i][k]-u_old[i][k]);
	deltaL+=deltaU;

	if(u[i][k]<err_max)u[i][k]=0.; 
	dist_u=std::max(abs(u[i][k]-u_old[i][k]),dist_u); // calculate max difference btw u_lod and new u membership 
	u_old[i][k]= u[i][k];  // update so that it will appear correctly in the normalization Zij_a

      } // cycle over i --------------------------------------------------
  }  // end cycle over k --------------------------------------------------
  d_u=dist_u;
  /*--------------- 1.  end cycle  over u_ik --------------------------*/
  /* ---------------------------------------------------------------------------
     --------------- 2.  Cycle  over v_ik -----------------------------------------
     ---------------------------------------------------------------------------*/
  for(int k=0;k<K;k++){
    // Calculate Z_u, is the same for all v_ik
    double Z_v=0.;
    for(int q=0;q<K;q++){ 
        double w_k=0.,Dv=0.;
        for(int a=0;a<L;a++)w_k+=w_old[q][k][a];
        for(int i=0;i<u_list.size();i++)Dv+=u[u_list[i]][q];
        Z_v+=w_k*Dv;  
    }
    // Cycle over i ----
     for(int z=0;z<v_list.size();z++)if(v_old[v_list[z]][k]>0.){     // Update only if v_ik >0 
      int i=v_list[z];  
      double deltaV=0.;   // to calculate deltaL, difference in Likelihood old vs new
      v[i][k]=0.;
      for(int a=0;a<L;a++){ // Cycle over layers a -----------
        in_edge_iterator eit, eend;   // Cycle over in-neighbors of i in layer a
        for (tie(eit, eend) = in_edges(i,A[a]); eit != eend; ++eit) {
          Vertex j=source(*eit,A[a]);  // IN-EDGE
          double rho_ijka=0.;
          double Zij_a=0.;
          for(int m=0;m<K;m++)for(int l=0;l<K;l++)Zij_a+=u[j][m]*v_old[i][l]*w_old[m][l][a];
          if(Zij_a==0)out_d<<"Zij_a in v "<<A[a][i].name<<"--> "<<A[a][j].name<<" "<<Zij_a<<" "<<a<<endl;
          if(Zij_a!=0.){
            for(int q=0;q<K;q++)rho_ijka+=u[j][q]*w_old[q][k][a];
            rho_ijka/=Zij_a;   //  v_old[i][k] will be multiplied only at the end!!!!
          }  // end if Zij_a
             
          v[i][k]+=rho_ijka;
    
        } //end cycle over in-j  in layer a
      } // end cycle over layer a - --- --- --- - - - - - - -- - - - - - - - - 
      deltaV=Z_v - v[i][k];
      if(Z_v!=0.)v[i][k]*= v_old[i][k]/Z_v;   // final update for v_ik
      else v[i][k]=0.;
    
      deltaV*=(v[i][k]-v_old[i][k]);
      deltaL+=deltaV;

      if(v[i][k]<err_max)v[i][k]=0.;
      dist_v=std::max(abs(v[i][k]-v_old[i][k]),dist_v);  // update max distance btw new and old v membership
      v_old[i][k]= v[i][k];
     } // cycle over i --------------------------------------------------
  }  // end cycle over k -------------------------------------------------- 
  d_v=dist_v;
  /*--------------- 2.  end cycle  over v_ik --------------------------*/
  /* ---------------------------------------------------------------------------
     --------------- 3.  Cycle  over w_ka -----------------------------------------
     ---------------------------------------------------------------------------*/
  
  for(int k=0;k<K;k++)for(int q=0;q<K;q++){
    // Calculate Z_kq, they are the same for all w_a
    double Z_kq=0.,Du=0.,Dv=0.;
    for(int i=0;i<v_list.size();i++)Dv+=v[v_list[i]][q];
    for(int i=0;i<u_list.size();i++)Du+=u[u_list[i]][k];
    Z_kq=Du*Dv;

    for(int a=0;a<L;a++)if(w_old[k][q][a]>0.){  // Update only if w_ka >0 
	w[k][q][a]=0.;
	double deltaW=0.;
	for(int i=0;i<N;i++){ // Cycle over i ----
	  edge_iterator eit, eend;   // Cycle over out-neighbors of i in layer a
	  double rho_ijka=0.;
	  for (tie(eit, eend) = out_edges(i,A[a]); eit != eend; ++eit) {
	    Vertex j=target(*eit,A[a]); 
	    double Zij_a=0.;
	    for(int m=0;m<K;m++)for(int l=0;l<K;l++)Zij_a+=u[i][m]*v_old[j][l]*w_old[m][l][a];
	    if(Zij_a==0)out_d<<"Zij_a (in updating w_ka) "<<A[a][i].name<<"--> "<<A[a][j].name<<" "<<Zij_a<<" "<<a<<endl;
	    if(Zij_a!=0.)rho_ijka+=v_old[j][q]/Zij_a;   //  w_old[k][a] will be multiplied only at the end!!!!
	  } //end cycle over out-j  of i
	  w[k][q][a]+=u[i][k]*rho_ijka;
	} // end cycle over i - --- --- --- - - - - - - -- - - - - - - - - 

	deltaW=Z_kq - w[k][q][a];
	if(Z_kq!=0.)w[k][q][a]*=w_old[k][q][a]/Z_kq;   // final update for w_ka
	else {w[k][q][a]=0.;
	  cout<<"Z_kq==0 "<<k<<" "<<q<<" "<<a<<endl;
	}   

	deltaW*=(w[k][q][a]-w_old[k][q][a]);
	deltaL+=deltaW;
	if(w[k][q][a]<err_max)w[k][q][a]=0.;
	if(k!=q)dist_w=std::max(abs(w[k][q][a]-w_old[k][q][a]),dist_w); // Update max distance old vs new w_kq
	w_old[k][q][a]=w[k][q][a];
      } // cycle over a --------------------------------------------------
    }  // end cycle over k --------------------------------------------------
  d_w=dist_w;
  
  return deltaL;
} // end update

// ---------------------- end update sequentially --------------------------------------------------------


//  -------   Calculate Likelihood  ------------------------------------------------------------------------------------
double Likelihood(vector< vector<double> > u, vector< vector<double> > v,vector<vector< vector<double> > > w, vector<Graph>  A){

  double l=0.;
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
        if(log_arg>0.)l+=log(log_arg);
      }// end cycle over j      
    }// end cycle over i
  }// end cycle over a
  return l;
}

//--------------------------------------------------

// 2-norm of a K-dimensional vector : can be used to check convergence of the parameters during the update iterations

double distance(vector< vector<double> > u,vector< vector<double> > u_old){
  // Input is the membership u[i][k]  ---> u[0].size() is the number of groups K
  double maxD=0.;
  for(int k=0;k<u[0].size();k++){
    double d=0.;
    for(int i=0;i<u.size();i++)d+=(u[i][k]-u_old[i][k])*(u[i][k]-u_old[i][k]);
    maxD=std::max(maxD,d);  
  }
  return std::sqrt(maxD);
}

double distance_w(vector< vector<vector<double> > > u,vector< vector<vector<double> > > u_old){
  // Input is the affinity matrix w[k][q][a]
  double maxD=0.;
  for(int a=0;a<L;a++){
    double d=0.;
    for(int k=0;k<K;k++){
      for(int i=0;i<K;i++)d+=(u[i][k][a]-u_old[i][k][a])*(u[i][k][a]-u_old[i][k][a]); 
    }
    maxD=std::max(maxD,d);  
  }
  return std::sqrt(maxD);
}
//--------------------------------------------------

/*  
    --------------  Input command-line parameters  -------------- 
    --------------  NEED boost/program_options     --------------
*/

namespace po = boost::program_options;

po::variables_map parse_command_line(int ac, char ** av)
{
  po::options_description desc("Usage: " + string(av[0]) + " <option> ... \n\twhere <option> is one or more of");
  desc.add_options()
    ("help", "produce help message")
    ("folder,f", po::value(&folder)->default_value("Elly/"), "set folder name")
    ("adj,a", po::value(&adj)->default_value("adjacency1.dat"), "set adjacency file name") 
    ("end_file,E", po::value(&end_file)->default_value(".dat"), "set end of file name")
    ("w_file,w", po::value(&w_file)->default_value("w.dat"), "set initialization file for w")
    ("L,l", po::value(&L)->default_value(12), "set number of layers")
    ("intialization,i", po::value(&initialization)->default_value(0), "set intialization flag")
    ("K,k", po::value(&K)->default_value(4), "set number of groups")
    ("N_real,r", po::value(&N_real)->default_value(1), "set number of realizations")
    ("N,N", po::value(&N)->default_value(100), "set number of nodes")
    ("maxit,t", po::value(&maxit)->default_value(10), "set maximum number of iterations")
    ("tolerance,e", po::value(&tolerance)->default_value(1.), "set convergence tolerance")
    ("err,g", po::value(&err)->default_value(0.00001), "setperturbation in w")
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


/* ----------------------------------------------------------------------
   ----------------                MAIN             ---------------------
   --------------------------------------------------------------------- */

int main(int ac, char** av)
{
  
  cout.setf(ios_base::fixed, ios_base::floatfield);
  po::variables_map vm = parse_command_line(ac, av);
  
  //Variables
  //vector<Graph> A(L);    // Commented out if defined before as a global variable instead
  //clock_t t=0.;//time to process one m (N_realizations)
  
  /*-------------  FILE I/O ------------------------------*/
  string file="data/"+folder;
  in1.open(("data/"+folder+adj).c_str());   // adjacency matrix input file
  cout<<"Infile = "<<"data/"+folder+adj<<endl;
  out_d.open((file+"debug_K"+to_string(K)+end_file).c_str()/*, ios::app*/); //  debug output file

  read_graph(in1,A);
  in1.close();

  N=(int)num_vertices(A[0]);
  cout<<"N= "<<N<<endl;
  vector<int> E(L,0);  // number of edges in each layer
  for(int a=0;a<L;a++){
    E[a]=(int)num_edges(A[a]);
    cout<<"E["<<a<<"] = "<<E[a]<<"  density= "<<100.*0.5*((double)E[a]/((double)(N*(N-1.))))<< endl;
  }
  //  out_graph("data/"+folder+"Graph/",A);  // uncomment if you want to output the graph

  vector<int> u_list=remove_zero_entries_u(A);  // Contains vertex indices of vertices having at least one out-edge
  vector<int> v_list=remove_zero_entries_v(A); // Contains vertex indices of vertices having at least one in-edge

/*--------------------------------------------------*/

/* ----------------RESULTS VARIABLES ----------------*/
  vector< vector<double> > u_f(N, vector<double>(K,0.)),v_f(N, vector<double>(K,0.));
  vector< vector< vector<double> > > w_f(K, vector< vector<double> >(K,vector<double>(L,0.)));   // Output results
  double maxL=-inf;

  // Cycle over different realizations ---> try to avoid local minima
  out1.open((file+"Likelihood_K"+to_string(K)+end_file).c_str());

  vector< vector<double> > u(N, vector<double>(K,0.)),v(N, vector<double>(K,0.));
  vector< vector< vector<double> > > w(K, vector< vector<double> >(K,vector<double>(L,0.)));   // Output results
  vector< vector<double> > u_old(N, vector<double>(K,0.)),v_old(N, vector<double>(K,0.));
  vector< vector< vector<double> > > w_old(K, vector< vector<double> >(K,vector<double>(L,0.)));   // Previous time step values 

  int v_length=v_list.size();  // Consider only the nodes that have at least one in-neighbor
  int u_length=u_list.size();  // Consider only the nodes that have at least one out-neighbor

  // Start cycle over realizations ------------------------------
  for(int r=0;r<N_real;r++){
	
    if(initialization==0){  // everything is random
      cout<<" Random initializations"<<endl;
      randomize_w(w);
      randomize_v_in_i(v,v_list);  
      randomize_v_in_i(u,u_list);  
    }
    if(initialization==1){ // everything is initialized
      cout<< " W, U and V are initialized using: "<< file+"v_K"+to_string(K)+w_file<<" ";
      cout<<file+"u_K"+to_string(K)+w_file<<" "; 
      cout<<file+"w_K"+to_string(K)+w_file<<endl; 
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
    if(initialization==2){  // only w is initialized
      cout<< "W is initializaed using: "<< file+"w_K"+to_string(K)+w_file<<endl; 
      in_w.open((file+"w_K"+to_string(K)+w_file).c_str());
      initialize_w(w,in_w);
      in_w.close();
      randomize_v_in_i(v,v_list);  
      randomize_v_in_i(u,u_list);
    }
    if(initialization==3){  // only u and v are initialized
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

    // Update old varaibles
    for(int i=0;i<v_length;i++)for(int k=0;k<K;k++) v_old[v_list[i]][k]=v[v_list[i]][k]; 
    for(int i=0;i<u_length;i++)for(int k=0;k<K;k++) u_old[u_list[i]][k]=u[u_list[i]][k]; 
    for(int k=0;k<K;k++)for(int q=0;q<K;q++)for(int a=0;a<L;a++)w_old[k][q][a]=w[k][q][a];  
  	
    // Convergence variables
    int coincide=0;  // Number of consecutive times error is lower than the threshold
    bool convergence=false; // flag checking for convergence
    int it=0; // Iteration time

    double l2=-inf;   // initialize Likelihood
    //double maxD=0.;   // max distance btw old vs new variables u,v and w
    double delta_u,delta_v,delta_w;

    /* ------------------- Single step iteration update ------------------*/
    while(!convergence and it<maxit){
      
      // out_d<<" it= "<<it<<endl; // output in debug file
      // Re-set variables: in case you want to reset them in every iteration
      /*
	for(int i=0;i<v_length;i++)for(int k=0;k<K;k++) v_old[v_list[i]][k]=v[v_list[i]][k]; 
	for(int i=0;i<u_length;i++)for(int k=0;k<K;k++) u_old[u_list[i]][k]=u[u_list[i]][k]; 
	for(int k=0;k<K;k++)for(int q=0;q<K;q++)for(int a=0;a<L;a++)w_old[k][q][a]=w[k][q][a];	
      */

      // Main EM update: updates membership and calculates max difference new vs old
      update_sequential(u,u_old,v,v_old,w,w_old,A,u_list,v_list,delta_u,delta_v,delta_w);

      // Check new vs old every 10 iterations
      if(it % 10 ==0){    
	double old_L=l2; 
	double maxD=0.;
	/* ----   Use the 2-norm of a vector  -----
	  double delta_u=distance(u,u_old);//cout<<" delta_u= "<<delta_u<<endl;
	  double delta_v=distance(v,v_old);//cout<<" delta_v= "<<delta_v<<endl;
	  double delta_w=distance_w(w,w_old);//cout<<" delta_w= "<<delta_w<<endl;
	*/
	maxD=std::max(maxD,delta_u);maxD=std::max(maxD,delta_v);std::max(maxD,delta_w);
	l2=Likelihood(u,v,w,A);   // Calculate likelihood
  
	out_d<<it<<" "<<maxD<<" "<<l2<<endl;   // output to debug file  
	if(/*maxD*/abs(l2-old_L)<tolerance)coincide++;  // can use maxD or abs(l2-old_L) to check convergence
	else coincide=0;
      }   // end if it multiple of 10

      out1<<it<< " "<<l2<<" "<<delta_u<<" "<<delta_v<<" "<<delta_w<< endl;    // output to likelihood file
     // cout<<it<< " "<<l2<<" "<<delta_u<<" "<<delta_v<<" "<<delta_w<< endl;    // output to likelihood file

      if(coincide==decision)convergence=true;
      if(l2<maxL and it>teq)convergence=true;
      it++;
    }// end while convergence

    l2=Likelihood(u,v,w,A);
    cout<<"r= "<<r<<" Final Likelihood= "<<l2<<endl;   // for this realization
    out_d<<"r= "<<r<<" Final Likelihood= "<<l2<<endl;  //  output to debug file
    out_d<<"----------------------"<<endl;

    // Keep track of the realization that achieves the max likelihood
    if(maxL<l2){
      maxL=l2;
      u_f=u;v_f=v;w_f=w;
    }

  } // end cycle over realizations

  maxL=(int)Likelihood(u_f,v_f,w_f,A); 
  cout<<" Final Likelihood= "<<maxL<<endl;  // over the all realizations
  out1.close();     // close likelihood file

  // Output results
  out1.open((file+"u_K"+to_string(K)+end_file).c_str());
  out2.open((file+"v_K"+to_string(K)+end_file).c_str());
  out3.open((file+"w_K"+to_string(K)+end_file).c_str());
  out1<<"# Max likelihood= "<<maxL<<" N_real="<<N_real<<endl;
  out2<<"# Max likelihood= "<<maxL<<" N_real="<<N_real<<endl;
  out3<<"# Max likelihood= "<<maxL<<" N_real="<<N_real<<endl;
  cout<< " Data saved in "<<file+"u_K"+to_string(K)+end_file;
  cout<< "  "<<file+"v_K"+to_string(K)+end_file;
  cout<< "  "<<file+"w_K"+to_string(K)+end_file<<endl;
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
  for(int a=0;a<L;a++){
    out3<<"a= "<<a<<endl;
    for(int k=0;k<K;k++){
      for(int q=0;q<K;q++){out3<<w_f[k][q][a]<<" ";cout<<w_f[k][q][a]<<" ";}
      out3<<endl; cout<<endl; 
    }
    out3<<endl; cout<<endl;
  }

  out1.close();out2.close();out3.close();
  out_d.close();
  // Output graph information --------------------
  cout<<"N= "<<N<<endl;
  for(int a=0;a<L;a++){
    E[a]=(int)num_edges(A[a]);
    cout<<"E["<<a<<"] = "<<E[a]<<"  density= "<<100.*0.5*((double)E[a]/((double)(N*(N-1.))))<< endl;
  }
  
return 1;
}


