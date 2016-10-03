
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
int teq=500;

unsigned int seed=0;
unsigned int rseed=0;

double inf=1e10; 
double err_max=0.00001; // max err in the converge() routine
double err=0.01;
int initialization=0;

/* ----------------I/O FILES ----------------*/
std::ofstream out1,out2,out3,out_d;
std::ifstream in1, in2 ,in_w;;
std::string folder="Elly";
std::string end_file="";
std::string adj="adjacency1.dat";
std::string w_file="w.dat";

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

void randomize_v_in_i(vector< vector<double> > & v, vector<int> v_list) {
  // Assign a random weight to each entry
  int list_length=v_list.size()	;
    for (int k=0; k<K; k++ )for(int i=0;i<list_length;i++){ 
    	int j=v_list[i];
	v[j][k] =real01();
      }
}

void randomize_w(vector<vector< vector<double> > > & w) {
  for(int i=0;i<L;i++) for (int k=0; k<K; k++ )for (int q=0; q<K; q++ ){
	if(k==q)w[k][q][i]=real01();  // for better convergence ...
	else w[k][q][i]=w[q][k][i]=err*real01();
      }
}

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
    for(int a=0;a<L;a++)for(int k=0;k<K;k++)for(int q=0;q<K;q++){
	  w[k][q][a]+=w[k][k][a]*err*real01();    
	}
}

void initialize_v(vector< vector<double> > & v, istream & file) {
    int a=-1;double max_entry=0.;
    for(string line; getline( file, line ); ){
      std::istringstream is( line );
      if(a>=0){
      int i;
      is>>i; // node index from 0 to N-1
      for(int k=0;k<K;k++){is>>v[a][k];max_entry=std::max(max_entry,v[a][k]);}
      }
      a++;
    }
    // Add noise (perturbation)
    for(int i=0;i<N;i++)for(int k=0;k<K;k++)v[i][k]+=max_entry*err*real01();
}
//-------------------------------------------------------------------------------
//--------------------------------------------------------------------

struct VertexProperties  {
  VertexProperties() : name("0"){}
  VertexProperties(string const & name) : name(name){}
  string name;
};
//--------------------------------------------------------------------
typedef adjacency_list<setS, vecS,/*directedS,undirectedS bidirectionalS*/undirectedS,VertexProperties
		       /*, EdgeProperties*/> GraphBase;

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
    } 
    out.close();
  }
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
  for(int i=0;i<N;i++)if(out_neigh(i,A)==0)for(int k=0;k<K;k++)u[i][k]=0.;
}

vector<int> remove_zero_entries_u(  vector<Graph>  A){
	vector<int> u;
	for(int i=0;i<N;i++)if(out_neigh(i,A)>0)u.push_back(i);else cout<<A[0][i].name<<endl;
	return u;	
}

vector<int> collect_zero_entries_u(  vector<Graph>  A){
  vector<int> v;
  for(int i=0;i<N;i++)if(out_neigh(i,A)==0)v.push_back(i);
  return v;
}


// ------  Update in sequence  -------------
double update_sequential(vector< vector<double> > & u,vector< vector<double> >  & u_old,vector<vector< vector<double> > > & w,vector< vector< vector<double> > > & w_old, vector<Graph> A, vector<int> u_list,double & d_u,double & d_w){
 
 double deltaL=0.;  // Delta Likelihood
 double dist_u=0.,dist_v=0.,dist_w=0.;
  /* ---------------------------------------------------------------------------
  --------------- 1.  Cycle over u_ik -----------------------------------------
  ---------------------------------------------------------------------------*/
  for(int k=0;k<K;k++){
    // Calculate Z_u, is the same for all u_ik
    double Z_u=0.;
    for(int q=0;q<K;q++){ 
      double w_k=0.,Du=0.;
      for(int a=0;a<L;a++)w_k+=w_old[k][q][a];
      for(int i=0;i<u_list.size();i++)Du+=u_old[u_list[i]][q];
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
          for(int q=0;q<K;q++){
            double Zij_a=0.;
            for(int q=0;q<K;q++)for(int l=0;l<K;l++)Zij_a+=u_old[i][q]*u_old[j][l]*w_old[q][l][a];            
            if(Zij_a==0.)out_d<<"Zij_a in u "<<A[a][i].name<<"--> "<<A[a][j].name<<" "<<Zij_a<<" "<<a<<endl;
            if(Zij_a!=0.)rho_ijka+=u_old[j][q]*w_old[k][q][a]/Zij_a;   //  u_old[i][k] will be multiplied only at the end!!!!
          }  // end cycle over q
             
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
      dist_u=std::max(abs(u[i][k]-u_old[i][k]),dist_u);
      u_old[i][k]= u[i][k];
    } // cycle over i --------------------------------------------------
  }  // end cycle over k --------------------------------------------------
  d_u=dist_u;
    /*--------------- 1.  end cycle  over u_ik --------------------------*/ 
  /* ---------------------------------------------------------------------------
  --------------- 3.  Cycle  over w_ka -----------------------------------------
  ---------------------------------------------------------------------------*/
  for(int k=0;k<K;k++)for(int q=0;q<K;q++){
    // Calculate Z_kq, they are the same for all w_a
    double Z_kq=0.;
    for(int i=0;i<N;i++)for(int j=0;j<N;j++) Z_kq+=u[i][k]*u[j][q];

    for(int a=0;a<L;a++)if(w_old[k][q][a]>0.){  // Update only if w_ka >0 
      w[k][q][a]=0.;
      double deltaW=0.;

      for(int i=0;i<N;i++){ // Cycle over i ----
        edge_iterator eit, eend;   // Cycle over out-neighbors of i in layer a
        double rho_ijka=0.;
        for (tie(eit, eend) = out_edges(i,A[a]); eit != eend; ++eit) {
          Vertex j=target(*eit,A[a]);
          double Zij_a=0.;
          for(int m=0;m<K;m++)for(int l=0;l<K;l++)Zij_a+=u_old[i][m]*u_old[j][l]*w_old[m][l][a];
          if(Zij_a==0)out_d<<"Zij_a (in updating w_ka) "<<A[a][i].name<<"--> "<<A[a][j].name<<" "<<Zij_a<<" "<<a<<endl;
          if(Zij_a!=0.)rho_ijka+=u_old[j][q]/Zij_a;   //  w_old[k][a] will be multiplied only at the end!!!!
        } //end cycle over out-j  of i
        w[k][q][a]+=u_old[i][k]*rho_ijka;
      } // end cycle over i - --- --- --- - - - - - - -- - - - - - - - - 
      deltaW=Z_kq - w[k][q][a];///w_old[k][a];
      if(Z_kq!=0.)w[k][q][a]*=w_old[k][q][a]/Z_kq;   // final update for w_ka
      else {w[k][q][a]=0.;
        cout<<"Z_kq==0 "<<k<<" "<<q<<" "<<a<<endl;
      }

      deltaW*=(w[k][q][a]-w_old[k][q][a]);
      deltaL+=deltaW;
      if(w[k][q][a]<err_max)w[k][q][a]=0.;
      if(k!=q)dist_w=std::max(abs(w[k][q][a]-w_old[k][q][a]),dist_w);
      w_old[k][q][a]=w[k][q][a];
    } // cycle over a --------------------------------------------------
  }  // end cycle over k --------------------------------------------------
  d_w=dist_w;

  return -deltaL;
} // end update
// ---------------------- end update sequentially


 //-------Calculate Likelihood ---------------------
double Likelihood(vector< vector<double> > u,vector<vector< vector<double> > > w, vector<Graph>  A){
  double l=0.;
  for(int a=0;a<L;a++){
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	double log_arg=0.;
        for(int k=0;k<K;k++)for(int q=0;q<K;q++){
	    double uvw=u[i][k]*u[j][q]*w[k][q][a];
	    l-=uvw;
	    if(edge(i, j, A[a]).second)log_arg+=uvw; // consider this term only if edge exists
        }// end cycle over k and q
        if(log_arg>0.)l+=log(log_arg);
      }// end cycle over j
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
    ("adj,a", po::value(&adj)->default_value("adjacency.dat"), "set adjacency file name") 
    ("end_file,E", po::value(&end_file)->default_value(".dat"), "set end of file name")
    ("w_file,w", po::value(&w_file)->default_value("w.dat"), "set initialization file for w")
    ("L,l", po::value(&L)->default_value(12), "set number of layers")
    ("intialization,i", po::value(&initialization)->default_value(0), "set intialization flag")
    ("K,k", po::value(&K)->default_value(4), "set number of groups")
    ("N_real,r", po::value(&N_real)->default_value(1), "set number of realizations")
    ("N,N", po::value(&N)->default_value(100), "set number of nodes")
    ("maxit,t", po::value(&maxit)->default_value(500), "set maximum number of iterations")
    ("tolerance,e", po::value(&tolerance)->default_value(1.), "set convergence tolerance")
    ("err,g", po::value(&err)->default_value(0.0001), "setperturbation in w")
    ("seed,s", po::value<unsigned>(), "sets instance seed")
    ("rseed,z", po::value<unsigned>(), "sets biases seed")
    ("decision,y", po::value<int>(&decision)->default_value(5), "program converges after this # repeats of the decision variables");

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

  //clock_t t=0.;//time to process one m (N_realizations)
  
  /*-------------  FILE I/O ------------------------------*/
  string file="data/"+folder;
  in1.open(("data/"+folder+adj).c_str());
  cout<<"Infile = "<<"data/"+folder+adj<<endl;
  //  out_d.open((file+"debug_K"+to_string(K)+"_i"+to_string(initialization)+".dat").c_str()/*, ios::app*/);

  read_graph(in1,A);
  in1.close();

  N=(int)num_vertices(A[0]);
  cout<<"N= "<<N<<endl;
  vector<int> E(L,0);
  for(int a=0;a<L;a++){
    E[a]=(int)num_edges(A[a]);
    cout<<"E["<<a<<"] = "<<E[a]<<"  density= "<<100.*0.5*((double)E[a]/((double)(N*(N-1.))))<< endl;
  }
  //  out_graph("data/"+folder+"Graph/",A);

  vector<int> u_list=remove_zero_entries_u(A);  // Contains vertex indices of vertices having  out-edges

  /*--------------------------------------------------*/

/* ----------------RESULTS VARIABLES ----------------*/
  vector< vector<double> > u_f(N, vector<double>(K,0.));
 vector< vector< vector<double> > > w_f(K, vector< vector<double> >(K,vector<double>(L,0.)));   // Output results
 double maxL=-inf;

 // out1.open((file+"Likelihood_K"+to_string(K)+"_i"+to_string(initialization)+end_file).c_str());
 
 vector< vector<double> > u(N, vector<double>(K,0.));
 vector< vector< vector<double> > > w(K, vector< vector<double> >(K,vector<double>(L,0.)));   // Output results
 vector< vector<double> > u_old(N, vector<double>(K,0.));
 vector< vector< vector<double> > > w_old(K, vector< vector<double> >(K,vector<double>(L,0.)));   // Previous time step values 
 int u_length=u_list.size();
// Cycle over different realizations ---> try to avoid local minima
 for(int r=0;r<N_real;r++){
 
   if(initialization==0){  // everything is random
     cout<<" Random initializations"<<endl;
     randomize_w(w);
     randomize_v_in_i(u,u_list);
   }
   if(initialization==1){ // everything is initialized
     cout<< " W and U are initialized using:  ";
     cout<<file+"u_K"+to_string(K)+w_file<<" "; 
     cout<<file+"w_K"+to_string(K)+w_file<<endl; 
     in_w.open((file+"w_K"+to_string(K)+w_file).c_str());
     initialize_w(w,in_w);
     in_w.close();
     in_w.open((file+"u_K"+to_string(K)+w_file).c_str());
     initialize_v(u,in_w);
     in_w.close();

     vector<int> u_to_be_removed=collect_zero_entries_u(A);
     for(int i=0;i<u_to_be_removed.size();i++)for(int k=0;k<K;k++)u[u_to_be_removed[i]][k]=0.;
   }
   if(initialization==2){  // only w is initialized
     cout<< "W is initializaed using: "<< file+"w_K"+to_string(K)+w_file<<endl; 
     in_w.open((file+"w_K"+to_string(K)+w_file).c_str());
     initialize_w(w,in_w);
     in_w.close();
     randomize_v_in_i(u,u_list);
   }
   if(initialization==3){  // only u and v are initialized
     cout<< " U is initialized using: ";
     cout<<file+"u_K"+to_string(K)+w_file<<endl;
     in_w.open((file+"u_K"+to_string(K)+w_file).c_str());
     initialize_v(u,in_w);
     in_w.close();
     vector<int> u_to_be_removed=collect_zero_entries_u(A);
     for(int i=0;i<u_to_be_removed.size();i++)for(int k=0;k<K;k++)u[u_to_be_removed[i]][k]=0.;
     randomize_w(w);
   }

   // Update Old variables
   for(int i=0;i<u_length;i++)for(int k=0;k<K;k++) u_old[u_list[i]][k]=u[u_list[i]][k]; 
   for(int k=0;k<K;k++)for(int q=0;q<K;q++)for(int a=0;a<L;a++)w_old[k][q][a]=w[k][q][a];  
      
  // Convergence variables
  int coincide=0;  // Number of consecutive times error is lower than the threshold
  bool convergence=false; 
  int it=0; // Iteration time

  double l2=-inf;
  // double maxD=0.;
  double delta_u,delta_v,delta_w;
/* ------------------- Single step iteration update ------------------*/
  while(!convergence and it<maxit){
    // Main EM update step
    update_sequential(u,u_old,w,w_old,A,u_list,delta_u,delta_w);
    // Check new vs old every 10 iterations 
    if(it % 10 ==0){
      double old_L=l2; 
      double maxD=0.;
   
      /* If you want to use 2-norm
	double delta_u=distance(u,u_old);//cout<<" delta_u= "<<delta_u<<endl;
	double delta_v=distance(v,v_old);//cout<<" delta_v= "<<delta_v<<endl;
	double delta_w=distance(w,w_old);//cout<<" delta_w= "<<delta_w<<endl;
      */
      maxD=std::max(maxD,delta_u);std::max(maxD,delta_w);

      l2=Likelihood(u,w,A);
      //      out_d<<it<<" "<<maxD<<" "<<l2<<" "<<l2<<endl;

      if(/*maxD*/abs(l2-old_L)<tolerance)coincide++;  // can use maxD or difference in likelihood
      else coincide=0;
    }

    //    out1<<it<< " "<<l2<<" "<<delta_u<<" "<<delta_w<< endl;

  if(coincide==decision)convergence=true;
  if(l2<maxL and it>teq)convergence=true;
  it++;
  }// end while convergence

  l2=Likelihood(u,w,A);
  cout<<"r= "<<r<<" Final Likelihood= "<<l2<<endl;
  //out_d<<"r= "<<r<<" Final Likelihood= "<<l2<<endl;
  //  out_d<<"----------------------"<<endl;

  // Keep track of the realization that achieves the max likelihood
  if(maxL<l2){
    maxL=l2;
    u_f=u;w_f=w;
  }

 } // end cycle over realizations
 maxL=(int)Likelihood(u_f,w_f,A); 
 cout<<" Final Likelihood= "<<maxL<<endl;
 // out1.close();
 // Output results
 out1.open((file+"u_K"+to_string(K)+"_i"+to_string(initialization)+end_file).c_str());
 out3.open((file+"w_K"+to_string(K)+"_i"+to_string(initialization)+end_file).c_str());
 out1<<"# Max likelihood= "<<maxL<<" N_real="<<N_real<<endl;
 out3<<"# Max likelihood= "<<maxL<<" N_real="<<N_real<<endl;
 cout<< " Data saved in "<<file+"u_K"+to_string(K)+"_i"+to_string(initialization)+end_file;
 cout<< "  "<<file+"w_K"+to_string(K)+"_i"+to_string(initialization)+end_file<<endl;
 for(int i=0;i<N;i++){
  out1<<A[0][i].name<<" ";
  for(int k=0;k<K;k++){
  out1<<u_f[i][k]<<" ";
  } 
  out1<<endl;
 } 

 for(int a=0;a<L;a++){
  out3<<"a= "<<a<<endl;
  for(int k=0;k<K;k++){
    for(int q=0;q<K;q++)out3<<w_f[k][q][a]<<" ";
    out3<<endl;  
  }
  out3<<endl;
 }

 out1.close();out3.close();
 // out_d.close();

 cout<<"N= "<<N<<endl;
 
  for(int a=0;a<L;a++){
    E[a]=(int)num_edges(A[a]);
    cout<<"E["<<a<<"] = "<<E[a]<<"  density= "<<100.*0.5*((double)E[a]/((double)(N*(N-1.))))<< endl;
  }

  return 1;
}


