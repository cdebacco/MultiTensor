

#include <boost/program_options.hpp>

#include "main.hpp"
#include "cycle_over_realizations.hpp"

#include <iostream>
#include <ctime>
#include <utility>

using namespace boost;
using namespace std;
namespace po = boost::program_options;

int read_graph(string folder,string adj, vector<Graph> & A);
void out_graph(string folder, vector<Graph> A);
vector<int> remove_zero_entries_v(  vector<Graph>  A);
vector<int> remove_zero_entries_u(  vector<Graph>  A);
void print_graph_stat(std::vector<Graph> );

/*  
    --------------  Input command-line parameters  -------------- 
    --------------  NEED boost/program_options     --------------
*/
po::variables_map parse_command_line(int ac, char ** av)
{
  po::options_description desc("Usage: " + string(av[0]) + " <option> ... \n\twhere <option> is one or more of");
  desc.add_options()
    ("help", "produce help message")
    ("folder,f", po::value(&folder)->default_value(""), "set folder name")
    ("adj,a", po::value(&adj)->default_value("adjacency.dat"), "set adjacency file name") 
    ("end_file,E", po::value(&end_file)->default_value(".dat"), "set end of file name")
    ("w_file,w", po::value(&w_file)->default_value("w.dat"), "set initialization file for w")
    ("L,l", po::value(&L)->default_value(4), "set number of layers")
    ("intialization,i", po::value(&initialization)->default_value(0), "set intialization flag")
    ("K,k", po::value(&K)->default_value(5), "set number of groups")
    ("N_real,r", po::value(&N_real)->default_value(1), "set number of realizations")
    ("maxit,t", po::value(&maxit)->default_value(500), "set maximum number of iterations")
    ("tolerance,e", po::value(&tolerance)->default_value(0.1), "set convergence tolerance")
    ("err,g", po::value(&err)->default_value(0.1), "set perturbation in w")
    ("output,o", po::value(&out_adjacency)->default_value(false), "flag to output single-layer adjacencies")
    ("assortative,A", po::value(&assortative)->default_value(false), "flag to restrict to assortative model")
    ("rseed,z", po::value<unsigned>(), "sets random number seed")
    ("decision,y", po::value<int>(&decision)->default_value(10), "program converges after this # repeats of the decision variables");

  po::variables_map vm;
  po::store(po::parse_command_line(ac, av, desc), vm);
  po::notify(vm);

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
 
  vector<Graph> A(L);  // a vector of graphs, one for each of the L layers; local variable 
  string file="../data/"+folder;

  N=read_graph(folder,adj,A);   // number of nodes
  print_graph_stat(A);
  
  if(out_adjacency) out_graph(folder,A);  // uncomment if you want to output the graph

  vector<int> u_list=remove_zero_entries_u(A);  // Contains vertex indices of vertices having at least one out-edge
  vector<int> v_list=remove_zero_entries_v(A); // Contains vertex indices of vertices having at least one in-edge

/* ----------------RESULTS VARIABLES ----------------*/
  vector< vector<double> > u_f(N, vector<double>(K,0.)),v_f(N, vector<double>(K,0.));
  vector< vector<double> > u(N, vector<double>(K,0.)),v(N, vector<double>(K,0.));
  vector< vector<double> > u_old(N, vector<double>(K,0.)),v_old(N, vector<double>(K,0.));

  time_t tstart, tend; 
  tstart = time(0);

  iterate(assortative,N_real,err, inf,initialization,maxit,tolerance,decision,file, end_file,w_file,A,u_list,v_list,u,v,u_f,v_f,u_old,v_old);

  print_graph_stat(A);

  tend = time(0); 
  cout << "It took "<< difftime(tend, tstart) <<" second(s). "<< endl;

return 1;
}


