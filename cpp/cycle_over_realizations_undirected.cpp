#include "cycle_over_realizations_undirected.hpp"

extern double err_max,inf;

void cycle_over_realizations_undirected(int N_real,double err, double inf,int initialization,int maxit, int tolerance,int decision,string file, string end_file,string w_file,vector<Graph_undirected> & A, vector<int> & u_list,vector< vector<double> > & u, vector< vector<double> > & u_f,vector< vector<double> > & u_old,vector< vector< vector<double> > >  & w, vector< vector< vector<double> > >& w_f, vector< vector< vector<double> > > & w_old) 
{
  double maxL=-inf;;  
  for(int r=0;r<N_real;r++){
    
      initialize_undirected(initialization,u_list,u,w,file,w_file,A);
      update_old_variables_undirected(u_old,u,w_old,w,u_list);

      // Convergence variables
      int coincide=0;  // Number of consecutive times error is lower than the threshold
      bool convergence=false; // flag checking for convergence
      int it=0;
      double l2=-inf;   // initialize Likelihood
      double delta_u,delta_w;

      /* ------------------- Single step iteration update ------------------*/
      while(!convergence and it<maxit){
        // Main EM update: updates membership and calculates max difference new vs old
        update_em_undirected(err_max,u,u_old,w,w_old,A,u_list,delta_u,delta_w);

        check_for_convergence_undirected(A,it,l2,tolerance,coincide,decision,convergence,u,w);
      }  // end while   

      cout<<"r="<<r<<" final Likelihood= "<<l2<<" iterations:"<<it<<endl;  // over the all realizations
      update_optimal_parameters_undirected(maxL,l2,u_f,u,w_f,w);
      cout<<"Max L= "<<maxL<<endl;
    } // end cycle over realizations

    maxL=(int)Likelihood_undirected(u_f,w_f,A); 

    cout<<" Final Likelihood= "<<maxL<<endl;  // over the all realizations
    output_results_undirected(A,file, end_file,N_real,maxL,u_f, w_f);

}

// Assortative version 
void cycle_over_realizations_undirected(int N_real,double err, double inf,int initialization,int maxit, int tolerance,int decision,string file, string end_file,string w_file,vector<Graph_undirected> & A, vector<int> & u_list,vector< vector<double> > & u, vector< vector<double> > & u_f,vector< vector<double> > & u_old,vector< vector<double> >  & w,vector<vector<double> > & w_f, vector< vector<double> > & w_old) 
{
  double maxL=-inf;;  
  for(int r=0;r<N_real;r++){
    
      initialize_undirected(initialization,u_list,u,w,file,w_file,A);
      update_old_variables_undirected(u_old,u,w_old,w,u_list);

      // Convergence variables
      int coincide=0;  // Number of consecutive times error is lower than the threshold
      bool convergence=false; // flag checking for convergence
      int it=0;
      double l2=-inf;   // initialize Likelihood
      double delta_u,delta_w;

      /* ------------------- Single step iteration update ------------------*/
      while(!convergence and it<maxit){
        // Main EM update: updates membership and calculates max difference new vs old
        update_em_undirected(err_max,u,u_old,w,w_old,A,u_list,delta_u,delta_w);

        check_for_convergence_undirected(A,it,l2,tolerance,coincide,decision,convergence,u,w);
      }  // end while   

      cout<<"r="<<r<<" final Likelihood= "<<l2<<" iterations:"<<it<<endl;  // over the all realizations
      update_optimal_parameters_undirected(maxL,l2,u_f,u,w_f,w);

    } // end cycle over realizations

    maxL=(int)Likelihood_undirected(u_f,w_f,A); 

    cout<<" Final Likelihood= "<<maxL<<endl;  // over the all realizations
    output_results_undirected(A,file, end_file,N_real,maxL,u_f, w_f);

}

// Selects btw full and restricted (Diagonal=Assortative for undirected networks) models, undirected network case: u=v
void iterate_undirected(bool assortative,int N_real,double err, double inf,int initialization,int maxit, int tolerance,int decision,string file, string end_file,string w_file,vector<Graph_undirected> & A, vector<int> & u_list,vector< vector<double> > & u, vector< vector<double> > & u_f,vector< vector<double> > & u_old)
{
    if(assortative){
      vector< vector<double> > w(K, vector<double>(L));
      vector< vector<double> > w_old(K, vector<double>(L));
      vector< vector<double> > w_f(K, vector<double>(L));

      cycle_over_realizations_undirected(N_real,err, inf,initialization,maxit,tolerance,decision,file, end_file,w_file,A,u_list,u,u_f,u_old,w,w_f,w_old);

    }
    else{
      vector< vector< vector<double> > > w(K, vector< vector<double> >(K,vector<double>(L,0.)));
      vector< vector< vector<double> > > w_old(K, vector< vector<double> >(K,vector<double>(L,0.)));
      vector< vector< vector<double> > > w_f(K, vector< vector<double> >(K,vector<double>(L,0.)));
      
      cycle_over_realizations_undirected(N_real,err, inf,initialization,maxit,tolerance,decision,file, end_file,w_file,A,u_list,u,u_f,u_old,w,w_f,w_old);

    }
 } 
