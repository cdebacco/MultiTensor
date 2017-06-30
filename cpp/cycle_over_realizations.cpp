

#include "cycle_over_realizations.hpp"

extern double err_max,inf;

void cycle_over_realizations(int N_real,double err, double inf,int initialization,int maxit, int tolerance,int decision,string file, string end_file,string w_file,vector<Graph> & A, vector<int> & u_list,vector<int> & v_list,vector< vector<double> > & u,vector< vector<double> > & v, vector< vector<double> > & u_f,vector< vector<double> > & v_f,vector< vector<double> > & u_old,vector< vector<double> > & v_old,vector< vector< vector<double> > > & w,vector< vector< vector<double> > > & w_f,vector< vector< vector<double> > > & w_old)
{
  double maxL=-inf;; 
  for(int r=0;r<N_real;r++){
  	
      initialize(initialization,u_list,v_list,u,v,w,file,w_file,A);
    	update_old_variables(u_old,v_old,u,v,w_old,w,u_list,v_list);

      // Convergence variables
      int coincide=0;  // Number of consecutive times error is lower than the threshold
      bool convergence=false; // flag checking for convergence
      int it=0;
      double l2=-inf;   // initialize Likelihood
      double delta_u,delta_v,delta_w;

      /* ------------------- Single step iteration update ------------------*/
      while(!convergence and it<maxit){
        // Main EM update: updates membership and calculates max difference new vs old
        update_em(err_max,u,u_old,v,v_old,w,w_old,A,u_list,v_list,delta_u,delta_v,delta_w);

        check_for_convergence(A,it,l2,tolerance,coincide,decision,convergence,u,v,w);
      }  // end while   
      cout<<"r="<<r<<" final Likelihood= "<<l2<<" iterations:"<<it<<endl;  // over the all realizations
      
      update_optimal_parameters(maxL,l2,u_f,v_f,u,v,w_f,w);
      cout<<"MaxL= "<<maxL<<endl;

    } // end cycle over realizations

    maxL=(int)Likelihood(u_f,v_f,w_f,A); 

    cout<<" Final Likelihood= "<<maxL<<endl;  // over the all realizations
    output_results(A,file, end_file,N_real,maxL,u_f,v_f, w_f);

}

void cycle_over_realizations(int N_real,double err, double inf,int initialization,int maxit, int tolerance,int decision,string file, string end_file,string w_file,vector<Graph> & A, vector<int> & u_list,vector<int> & v_list,vector< vector<double> > & u,vector< vector<double> > & v, vector< vector<double> > & u_f,vector< vector<double> > & v_f,vector< vector<double> > & u_old,vector< vector<double> > & v_old,vector< vector<double> >  & w,vector<vector<double> > & w_f, vector< vector<double> > & w_old) 
{
  double maxL=-inf;;  
  for(int r=0;r<N_real;r++){
    
      initialize(initialization,u_list,v_list,u,v,w,file,w_file,A);
      update_old_variables(u_old,v_old,u,v,w_old,w,u_list,v_list);

      // Convergence variables
      int coincide=0;  // Number of consecutive times error is lower than the threshold
      bool convergence=false; // flag checking for convergence
      int it=0;
      double l2=-inf;   // initialize Likelihood
      double delta_u,delta_v,delta_w;

      /* ------------------- Single step iteration update ------------------*/
      while(!convergence and it<maxit){
        // Main EM update: updates membership and calculates max difference new vs old
        update_em(err_max,u,u_old,v,v_old,w,w_old,A,u_list,v_list,delta_u,delta_v,delta_w);

        check_for_convergence(A,it,l2,tolerance,coincide,decision,convergence,u,v,w);
      }  // end while   

      cout<<"r="<<r<<" final Likelihood= "<<l2<<" iterations:"<<it<<endl;  // over the all realizations
      update_optimal_parameters(maxL,l2,u_f,v_f,u,v,w_f,w);


    } // end cycle over realizations

    maxL=(int)Likelihood(u_f,v_f,w_f,A); 

    cout<<" Final Likelihood= "<<maxL<<endl;  // over the all realizations
    output_results(A,file, end_file,N_real,maxL,u_f,v_f, w_f);

}



// Selects btw full and restricted (Diagonal) models
void iterate(bool Diagonal,int N_real,double err, double inf,int initialization,int maxit, int tolerance,int decision,string file, string end_file,string w_file,vector<Graph> & A, vector<int> & u_list,vector<int> & v_list,vector< vector<double> > & u,vector< vector<double> > & v, vector< vector<double> > & u_f,vector< vector<double> > & v_f,vector< vector<double> > & u_old,vector< vector<double> > & v_old)
{
    if(Diagonal){
      vector< vector<double> > w(K, vector<double>(L));
      vector< vector<double> > w_old(K, vector<double>(L));
      vector< vector<double> > w_f(K, vector<double>(L));

      cycle_over_realizations(N_real,err, inf,initialization,maxit,tolerance,decision,file, end_file,w_file,A,u_list,v_list,u,v,u_f,v_f,u_old,v_old,w,w_f,w_old);

    }
    else{
      vector< vector< vector<double> > > w(K, vector< vector<double> >(K,vector<double>(L,0.)));
      vector< vector< vector<double> > > w_old(K, vector< vector<double> >(K,vector<double>(L,0.)));
      vector< vector< vector<double> > > w_f(K, vector< vector<double> >(K,vector<double>(L,0.)));
      
      cycle_over_realizations(N_real,err, inf,initialization,maxit,tolerance,decision,file, end_file,w_file,A,u_list,v_list,u,v,u_f,v_f,u_old,v_old,w,w_f,w_old);

    }
 } 

