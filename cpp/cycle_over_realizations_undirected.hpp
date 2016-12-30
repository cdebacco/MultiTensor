
#ifndef CORU_H
#define CORU_H

#include <iostream>
#include "mlg.hpp"

extern int L;
extern int K;
extern int N;

// Undirected version u==v
double Likelihood_undirected(vector< vector<double> > & u, vector<vector< vector<double> > > & w, vector<Graph_undirected> & A);
double Likelihood_undirected(vector< vector<double> > & u, vector<vector<double> > & w, vector<Graph_undirected> & A);
void initialize_undirected(int initialization,  vector<int> u_list, vector< vector<double> > & u, vector< vector< vector<double> > > & w, string file, string w_file,vector<Graph_undirected> A);
void initialize_undirected(int initialization,  vector<int> u_list,  vector< vector<double> > & u, vector< vector< double> > & w, string file, string w_file,vector<Graph_undirected> A);
void update_old_variables_undirected(  vector< vector<double> > & u_old,   vector< vector<double> > u, vector< vector< vector<double> > > & w_old,  vector< vector< vector<double> > > w, vector<int> u_list);
void update_old_variables_undirected(  vector< vector<double> > & u_old,  vector< vector<double> > u,  vector< vector<double> > & w_old, vector< vector<double> > w, vector<int> u_list);
void check_for_convergence_undirected(vector<Graph_undirected> A,int & it,double & l2,double tolerance,int & coincide, int decision, bool & convergence,vector< vector<double> > u, vector< vector< vector<double> > > w);
void check_for_convergence_undirected(vector<Graph_undirected> A,int & it,double & l2,double tolerance,int & coincide, int decision, bool & convergence,vector< vector<double> > u, vector< vector<double> > w);
void update_em_undirected(double err_max, vector< vector<double> > & u,vector< vector<double> > & u_old,vector<vector< vector<double> > > & w,vector< vector< vector<double> > > & w_old, vector<Graph_undirected> & A, vector<int> & u_list,double & d_u,double & d_w);
void update_em_undirected(double err_max, vector< vector<double> > & u,vector< vector<double> > & u_old,vector< vector<double> >  & w, vector< vector<double> >  & w_old, vector<Graph_undirected> & A, vector<int> & u_list,double & d_u,double & d_w);
void update_optimal_parameters_undirected(double & maxL, double l2,vector< vector<double> > & u_f,vector< vector<double> > u, vector< vector< vector<double> > > & w_f, vector< vector< vector<double> > > w);
void update_optimal_parameters_undirected(double & maxL, double l2,vector< vector<double> > & u_f,vector< vector<double> > u, vector< vector<double> > & w_f,  vector< vector<double> > w);
void output_results_undirected(std::vector<Graph_undirected> A,string file,string end_file,int N_real,double maxL,vector< vector<double> > u_f,vector< vector< vector<double> > > w_f);
void output_results_undirected(std::vector<Graph_undirected> A,string file,string end_file,int N_real,double maxL,vector< vector<double> > u_f,vector< vector<double> > w_f);
void cycle_over_realizations_undirected(int N_real,double err, double inf,int initialization,int maxit, int tolerance,int decision,string file, string end_file,string w_file,vector<Graph_undirected> & A, vector<int> & u_list,vector< vector<double> > & u, vector< vector<double> > & u_f,vector< vector<double> > & u_old,vector< vector< vector<double> > > & w,vector< vector< vector<double> > > & w_f,vector< vector< vector<double> > > & w_old);
void cycle_over_realizations_undirected(int N_real,double err, double inf,int initialization,int maxit, int tolerance,int decision,string file, string end_file,string w_file,vector<Graph_undirected> & A, vector<int> & u_list,vector< vector<double> > & u, vector< vector<double> > & u_f,vector< vector<double> > & u_old,vector< vector<double> >  & w,vector<vector<double> > & w_f, vector< vector<double> > & w_old); 
void iterate_undirected(bool assortative,int N_real,double err, double inf,int initialization,int maxit, int tolerance,int decision,string file, string end_file,string w_file,vector<Graph_undirected> & A, vector<int> & u_list,vector< vector<double> > & u, vector< vector<double> > & u_f,vector< vector<double> > & u_old);

#endif