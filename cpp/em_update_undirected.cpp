/* -----------------------------------------------------------------------------------------
   -------------------- Main routine for the EM algorithm ----------------------------------
   ----------------------------------------------------------------------------------------- 
  There functions are overloaded: the first one is for the full tensor model, the second one is for the assortative restricted version.
*/

#include "mlg.hpp"
#include <iostream>
#include <vector>

using namespace std;
using namespace boost;
extern int L;
extern int K;
extern int N;

double calculate_Zu(int k,int L,vector< vector<double> > & u_old,vector<int> & u_list,vector< vector< vector<double> > > & w_old){
  // Calculate Z_u, is the same for all u_ik
  int K=(int)w_old.size();
  double Z_u=0.;
  for(int q=0;q<K;q++){ 
      double w_k=0.,Du=0.;
      for(int a=0;a<L;a++)w_k+=w_old[k][q][a];
      for(int i=0;i<(int)u_list.size();i++)Du+=u_old[u_list[i]][q];
      Z_u+=w_k*Du;  
  }
  return Z_u;  
}

// Assortative version
double calculate_Zu(int k,int L,vector< vector<double> > & u_old,vector<int> & u_list,vector<vector<double> >  & w_old){
  // Calculate Z_u, is the same for all u_ik
  double Z_u=0.;
   
  double w_k=0.,Du=0.;
  for(int a=0;a<L;a++)w_k+=w_old[k][a];
  for(int i=0;i<(int)u_list.size();i++)Du+=u_old[u_list[i]][k];
  Z_u+=w_k*Du;  

  return Z_u;  
}


double calculate_Zkq(int q,int k,vector< vector<double> > & u,vector<int> & u_list){
  double Z_kq=0.,Duk=0.,Duq=0.;
  for(int i=0;i<u_list.size();i++){
      Duk+=u[u_list[i]][k];
      Duq+=u[u_list[i]][q];
    }
  Z_kq=Duk*Duq;
  return Z_kq;
}

// Assortative version
double calculate_Zk(int k,vector< vector<double> > & u,vector<int> & u_list){
  double Z_k=0.,Du=0.;
  for(int i=0;i<u_list.size();i++)Du+=u[u_list[i]][k];
  Z_k=Du*Du;
  return Z_k;
}


double calculate_rho_ijkq(int i,int k,Vertex j,int a,vector< vector<double> > & u_old,vector< vector< vector<double> > > & w_old){
  
  double rho_ijka=0.;
  double Zij_a=0.;
  int K=(int)w_old.size();

  for(int m=0;m<K;m++)for(int l=0;l<K;l++)Zij_a+=u_old[i][m]*u_old[j][l]*w_old[m][l][a];
  if(Zij_a!=0.){
    for(int q=0;q<K;q++)rho_ijka+=u_old[j][q]*w_old[k][q][a];
    rho_ijka/=Zij_a;   //  u_old[i][k] will be multiplied only at the end!!!!
  }  // end if Zij_a
  return rho_ijka;
}

// Assortative version
double calculate_rho_ijkq(int i,int k,Vertex j,int a,vector< vector<double> > & u_old,vector< vector<double> > & w_old){
  
  double rho_ijka=0.;
  double Zij_a=0.;
  int K=(int)w_old.size();

  for(int m=0;m<K;m++)Zij_a+=u_old[i][m]*u_old[j][m]*w_old[m][a];
  if(Zij_a!=0.){
    rho_ijka+=u_old[j][k]*w_old[k][a];
    rho_ijka/=Zij_a;   //  u_old[i][k] will be multiplied only at the end!!!!
  }  // end if Zij_a
  return rho_ijka;
}

double calculate_rho_jiqk(int i, int k,Vertex j,int a,vector< vector<double> > & u,vector< vector<double> > & v_old,vector< vector< vector<double> > > & w_old){
  
  double rho_jiqk=0.;
  double Zij_a=0.;
  int K=(int)w_old.size();

  for(int m=0;m<K;m++)for(int l=0;l<K;l++)Zij_a+=u[j][m]*v_old[i][l]*w_old[m][l][a];
  if(Zij_a!=0.){
    for(int q=0;q<K;q++)rho_jiqk+=u[j][q]*w_old[q][k][a];
    rho_jiqk/=Zij_a;   //  u_old[i][k] will be multiplied only at the end!!!!
  }  // end if Zij_a
  return rho_jiqk;
}


double calculate_rho_jiqk(int i, int k,Vertex j,int a,vector< vector<double> > & u,vector< vector<double> > & v_old,vector< vector< double> > & w_old){
  
  double rho_jiqk=0.;
  double Zij_a=0.;
  int K=(int)w_old.size();

  for(int m=0;m<K;m++)Zij_a+=u[j][m]*v_old[i][m]*w_old[m][a];
  if(Zij_a!=0.){
    rho_jiqk+=u[j][k]*w_old[k][a];
    rho_jiqk/=Zij_a;   //  u_old[i][k] will be multiplied only at the end!!!!
  }  // end if Zij_a
  return rho_jiqk;
}

double calculate_rho_w(int i,int a,int q, vector<Graph_undirected> & A,vector< vector<double> > & u,vector< vector< vector<double> > > & w_old){
  
  double rho_w=0.;
  //int K=(int)w_old.size();

  g_edge_iterator_und eit, eend;   // Cycle over out-neighbors of i in layer a
  for (tie(eit, eend) = out_edges(i,A[a]); eit != eend; ++eit) {
    Vertex j=target(*eit,A[a]); 
    double Zij_a=0.;
    for(int m=0;m<K;m++)for(int l=0;l<K;l++)Zij_a+=u[i][m]*u[j][l]*w_old[m][l][a];
    if(Zij_a!=0.)rho_w+=u[j][q]/Zij_a;   //  w_old[k][a] will be multiplied only at the end!!!!
  } //end cycle over out-j  of i

  return rho_w;
}

// Assortative version
double calculate_rho_w(int i,int a,int q,vector<Graph_undirected> & A,vector< vector<double> > & u,vector< vector< double> > & w_old){
  
  double rho_w=0.;
  //int K=(int)w_old.size();

  g_edge_iterator_und eit, eend;   // Cycle over out-neighbors of i in layer a
  for (tie(eit, eend) = out_edges(i,A[a]); eit != eend; ++eit) {
    Vertex j=target(*eit,A[a]); 
    double Zij_a=0.;
    for(int m=0;m<K;m++)Zij_a+=u[i][m]*u[j][m]*w_old[m][a];
    if(Zij_a!=0.)rho_w+=u[j][q]/Zij_a;   //  w_old[k][a] will be multiplied only at the end!!!!
  } //end cycle over out-j  of i

  return rho_w;
}

double calculate_numerator_u(int i,int k, vector<Graph_undirected> & A,vector< vector<double> > & u_old,vector< vector< vector<double> > > & w_old){
  int L=(int)A.size();
  double u_ik=0.;
  for(int a=0;a<L;a++){ // Cycle over layers a -----------
    g_edge_iterator_und eit, eend;   // Cycle over out-neighbors of i in layer a
    for (tie(eit, eend) = out_edges(i,A[a]); eit != eend; ++eit) {
      Vertex j=target(*eit,A[a]);   // OUT-EDGE                         
      u_ik+=calculate_rho_ijkq(i,k,j,a,u_old,w_old);
    } //end cycle over out-j  in layer a
  } // end cycle over layer a - --- --- --- - - - - - - -- - - - - - - - - 
  return u_ik;
}

// Assortative version
double calculate_numerator_u(int i,int k, vector<Graph_undirected> & A,vector< vector<double> > & u_old,vector<vector<double> > & w_old){
  int L=(int)A.size();
  double u_ik=0.;
  for(int a=0;a<L;a++){ // Cycle over layers a -----------
    g_edge_iterator_und eit, eend;   // Cycle over out-neighbors of i in layer a
    for (tie(eit, eend) = out_edges(i,A[a]); eit != eend; ++eit) {
      Vertex j=target(*eit,A[a]);   // OUT-EDGE                         
      u_ik+=calculate_rho_ijkq(i,k,j,a,u_old,w_old);
    } //end cycle over out-j  in layer a
  } // end cycle over layer a - --- --- --- - - - - - - -- - - - - - - - - 
  return u_ik;
}



double calculate_numerator_w( int a,int k,int q,vector<Graph_undirected> & A,vector< vector<double> > & u,vector< vector< vector<double> > > & w_old){
  double w_kqa=0.;
  //int N=(int)num_vertices(A[0]);
  for(int i=0;i<N;i++)w_kqa+=u[i][k]*calculate_rho_w(i, a,q,A, u,w_old);
  return w_kqa ;
}

// Assortative version
double calculate_numerator_w( int a,int k,vector<Graph_undirected> & A,vector< vector<double> > & u,vector< vector<double> > & w_old){
  double w_ka=0.;
  //int N=(int)num_vertices(A[0]);
  for(int i=0;i<N;i++)w_ka+=u[i][k]*calculate_rho_w(i, a,k,A, u,w_old);
  return w_ka ;
}


double calculate_uik(double err_max,vector<Graph_undirected> & A, int i,int k,double Z_u,vector< vector<double> > & u_old,vector< vector< vector<double> > > & w_old){ 
  double u_ik=(u_old[i][k]/Z_u)*calculate_numerator_u(i,k,A,u_old,w_old);
  if(u_ik<err_max)u_ik=0.; 
  return u_ik;
}

// Assortative version
double calculate_uik(double err_max,vector<Graph_undirected> & A, int i,int k,double Z_u,vector< vector<double> > & u_old,vector< vector<double> > & w_old){ 
  double u_ik=(u_old[i][k]/Z_u)*calculate_numerator_u(i,k,A,u_old,w_old);
  if(u_ik<err_max)u_ik=0.; 
  return u_ik;
}


double calculate_wkqa(int a,int k, int q,double err_max,vector<Graph_undirected> & A,double Z_kq,vector< vector<double> > & u,vector< vector< vector<double> > > & w_old ){
  double w_kqaw=(w_old[k][q][a]/Z_kq)*calculate_numerator_w(a,k,q, A, u,w_old);   // final update for w_ka
  if(w_kqaw<err_max)w_kqaw=0.;
  return w_kqaw;
}

// Assortative version
double calculate_wka(int a,int k, double err_max,vector<Graph_undirected> & A,double Z_k,vector< vector<double> > & u,vector< vector<double> > & w_old ){
  double w_kaw=(w_old[k][a]/Z_k)*calculate_numerator_w(a,k, A, u,w_old);   // final update for w_ka
  if(w_kaw<err_max)w_kaw=0.;
  return w_kaw;
}

double update_u(double err_max, vector<Graph_undirected> & A, vector< vector<double> > & u,vector< vector<double> > & u_old,vector< vector< vector<double> > > & w_old,vector<int> & u_list){
   
  int L=(int)A.size();
  int K=(int)w_old.size();
  double dist_u=0.;

  for(int k=0;k<K;k++){
    double Z_u=calculate_Zu(k,L,u_old,u_list,w_old);
    if(Z_u>=0.){     // Cycle over i ----0
      for(int z=0;z<u_list.size();z++)if(u_old[u_list[z]][k]>0.){     // Update only if u_ik >0   
      int i=u_list[z];  
      u[i][k]=calculate_uik(err_max,A,i,k,Z_u,u_old,w_old);
      dist_u=std::max(abs(u[i][k]-u_old[i][k]),dist_u); // calculate max difference btw u_lod and new u membership 
     // u_old[i][k]= u[i][k];  // update so that it will appear correctly in the normalization Zij_a
      } // cycle over i --------------------------------------------------
    }
  }  // end cycle over k --------------------------------------------------
  for(int i=0;i<u.size();i++)for(int k=0;k<K;k++)u_old[i][k]= u[i][k];
  return dist_u;
}

// Assortative version
double update_u(double err_max, vector<Graph_undirected> & A, vector< vector<double> > & u,vector< vector<double> > & u_old,vector< vector<double> > & w_old,vector<int> & u_list){
   
  int L=(int)A.size();
  int K=(int)w_old.size();
  double dist_u=0.;

  for(int k=0;k<K;k++){
    double Z_u=calculate_Zu(k,L,u_old,u_list,w_old);
    if(Z_u>=0.){     // Cycle over i ----0
      for(int z=0;z<u_list.size();z++)if(u_old[u_list[z]][k]>0.){     // Update only if u_ik >0   
      int i=u_list[z];  
      u[i][k]=calculate_uik(err_max,A,i,k,Z_u,u_old,w_old);
      dist_u=std::max(abs(u[i][k]-u_old[i][k]),dist_u); // calculate max difference btw u_lod and new u membership 
      // u_old[i][k]= u[i][k];  // update so that it will appear correctly in the normalization Zij_a
      } // cycle over i --------------------------------------------------
    }
  }  // end cycle over k --------------------------------------------------
  for(int i=0;i<u.size();i++)for(int k=0;k<K;k++)u_old[i][k]= u[i][k];
  return dist_u;
}

double update_w(double err_max,vector<Graph_undirected> & A,vector< vector< vector<double> > > & w,vector< vector<double> > & u, vector< vector< vector<double> > > & w_old,vector<int> & u_list){

  //int L=(int)A.size();
  //int K=(int)w_old.size();
  double dist_w=0.;

  for(int k=0;k<K;k++)for(int q=0;q<K;q++){
    double Z_kq=calculate_Zkq(q,k,u,u_list);
    if(Z_kq!=0.){
      for(int a=0;a<L;a++)if(w_old[k][q][a]>0.){  // Update only if w_ka >0 
        w[k][q][a]=calculate_wkqa(a,k,q,err_max,A,Z_kq,u,w_old);
        if(w[k][q][a]<err_max)w[k][q][a]=0.;
        dist_w=std::max(abs(w[k][q][a]-w_old[k][q][a]),dist_w); // Update max distance old vs new w_kq
      //  w_old[k][q][a]=w[k][q][a];
      } // cycle over a --------------------------------------------------
    } 
  }  // end cycle over k --------------------------------------------------
  for(int k=0;k<K;k++)for(int q=0;q<K;q++)for(int a=0;a<L;a++)w_old[k][q][a]=w[k][q][a];
  return dist_w;
}

// Assortative version
double update_w(double err_max,vector<Graph_undirected> & A,vector< vector<double> > & w,vector< vector<double> > & u, vector< vector<double> > & w_old,vector<int> & u_list){

  //int L=(int)A.size();
  //int K=(int)w_old.size();
  double dist_w=0.;

  for(int k=0;k<K;k++){
    double Z_k=calculate_Zk(k,u,u_list);
    if(Z_k!=0.){
      for(int a=0;a<L;a++)if(w_old[k][a]>0.){  // Update only if w_ka >0 
        w[k][a]=calculate_wka(a,k,err_max,A,Z_k,u,w_old);
        if(w[k][a]<err_max)w[k][a]=0.;
        dist_w=std::max(abs(w[k][a]-w_old[k][a]),dist_w); // Update max distance old vs new w_kq
        // w_old[k][a]=w[k][a];
      } // cycle over a --------------------------------------------------
    } 
  }  // end cycle over k --------------------------------------------------
  for(int k=0;k<K;k++)for(int q=0;q<K;q++)for(int a=0;a<L;a++)w_old[k][a]=w[k][a];

  return dist_w;
}

// ------  Update in sequence  -------------
void update_em_undirected(double err_max, vector< vector<double> > & u,vector< vector<double> > & u_old,vector<vector< vector<double> > > & w,vector< vector< vector<double> > > & w_old, vector<Graph_undirected> & A, vector<int> & u_list,double & d_u,double & d_w){

  d_u=update_u(err_max,A,u,u_old,w_old,u_list);
  d_w=update_w(err_max,A,w,u,w_old,u_list);

} // end update

// Assortative version
void update_em_undirected(double err_max, vector< vector<double> > & u,vector< vector<double> > & u_old,vector<vector<double> > & w,vector< vector<double> > & w_old, vector<Graph_undirected> & A, vector<int> & u_list,double & d_u,double & d_w){

  d_u=update_u(err_max,A,u,u_old,w_old,u_list);
  d_w=update_w(err_max,A,w,u,w_old,u_list);

} // end update


