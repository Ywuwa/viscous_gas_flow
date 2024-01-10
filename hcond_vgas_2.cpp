#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "input.h"
#include "func.h"
#include "Eigen/Dense"
#include <Eigen/Sparse>

using vector_t = std::vector<double>;
using eigen_matrix_t = Eigen::SparseMatrix<double>;
using eigen_triplet_t = Eigen::Triplet<double>;
using eigen_vector_t = Eigen::VectorXd;
//using eigen_precond_t = Eigen::IncompleteLUT<double>;
using eigen_precond_t = Eigen::DiagonalPreconditioner<double>;
//using eigen_precond_t = Eigen::LeastSquareDiagonalPreconditioner<double>;
//using eigen_precond_t = Eigen::IdentityPreconditioner<double>;
using eigen_solver_t = Eigen::BiCGSTAB<eigen_matrix_t, eigen_precond_t>;
//using eigen_solver_t = Eigen::BiCGSTAB<eigen_matrix_t>;

int Sxema_hcond_vgas_2(double mu, int N,int M1, int M2, int Dim,
                      Str_param_gas *str_param_gas, Str_setka *str_setka, P_she *p_she,
                      int *st, double *X, double *Y, int *M0R, int *M0L,
                      double *n_c, double *n_l, double *n_w, int ij, int K, int cas, double *delta_mas)
{
printf("hcond_vgas_2 \n");

double  h1, h2, tau, R, kappa, cv, eps;
tau = str_param_gas->T/N;//p_she->tau;
h1 = str_param_gas->X1/M1;//p_she->h1;
h2 = str_param_gas->X2/M2;//p_she->h2;
R = str_param_gas->R;
kappa = str_param_gas->kappa;
cv = str_param_gas->cv;
eps = str_param_gas->eps;

double EX2 = str_param_gas->EX2;
// assigning functions
double (*cas_u1)(double , double , double);
double (*cas_u2)(double , double , double);
double (*cas_g)(double , double , double);
double (*cas_theta)(double , double , double);

double (*cas_fu1)(double , double , double, double , double , double);
double (*cas_fu2)(double , double , double, double , double , double);
double (*cas_fg)(double , double , double);
double (*cas_ft)(double , double , double, double , double , double , double );

if (cas<2)// for DEBUGGING TEST 1-2 functions u1, u2, g and theta are used
    {
        cas_u1=u1;
        cas_u2=u2;
        cas_g=g;
        cas_theta=theta;

        cas_fu1=fu1; //pravie chasti
        cas_fu2=fu2;
        cas_fg=fg;
        cas_ft=ft;
    }

//declare vectors
eigen_vector_t VEC (3*Dim);
eigen_vector_t T (Dim);
//VEC consists of g, v1, v2: VEC = (G_00, V1_00,V2_00, G_10, V1_10, V2_10, ..., G_M10, V1_M10, V2_M10, G_01, V1_01, V2_01, G_11, V1_11, V2_11,...)
//knots' coord
int m00, mR0, m0R, mL0, m0L;
//initial conditions---------------------------------------
for (int i=0; i<Dim; i++)
{
    VEC[3*i] = cas_g(0,X[i],Y[i]);
    VEC[3*i+1] = cas_u1(0,X[i],Y[i]);
    VEC[3*i+2] = cas_u2(0,X[i],Y[i]);
    T[i] = cas_theta(0,X[i],Y[i]);
}

//const
double th16=tau/(6*h1);
double th26=tau/(6*h2);
double th14=tau/(4*h1);
double th24=tau/(4*h2);
double th12=tau/(2*h1);
double th22=tau/(2*h2);
double t4h1_3=4*tau/(3*h1*h1);
double t4h2_3=4*tau/(3*h2*h2);
double th1_2=tau/(h1*h1);
double th2_2=tau/(h2*h2);

char txt[] = ".txt";
char file[20];
sprintf(file, "%d", ij);
strcat(file,txt);
FILE *OUT;
OUT = fopen(file, "w");
fprintf(OUT, "1 ");

double koef=1.;
double mu_C=0, kap_C=0;
double t_hat=tau;
//function values in the grid konts
double g00=0, gR0=0, g0R=0, gL0=0, g0L=0;
double t00=0, tR0=0, t0R=0, tL0=0, t0L=0;
double v100=0, v1R0=0, v10R=0, v1L0=0, v10L=0, v1RR=0, v1LL=0, v1RL=0, v1LR=0;
double v200=0, v2R0=0, v20R=0, v2L0=0, v20L=0, v2RR=0, v2LL=0, v2RL=0, v2LR=0;
double g2R0=0, g3R0=0, v12R0=0, v13R0=0, v22R0=0, v23R0=0; // extra knots values for border equations
//matrix coefficients for corresponding components
double ag_G00=0, ag_GR0=0, ag_GL0=0, ag_G0R=0, ag_G0L=0, ag_V100=0, ag_V200=0, ag_V1R0=0, ag_V1L0=0, ag_V20R=0, ag_V20L=0;
double av1_GR0=0, av1_GL0=0, av1_V100=0, av1_V1R0=0, av1_V1L0=0, av1_V10R=0, av1_V10L=0, av1_V20R=0, av1_V20L=0;
double av2_G0R=0, av2_G0L=0, av2_V200=0, av2_V2R0=0, av2_V2L0=0, av2_V20R=0, av2_V20L=0, av2_V1R0=0, av2_V1L0=0;
double a_T00=0, a_TR0=0, a_TL0=0, a_T0R=0, a_T0L=0;

//solve=========================================================================================================================================
for (int n=0; n<N; n++) //in the first acceptable version coeficients of V*00 was taken as 0 at all borders; in the second version try not to do it
{
    //declare matrix and right part
eigen_matrix_t matrix (3*Dim, 3*Dim);
eigen_matrix_t matrixt (Dim, Dim);
eigen_vector_t B (3*Dim);
eigen_vector_t Bt (Dim);

std::vector<eigen_triplet_t> triplets;
std::vector<eigen_triplet_t> triplets_t;
// firstly solve: matrix*VEC = B
// then solve: matrixt*T = Bt
  for (int i=0; i<Dim; i++)
  {
    if (mu*exp(-VEC[3*i])>mu_C) mu_C=mu*exp(-VEC[3*i]);
    //if (kappa*exp(-VEC[3*i])>kap_C) kap_C=kappa*exp(-VEC[3*i]);
  }

    R*=0.001;
  /*for (int i=0; i<Dim; i++)
  {
    //VEC[3*i] = cas_g(t_hat,X[i],Y[i]);
    //VEC[3*i+1] = cas_u1(t_hat,X[i],Y[i]);
    //VEC[3*i+2] = cas_u2(t_hat,X[i],Y[i]);
    T[i] = cas_theta(t_hat,X[i],Y[i]);
  }*/ // solve matrix*VEC=B
  for (int i=0; i<Dim; i++)
  {
    //fill matrix
    if (st[i]==0) // inner knots
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*mR0]; g0R=VEC[3*m0R]; gL0=VEC[3*mL0]; g0L=VEC[3*m0L];
        t00=T[m00]; tR0=T[mR0]; tL0=T[mL0]; t0R=T[m0R]; t0L=T[m0L];
        v100=VEC[3*m00+1]; v1R0=VEC[3*mR0+1]; v10R=VEC[3*m0R+1]; v1L0=VEC[3*mL0+1]; v10L=VEC[3*m0L+1]; v1RR=VEC[3*m0R+4]; v1LL=VEC[3*m0L-2]; v1RL=VEC[3*m0L+4]; v1LR=VEC[3*m0R-2];
        v200=VEC[3*m00+2]; v2R0=VEC[3*mR0+2]; v20R=VEC[3*m0R+2]; v2L0=VEC[3*mL0+2]; v20L=VEC[3*m0L+2]; v2RR=VEC[3*m0R+5]; v2LL=VEC[3*m0L-1]; v2RL=VEC[3*m0L+5]; v2LR=VEC[3*m0R-1];

        ag_G00 = 1.;
        ag_GL0 = -th14*(v100+v1L0);
        ag_GR0 = th14*(v100+v1R0);
        ag_G0L = -th24*(v200+v20L);
        ag_G0R = th24*(v200+v20R);
        ag_V1L0 = -th12;
        ag_V20L = -th22;
        ag_V1R0 = th12;
        ag_V20R = th22;

        av1_V100 = 1. + koef*2*mu_C*t4h1_3 + koef*2*mu_C*th2_2;
        av1_V1R0 = th16*(v100 + v1R0) - koef*mu_C*t4h1_3;
        av1_V1L0 = -th16*(v100 + v1L0) - koef*mu_C*t4h1_3;
        av1_V10R = th24*(v200 + v20R) - koef*mu_C*th2_2;
        av1_V10L = -th24*(v200 + v20L) - koef*mu_C*th2_2;
        av1_GL0 = -R*th12*t00;
        av1_GR0 = R*th12*t00;

        av2_V200 = 1. + koef*2*mu_C*th1_2 + koef*2*mu_C*t4h2_3;
        av2_V20R = th26*(v200 + v20R) - koef*mu_C*t4h2_3;
        av2_V20L = -th26*(v200 + v20L) - koef*mu_C*t4h2_3;
        av2_V2R0 = th14*(v100 + v1R0) - koef*mu_C*th1_2;
        av2_V2L0 = -th14*(v100 + v1L0) - koef*mu_C*th1_2;
        av2_G0L = -R*th22*t00;
        av2_G0R = R*th22*t00;

        B[3*m00] = g00 + g00*(th14*(v1R0-v1L0) + th24*(v20R-v20L)) + tau*cas_fg(t_hat, X[i],Y[i]);

        B[3*m00+1] = v100 - th12*R*(tR0-tL0) + th24*v100*(v20R-v20L)
                    + koef*(mu*exp(-g00)-mu_C)*(t4h1_3*(v1R0-2*v100+v1L0) + th2_2*(v10R-2*v100+v10L))
                    + koef*tau*mu*exp(-g00)*(v2RR-v2RL-v2LR+v2LL)/(12*h1*h2)
                    + tau*cas_fu1(t_hat, X[i],Y[i], R,mu, koef);

        B[3*m00+2] = v200 - th22*R*(t0R-t0L) + th14*v200*(v1R0-v1L0)
                    + koef*(mu*exp(-g00)-mu_C)*(t4h2_3*(v20R-2*v200+v20L) + th1_2*(v2R0-2*v200+v2L0))
                    + koef*tau*mu*exp(-g00)*(v1RR-v1RL-v1LR+v1LL)/(12*h1*h2)
                    + tau*cas_fu2(t_hat, X[i],Y[i], R,mu, koef);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*mL0,ag_GL0);
        triplets.emplace_back(3*m00,3*mR0,ag_GR0);
        triplets.emplace_back(3*m00,3*m0L,ag_G0L);
        triplets.emplace_back(3*m00,3*m0R,ag_G0R);
        triplets.emplace_back(3*m00,3*mL0+1,ag_V1L0);
        triplets.emplace_back(3*m00,3*mR0+1,ag_V1R0);
        triplets.emplace_back(3*m00,3*m0L+2,ag_V20L);
        triplets.emplace_back(3*m00,3*m0R+2,ag_V20R);
        //B[3*m00] = cas_g(t_hat,X[i],Y[i]);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+1,3*mR0+1,av1_V1R0);
        triplets.emplace_back(3*m00+1,3*mL0+1,av1_V1L0);
        triplets.emplace_back(3*m00+1,3*m0R+1,av1_V10R);
        triplets.emplace_back(3*m00+1,3*m0L+1,av1_V10L);
        triplets.emplace_back(3*m00+1,3*mL0,av1_GL0);
        triplets.emplace_back(3*m00+1,3*mR0,av1_GR0);

        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);
        triplets.emplace_back(3*m00+2,3*mR0+2,av2_V2R0);
        triplets.emplace_back(3*m00+2,3*mL0+2,av2_V2L0);
        triplets.emplace_back(3*m00+2,3*m0R+2,av2_V20R);
        triplets.emplace_back(3*m00+2,3*m0L+2,av2_V20L);
        triplets.emplace_back(3*m00+2,3*m0L,av2_G0L);
        triplets.emplace_back(3*m00+2,3*m0R,av2_G0R);
        /*av1_V100=1;
        av2_V200=1;
        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = cas_u2(t_hat,X[i],Y[i]);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);*/

    }
    else if (st[i]==1) // left border
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*mR0];  g2R0=VEC[3*mR0+3]; g3R0=VEC[3*mR0+6];
        v100=VEC[3*m00+1]; v1R0=VEC[3*mR0+1]; v12R0=VEC[3*mR0+4]; v13R0=VEC[3*mR0+7];
        v100=0;
        ag_G00 = 1. - th12*v100;
        ag_GR0 = th12*v1R0;
        ag_V100=-2*th12;
        ag_V1R0=2*th12;

        av1_V100=1;
        av2_V200=1;

        B[3*m00] = g00 + th12*g00*(v1R0-v100) + th12*((g2R0*v12R0-2*gR0*v1R0+g00*v100) - 0.5*(g3R0*v13R0-2*g2R0*v12R0+gR0*v1R0))
                                            + th12*(2-g00)*(v12R0-2*v1R0+v100 - 0.5*(v13R0-2*v12R0+v1R0))
                                            + tau*cas_fg(t_hat, X[i],Y[i]);
        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = cas_u2(t_hat,X[i],Y[i]);
        //ag_G00 = 1;
        //B[3*m00] = cas_g(t_hat,X[i],Y[i]);
        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*mR0,ag_GR0);
        triplets.emplace_back(3*m00,3*m00+1,ag_V100);
        triplets.emplace_back(3*m00,3*mR0+1,ag_V1R0);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);

    }
    else if (st[i]==2) // right border
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*mL0]; g2R0=VEC[3*mL0-3]; g3R0=VEC[3*mL0-6];
        v100=VEC[3*m00+1]; v1R0=VEC[3*mL0+1]; v12R0=VEC[3*mL0-2]; v13R0=VEC[3*mL0-5]; //here R means L in order to save memory
        v100=0;
        ag_G00 = 1. + th12*v100;
        ag_GL0 = -th12*v1R0;
        ag_V100=2*th12;
        ag_V1L0=-2*th12;

        B[3*m00] = g00 + th12*g00*(v100-v1R0) - th12*((g00*v100-2*gR0*v1R0+g2R0*v12R0) - 0.5*(gR0*v1R0-2*g2R0*v12R0+g3R0*v13R0))
                                            - th12*(2-g00)*(v100-2*v1R0+v12R0 - 0.5*(v1R0-2*v12R0+v13R0))
                                            + tau*cas_fg(t_hat, X[i],Y[i]);
        av1_V100=1;
        av2_V200=1;
        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = cas_u2(t_hat,X[i],Y[i]);
        //ag_G00 = 1;
        //B[3*m00] = cas_g(t_hat,X[i],Y[i]);
        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*mL0,ag_GL0);
        triplets.emplace_back(3*m00,3*m00+1,ag_V100);
        triplets.emplace_back(3*m00,3*mL0+1,ag_V1L0);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);

    }
    else if (st[i]==3) // down border
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*m0R];  g2R0=VEC[3*m0R+3*(m0R-m00)]; g3R0=VEC[3*m0R+6*(m0R-m00)]; //here R0 means 0R in order to save memory
        v200=VEC[3*m00+2]; v2R0=VEC[3*m0R+2]; v22R0=VEC[3*m0R+3*(m0R-m00)+2]; v23R0=VEC[3*m0R+6*(m0R-m00)+2];
        v200=0;
        ag_G00 = 1. - th22*v200;
        ag_G0R = th22*v2R0;
        ag_V200=-2*th22;
        ag_V20R=2*th22;

        B[3*m00] = g00 + th22*g00*(v2R0-v200) + th22*((g2R0*v22R0-2*gR0*v2R0+g00*v200) - 0.5*(g3R0*v23R0-2*g2R0*v22R0+gR0*v2R0))
                                            + th22*(2-g00)*(v22R0-2*v2R0+v200 - 0.5*(v23R0-2*v22R0+v2R0))
                                            + tau*cas_fg(t_hat, X[i],Y[i]);
        av1_V100=1;
        av2_V200=1;
        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = cas_u2(t_hat,X[i],Y[i]);
        //ag_G00 = 1;
        //B[3*m00] = cas_g(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*m0R,ag_G0R);
        triplets.emplace_back(3*m00,3*m00+2,ag_V200);
        triplets.emplace_back(3*m00,3*m0R+2,ag_V20R);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);
    }
    else if (st[i]==4) // up border
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*m0L];  g2R0=VEC[3*m0L-3*(m00-m0L)]; g3R0=VEC[3*m0L-6*(m00-m0L)]; //here R0 means 0R in order to save memory
        v200=VEC[3*m00+2]; v2R0=VEC[3*m0L+2]; v22R0=VEC[3*m0L-3*(m00-m0L)+2]; v23R0=VEC[3*m0L-6*(m00-m0L)+2]; // also here R means L in order to save memory
        v200=0;
        ag_G00 = 1. + th22*v200;
        ag_G0L = -th22*v2R0;
        ag_V200=2*th22;
        ag_V20L=-2*th22;

        B[3*m00] = g00 + th22*g00*(v200-v2R0) - th22*((g2R0*v22R0-2*gR0*v2R0+g00*v200) - 0.5*(g3R0*v23R0-2*g2R0*v22R0+gR0*v2R0))
                                            - th22*(2-g00)*(v22R0-2*v2R0+v200 - 0.5*(v23R0-2*v22R0+v2R0))
                                            + tau*cas_fg(t_hat, X[i],Y[i]);
        av1_V100=1;
        av2_V200=1;
        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = cas_u2(t_hat,X[i],Y[i]);
        //ag_G00 = 1;
        //B[3*m00] = cas_g(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*m0L,ag_G0L);
        triplets.emplace_back(3*m00,3*m00+2,ag_V200);
        triplets.emplace_back(3*m00,3*m0L+2,ag_V20L);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);
    }
    else if (st[i]==5) // left down corner (test)
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00];

        ag_G00 = 1;
        av1_V100=1;
        av2_V200=1;

        B[3*m00] = g00 + tau*cas_fg(t_hat, X[i],Y[i]);
        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = cas_u2(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);

        //B[3*m00] = cas_g(t_hat,X[i],Y[i]);
    }
    else if (st[i]==6) // right down corner (test)
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00];

        ag_G00 = 1;
        av1_V100=1;
        av2_V200=1;

        B[3*m00] = g00 + tau*cas_fg(t_hat, X[i],Y[i]);
        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = cas_u2(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);

        //B[3*m00] = cas_g(t_hat,X[i],Y[i]);
    }
    else if (st[i]==7) // left up corner (test)
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00];

        ag_G00 = 1;
        av1_V100=1;
        av2_V200=1;

        B[3*m00] = g00 + tau*cas_fg(t_hat, X[i],Y[i]);
        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = cas_u2(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);

        //B[3*m00] = cas_g(t_hat,X[i],Y[i]);
    }
    else if (st[i]==8) // right up corner (test)
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00];

        ag_G00 = 1;
        av1_V100=1;
        av2_V200=1;

        B[3*m00] = g00 + tau*cas_fg(t_hat, X[i],Y[i]);

        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = cas_u2(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);

        //B[3*m00] = cas_g(t_hat,X[i],Y[i]);
    }
    else if (st[i]==11) // inner left up corner
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*mR0]; g0R=VEC[3*m0R]; gL0=VEC[3*mL0]; g0L=VEC[3*m0L];
        t00=T[m00]; tR0=T[mR0]; tL0=T[mL0]; t0R=T[m0R]; t0L=T[m0L];
        v100=VEC[3*m00+1]; v1R0=VEC[3*mR0+1]; v10R=VEC[3*m0R+1]; v1L0=VEC[3*mL0+1]; v10L=VEC[3*m0L+1]; v1RR=VEC[3*m0R+4]; v1LL=VEC[3*m0L-2]; v1RL=VEC[3*m0L+4]; v1LR=0;
        v200=VEC[3*m00+2]; v2R0=VEC[3*mR0+2]; v20R=VEC[3*m0R+2]; v2L0=VEC[3*mL0+2]; v20L=VEC[3*m0L+2]; v2RR=VEC[3*m0R+5]; v2LL=VEC[3*m0L-1]; v2RL=VEC[3*m0L+5]; v2LR=0;

        ag_G00 = 1.;
        ag_GL0 = -th14*(v100+v1L0);
        ag_GR0 = th14*(v100+v1R0);
        ag_G0L = -th24*(v200+v20L);
        ag_G0R = th24*(v200+v20R);
        ag_V1L0 = -th12;
        ag_V20L = -th22;
        ag_V1R0 = th12;
        ag_V20R = th22;

        av1_V100 = 1. + koef*2*mu_C*t4h1_3 + koef*2*mu_C*th2_2;
        av1_V1R0 = th16*(v100 + v1R0) - koef*mu_C*t4h1_3;
        av1_V1L0 = -th16*(v100 + v1L0) - koef*mu_C*t4h1_3;
        av1_V10R = th24*(v200 + v20R) - koef*mu_C*th2_2;
        av1_V10L = -th24*(v200 + v20L) - koef*mu_C*th2_2;
        av1_GL0 = -R*th12*t00;
        av1_GR0 = R*th12*t00;

        av2_V200 = 1. + koef*2*mu_C*th1_2 + koef*2*mu_C*t4h2_3;
        av2_V20R = th26*(v200 + v20R) - koef*mu_C*t4h2_3;
        av2_V20L = -th26*(v200 + v20L) - koef*mu_C*t4h2_3;
        av2_V2R0 = th14*(v100 + v1R0) - koef*mu_C*th1_2;
        av2_V2L0 = -th14*(v100 + v1L0) - koef*mu_C*th1_2;
        av2_G0L = -R*th22*t00;
        av2_G0R = R*th22*t00;

        B[3*m00] = g00 + g00*(th14*(v1R0-v1L0) + th24*(v20R-v20L)) + tau*cas_fg(t_hat, X[i],Y[i]);

        B[3*m00+1] = v100 - th12*R*(tR0-tL0) + th24*v100*(v20R-v20L)
                    + koef*(mu*exp(-g00)-mu_C)*(t4h1_3*(v1R0-2*v100+v1L0) + th2_2*(v10R-2*v100+v10L))
                    + koef*tau*mu*exp(-g00)*(v2RR-v2RL-v2LR+v2LL)/(12*h1*h2)
                    + tau*cas_fu1(t_hat, X[i],Y[i], R,mu, koef);

        B[3*m00+2] = v200 - th22*R*(t0R-t0L) + th14*v200*(v1R0-v1L0)
                    + koef*(mu*exp(-g00)-mu_C)*(t4h2_3*(v20R-2*v200+v20L) + th1_2*(v2R0-2*v200+v2L0))
                    + koef*tau*mu*exp(-g00)*(v1RR-v1RL-v1LR+v1LL)/(12*h1*h2)
                    + tau*cas_fu2(t_hat, X[i],Y[i], R,mu, koef);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*mL0,ag_GL0);
        triplets.emplace_back(3*m00,3*mR0,ag_GR0);
        triplets.emplace_back(3*m00,3*m0L,ag_G0L);
        triplets.emplace_back(3*m00,3*m0R,ag_G0R);
        triplets.emplace_back(3*m00,3*mL0+1,ag_V1L0);
        triplets.emplace_back(3*m00,3*mR0+1,ag_V1R0);
        triplets.emplace_back(3*m00,3*m0L+2,ag_V20L);
        triplets.emplace_back(3*m00,3*m0R+2,ag_V20R);
        //B[3*m00] = cas_g(t_hat,X[i],Y[i]);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+1,3*mR0+1,av1_V1R0);
        triplets.emplace_back(3*m00+1,3*mL0+1,av1_V1L0);
        triplets.emplace_back(3*m00+1,3*m0R+1,av1_V10R);
        triplets.emplace_back(3*m00+1,3*m0L+1,av1_V10L);
        triplets.emplace_back(3*m00+1,3*mL0,av1_GL0);
        triplets.emplace_back(3*m00+1,3*mR0,av1_GR0);

        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);
        triplets.emplace_back(3*m00+2,3*mR0+2,av2_V2R0);
        triplets.emplace_back(3*m00+2,3*mL0+2,av2_V2L0);
        triplets.emplace_back(3*m00+2,3*m0R+2,av2_V20R);
        triplets.emplace_back(3*m00+2,3*m0L+2,av2_V20L);
        triplets.emplace_back(3*m00+2,3*m0L,av2_G0L);
        triplets.emplace_back(3*m00+2,3*m0R,av2_G0R);
        /*av1_V100=1;
        av2_V200=1;
        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = cas_u2(t_hat,X[i],Y[i]);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);*/

    }
  }
  //build matrix, solve with eigen
  matrix.setFromTriplets (triplets.begin(), triplets.end());
  eigen_solver_t solver(matrix);

  if (solver.info() != Eigen::Success)
  {
      fprintf(OUT, "Can not build preconditioner");
      return -1;
  }
  eigen_vector_t ress = solver.solveWithGuess(B, VEC);
  if (solver.info() != Eigen::Success)
  {
      fprintf(OUT, "Failed to solve the system with Eigen, n=%d", n);
      n_c[4*ij]=0; n_c[4*ij+1]=0; n_c[4*ij+2]=0; n_c[4*ij+3]=0;
      n_l[4*ij]=0; n_l[4*ij+1]=0; n_l[4*ij+2]=0; n_l[4*ij+3]=0;
      return 0;
  }
  VEC=ress;
  /*for (int i=0; i<Dim; i++)
  {
    VEC[3*i] = cas_g(t_hat,X[i],Y[i]);
    VEC[3*i+1] = cas_u1(t_hat,X[i],Y[i]);
    VEC[3*i+2] = cas_u2(t_hat,X[i],Y[i]);
    //T[i] = cas_theta(t_hat,X[i],Y[i]);
  }*/
  R*=1000;
  for (int i=0; i<Dim; i++)
  {
    if (kappa*exp(-VEC[3*i])>kap_C) kap_C=kappa*exp(-VEC[3*i]);
  }
  // solve matrixt*T=Bt
  for (int i=0; i<Dim; i++)
  {
    //T[i] = cas_theta(t_hat,X[i],Y[i]);
    // fill matrixt
    if (st[i]==0) // inner knots
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        //m0R=i+M1+1;//M0R[i];
        //m0L=i-M1-1;//M0L[i];
        m0R=M0R[i];
        m0L=M0L[i];

        g00=VEC[3*m00];
        t00=T[m00]; tR0=T[mR0]; tL0=T[mL0]; t0R=T[m0R]; t0L=T[m0L];
        v100=VEC[3*m00+1]; v1R0=VEC[3*mR0+1]; v10R=VEC[3*m0R+1]; v1L0=VEC[3*mL0+1]; v10L=VEC[3*m0L+1]; // try to take 3*i+3*r+1 instead of 3*i+3*r+4 & 3*i-3*l+1 instead of 3*i-3*l-2
        v200=VEC[3*m00+2]; v2R0=VEC[3*mR0+2]; v20R=VEC[3*m0R+2]; v2L0=VEC[3*mL0+2]; v20L=VEC[3*m0L+2];

        a_T00 = cv + 2*kap_C*th1_2 + 2*kap_C*th2_2;
        a_TR0 = cv*th14*(v100+v1R0) - kap_C*th1_2;
        a_TL0 = -cv*th14*(v100+v1L0) - kap_C*th1_2;
        a_T0R = cv*th24*(v200+v20R) - kap_C*th2_2;
        a_T0L = -cv*th24*(v200+v20L) - kap_C*th2_2;

        Bt[m00] = cv*t00 + cv*t00*(th14*(v1R0-v1L0)+th24*(v20R-v20L)) + (kappa*exp(-g00)-kap_C)*(th1_2*(tR0-2*t00+tL0)+th2_2*(t0R-2*t00+t0L))
                                                       - R*t00*(th12*(v1R0-v1L0)+th22*(v20R-v20L))
                                                       - (1./6)*mu*exp(-g00)*(th1_2*(v1R0-v1L0)*(v1R0-v1L0)+th2_2*(v20R-v20L)*(v20R-v20L))
                                                       - (tau/(3*h1*h2))*mu*exp(-g00)*(v1R0-v1L0)*(v20R-v20L)
                                                       + 0.5*mu*exp(-g00) * (th2_2*(v20R-v20L)*(v20R-v20L) + th1_2*(v1R0-v1L0)*(v1R0-v1L0))
                                                       + 0.25*tau*mu*exp(-g00) * ((v2R0-v2L0)/h1+(v10R-v10L)/h2) * ((v2R0-v2L0)/h1+(v10R-v10L)/h2)
                                                       + tau*cas_ft(t_hat,X[i],Y[i],R,mu,cv,kappa);
        triplets_t.emplace_back(m00,m00,a_T00);
        triplets_t.emplace_back(m00,mR0,a_TR0);
        triplets_t.emplace_back(m00,mL0,a_TL0);
        triplets_t.emplace_back(m00,m0R,a_T0R);
        triplets_t.emplace_back(m00,m0L,a_T0L);
    }
    else if (st[i]==1) // left border
    {
        m00=i;
        a_T00=1;
        Bt[m00] = cas_theta(t_hat,X[i],Y[i]);
        triplets_t.emplace_back(m00,m00,a_T00); // matrixt is Dim*Dim, T is Dim, i=0,..Dim-1; first coord is matrixt's string, second coord is matrixt's row
    }
    else if (st[i]==2) // right border
    {
        m00=i;
        a_T00=1;
        Bt[m00] = cas_theta(t_hat,X[i],Y[i]);
        triplets_t.emplace_back(m00,m00,a_T00);
    }
    else if (st[i]==3) // down border
    {
        m00=i;
        a_T00=1;
        Bt[m00] = cas_theta(t_hat,X[i],Y[i]);
        triplets_t.emplace_back(m00,m00,a_T00);
    }
    else if (st[i]==4) // up border
    {
        m00=i;
        a_T00=1;
        Bt[m00] = cas_theta(t_hat,X[i],Y[i]);
        triplets_t.emplace_back(m00,m00,a_T00);
    }
    else if (st[i]==5) // left down corner (test)
    {
        m00=i;
        a_T00=1;
        Bt[m00] = cas_theta(t_hat,X[i],Y[i]);
        triplets_t.emplace_back(m00,m00,a_T00);
    }
    else if (st[i]==6) // right down corner (test)
    {
        m00=i;
        a_T00=1;
        Bt[i] = cas_theta(t_hat,X[i],Y[i]);
        triplets_t.emplace_back(m00,m00,a_T00);
    }
    else if (st[i]==7) // left up corner (test)
    {
        m00=i;
        a_T00=1;
        Bt[m00] = cas_theta(t_hat,X[i],Y[i]);
        triplets_t.emplace_back(m00,m00,a_T00);
    }
    else if (st[i]==8) // right up corner (test)
    {
        m00=i;
        a_T00=1;
        Bt[m00] = cas_theta(t_hat,X[i],Y[i]);
        triplets_t.emplace_back(m00,m00,a_T00);
    }
    else if (st[i]==11) // inner left up corner
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00];
        t00=T[m00]; tR0=T[mR0]; t0R=T[m0R]; tL0=T[mL0]; t0L=T[m0L];
        v100=VEC[3*m00+1]; v1R0=VEC[3*mR0+1]; v10R=VEC[3*m0R+1]; v1L0=VEC[3*mL0+1]; v10L=VEC[3*m0L+1]; // try to take 3*i+3*r+1 instead of 3*i+3*r+4 & 3*i-3*l+1 instead of 3*i-3*l-2
        v200=VEC[3*m00+2]; v2R0=VEC[3*mR0+2]; v20R=VEC[3*m0R+2]; v2L0=VEC[3*mL0+2]; v20L=VEC[3*m0L+2];

        a_T00 = cv + 2*kap_C*th1_2 + 2*kap_C*th2_2;
        a_TR0 = cv*th14*(v100+v1R0) - kap_C*th1_2;
        a_TL0 = -cv*th14*(v100+v1L0) - kap_C*th1_2;
        a_T0R = cv*th24*(v200+v20R) - kap_C*th2_2;
        a_T0L = -cv*th24*(v200+v20L) - kap_C*th2_2;

        Bt[m00] = cv*t00 + cv*t00*(th14*(v1R0-v1L0)+th24*(v20R-v20L)) + (kappa*exp(-g00)-kap_C)*(th1_2*(tR0-2*t00+tL0)+th2_2*(t0R-2*t00+t0L))
                                                       - R*t00*(th12*(v1R0-v1L0)+th22*(v20R-v20L))
                                                       - (1./6)*mu*exp(-g00)*(th1_2*(v1R0-v1L0)*(v1R0-v1L0)+th2_2*(v20R-v20L)*(v20R-v20L))
                                                       - (tau/(3*h1*h2))*mu*exp(-g00)*(v1R0-v1L0)*(v20R-v20L)
                                                       + 0.5*mu*exp(-g00) * (th2_2*(v20R-v20L)*(v20R-v20L) + th1_2*(v1R0-v1L0)*(v1R0-v1L0))
                                                       + 0.25*tau*mu*exp(-g00) * ((v2R0-v2L0)/h1+(v10R-v10L)/h2) * ((v2R0-v2L0)/h1+(v10R-v10L)/h2)
                                                       + tau*cas_ft(t_hat,X[i],Y[i],R,mu,cv,kappa);
        triplets_t.emplace_back(m00,m00,a_T00);
        triplets_t.emplace_back(m00,mR0,a_TR0);
        triplets_t.emplace_back(m00,mL0,a_TL0);
        triplets_t.emplace_back(m00,m0R,a_T0R);
        triplets_t.emplace_back(m00,m0L,a_T0L);
    }
  }
  //build matrixt, solve with eigen
  matrixt.setFromTriplets (triplets_t.begin(), triplets_t.end());
  eigen_solver_t solvert(matrixt);

  if (solvert.info() != Eigen::Success)
  {
      fprintf(OUT, "Can not build preconditioner");
      return -1;
  }
  eigen_vector_t res (Dim);
  res = solvert.solveWithGuess(Bt, T);
  if (solvert.info() != Eigen::Success)
  {
      fprintf(OUT, "Failed to solve the system with Eigen, n=%d", n);
      n_c[4*ij]=0; n_c[4*ij+1]=0; n_c[4*ij+2]=0; n_c[4*ij+3]=0;
      n_l[4*ij]=0; n_l[4*ij+1]=0; n_l[4*ij+2]=0; n_l[4*ij+3]=0;
      return 0;
  }
  T=res;
  t_hat+=tau; mu_C=0; kap_C=0;
}
// solve end =========================================================================================================
fprintf(OUT, "1");
for (int i=0; i<Dim; i++)
  {
    VEC[3*i] -= cas_g(str_param_gas->T,X[i],Y[i]);
    VEC[3*i+1] -= cas_u1(str_param_gas->T,X[i],Y[i]);
    VEC[3*i+2] -= cas_u2(str_param_gas->T,X[i],Y[i]);
    T[i] -= cas_theta(str_param_gas->T,X[i],Y[i]);
  }
//normi
double maxg=0, maxu1=0, maxu2=0, maxt=0;
double sumu1=0, sumu2=0, sumg=0, sumt=0;
for (int i=0; i<Dim; i++)
{
    if (maxg<fabs(VEC[3*i])) maxg=fabs(VEC[3*i]);
    if (maxu1<fabs(VEC[3*i+1])) maxu1=fabs(VEC[3*i+1]);
    if (maxu2<fabs(VEC[3*i+2])) maxu2=fabs(VEC[3*i+2]);
    if (maxt<fabs(T[i])) {maxt=fabs(T[i]);}// fprintf(OUT, " %d/%d ", i,st[i]);}
    sumg+=(VEC[3*i])*(VEC[3*i]);
    sumu1+=(VEC[3*i+1])*(VEC[3*i+1]);
    sumu2+=(VEC[3*i+2])*(VEC[3*i+2]);
    sumt+=(T[i])*(T[i]);
}
n_c[4*ij]=maxu1; n_c[4*ij+1]=maxu2; n_c[4*ij+2]=maxg; n_c[4*ij+3]=maxt;
n_l[4*ij]=sqrt(h1*h2*sumu1); n_l[4*ij+1]=sqrt(h1*h2*sumu2); n_l[4*ij+2]=sqrt(h1*h2*sumg); n_l[4*ij+3]=sqrt(h1*h2*sumt);
fclose(OUT);
return 0;
}
//====================================================================================================================
int Sxema_hcond_vgas_2_test(double mu, int N,int M1, int M2, int Dim,
                      Str_param_gas *str_param_gas, Str_setka *str_setka, P_she *p_she,
                      int *st, double *X, double *Y, int *M0R, int *M0L,
                      double *n_c, double *n_l, double *n_w, int ij, int K, int cas, double *delta_mas)
{
    printf("hcond_vgas_2_test \n");

double  h1, h2, tau, R, kappa, cv, eps;
tau = str_param_gas->T/N;//p_she->tau;
h1 = str_param_gas->X1/M1;//p_she->h1;
h2 = str_param_gas->X2/M2;//p_she->h2;
R = str_param_gas->R;
kappa = str_param_gas->kappa;
cv = str_param_gas->cv;
eps = str_param_gas->eps;

double EX2 = str_param_gas->EX2;
// assigning functions
double (*cas_u1)(double , double , double);
double (*cas_u2)(double , double , double);
double (*cas_g)(double , double , double);
double (*cas_theta)(double , double , double);

double (*cas_fu1)(double , double , double, double , double , double );
double (*cas_fu2)(double , double , double, double , double , double );
double (*cas_fg)(double , double , double);
double (*cas_ft)(double , double , double, double , double , double , double );

if (cas<2)// for DEBUGGING TEST 1-2 functions u1, u2, g and theta are used
    {
        cas_u1=u1;
        cas_u2=u2;
        cas_g=g;
        cas_theta=theta;

        cas_fu1=fu12; //pravie chasti
        cas_fu2=fu22;
        cas_fg=fg;
        cas_ft=ft;
    }

//declare vectors
eigen_vector_t VEC (3*Dim);
//VEC consists of g, v1, v2: VEC = (G_00, V1_00,V2_00, G_10, V1_10, V2_10, ..., G_M10, V1_M10, V2_M10, G_01, V1_01, V2_01, G_11, V1_11, V2_11,...)
//knots' coord
int m00, mR0, m0R, mL0, m0L;
//initial conditions---------------------------------------
for (int i=0; i<Dim; i++)
{
    VEC[3*i] = cas_g(0,X[i],Y[i]);
    VEC[3*i+1] = cas_u1(0,X[i],Y[i]);
    VEC[3*i+2] = cas_u2(0,X[i],Y[i]);
}

//const
double th16=tau/(6*h1);
double th26=tau/(6*h2);
double th14=tau/(4*h1);
double th24=tau/(4*h2);
double th12=tau/(2*h1);
double th22=tau/(2*h2);
double t4h1_3=4*tau/(3*h1*h1);
double t4h2_3=4*tau/(3*h2*h2);
double th1_2=tau/(h1*h1);
double th2_2=tau/(h2*h2);

char txt[] = ".txt";
char file[20];
sprintf(file, "%d", ij);
strcat(file,txt);
FILE *OUT;
OUT = fopen(file, "w");
fprintf(OUT, "1 ");

double koef=1.;
double mu_C=0, kap_C=0;
double t_hat=tau;
//function values in the grid konts
double g00=0, gR0=0, g0R=0, gL0=0, g0L=0;
double t00=0, tR0=0, t0R=0, tL0=0, t0L=0;
double v100=0, v1R0=0, v10R=0, v1L0=0, v10L=0, v1RR=0, v1LL=0, v1RL=0, v1LR=0;
double v200=0, v2R0=0, v20R=0, v2L0=0, v20L=0, v2RR=0, v2LL=0, v2RL=0, v2LR=0;
double g2R0=0, g3R0=0, v12R0=0, v13R0=0, v22R0=0, v23R0=0; // extra knots values for border equations
//matrix coefficients for corresponding components
double ag_G00=0, ag_GR0=0, ag_GL0=0, ag_G0R=0, ag_G0L=0, ag_V100=0, ag_V200=0, ag_V1R0=0, ag_V1L0=0, ag_V20R=0, ag_V20L=0;
double av1_GR0=0, av1_GL0=0, av1_V100=0, av1_V1R0=0, av1_V1L0=0, av1_V10R=0, av1_V10L=0, av1_V20R=0, av1_V20L=0;
double av2_G0R=0, av2_G0L=0, av2_V200=0, av2_V2R0=0, av2_V2L0=0, av2_V20R=0, av2_V20L=0, av2_V1R0=0, av2_V1L0=0;
double a_T00=0, a_TR0=0, a_TL0=0, a_T0R=0, a_T0L=0;

//solve=========================================================================================================================================
for (int n=0; n<N; n++)
{
    //declare matrix and right part
eigen_matrix_t matrix (3*Dim, 3*Dim);
eigen_vector_t B (3*Dim);

std::vector<eigen_triplet_t> triplets;
// firstly solve: matrix*VEC = B

  /*for (int i=0; i<Dim; i++)
  {
    //VEC[3*i] = cas_g(t_hat,X[i],Y[i]);
    //VEC[3*i+1] = cas_u1(t_hat,X[i],Y[i]);
    //VEC[3*i+2] = cas_u2(t_hat,X[i],Y[i]);
  }*/
  for (int i=0; i<Dim; i++)
  {
    if (mu*exp(-VEC[3*i])>mu_C) mu_C=mu*exp(-VEC[3*i]);
  }
  for (int i=0; i<Dim; i++)
  {
    //fill matrix
    if (st[i]==0) // inner knots
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*mR0]; g0R=VEC[3*m0R]; gL0=VEC[3*mL0]; g0L=VEC[3*m0L];
        v100=VEC[3*m00+1]; v1R0=VEC[3*mR0+1]; v10R=VEC[3*m0R+1]; v1L0=VEC[3*mL0+1]; v10L=VEC[3*m0L+1]; v1RR=VEC[3*m0R+4]; v1LL=VEC[3*m0L-2]; v1RL=VEC[3*m0L+4]; v1LR=VEC[3*m0R-2];
        v200=VEC[3*m00+2]; v2R0=VEC[3*mR0+2]; v20R=VEC[3*m0R+2]; v2L0=VEC[3*mL0+2]; v20L=VEC[3*m0L+2]; v2RR=VEC[3*m0R+5]; v2LL=VEC[3*m0L-1]; v2RL=VEC[3*m0L+5]; v2LR=VEC[3*m0R-1];

        ag_G00 = 1.;
        ag_GL0 = -th14*(v100+v1L0);
        ag_GR0 = th14*(v100+v1R0);
        ag_G0L = -th24*(v200+v20L);
        ag_G0R = th24*(v200+v20R);
        ag_V1L0 = -th12;
        ag_V20L = -th22;
        ag_V1R0 = th12;
        ag_V20R = th22;

        av1_V100 = 1. + koef*2*mu_C*t4h1_3 + koef*2*mu_C*th2_2;
        av1_V1R0 = th16*(v100 + v1R0) - koef*mu_C*t4h1_3;
        av1_V1L0 = -th16*(v100 + v1L0) - koef*mu_C*t4h1_3;
        av1_V10R = th24*(v200 + v20R) - koef*mu_C*th2_2;
        av1_V10L = -th24*(v200 + v20L) - koef*mu_C*th2_2;
        av1_GL0 = -10*th12;
        av1_GR0 = 10*th12;

        av2_V200 = 1. + koef*2*mu_C*th1_2 + koef*2*mu_C*t4h2_3;
        av2_V20R = th26*(v200 + v20R) - koef*mu_C*t4h2_3;
        av2_V20L = -th26*(v200 + v20L) - koef*mu_C*t4h2_3;
        av2_V2R0 = th14*(v100 + v1R0) - koef*mu_C*th1_2;
        av2_V2L0 = -th14*(v100 + v1L0) - koef*mu_C*th1_2;
        av2_G0L = -10*th22;
        av2_G0R = 10*th22;

        B[3*m00] = g00 + g00*(th14*(v1R0-v1L0) + th24*(v20R-v20L)) + tau*cas_fg(t_hat, X[i],Y[i]);

        B[3*m00+1] = v100 + th24*v100*(v20R-v20L)
                    + koef*(mu*exp(-g00)-mu_C)*(t4h1_3*(v1R0-2*v100+v1L0) + th2_2*(v10R-2*v100+v10L))
                    + koef*tau*mu*exp(-g00)*(v2RR-v2RL-v2LR+v2LL)/(12*h1*h2)
                    + tau*cas_fu1(t_hat, X[i],Y[i], R,mu, koef);

        B[3*m00+2] = v200 + th14*v200*(v1R0-v1L0)
                    + koef*(mu*exp(-g00)-mu_C)*(t4h2_3*(v20R-2*v200+v20L) + th1_2*(v2R0-2*v200+v2L0))
                    + koef*tau*mu*exp(-g00)*(v1RR-v1RL-v1LR+v1LL)/(12*h1*h2)
                    + tau*cas_fu2(t_hat, X[i],Y[i], R,mu, koef);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*mL0,ag_GL0);
        triplets.emplace_back(3*m00,3*mR0,ag_GR0);
        triplets.emplace_back(3*m00,3*m0L,ag_G0L);
        triplets.emplace_back(3*m00,3*m0R,ag_G0R);
        triplets.emplace_back(3*m00,3*mL0+1,ag_V1L0);
        triplets.emplace_back(3*m00,3*mR0+1,ag_V1R0);
        triplets.emplace_back(3*m00,3*m0L+2,ag_V20L);
        triplets.emplace_back(3*m00,3*m0R+2,ag_V20R);
        //B[3*m00] = cas_g(t_hat,X[i],Y[i]);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+1,3*mR0+1,av1_V1R0);
        triplets.emplace_back(3*m00+1,3*mL0+1,av1_V1L0);
        triplets.emplace_back(3*m00+1,3*m0R+1,av1_V10R);
        triplets.emplace_back(3*m00+1,3*m0L+1,av1_V10L);
        triplets.emplace_back(3*m00+1,3*mL0,av1_GL0);
        triplets.emplace_back(3*m00+1,3*mR0,av1_GR0);

        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);
        triplets.emplace_back(3*m00+2,3*mR0+2,av2_V2R0);
        triplets.emplace_back(3*m00+2,3*mL0+2,av2_V2L0);
        triplets.emplace_back(3*m00+2,3*m0R+2,av2_V20R);
        triplets.emplace_back(3*m00+2,3*m0L+2,av2_V20L);
        triplets.emplace_back(3*m00+2,3*m0L,av2_G0L);
        triplets.emplace_back(3*m00+2,3*m0R,av2_G0R);
        /*av1_V100=1;
        av2_V200=1;
        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = cas_u2(t_hat,X[i],Y[i]);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);*/

    }
    else if (st[i]==1) // left border
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*mR0];  g2R0=VEC[3*mR0+3]; g3R0=VEC[3*mR0+6];
        v100=VEC[3*m00+1]; v1R0=VEC[3*mR0+1]; v12R0=VEC[3*mR0+4]; v13R0=VEC[3*mR0+7];
        v100=0;
        ag_G00 = 1. - th12*v100;
        ag_GR0 = th12*v1R0;
        //ag_V100=-2*th12;
        ag_V1R0=2*th12;

        av1_V100=1;
        av2_V200=1;

        B[3*m00] = g00 + th12*g00*(v1R0-v100) + th12*((g2R0*v12R0-2*gR0*v1R0+g00*v100) - 0.5*(g3R0*v13R0-2*g2R0*v12R0+gR0*v1R0))
                                            + th12*(2-g00)*(v12R0-2*v1R0+v100 - 0.5*(v13R0-2*v12R0+v1R0))
                                            + tau*cas_fg(t_hat, X[i],Y[i]);
        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = cas_u2(t_hat,X[i],Y[i]);
        //ag_G00 = 1;
        //B[3*m00] = cas_g(t_hat,X[i],Y[i]);
        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*mR0,ag_GR0);
        //triplets.emplace_back(3*m00,3*m00+1,ag_V100);
        triplets.emplace_back(3*m00,3*mR0+1,ag_V1R0);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);

    }
    else if (st[i]==2) // right border
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*mL0]; g2R0=VEC[3*mL0-3]; g3R0=VEC[3*mL0-6];
        v100=VEC[3*m00+1]; v1R0=VEC[3*mL0+1]; v12R0=VEC[3*mL0-2]; v13R0=VEC[3*mL0-5]; //here R means L in order to save memory
        v100=0;
        ag_G00 = 1. + th12*v100;
        ag_GL0 = -th12*v1R0;
        //ag_V100=2*th12;
        ag_V1L0=-2*th12;

        B[3*m00] = g00 + th12*g00*(v100-v1R0) - th12*((g00*v100-2*gR0*v1R0+g2R0*v12R0) - 0.5*(gR0*v1R0-2*g2R0*v12R0+g3R0*v13R0))
                                            - th12*(2-g00)*(v100-2*v1R0+v12R0 - 0.5*(v1R0-2*v12R0+v13R0))
                                            + tau*cas_fg(t_hat, X[i],Y[i]);
        av1_V100=1;
        av2_V200=1;
        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = cas_u2(t_hat,X[i],Y[i]);
        //ag_G00 = 1;
        //B[3*m00] = cas_g(t_hat,X[i],Y[i]);
        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*mL0,ag_GL0);
        //triplets.emplace_back(3*m00,3*m00+1,ag_V100);
        triplets.emplace_back(3*m00,3*mL0+1,ag_V1L0);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);

    }
    else if (st[i]==3) // down border
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*m0R];  g2R0=VEC[3*m0R+3*(m0R-m00)]; g3R0=VEC[3*m0R+6*(m0R-m00)]; //here R0 means 0R in order to save memory
        v200=VEC[3*m00+2]; v2R0=VEC[3*m0R+2]; v22R0=VEC[3*m0R+3*(m0R-m00)+2]; v23R0=VEC[3*m0R+6*(m0R-m00)+2];
        v200=0;
        ag_G00 = 1. - th22*v200;
        ag_G0R = th22*v2R0;
        //ag_V200=-2*th22;
        ag_V20R=2*th22;

        B[3*m00] = g00 + th22*g00*(v2R0-v200) + th22*((g2R0*v22R0-2*gR0*v2R0+g00*v200) - 0.5*(g3R0*v23R0-2*g2R0*v22R0+gR0*v2R0))
                                            + th22*(2-g00)*(v22R0-2*v2R0+v200 - 0.5*(v23R0-2*v22R0+v2R0))
                                            + tau*cas_fg(t_hat, X[i],Y[i]);
        av1_V100=1;
        av2_V200=1;
        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = cas_u2(t_hat,X[i],Y[i]);
        //ag_G00 = 1;
        //B[3*m00] = cas_g(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*m0R,ag_G0R);
        //triplets.emplace_back(3*m00,3*m00+2,ag_V200);
        triplets.emplace_back(3*m00,3*m0R+2,ag_V20R);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);
    }
    else if (st[i]==4) // up border
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*m0L];  g2R0=VEC[3*m0L-3*(m00-m0L)]; g3R0=VEC[3*m0L-6*(m00-m0L)]; //here R0 means 0R in order to save memory
        v200=VEC[3*m00+2]; v2R0=VEC[3*m0L+2]; v22R0=VEC[3*m0L-3*(m00-m0L)+2]; v23R0=VEC[3*m0L-6*(m00-m0L)+2]; // also here R means L in order to save memory
        v200=0;
        ag_G00 = 1. + th22*v200;
        ag_G0L = -th22*v2R0;
        //ag_V200=2*th22;
        ag_V20L=-2*th22;

        B[3*m00] = g00 + th22*g00*(v200-v2R0) - th22*((g2R0*v22R0-2*gR0*v2R0+g00*v200) - 0.5*(g3R0*v23R0-2*g2R0*v22R0+gR0*v2R0))
                                            - th22*(2-g00)*(v22R0-2*v2R0+v200 - 0.5*(v23R0-2*v22R0+v2R0))
                                            + tau*cas_fg(t_hat, X[i],Y[i]);
        av1_V100=1;
        av2_V200=1;
        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = cas_u2(t_hat,X[i],Y[i]);
        //ag_G00 = 1;
        //B[3*m00] = cas_g(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*m0L,ag_G0L);
        //triplets.emplace_back(3*m00,3*m00+2,ag_V200);
        triplets.emplace_back(3*m00,3*m0L+2,ag_V20L);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);
    }
    else if (st[i]==5) // left down corner (test)
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00];

        ag_G00 = 1;
        av1_V100=1;
        av2_V200=1;

        B[3*m00] = g00 + tau*cas_fg(t_hat, X[i],Y[i]);
        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = cas_u2(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);

        //B[3*m00] = cas_g(t_hat,X[i],Y[i]);
    }
    else if (st[i]==6) // right down corner (test)
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00];

        ag_G00 = 1;
        av1_V100=1;
        av2_V200=1;

        B[3*m00] = g00 + tau*cas_fg(t_hat, X[i],Y[i]);
        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = cas_u2(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);

        //B[3*m00] = cas_g(t_hat,X[i],Y[i]);
    }
    else if (st[i]==7) // left up corner (test)
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00];

        ag_G00 = 1;
        av1_V100=1;
        av2_V200=1;

        B[3*m00] = g00 + tau*cas_fg(t_hat, X[i],Y[i]);
        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = cas_u2(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);

        //B[3*m00] = cas_g(t_hat,X[i],Y[i]);
    }
    else if (st[i]==8) // right up corner (test)
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00];

        ag_G00 = 1;
        av1_V100=1;
        av2_V200=1;

        B[3*m00] = g00 + tau*cas_fg(t_hat, X[i],Y[i]);

        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = cas_u2(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);

        //B[3*m00] = cas_g(t_hat,X[i],Y[i]);
    }
    else if (st[i]==11) // inner left up corner
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*mR0]; g0R=VEC[3*m0R]; gL0=VEC[3*mL0]; g0L=VEC[3*m0L];
        v100=VEC[3*m00+1]; v1R0=VEC[3*mR0+1]; v10R=VEC[3*m0R+1]; v1L0=VEC[3*mL0+1]; v10L=VEC[3*m0L+1]; v1RR=VEC[3*m0R+4]; v1LL=VEC[3*m0L-2]; v1RL=VEC[3*m0L+4]; v1LR=0;
        v200=VEC[3*m00+2]; v2R0=VEC[3*mR0+2]; v20R=VEC[3*m0R+2]; v2L0=VEC[3*mL0+2]; v20L=VEC[3*m0L+2]; v2RR=VEC[3*m0R+5]; v2LL=VEC[3*m0L-1]; v2RL=VEC[3*m0L+5]; v2LR=0;

        ag_G00 = 1.;
        ag_GL0 = -th14*(v100+v1L0);
        ag_GR0 = th14*(v100+v1R0);
        ag_G0L = -th24*(v200+v20L);
        ag_G0R = th24*(v200+v20R);
        ag_V1L0 = -th12;
        ag_V20L = -th22;
        ag_V1R0 = th12;
        ag_V20R = th22;

        av1_V100 = 1. + koef*2*mu_C*t4h1_3 + koef*2*mu_C*th2_2;
        av1_V1R0 = th16*(v100 + v1R0) - koef*mu_C*t4h1_3;
        av1_V1L0 = -th16*(v100 + v1L0) - koef*mu_C*t4h1_3;
        av1_V10R = th24*(v200 + v20R) - koef*mu_C*th2_2;
        av1_V10L = -th24*(v200 + v20L) - koef*mu_C*th2_2;
        av1_GL0 = -R*th12*t00;
        av1_GR0 = R*th12*t00;

        av2_V200 = 1. + koef*2*mu_C*th1_2 + koef*2*mu_C*t4h2_3;
        av2_V20R = th26*(v200 + v20R) - koef*mu_C*t4h2_3;
        av2_V20L = -th26*(v200 + v20L) - koef*mu_C*t4h2_3;
        av2_V2R0 = th14*(v100 + v1R0) - koef*mu_C*th1_2;
        av2_V2L0 = -th14*(v100 + v1L0) - koef*mu_C*th1_2;
        av2_G0L = -R*th22*t00;
        av2_G0R = R*th22*t00;

        B[3*m00] = g00 + g00*(th14*(v1R0-v1L0) + th24*(v20R-v20L)) + tau*cas_fg(t_hat, X[i],Y[i]);

        B[3*m00+1] = v100 - th12*R*(tR0-tL0) + th24*v100*(v20R-v20L)
                    + koef*(mu*exp(-g00)-mu_C)*(t4h1_3*(v1R0-2*v100+v1L0) + th2_2*(v10R-2*v100+v10L))
                    + koef*tau*mu*exp(-g00)*(v2RR-v2RL-v2LR+v2LL)/(12*h1*h2)
                    + tau*cas_fu1(t_hat, X[i],Y[i], R,mu, koef);

        B[3*m00+2] = v200 - th22*R*(t0R-t0L) + th14*v200*(v1R0-v1L0)
                    + koef*(mu*exp(-g00)-mu_C)*(t4h2_3*(v20R-2*v200+v20L) + th1_2*(v2R0-2*v200+v2L0))
                    + koef*tau*mu*exp(-g00)*(v1RR-v1RL-v1LR+v1LL)/(12*h1*h2)
                    + tau*cas_fu2(t_hat, X[i],Y[i], R,mu, koef);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*mL0,ag_GL0);
        triplets.emplace_back(3*m00,3*mR0,ag_GR0);
        triplets.emplace_back(3*m00,3*m0L,ag_G0L);
        triplets.emplace_back(3*m00,3*m0R,ag_G0R);
        triplets.emplace_back(3*m00,3*mL0+1,ag_V1L0);
        triplets.emplace_back(3*m00,3*mR0+1,ag_V1R0);
        triplets.emplace_back(3*m00,3*m0L+2,ag_V20L);
        triplets.emplace_back(3*m00,3*m0R+2,ag_V20R);
        //B[3*m00] = cas_g(t_hat,X[i],Y[i]);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+1,3*mR0+1,av1_V1R0);
        triplets.emplace_back(3*m00+1,3*mL0+1,av1_V1L0);
        triplets.emplace_back(3*m00+1,3*m0R+1,av1_V10R);
        triplets.emplace_back(3*m00+1,3*m0L+1,av1_V10L);
        triplets.emplace_back(3*m00+1,3*mL0,av1_GL0);
        triplets.emplace_back(3*m00+1,3*mR0,av1_GR0);

        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);
        triplets.emplace_back(3*m00+2,3*mR0+2,av2_V2R0);
        triplets.emplace_back(3*m00+2,3*mL0+2,av2_V2L0);
        triplets.emplace_back(3*m00+2,3*m0R+2,av2_V20R);
        triplets.emplace_back(3*m00+2,3*m0L+2,av2_V20L);
        triplets.emplace_back(3*m00+2,3*m0L,av2_G0L);
        triplets.emplace_back(3*m00+2,3*m0R,av2_G0R);
        /*av1_V100=1;
        av2_V200=1;
        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = cas_u2(t_hat,X[i],Y[i]);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);*/

    }
  }
  //build matrix, solve with eigen
  matrix.setFromTriplets (triplets.begin(), triplets.end());
  eigen_solver_t solver(matrix);

  if (solver.info() != Eigen::Success)
  {
      fprintf(OUT, "Can not build preconditioner");
      return -1;
  }
  eigen_vector_t ress = solver.solveWithGuess(B, VEC);
  if (solver.info() != Eigen::Success)
  {
      fprintf(OUT, "Failed to solve the system with Eigen, n=%d", n);
      return -1;
  }
  VEC=ress;
  t_hat+=tau; mu_C=0; kap_C=0;
}
// solve end =========================================================================================================
fprintf(OUT, "1");
for (int i=0; i<Dim; i++)
  {
    VEC[3*i] -= cas_g(str_param_gas->T,X[i],Y[i]);
    VEC[3*i+1] -= cas_u1(str_param_gas->T,X[i],Y[i]);
    VEC[3*i+2] -= cas_u2(str_param_gas->T,X[i],Y[i]);
  }
//normi
double maxg=0, maxu1=0, maxu2=0, maxt=0;
double sumu1=0, sumu2=0, sumg=0, sumt=0;
for (int i=0; i<Dim; i++)
{
    if (maxg<fabs(VEC[3*i])) maxg=fabs(VEC[3*i]);
    if (maxu1<fabs(VEC[3*i+1])) maxu1=fabs(VEC[3*i+1]);
    if (maxu2<fabs(VEC[3*i+2])) maxu2=fabs(VEC[3*i+2]);

    sumg+=(VEC[3*i])*(VEC[3*i]);
    sumu1+=(VEC[3*i+1])*(VEC[3*i+1]);
    sumu2+=(VEC[3*i+2])*(VEC[3*i+2]);
}
n_c[4*ij]=maxu1; n_c[4*ij+1]=maxu2; n_c[4*ij+2]=maxg; n_c[4*ij+3]=0;
n_l[4*ij]=sqrt(h1*h2*sumu1); n_l[4*ij+1]=sqrt(h1*h2*sumu2); n_l[4*ij+2]=sqrt(h1*h2*sumg); n_l[4*ij+3]=0;
fclose(OUT);
return 0;
}


int Sxema_hcond_vgas_2_nested_grid(double mu, int N,int M1, int M2, int Dim,
                      Str_param_gas *str_param_gas, Str_setka *str_setka, P_she *p_she,
                      int *st, double *X, double *Y, int *M0R, int *M0L,
                      double *n_c, double *n_l, double *n_w, int ij, int K, int cas, int k, double *VGth, double *Tth, double *Xth, double *Yth)
{
    printf("hcond_vgas_2 \n");

double  h1, h2, tau, R, kappa, cv, eps;
tau = str_param_gas->T/N;//p_she->tau;
h1 = str_param_gas->X1/M1;//p_she->h1;
h2 = str_param_gas->X2/M2;//p_she->h2;
R = str_param_gas->R;
kappa = str_param_gas->kappa;
cv = str_param_gas->cv;
eps = str_param_gas->eps;

double EX2 = str_param_gas->EX2;
// assigning functions
double (*cas_u1)(double , double , double);
double (*cas_u2)(double , double , double);
double (*cas_g)(double , double , double);
double (*cas_theta)(double , double , double);

double (*cas_fu1)(double , double , double, double , double , double);
double (*cas_fu2)(double , double , double, double , double , double);
double (*cas_fg)(double , double , double);
double (*cas_ft)(double , double , double, double , double , double , double );

if (cas>1)// for DEBUGGING TEST 1-2 functions u1, u2, g and theta are used
    {
        cas_u1=u1;
        cas_u2=u2;
        cas_g=g;
        cas_theta=theta;

        cas_fu1=fu1; //pravie chasti
        cas_fu2=fu2;
        cas_fg=fg;
        cas_ft=ft;
    }

//declare vectors
eigen_vector_t VEC (3*Dim);
eigen_vector_t T (Dim);
//VEC consists of g, v1, v2: VEC = (G_00, V1_00,V2_00, G_10, V1_10, V2_10, ..., G_M10, V1_M10, V2_M10, G_01, V1_01, V2_01, G_11, V1_11, V2_11,...)
//knots' coord
int m00, mR0, m0R, mL0, m0L;
//initial conditions---------------------------------------
for (int i=0; i<Dim; i++)
{
    VEC[3*i] = cas_g(0,X[i],Y[i]);
    VEC[3*i+1] = cas_u1(0,X[i],Y[i]);
    VEC[3*i+2] = cas_u2(0,X[i],Y[i]);
    T[i] = cas_theta(0,X[i],Y[i]);
}

//const
double th16=tau/(6*h1);
double th26=tau/(6*h2);
double th14=tau/(4*h1);
double th24=tau/(4*h2);
double th12=tau/(2*h1);
double th22=tau/(2*h2);
double t4h1_3=4*tau/(3*h1*h1);
double t4h2_3=4*tau/(3*h2*h2);
double th1_2=tau/(h1*h1);
double th2_2=tau/(h2*h2);

char txt[] = ".txt";
char file[20];
sprintf(file, "%d", ij);
strcat(file,txt);
FILE *OUT;
OUT = fopen(file, "w");
fprintf(OUT, "1 ");

double koef=1.;
double mu_C=0, kap_C=0;
double t_hat=tau;
//function values in the grid konts
double g00=0, gR0=0, g0R=0, gL0=0, g0L=0;
double t00=0, tR0=0, t0R=0, tL0=0, t0L=0;
double v100=0, v1R0=0, v10R=0, v1L0=0, v10L=0, v1RR=0, v1LL=0, v1RL=0, v1LR=0;
double v200=0, v2R0=0, v20R=0, v2L0=0, v20L=0, v2RR=0, v2LL=0, v2RL=0, v2LR=0;
double g2R0=0, g3R0=0, v12R0=0, v13R0=0, v22R0=0, v23R0=0; // extra knots values for border equations
//matrix coefficients for corresponding components
double ag_G00=0, ag_GR0=0, ag_GL0=0, ag_G0R=0, ag_G0L=0, ag_V100=0, ag_V200=0, ag_V1R0=0, ag_V1L0=0, ag_V20R=0, ag_V20L=0;
double av1_GR0=0, av1_GL0=0, av1_V100=0, av1_V1R0=0, av1_V1L0=0, av1_V10R=0, av1_V10L=0, av1_V20R=0, av1_V20L=0;
double av2_G0R=0, av2_G0L=0, av2_V200=0, av2_V2R0=0, av2_V2L0=0, av2_V20R=0, av2_V20L=0, av2_V1R0=0, av2_V1L0=0;
double a_T00=0, a_TR0=0, a_TL0=0, a_T0R=0, a_T0L=0;

//solve=========================================================================================================================================
for (int n=0; n<N; n++)
{
    //declare matrix and right part
eigen_matrix_t matrix (3*Dim, 3*Dim);
eigen_matrix_t matrixt (Dim, Dim);
eigen_vector_t B (3*Dim);
eigen_vector_t Bt (Dim);

std::vector<eigen_triplet_t> triplets;
std::vector<eigen_triplet_t> triplets_t;
// firstly solve: matrix*VEC = B
// then solve: matrixt*T = Bt
  for (int i=0; i<Dim; i++)
  {
    if (mu*exp(-VEC[3*i])>mu_C) mu_C=mu*exp(-VEC[3*i]);
    //if (kappa*exp(-VEC[3*i])>kap_C) kap_C=kappa*exp(-VEC[3*i]);
  }

    R*=0.001;
  /*for (int i=0; i<Dim; i++)
  {
    //VEC[3*i] = cas_g(t_hat,X[i],Y[i]);
    //VEC[3*i+1] = cas_u1(t_hat,X[i],Y[i]);
    //VEC[3*i+2] = cas_u2(t_hat,X[i],Y[i]);
    T[i] = cas_theta(t_hat,X[i],Y[i]);
  }*/
  for (int i=0; i<Dim; i++)
  {
    //fill matrix
    if (st[i]==0) // inner knots
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*mR0]; g0R=VEC[3*m0R]; gL0=VEC[3*mL0]; g0L=VEC[3*m0L];
        t00=T[m00]; tR0=T[mR0]; tL0=T[mL0]; t0R=T[m0R]; t0L=T[m0L];
        v100=VEC[3*m00+1]; v1R0=VEC[3*mR0+1]; v10R=VEC[3*m0R+1]; v1L0=VEC[3*mL0+1]; v10L=VEC[3*m0L+1]; v1RR=VEC[3*m0R+4]; v1LL=VEC[3*m0L-2]; v1RL=VEC[3*m0L+4]; v1LR=VEC[3*m0R-2];
        v200=VEC[3*m00+2]; v2R0=VEC[3*mR0+2]; v20R=VEC[3*m0R+2]; v2L0=VEC[3*mL0+2]; v20L=VEC[3*m0L+2]; v2RR=VEC[3*m0R+5]; v2LL=VEC[3*m0L-1]; v2RL=VEC[3*m0L+5]; v2LR=VEC[3*m0R-1];

        ag_G00 = 1.;
        ag_GL0 = -th14*(v100+v1L0);
        ag_GR0 = th14*(v100+v1R0);
        ag_G0L = -th24*(v200+v20L);
        ag_G0R = th24*(v200+v20R);
        ag_V1L0 = -th12;
        ag_V20L = -th22;
        ag_V1R0 = th12;
        ag_V20R = th22;

        av1_V100 = 1. + koef*2*mu_C*t4h1_3 + koef*2*mu_C*th2_2;
        av1_V1R0 = th16*(v100 + v1R0) - koef*mu_C*t4h1_3;
        av1_V1L0 = -th16*(v100 + v1L0) - koef*mu_C*t4h1_3;
        av1_V10R = th24*(v200 + v20R) - koef*mu_C*th2_2;
        av1_V10L = -th24*(v200 + v20L) - koef*mu_C*th2_2;
        av1_GL0 = -R*th12*t00;
        av1_GR0 = R*th12*t00;

        av2_V200 = 1. + koef*2*mu_C*th1_2 + koef*2*mu_C*t4h2_3;
        av2_V20R = th26*(v200 + v20R) - koef*mu_C*t4h2_3;
        av2_V20L = -th26*(v200 + v20L) - koef*mu_C*t4h2_3;
        av2_V2R0 = th14*(v100 + v1R0) - koef*mu_C*th1_2;
        av2_V2L0 = -th14*(v100 + v1L0) - koef*mu_C*th1_2;
        av2_G0L = -R*th22*t00;
        av2_G0R = R*th22*t00;

        B[3*m00] = g00 + g00*(th14*(v1R0-v1L0) + th24*(v20R-v20L)) + tau*cas_fg(t_hat, X[i],Y[i]);

        B[3*m00+1] = v100 - th12*R*(tR0-tL0) + th24*v100*(v20R-v20L)
                    + koef*(mu*exp(-g00)-mu_C)*(t4h1_3*(v1R0-2*v100+v1L0) + th2_2*(v10R-2*v100+v10L))
                    + koef*tau*mu*exp(-g00)*(v2RR-v2RL-v2LR+v2LL)/(12*h1*h2)
                    + tau*cas_fu1(t_hat, X[i],Y[i], R,mu, koef);

        B[3*m00+2] = v200 - th22*R*(t0R-t0L) + th14*v200*(v1R0-v1L0)
                    + koef*(mu*exp(-g00)-mu_C)*(t4h2_3*(v20R-2*v200+v20L) + th1_2*(v2R0-2*v200+v2L0))
                    + koef*tau*mu*exp(-g00)*(v1RR-v1RL-v1LR+v1LL)/(12*h1*h2)
                    + tau*cas_fu2(t_hat, X[i],Y[i], R,mu, koef);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*mL0,ag_GL0);
        triplets.emplace_back(3*m00,3*mR0,ag_GR0);
        triplets.emplace_back(3*m00,3*m0L,ag_G0L);
        triplets.emplace_back(3*m00,3*m0R,ag_G0R);
        triplets.emplace_back(3*m00,3*mL0+1,ag_V1L0);
        triplets.emplace_back(3*m00,3*mR0+1,ag_V1R0);
        triplets.emplace_back(3*m00,3*m0L+2,ag_V20L);
        triplets.emplace_back(3*m00,3*m0R+2,ag_V20R);
        //B[3*m00] = cas_g(t_hat,X[i],Y[i]);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+1,3*mR0+1,av1_V1R0);
        triplets.emplace_back(3*m00+1,3*mL0+1,av1_V1L0);
        triplets.emplace_back(3*m00+1,3*m0R+1,av1_V10R);
        triplets.emplace_back(3*m00+1,3*m0L+1,av1_V10L);
        triplets.emplace_back(3*m00+1,3*mL0,av1_GL0);
        triplets.emplace_back(3*m00+1,3*mR0,av1_GR0);

        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);
        triplets.emplace_back(3*m00+2,3*mR0+2,av2_V2R0);
        triplets.emplace_back(3*m00+2,3*mL0+2,av2_V2L0);
        triplets.emplace_back(3*m00+2,3*m0R+2,av2_V20R);
        triplets.emplace_back(3*m00+2,3*m0L+2,av2_V20L);
        triplets.emplace_back(3*m00+2,3*m0L,av2_G0L);
        triplets.emplace_back(3*m00+2,3*m0R,av2_G0R);
        /*av1_V100=1;
        av2_V200=1;
        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = cas_u2(t_hat,X[i],Y[i]);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);*/

    }
    else if (st[i]==1) // left border
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*mR0];  g2R0=VEC[3*mR0+3]; g3R0=VEC[3*mR0+6];
        v100=VEC[3*m00+1]; v1R0=VEC[3*mR0+1]; v12R0=VEC[3*mR0+4]; v13R0=VEC[3*mR0+7];
        v100=0;
        ag_G00 = 1. - th12*v100;
        ag_GR0 = th12*v1R0;
        //ag_V100=-2*th12;
        ag_V1R0=2*th12;

        av1_V100=1;
        av2_V200=1;

        B[3*m00] = g00 + th12*g00*(v1R0-v100) + th12*((g2R0*v12R0-2*gR0*v1R0+g00*v100) - 0.5*(g3R0*v13R0-2*g2R0*v12R0+gR0*v1R0))
                                            + th12*(2-g00)*(v12R0-2*v1R0+v100 - 0.5*(v13R0-2*v12R0+v1R0))
                                            + tau*cas_fg(t_hat, X[i],Y[i]);
        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = cas_u2(t_hat,X[i],Y[i]);
        //ag_G00 = 1;
        //B[3*m00] = cas_g(t_hat,X[i],Y[i]);
        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*mR0,ag_GR0);
        //triplets.emplace_back(3*m00,3*m00+1,ag_V100);
        triplets.emplace_back(3*m00,3*mR0+1,ag_V1R0);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);

    }
    else if (st[i]==2) // right border
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*mL0]; g2R0=VEC[3*mL0-3]; g3R0=VEC[3*mL0-6];
        v100=VEC[3*m00+1]; v1R0=VEC[3*mL0+1]; v12R0=VEC[3*mL0-2]; v13R0=VEC[3*mL0-5]; //here R means L in order to save memory
        v100=0;
        ag_G00 = 1. + th12*v100;
        ag_GL0 = -th12*v1R0;
        //ag_V100=2*th12;
        ag_V1L0=-2*th12;

        B[3*m00] = g00 + th12*g00*(v100-v1R0) - th12*((g00*v100-2*gR0*v1R0+g2R0*v12R0) - 0.5*(gR0*v1R0-2*g2R0*v12R0+g3R0*v13R0))
                                            - th12*(2-g00)*(v100-2*v1R0+v12R0 - 0.5*(v1R0-2*v12R0+v13R0))
                                            + tau*cas_fg(t_hat, X[i],Y[i]);
        av1_V100=1;
        av2_V200=1;
        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = cas_u2(t_hat,X[i],Y[i]);
        //ag_G00 = 1;
        //B[3*m00] = cas_g(t_hat,X[i],Y[i]);
        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*mL0,ag_GL0);
        //triplets.emplace_back(3*m00,3*m00+1,ag_V100);
        triplets.emplace_back(3*m00,3*mL0+1,ag_V1L0);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);

    }
    else if (st[i]==3) // down border
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*m0R];  g2R0=VEC[3*m0R+3*(m0R-m00)]; g3R0=VEC[3*m0R+6*(m0R-m00)]; //here R0 means 0R in order to save memory
        v200=VEC[3*m00+2]; v2R0=VEC[3*m0R+2]; v22R0=VEC[3*m0R+3*(m0R-m00)+2]; v23R0=VEC[3*m0R+6*(m0R-m00)+2];
        v200=0;
        ag_G00 = 1. - th22*v200;
        ag_G0R = th22*v2R0;
        //ag_V200=-2*th22;
        ag_V20R=2*th22;

        B[3*m00] = g00 + th22*g00*(v2R0-v200) + th22*((g2R0*v22R0-2*gR0*v2R0+g00*v200) - 0.5*(g3R0*v23R0-2*g2R0*v22R0+gR0*v2R0))
                                            + th22*(2-g00)*(v22R0-2*v2R0+v200 - 0.5*(v23R0-2*v22R0+v2R0))
                                            + tau*cas_fg(t_hat, X[i],Y[i]);
        av1_V100=1;
        av2_V200=1;
        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = cas_u2(t_hat,X[i],Y[i]);
        //ag_G00 = 1;
        //B[3*m00] = cas_g(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*m0R,ag_G0R);
        //triplets.emplace_back(3*m00,3*m00+2,ag_V200);
        triplets.emplace_back(3*m00,3*m0R+2,ag_V20R);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);
    }
    else if (st[i]==4) // up border
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*m0L];  g2R0=VEC[3*m0L-3*(m00-m0L)]; g3R0=VEC[3*m0L-6*(m00-m0L)]; //here R0 means 0R in order to save memory
        v200=VEC[3*m00+2]; v2R0=VEC[3*m0L+2]; v22R0=VEC[3*m0L-3*(m00-m0L)+2]; v23R0=VEC[3*m0L-6*(m00-m0L)+2]; // also here R means L in order to save memory
        v200=0;
        ag_G00 = 1. + th22*v200;
        ag_G0L = -th22*v2R0;
        //ag_V200=2*th22;
        ag_V20L=-2*th22;

        B[3*m00] = g00 + th22*g00*(v200-v2R0) - th22*((g2R0*v22R0-2*gR0*v2R0+g00*v200) - 0.5*(g3R0*v23R0-2*g2R0*v22R0+gR0*v2R0))
                                            - th22*(2-g00)*(v22R0-2*v2R0+v200 - 0.5*(v23R0-2*v22R0+v2R0))
                                            + tau*cas_fg(t_hat, X[i],Y[i]);
        av1_V100=1;
        av2_V200=1;
        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = cas_u2(t_hat,X[i],Y[i]);
        //ag_G00 = 1;
        //B[3*m00] = cas_g(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*m0L,ag_G0L);
        //triplets.emplace_back(3*m00,3*m00+2,ag_V200);
        triplets.emplace_back(3*m00,3*m0L+2,ag_V20L);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);
    }
    else if (st[i]==5) // left down corner (test)
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00];

        ag_G00 = 1;
        av1_V100=1;
        av2_V200=1;

        B[3*m00] = g00 + tau*cas_fg(t_hat, X[i],Y[i]);
        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = cas_u2(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);

        //B[3*m00] = cas_g(t_hat,X[i],Y[i]);
    }
    else if (st[i]==6) // right down corner (test)
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00];

        ag_G00 = 1;
        av1_V100=1;
        av2_V200=1;

        B[3*m00] = g00 + tau*cas_fg(t_hat, X[i],Y[i]);
        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = cas_u2(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);

        //B[3*m00] = cas_g(t_hat,X[i],Y[i]);
    }
    else if (st[i]==7) // left up corner (test)
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00];

        ag_G00 = 1;
        av1_V100=1;
        av2_V200=1;

        B[3*m00] = g00 + tau*cas_fg(t_hat, X[i],Y[i]);
        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = cas_u2(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);

        //B[3*m00] = cas_g(t_hat,X[i],Y[i]);
    }
    else if (st[i]==8) // right up corner (test)
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00];

        ag_G00 = 1;
        av1_V100=1;
        av2_V200=1;

        B[3*m00] = g00 + tau*cas_fg(t_hat, X[i],Y[i]);

        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = cas_u2(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);

        //B[3*m00] = cas_g(t_hat,X[i],Y[i]);
    }
    else if (st[i]==11) // inner left up corner
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*mR0]; g0R=VEC[3*m0R]; gL0=VEC[3*mL0]; g0L=VEC[3*m0L];
        t00=T[m00]; tR0=T[mR0]; tL0=T[mL0]; t0R=T[m0R]; t0L=T[m0L];
        v100=VEC[3*m00+1]; v1R0=VEC[3*mR0+1]; v10R=VEC[3*m0R+1]; v1L0=VEC[3*mL0+1]; v10L=VEC[3*m0L+1]; v1RR=VEC[3*m0R+4]; v1LL=VEC[3*m0L-2]; v1RL=VEC[3*m0L+4]; v1LR=0;
        v200=VEC[3*m00+2]; v2R0=VEC[3*mR0+2]; v20R=VEC[3*m0R+2]; v2L0=VEC[3*mL0+2]; v20L=VEC[3*m0L+2]; v2RR=VEC[3*m0R+5]; v2LL=VEC[3*m0L-1]; v2RL=VEC[3*m0L+5]; v2LR=0;

        ag_G00 = 1.;
        ag_GL0 = -th14*(v100+v1L0);
        ag_GR0 = th14*(v100+v1R0);
        ag_G0L = -th24*(v200+v20L);
        ag_G0R = th24*(v200+v20R);
        ag_V1L0 = -th12;
        ag_V20L = -th22;
        ag_V1R0 = th12;
        ag_V20R = th22;

        av1_V100 = 1. + koef*2*mu_C*t4h1_3 + koef*2*mu_C*th2_2;
        av1_V1R0 = th16*(v100 + v1R0) - koef*mu_C*t4h1_3;
        av1_V1L0 = -th16*(v100 + v1L0) - koef*mu_C*t4h1_3;
        av1_V10R = th24*(v200 + v20R) - koef*mu_C*th2_2;
        av1_V10L = -th24*(v200 + v20L) - koef*mu_C*th2_2;
        av1_GL0 = -R*th12*t00;
        av1_GR0 = R*th12*t00;

        av2_V200 = 1. + koef*2*mu_C*th1_2 + koef*2*mu_C*t4h2_3;
        av2_V20R = th26*(v200 + v20R) - koef*mu_C*t4h2_3;
        av2_V20L = -th26*(v200 + v20L) - koef*mu_C*t4h2_3;
        av2_V2R0 = th14*(v100 + v1R0) - koef*mu_C*th1_2;
        av2_V2L0 = -th14*(v100 + v1L0) - koef*mu_C*th1_2;
        av2_G0L = -R*th22*t00;
        av2_G0R = R*th22*t00;

        B[3*m00] = g00 + g00*(th14*(v1R0-v1L0) + th24*(v20R-v20L)) + tau*cas_fg(t_hat, X[i],Y[i]);

        B[3*m00+1] = v100 - th12*R*(tR0-tL0) + th24*v100*(v20R-v20L)
                    + koef*(mu*exp(-g00)-mu_C)*(t4h1_3*(v1R0-2*v100+v1L0) + th2_2*(v10R-2*v100+v10L))
                    + koef*tau*mu*exp(-g00)*(v2RR-v2RL-v2LR+v2LL)/(12*h1*h2)
                    + tau*cas_fu1(t_hat, X[i],Y[i], R,mu, koef);

        B[3*m00+2] = v200 - th22*R*(t0R-t0L) + th14*v200*(v1R0-v1L0)
                    + koef*(mu*exp(-g00)-mu_C)*(t4h2_3*(v20R-2*v200+v20L) + th1_2*(v2R0-2*v200+v2L0))
                    + koef*tau*mu*exp(-g00)*(v1RR-v1RL-v1LR+v1LL)/(12*h1*h2)
                    + tau*cas_fu2(t_hat, X[i],Y[i], R,mu, koef);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*mL0,ag_GL0);
        triplets.emplace_back(3*m00,3*mR0,ag_GR0);
        triplets.emplace_back(3*m00,3*m0L,ag_G0L);
        triplets.emplace_back(3*m00,3*m0R,ag_G0R);
        triplets.emplace_back(3*m00,3*mL0+1,ag_V1L0);
        triplets.emplace_back(3*m00,3*mR0+1,ag_V1R0);
        triplets.emplace_back(3*m00,3*m0L+2,ag_V20L);
        triplets.emplace_back(3*m00,3*m0R+2,ag_V20R);
        //B[3*m00] = cas_g(t_hat,X[i],Y[i]);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+1,3*mR0+1,av1_V1R0);
        triplets.emplace_back(3*m00+1,3*mL0+1,av1_V1L0);
        triplets.emplace_back(3*m00+1,3*m0R+1,av1_V10R);
        triplets.emplace_back(3*m00+1,3*m0L+1,av1_V10L);
        triplets.emplace_back(3*m00+1,3*mL0,av1_GL0);
        triplets.emplace_back(3*m00+1,3*mR0,av1_GR0);

        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);
        triplets.emplace_back(3*m00+2,3*mR0+2,av2_V2R0);
        triplets.emplace_back(3*m00+2,3*mL0+2,av2_V2L0);
        triplets.emplace_back(3*m00+2,3*m0R+2,av2_V20R);
        triplets.emplace_back(3*m00+2,3*m0L+2,av2_V20L);
        triplets.emplace_back(3*m00+2,3*m0L,av2_G0L);
        triplets.emplace_back(3*m00+2,3*m0R,av2_G0R);
        /*av1_V100=1;
        av2_V200=1;
        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = cas_u2(t_hat,X[i],Y[i]);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);*/

    }
  }
  //build matrix, solve with eigen
  matrix.setFromTriplets (triplets.begin(), triplets.end());
  eigen_solver_t solver(matrix);

  if (solver.info() != Eigen::Success)
  {
      fprintf(OUT, "Can not build preconditioner");
      return -1;
  }
  eigen_vector_t ress = solver.solveWithGuess(B, VEC);
  if (solver.info() != Eigen::Success)
  {
      fprintf(OUT, "Failed to solve the system with Eigen, n=%d", n);
      n_c[4*ij]=0; n_c[4*ij+1]=0; n_c[4*ij+2]=0; n_c[4*ij+3]=0;
      n_l[4*ij]=0; n_l[4*ij+1]=0; n_l[4*ij+2]=0; n_l[4*ij+3]=0;
      return 0;
  }
  VEC=ress;
  /*for (int i=0; i<Dim; i++)
  {
    VEC[3*i] = cas_g(t_hat,X[i],Y[i]);
    VEC[3*i+1] = cas_u1(t_hat,X[i],Y[i]);
    VEC[3*i+2] = cas_u2(t_hat,X[i],Y[i]);
    //T[i] = cas_theta(t_hat,X[i],Y[i]);
  }*/
  R*=1000;
  for (int i=0; i<Dim; i++)
  {
    if (kappa*exp(-VEC[3*i])>kap_C) kap_C=kappa*exp(-VEC[3*i]);
  }
  // solve matrixt*T=Bt
  for (int i=0; i<Dim; i++)
  {
    //T[i] = cas_theta(t_hat,X[i],Y[i]);
    // fill matrixt
    if (st[i]==0) // inner knots
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        //m0R=i+M1+1;//M0R[i];
        //m0L=i-M1-1;//M0L[i];
        m0R=M0R[i];
        m0L=M0L[i];

        g00=VEC[3*m00];
        t00=T[m00]; tR0=T[mR0]; tL0=T[mL0]; t0R=T[m0R]; t0L=T[m0L];
        v100=VEC[3*m00+1]; v1R0=VEC[3*mR0+1]; v10R=VEC[3*m0R+1]; v1L0=VEC[3*mL0+1]; v10L=VEC[3*m0L+1]; // try to take 3*i+3*r+1 instead of 3*i+3*r+4 & 3*i-3*l+1 instead of 3*i-3*l-2
        v200=VEC[3*m00+2]; v2R0=VEC[3*mR0+2]; v20R=VEC[3*m0R+2]; v2L0=VEC[3*mL0+2]; v20L=VEC[3*m0L+2];

        a_T00 = cv + 2*kap_C*th1_2 + 2*kap_C*th2_2;
        a_TR0 = cv*th14*(v100+v1R0) - kap_C*th1_2;
        a_TL0 = -cv*th14*(v100+v1L0) - kap_C*th1_2;
        a_T0R = cv*th24*(v200+v20R) - kap_C*th2_2;
        a_T0L = -cv*th24*(v200+v20L) - kap_C*th2_2;

        Bt[m00] = cv*t00 + cv*t00*(th14*(v1R0-v1L0)+th24*(v20R-v20L)) + (kappa*exp(-g00)-kap_C)*(th1_2*(tR0-2*t00+tL0)+th2_2*(t0R-2*t00+t0L))
                                                       - R*t00*(th12*(v1R0-v1L0)+th22*(v20R-v20L))
                                                       - (1./6)*mu*exp(-g00)*(th1_2*(v1R0-v1L0)*(v1R0-v1L0)+th2_2*(v20R-v20L)*(v20R-v20L))
                                                       - (tau/(3*h1*h2))*mu*exp(-g00)*(v1R0-v1L0)*(v20R-v20L)
                                                       + 0.5*mu*exp(-g00) * (th2_2*(v20R-v20L)*(v20R-v20L) + th1_2*(v1R0-v1L0)*(v1R0-v1L0))
                                                       + 0.25*tau*mu*exp(-g00) * ((v2R0-v2L0)/h1+(v10R-v10L)/h2) * ((v2R0-v2L0)/h1+(v10R-v10L)/h2)
                                                       + tau*cas_ft(t_hat,X[i],Y[i],R,mu,cv,kappa);
        triplets_t.emplace_back(m00,m00,a_T00);
        triplets_t.emplace_back(m00,mR0,a_TR0);
        triplets_t.emplace_back(m00,mL0,a_TL0);
        triplets_t.emplace_back(m00,m0R,a_T0R);
        triplets_t.emplace_back(m00,m0L,a_T0L);
    }
    else if (st[i]==1) // left border
    {
        m00=i;
        a_T00=1;
        Bt[m00] = cas_theta(t_hat,X[i],Y[i]);
        triplets_t.emplace_back(m00,m00,a_T00); // matrixt is Dim*Dim, T is Dim, i=0,..Dim-1; first coord is matrixt's string, second coord is matrixt's row
    }
    else if (st[i]==2) // right border
    {
        m00=i;
        a_T00=1;
        Bt[m00] = cas_theta(t_hat,X[i],Y[i]);
        triplets_t.emplace_back(m00,m00,a_T00);
    }
    else if (st[i]==3) // down border
    {
        m00=i;
        a_T00=1;
        Bt[m00] = cas_theta(t_hat,X[i],Y[i]);
        triplets_t.emplace_back(m00,m00,a_T00);
    }
    else if (st[i]==4) // up border
    {
        m00=i;
        a_T00=1;
        Bt[m00] = cas_theta(t_hat,X[i],Y[i]);
        triplets_t.emplace_back(m00,m00,a_T00);
    }
    else if (st[i]==5) // left down corner (test)
    {
        m00=i;
        a_T00=1;
        Bt[m00] = cas_theta(t_hat,X[i],Y[i]);
        triplets_t.emplace_back(m00,m00,a_T00);
    }
    else if (st[i]==6) // right down corner (test)
    {
        m00=i;
        a_T00=1;
        Bt[i] = cas_theta(t_hat,X[i],Y[i]);
        triplets_t.emplace_back(m00,m00,a_T00);
    }
    else if (st[i]==7) // left up corner (test)
    {
        m00=i;
        a_T00=1;
        Bt[m00] = cas_theta(t_hat,X[i],Y[i]);
        triplets_t.emplace_back(m00,m00,a_T00);
    }
    else if (st[i]==8) // right up corner (test)
    {
        m00=i;
        a_T00=1;
        Bt[m00] = cas_theta(t_hat,X[i],Y[i]);
        triplets_t.emplace_back(m00,m00,a_T00);
    }
    else if (st[i]==11) // inner left up corner
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00];
        t00=T[m00]; tR0=T[mR0]; t0R=T[m0R]; tL0=T[mL0]; t0L=T[m0L];
        v100=VEC[3*m00+1]; v1R0=VEC[3*mR0+1]; v10R=VEC[3*m0R+1]; v1L0=VEC[3*mL0+1]; v10L=VEC[3*m0L+1]; // try to take 3*i+3*r+1 instead of 3*i+3*r+4 & 3*i-3*l+1 instead of 3*i-3*l-2
        v200=VEC[3*m00+2]; v2R0=VEC[3*mR0+2]; v20R=VEC[3*m0R+2]; v2L0=VEC[3*mL0+2]; v20L=VEC[3*m0L+2];

        a_T00 = cv + 2*kap_C*th1_2 + 2*kap_C*th2_2;
        a_TR0 = cv*th14*(v100+v1R0) - kap_C*th1_2;
        a_TL0 = -cv*th14*(v100+v1L0) - kap_C*th1_2;
        a_T0R = cv*th24*(v200+v20R) - kap_C*th2_2;
        a_T0L = -cv*th24*(v200+v20L) - kap_C*th2_2;

        Bt[m00] = cv*t00 + cv*t00*(th14*(v1R0-v1L0)+th24*(v20R-v20L)) + (kappa*exp(-g00)-kap_C)*(th1_2*(tR0-2*t00+tL0)+th2_2*(t0R-2*t00+t0L))
                                                       - R*t00*(th12*(v1R0-v1L0)+th22*(v20R-v20L))
                                                       - (1./6)*mu*exp(-g00)*(th1_2*(v1R0-v1L0)*(v1R0-v1L0)+th2_2*(v20R-v20L)*(v20R-v20L))
                                                       - (tau/(3*h1*h2))*mu*exp(-g00)*(v1R0-v1L0)*(v20R-v20L)
                                                       + 0.5*mu*exp(-g00) * (th2_2*(v20R-v20L)*(v20R-v20L) + th1_2*(v1R0-v1L0)*(v1R0-v1L0))
                                                       + 0.5*tau*mu*exp(-g00) * ((v2R0-v2L0)/h1+(v10R-v10L)/h2) * ((v2R0-v2L0)/h1+(v10R-v10L)/h2)
                                                       + tau*cas_ft(t_hat,X[i],Y[i],R,mu,cv,kappa);
        triplets_t.emplace_back(m00,m00,a_T00);
        triplets_t.emplace_back(m00,mR0,a_TR0);
        triplets_t.emplace_back(m00,mL0,a_TL0);
        triplets_t.emplace_back(m00,m0R,a_T0R);
        triplets_t.emplace_back(m00,m0L,a_T0L);
    }
  }
  //build matrixt, solve with eigen
  matrixt.setFromTriplets (triplets_t.begin(), triplets_t.end());
  eigen_solver_t solvert(matrixt);

  if (solvert.info() != Eigen::Success)
  {
      fprintf(OUT, "Can not build preconditioner");
      return -1;
  }
  eigen_vector_t res (Dim);
  res = solvert.solveWithGuess(Bt, T);
  if (solvert.info() != Eigen::Success)
  {
      fprintf(OUT, "Failed to solve the system with Eigen, n=%d", n);
      n_c[4*ij]=0; n_c[4*ij+1]=0; n_c[4*ij+2]=0; n_c[4*ij+3]=0;
      n_l[4*ij]=0; n_l[4*ij+1]=0; n_l[4*ij+2]=0; n_l[4*ij+3]=0;
      return 0;
  }
  T=res;
  t_hat+=tau; mu_C=0; kap_C=0;
}
// solve end =========================================================================================================
fprintf(OUT, "1");
int dim=0;
/*for (int i=0; i<Dim; i++)
{
    if (fabs(X[i]-Xth[dim])<h1/2.) dim++;
}*/
double maxg=0, maxu1=0, maxu2=0, maxt=0;
double sumu1=0, sumu2=0, sumg=0, sumt=0;
if (k==0)
{
  for (int i=0; i<Dim; i++)
  {
    Xth[i]=X[i];
    Yth[i]=Y[i];
    VGth[3*i]=VEC[3*i];
    VGth[3*i+1]=VEC[3*i+1];
    VGth[3*i+2]=VEC[3*i+2];
    Tth[i]=T[i];
    VEC[3*i] -= cas_g(str_param_gas->T,X[i],Y[i]);
    VEC[3*i+1] -= cas_u1(str_param_gas->T,X[i],Y[i]);
    VEC[3*i+2] -= cas_u2(str_param_gas->T,X[i],Y[i]);
    T[i] -= cas_theta(str_param_gas->T,X[i],Y[i]);
  }
  dim=0;
  for (int i=0; i<Dim; i++)
  {
    if (maxg<fabs(VEC[3*i])) maxg=fabs(VEC[3*i]);
    if (maxu1<fabs(VEC[3*i+1])) maxu1=fabs(VEC[3*i+1]);
    if (maxu2<fabs(VEC[3*i+2])) maxu2=fabs(VEC[3*i+2]);
    if (maxt<fabs(T[i])) {maxt=fabs(T[i]);}// fprintf(OUT, " %d/%d ", i,st[i]);}
    sumg+=(VEC[3*i])*(VEC[3*i]);
    sumu1+=(VEC[3*i+1])*(VEC[3*i+1]);
    sumu2+=(VEC[3*i+2])*(VEC[3*i+2]);
    sumt+=(T[i])*(T[i]);
  }
}
else
{
  for (int i=0; i<Dim; i++)
  {
    if (fabs(X[i]-Xth[dim])<h1*h1*h1*h1 && fabs(Y[i]-Yth[dim])<h2*h2*h2*h2)
    {
      VEC[3*i] -= VGth[3*dim];
      VEC[3*i+1] -= VGth[3*dim+1];
      VEC[3*i+2] -= VGth[3*dim+2];
      T[i] -= Tth[dim];
      dim++;
    }
  }
  dim=0;
  for (int i=0; i<Dim; i++)
  {
    if (fabs(X[i]-Xth[dim])<h1*h1*h1*h1 && fabs(Y[i]-Yth[dim])<h2*h2*h2*h2)
    {
      if (maxg<fabs(VEC[3*i])) maxg=fabs(VEC[3*i]);
      if (maxu1<fabs(VEC[3*i+1])) maxu1=fabs(VEC[3*i+1]);
      if (maxu2<fabs(VEC[3*i+2])) maxu2=fabs(VEC[3*i+2]);
      if (maxt<fabs(T[i])) {maxt=fabs(T[i]);}
      sumg+=(VEC[3*i])*(VEC[3*i]);
      sumu1+=(VEC[3*i+1])*(VEC[3*i+1]);
      sumu2+=(VEC[3*i+2])*(VEC[3*i+2]);
      sumt+=(T[i])*(T[i]);
      dim++;
    }
  }
}
fprintf(OUT, "2");
if (dim==0) dim=1./(h1*h2);
n_c[4*ij]=maxu1; n_c[4*ij+1]=maxu2; n_c[4*ij+2]=maxg; n_c[4*ij+3]=maxt;
n_l[4*ij]=sqrt(sumu1/dim); n_l[4*ij+1]=sqrt(sumu2/dim); n_l[4*ij+2]=sqrt(sumg/dim); n_l[4*ij+3]=sqrt(sumt/dim);
fclose(OUT);
return 0;
}

int Sxema_hcond_vgas_2_flowp(double mu, int N,int M1, int M2, int Dim,
                      Str_param_gas *str_param_gas, Str_setka *str_setka, P_she *p_she,
                      int *st, double *X, double *Y, int *M0R, int *M0L,
                      double *n_c, double *n_l, double *n_w, int ij, int K, int cas, double *t_len, double tv, double tg, double tt)
{
printf("hcond_vgas_2_flowp \n");

double  h1, h2, tau, R, kappa, cv, eps;
tau = str_param_gas->T/N;//p_she->tau;
h1 = str_param_gas->X1/M1;//p_she->h1;
h2 = str_param_gas->X2/M2;//p_she->h2;
R = str_param_gas->R;
kappa = str_param_gas->kappa;
cv = str_param_gas->cv;
eps = str_param_gas->eps;

double EX2 = str_param_gas->EX2;
// assigning functions
double (*cas_u1)(double , double , double, double , double );
double (*cas_u2)(double , double , double );
double (*cas_g)(double , double , double, double , double );
double (*cas_theta)(double , double , double, double , double );

if (cas==4)// for DEBUGGING TEST 1-2 functions u1, u2, g and theta are used
    {
        cas_u1=u1_1;
        cas_u2=u2_1;
        cas_g=g_1;
        cas_theta=theta_1;
    }

//declare vectors
eigen_vector_t VEC (3*Dim);
eigen_vector_t T (Dim);
//VEC consists of g, v1, v2: VEC = (G_00, V1_00,V2_00, G_10, V1_10, V2_10, ..., G_M10, V1_M10, V2_M10, G_01, V1_01, V2_01, G_11, V1_11, V2_11,...)
//knots' coord
int m00, mR0, m0R, mL0, m0L;
//initial conditions---------------------------------------
for (int i=0; i<Dim; i++)
{
    VEC[3*i] = cas_g(0,X[i],Y[i],h1,tg);
    VEC[3*i+1] = cas_u1(0,X[i],Y[i],h1,tv);
    VEC[3*i+2] = cas_u2(0,X[i],Y[i]);
    T[i] = cas_theta(0,X[i],Y[i],h1,tt);
}

//const
double th16=tau/(6*h1);
double th26=tau/(6*h2);
double th14=tau/(4*h1);
double th24=tau/(4*h2);
double th12=tau/(2*h1);
double th22=tau/(2*h2);
double t4h1_3=4*tau/(3*h1*h1);
double t4h2_3=4*tau/(3*h2*h2);
double th1_2=tau/(h1*h1);
double th2_2=tau/(h2*h2);

char txt[] = ".txt";
char file[20];
sprintf(file, "%d", ij);
strcat(file,txt);
FILE *OUT;
OUT = fopen(file, "w");
fprintf(OUT, "1 ");
FILE *OUT1;
FILE *OUT2;
FILE *OUT3;
double koef=1;
double mu_C=0, kap_C=0;
double t_hat=tau;
//function values in the grid konts
double g00=0, gR0=0, g0R=0, gL0=0, g0L=0;
double t00=0, tR0=0, t0R=0, tL0=0, t0L=0;
double v100=0, v1R0=0, v10R=0, v1L0=0, v10L=0, v1RR=0, v1LL=0, v1RL=0, v1LR=0;
double v200=0, v2R0=0, v20R=0, v2L0=0, v20L=0, v2RR=0, v2LL=0, v2RL=0, v2LR=0;
double g2R0=0, g3R0=0, v12R0=0, v13R0=0, v22R0=0, v23R0=0; // extra knots values for border equations
//matrix coefficients for corresponding components
double ag_G00=0, ag_GR0=0, ag_GL0=0, ag_G0R=0, ag_G0L=0, ag_V100=0, ag_V200=0, ag_V1R0=0, ag_V1L0=0, ag_V20R=0, ag_V20L=0;
double av1_GR0=0, av1_GL0=0, av1_V100=0, av1_V1R0=0, av1_V1L0=0, av1_V10R=0, av1_V10L=0, av1_V20R=0, av1_V20L=0;
double av2_G0R=0, av2_G0L=0, av2_V200=0, av2_V2R0=0, av2_V2L0=0, av2_V20R=0, av2_V20L=0, av2_V1R0=0, av2_V1L0=0;
double a_T00=0, a_TR0=0, a_TL0=0, a_T0R=0, a_T0L=0;

//solve=========================================================================================================================================
int n=0; double dif=-1;
double maxg0=-1, maxu10=-1, maxu20=-1, maxt0=-1;
double maxg=-1, maxu1=-1, maxu2=-1, maxt=-1;
double sumu1=0, sumu2=0, sumg=0, sumt=0; double eta=0.0;
eigen_vector_t ress_st (3*Dim);
eigen_vector_t res_st (Dim);
while (fabs(dif)>eps && n<3003)
{
    //declare matrix and right part
eigen_matrix_t matrix (3*Dim, 3*Dim);
eigen_matrix_t matrixt (Dim, Dim);
eigen_vector_t B (3*Dim);
eigen_vector_t Bt (Dim);

std::vector<eigen_triplet_t> triplets;
std::vector<eigen_triplet_t> triplets_t;
// firstly solve: matrix*VEC = B
// then solve: matrixt*T = Bt
  for (int i=0; i<Dim; i++)
  {
    if (mu*exp(-VEC[3*i])>mu_C) mu_C=mu*exp(-VEC[3*i]);
  }

    R*=0.001;
    // solve matrix*VEC=B
  for (int i=0; i<Dim; i++)
  {
    //fill matrix
    if (st[i]==0) // inner knots
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*mR0]; g0R=VEC[3*m0R]; gL0=VEC[3*mL0]; g0L=VEC[3*m0L];
        t00=T[m00]; tR0=T[mR0]; tL0=T[mL0]; t0R=T[m0R]; t0L=T[m0L];
        v100=VEC[3*m00+1]; v1R0=VEC[3*mR0+1]; v10R=VEC[3*m0R+1]; v1L0=VEC[3*mL0+1]; v10L=VEC[3*m0L+1]; v1RR=VEC[3*m0R+4]; v1LL=VEC[3*m0L-2]; v1RL=VEC[3*m0L+4]; v1LR=VEC[3*m0R-2];
        v200=VEC[3*m00+2]; v2R0=VEC[3*mR0+2]; v20R=VEC[3*m0R+2]; v2L0=VEC[3*mL0+2]; v20L=VEC[3*m0L+2]; v2RR=VEC[3*m0R+5]; v2LL=VEC[3*m0L-1]; v2RL=VEC[3*m0L+5]; v2LR=VEC[3*m0R-1];

        ag_G00 = 1. + 2*eta*tau*th1_2 + 2*eta*tau*th2_2;
        ag_GL0 = -th14*(v100+v1L0) - eta*tau*th1_2;
        ag_GR0 = th14*(v100+v1R0) - eta*tau*th1_2;
        ag_G0L = -th24*(v200+v20L) - eta*tau*th2_2;
        ag_G0R = th24*(v200+v20R) - eta*tau*th2_2;
        ag_V1L0 = -th12;
        ag_V20L = -th22;
        ag_V1R0 = th12;
        ag_V20R = th22;

        av1_V100 = 1. + koef*2*mu_C*t4h1_3 + koef*2*mu_C*th2_2;
        av1_V1R0 = th16*(v100 + v1R0) - koef*mu_C*t4h1_3;
        av1_V1L0 = -th16*(v100 + v1L0) - koef*mu_C*t4h1_3;
        av1_V10R = th24*(v200 + v20R) - koef*mu_C*th2_2;
        av1_V10L = -th24*(v200 + v20L) - koef*mu_C*th2_2;
        av1_GL0 = -R*th12*t00;
        av1_GR0 = R*th12*t00;

        av2_V200 = 1. + koef*2*mu_C*th1_2 + koef*2*mu_C*t4h2_3;
        av2_V20R = th26*(v200 + v20R) - koef*mu_C*t4h2_3;
        av2_V20L = -th26*(v200 + v20L) - koef*mu_C*t4h2_3;
        av2_V2R0 = th14*(v100 + v1R0) - koef*mu_C*th1_2;
        av2_V2L0 = -th14*(v100 + v1L0) - koef*mu_C*th1_2;
        av2_G0L = -R*th22*t00;
        av2_G0R = R*th22*t00;

        B[3*m00] = g00 + g00*(th14*(v1R0-v1L0) + th24*(v20R-v20L));// + tau*cas_fg(t_hat, X[i],Y[i]);

        B[3*m00+1] = v100 - th12*R*(tR0-tL0) + th24*v100*(v20R-v20L)
                    + koef*(mu*exp(-g00)-mu_C)*(t4h1_3*(v1R0-2*v100+v1L0) + th2_2*(v10R-2*v100+v10L))
                    + koef*tau*mu*exp(-g00)*(v2RR-v2RL-v2LR+v2LL)/(12*h1*h2)
                    ;//+ tau*cas_fu1(t_hat, X[i],Y[i], R,mu, koef);

        B[3*m00+2] = v200 - th22*R*(t0R-t0L) + th14*v200*(v1R0-v1L0)
                    + koef*(mu*exp(-g00)-mu_C)*(t4h2_3*(v20R-2*v200+v20L) + th1_2*(v2R0-2*v200+v2L0))
                    + koef*tau*mu*exp(-g00)*(v1RR-v1RL-v1LR+v1LL)/(12*h1*h2)
                    ;//+ tau*cas_fu2(t_hat, X[i],Y[i], R,mu, koef);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*mL0,ag_GL0);
        triplets.emplace_back(3*m00,3*mR0,ag_GR0);
        triplets.emplace_back(3*m00,3*m0L,ag_G0L);
        triplets.emplace_back(3*m00,3*m0R,ag_G0R);
        triplets.emplace_back(3*m00,3*mL0+1,ag_V1L0);
        triplets.emplace_back(3*m00,3*mR0+1,ag_V1R0);
        triplets.emplace_back(3*m00,3*m0L+2,ag_V20L);
        triplets.emplace_back(3*m00,3*m0R+2,ag_V20R);

        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+1,3*mR0+1,av1_V1R0);
        triplets.emplace_back(3*m00+1,3*mL0+1,av1_V1L0);
        triplets.emplace_back(3*m00+1,3*m0R+1,av1_V10R);
        triplets.emplace_back(3*m00+1,3*m0L+1,av1_V10L);
        triplets.emplace_back(3*m00+1,3*mL0,av1_GL0);
        triplets.emplace_back(3*m00+1,3*mR0,av1_GR0);

        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);
        triplets.emplace_back(3*m00+2,3*mR0+2,av2_V2R0);
        triplets.emplace_back(3*m00+2,3*mL0+2,av2_V2L0);
        triplets.emplace_back(3*m00+2,3*m0R+2,av2_V20R);
        triplets.emplace_back(3*m00+2,3*m0L+2,av2_V20L);
        triplets.emplace_back(3*m00+2,3*m0L,av2_G0L);
        triplets.emplace_back(3*m00+2,3*m0R,av2_G0R);
    }
    else if (st[i]==1) // left border flow start
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*mR0];  g2R0=VEC[3*mR0+3]; g3R0=VEC[3*mR0+6];
        v100=VEC[3*m00+1]; v1R0=VEC[3*mR0+1]; v12R0=VEC[3*mR0+4]; v13R0=VEC[3*mR0+7];
        //v100=0;
        ag_G00 = 1.;// - th12*v100;
        ag_GR0 = 0;//th12*v1R0;
        ag_V1R0=0;//2*th12;

        av1_V100=1;
        av2_V200=1;

        B[3*m00] = cas_g(t_hat,X[i],Y[i],h1,tg);//g00 + th12*g00*(v1R0-v100) + th12*((g2R0*v12R0-2*gR0*v1R0+g00*v100) - 0.5*(g3R0*v13R0-2*g2R0*v12R0+gR0*v1R0))
                                            //+ th12*(2-g00)*(v12R0-2*v1R0+v100 - 0.5*(v13R0-2*v12R0+v1R0));
        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i],h1,tv);
        B[3*m00+2] = 0;//cas_u2(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*mR0,ag_GR0);
        triplets.emplace_back(3*m00,3*mR0+1,ag_V1R0);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);

    }
    else if (st[i]==2) // left border no flow
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*mR0];  g2R0=VEC[3*mR0+3]; g3R0=VEC[3*mR0+6];
        v100=VEC[3*m00+1]; v1R0=VEC[3*mR0+1]; v12R0=VEC[3*mR0+4]; v13R0=VEC[3*mR0+7];
        v100=0;
        ag_G00 = 1. - th12*v100 + eta*tau*th1_2;
        ag_GR0 = th12*v1R0 -eta*tau*th1_2;
        ag_V100=-2*th12;
        ag_V1R0=2*th12;

        av1_V100=1;
        av2_V200=1;

        B[3*m00] = g00 + th12*g00*(v1R0-v100) + th12*((g2R0*v12R0-2*gR0*v1R0+g00*v100) - 0.5*(g3R0*v13R0-2*g2R0*v12R0+gR0*v1R0))
                                            + th12*(2-g00)*(v12R0-2*v1R0+v100 - 0.5*(v13R0-2*v12R0+v1R0));
        B[3*m00+1] = 0;//cas_u1(t_hat,X[i],Y[i],h1,tv);
        B[3*m00+2] = 0;//cas_u2(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*mR0,ag_GR0);
        triplets.emplace_back(3*m00,3*mR0+1,ag_V1R0);
        triplets.emplace_back(3*m00,3*m00+1,ag_V100);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);

    }
    else if (st[i]==3) // right border
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*mL0]; g2R0=VEC[3*mL0-3]; g3R0=VEC[3*mL0-6];
        v100=VEC[3*m00+1]; v1R0=VEC[3*mL0+1]; v12R0=VEC[3*mL0-2]; v13R0=VEC[3*mL0-5]; //here R means L in order to save memory
        //v100=0;
        ag_G00 = 1. + 2*th12*v100 - eta*tau*th1_2;
        ag_GL0 = -2*th12*v100 + eta*tau*th1_2;
        //ag_V1L0=-2*th12;

        B[3*m00] = g00;// + th12*g00*(v100-v1R0) - th12*((g00*v100-2*gR0*v1R0+g2R0*v12R0) - 0.5*(gR0*v1R0-2*g2R0*v12R0+g3R0*v13R0))
                                            //- th12*(2-g00)*(v100-2*v1R0+v12R0 - 0.5*(v1R0-2*v12R0+v13R0));
        av1_V100=1;
        av1_V1R0=-1;
        av2_V200=1;
        B[3*m00+1] = 0;//cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = 0;//cas_u2(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*mL0,ag_GL0);
        //triplets.emplace_back(3*m00,3*mL0+1,ag_V1L0);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+1,3*mL0+1,av1_V1R0);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);

    }
    else if (st[i]==4) // down border
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*m0R];  g2R0=VEC[3*m0R+3*(m0R-m00)]; g3R0=VEC[3*m0R+6*(m0R-m00)]; //here R0 means 0R in order to save memory
        v200=VEC[3*m00+2]; v2R0=VEC[3*m0R+2]; v22R0=VEC[3*m0R+3*(m0R-m00)+2]; v23R0=VEC[3*m0R+6*(m0R-m00)+2];
        v200=0;
        ag_G00 = 1. - th22*v200 + eta*tau*th2_2;
        ag_G0R = th22*v2R0 - eta*tau*th2_2;
        ag_V200=-2*th22;
        ag_V20R=2*th22;

        B[3*m00] = g00 + th22*g00*(v2R0-v200) + th22*((g2R0*v22R0-2*gR0*v2R0+g00*v200) - 0.5*(g3R0*v23R0-2*g2R0*v22R0+gR0*v2R0))
                                            + th22*(2-g00)*(v22R0-2*v2R0+v200 - 0.5*(v23R0-2*v22R0+v2R0));
        av1_V100=1;
        //av1_V1R0=1;
        av2_V200=1;
        B[3*m00+1] = v100;//cas_u1(t_hat,X[i],Y[i],h1,tv);
        B[3*m00+2] = 0;//cas_u2(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*m0R,ag_G0R);
        triplets.emplace_back(3*m00,3*m00+2,ag_V200);
        triplets.emplace_back(3*m00,3*m0R+2,ag_V20R);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        //triplets.emplace_back(3*m00+1,3*m0R+1,av1_V1R0);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);
    }
    else if (st[i]==5) // up border
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*m0L];  g2R0=VEC[3*m0L-3*(m00-m0L)]; g3R0=VEC[3*m0L-6*(m00-m0L)]; //here R0 means 0R in order to save memory
        v200=VEC[3*m00+2]; v2R0=VEC[3*m0L+2]; v22R0=VEC[3*m0L-3*(m00-m0L)+2]; v23R0=VEC[3*m0L-6*(m00-m0L)+2]; // also here R means L in order to save memory
        v200=0;
        ag_G00 = 1. + th22*v200 - eta*tau*th2_2;
        ag_G0L = -th22*v2R0 + eta*tau*th2_2;
        ag_V200=2*th22;
        ag_V20L=-2*th22;

        B[3*m00] = g00 + th22*g00*(v200-v2R0) - th22*((g2R0*v22R0-2*gR0*v2R0+g00*v200) - 0.5*(g3R0*v23R0-2*g2R0*v22R0+gR0*v2R0))
                                            - th22*(2-g00)*(v22R0-2*v2R0+v200 - 0.5*(v23R0-2*v22R0+v2R0));
        av1_V100=1;
        //av1_V1R0=-1;
        av2_V200=1;
        B[3*m00+1] = 0;//cas_u1(t_hat,X[i],Y[i],h1,tv);
        B[3*m00+2] = 0;//cas_u2(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*m0L,ag_G0L);
        triplets.emplace_back(3*m00,3*m00+2,ag_V200);
        triplets.emplace_back(3*m00,3*m0L+2,ag_V20L);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        //triplets.emplace_back(3*m00+1,3*m0L+1,av1_V1R0);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);
    }
    else if (st[i]==6) // left down corner (test) flow start
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00];

        ag_G00 = 1;
        av1_V100=1;
        av2_V200=1;

        B[3*m00] = g00;
        B[3*m00+1] = 0;//cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = 0;//cas_u2(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);
    }
    else if (st[i]==7) // right down corner (test)
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*m0R];  g2R0=VEC[3*m0R+3*(m0R-m00)]; g3R0=VEC[3*m0R+6*(m0R-m00)]; //here R0 means 0R in order to save memory
        v200=VEC[3*m00+2]; v2R0=VEC[3*m0R+2]; v22R0=VEC[3*m0R+3*(m0R-m00)+2]; v23R0=VEC[3*m0R+6*(m0R-m00)+2];
        v200=0;
        ag_G00 = 1. - th22*v200;
        ag_G0R = th22*v200;
        //ag_V200=-2*th22;
        //ag_V20R=2*th22;

        B[3*m00] =  gR0;//g00+ th22*g00*(v2R0-v200) + th22*((g2R0*v22R0-2*gR0*v2R0+g00*v200) - 0.5*(g3R0*v23R0-2*g2R0*v22R0+gR0*v2R0))
                          //                  + th22*(2-g00)*(v22R0-2*v2R0+v200 - 0.5*(v23R0-2*v22R0+v2R0));
        av1_V100=1;
        //av1_V1R0=-1;
        av2_V200=1;
        B[3*m00+1] = 0;//cas_u1(t_hat,X[i],Y[i],h1,tv);
        B[3*m00+2] = 0;//cas_u2(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*m0R,ag_G0R);
        triplets.emplace_back(3*m00,3*m00+2,ag_V200);
        triplets.emplace_back(3*m00,3*m0R+2,ag_V20R);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        //triplets.emplace_back(3*m00+1,3*mL0+1,av1_V1R0);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);
    }
    else if (st[i]==8) // left up corner (test) flow start
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00];

        ag_G00 = 1;
        av1_V100=1;
        av2_V200=1;

        B[3*m00] = g00;
        B[3*m00+1] = 0;//cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = 0;//cas_u2(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);
    }
    else if (st[i]==9) // left up corner (test) no flow
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*mR0];  g2R0=VEC[3*mR0+3]; g3R0=VEC[3*mR0+6];
        v100=VEC[3*m00+1]; v1R0=VEC[3*mR0+1]; v12R0=VEC[3*mR0+4]; v13R0=VEC[3*mR0+7];
        v100=0;
        ag_G00 = 1. - th12*v100;
        //ag_GR0 = th12*v1R0;
        //ag_V100=-2*th12;
        //ag_V1R0=2*th12;

        av1_V100=1;
        av2_V200=1;

        B[3*m00] = gR0;// + th12*g00*(v1R0-v100) + th12*((g2R0*v12R0-2*gR0*v1R0+g00*v100) - 0.5*(g3R0*v13R0-2*g2R0*v12R0+gR0*v1R0))
                         //                   + th12*(2-g00)*(v12R0-2*v1R0+v100 - 0.5*(v13R0-2*v12R0+v1R0));
        B[3*m00+1] = 0;//cas_u1(t_hat,X[i],Y[i],h1,tv);
        B[3*m00+2] = 0;//cas_u2(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        //triplets.emplace_back(3*m00,3*mR0,ag_GR0);
        //triplets.emplace_back(3*m00,3*mR0+1,ag_V1R0);
        //triplets.emplace_back(3*m00,3*m00+1,ag_V100);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);
    }
    else if (st[i]==10) // right up corner (test)
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*m0L];  g2R0=VEC[3*m0L-3*(m00-m0L)]; g3R0=VEC[3*m0L-6*(m00-m0L)]; //here R0 means 0R in order to save memory
        v200=VEC[3*m00+2]; v2R0=VEC[3*m0L+2]; v22R0=VEC[3*m0L-3*(m00-m0L)+2]; v23R0=VEC[3*m0L-6*(m00-m0L)+2]; // also here R means L in order to save memory
        v200=0;
        ag_G00 = 1. + th22*v200;
        ag_G0L = -th22*v200;
        //ag_V200=2*th22;
        //ag_V20L=-2*th22;

        B[3*m00] = gR0;//g00 + th22*g00*(v200-v2R0) - th22*((g2R0*v22R0-2*gR0*v2R0+g00*v200) - 0.5*(g3R0*v23R0-2*g2R0*v22R0+gR0*v2R0))
                                            //- th22*(2-g00)*(v22R0-2*v2R0+v200 - 0.5*(v23R0-2*v22R0+v2R0));
        av1_V100=1;
        //av1_V1R0=-1;
        av2_V200=1;
        B[3*m00+1] = 0;//cas_u1(t_hat,X[i],Y[i],h1,tv);
        B[3*m00+2] = 0;//cas_u2(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*m0L,ag_G0L);
        triplets.emplace_back(3*m00,3*m00+2,ag_V200);
        triplets.emplace_back(3*m00,3*m0L+2,ag_V20L);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        //triplets.emplace_back(3*m00+1,3*mL0+1,av1_V1R0);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);
    }
    else if (st[i]==11) // inner left up corner
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*mR0]; g0R=VEC[3*m0R]; gL0=VEC[3*mL0]; g0L=VEC[3*m0L];
        t00=T[m00]; tR0=T[mR0]; tL0=T[mL0]; t0R=T[m0R]; t0L=T[m0L];
        v100=VEC[3*m00+1]; v1R0=VEC[3*mR0+1]; v10R=VEC[3*m0R+1]; v1L0=VEC[3*mL0+1]; v10L=VEC[3*m0L+1]; v1RR=VEC[3*m0R+4]; v1LL=VEC[3*m0L-2]; v1RL=VEC[3*m0L+4]; v1LR=0;
        v200=VEC[3*m00+2]; v2R0=VEC[3*mR0+2]; v20R=VEC[3*m0R+2]; v2L0=VEC[3*mL0+2]; v20L=VEC[3*m0L+2]; v2RR=VEC[3*m0R+5]; v2LL=VEC[3*m0L-1]; v2RL=VEC[3*m0L+5]; v2LR=0;
        v100=0;v200=0;
        ag_G00 = 1.;
        ag_GL0 = -th14*(v100+v1L0);
        ag_GR0 = th14*(v100+v1R0);
        ag_G0L = -th24*(v200+v20L);
        ag_G0R = th24*(v200+v20R);
        ag_V1L0 = -th12;
        ag_V20L = -th22;
        ag_V1R0 = th12;
        ag_V20R = th22;

        /*av1_V100 = 1. + koef*2*mu_C*t4h1_3 + koef*2*mu_C*th2_2;
        av1_V1R0 = th16*(v100 + v1R0) - koef*mu_C*t4h1_3;
        av1_V1L0 = -th16*(v100 + v1L0) - koef*mu_C*t4h1_3;
        av1_V10R = th24*(v200 + v20R) - koef*mu_C*th2_2;
        av1_V10L = -th24*(v200 + v20L) - koef*mu_C*th2_2;
        av1_GL0 = -R*th12*t00;
        av1_GR0 = R*th12*t00;

        av2_V200 = 1. + koef*2*mu_C*th1_2 + koef*2*mu_C*t4h2_3;
        av2_V20R = th26*(v200 + v20R) - koef*mu_C*t4h2_3;
        av2_V20L = -th26*(v200 + v20L) - koef*mu_C*t4h2_3;
        av2_V2R0 = th14*(v100 + v1R0) - koef*mu_C*th1_2;
        av2_V2L0 = -th14*(v100 + v1L0) - koef*mu_C*th1_2;
        av2_G0L = -R*th22*t00;
        av2_G0R = R*th22*t00;*/

        B[3*m00] = g00 + g00*(th14*(v1R0-v1L0) + th24*(v20R-v20L));

        /*B[3*m00+1] = v100 - th12*R*(tR0-tL0) + th24*v100*(v20R-v20L)
                    + koef*(mu*exp(-g00)-mu_C)*(t4h1_3*(v1R0-2*v100+v1L0) + th2_2*(v10R-2*v100+v10L))
                    + koef*tau*mu*exp(-g00)*(v2RR-v2RL-v2LR+v2LL)/(12*h1*h2)
                    ;

        B[3*m00+2] = v200 - th22*R*(t0R-t0L) + th14*v200*(v1R0-v1L0)
                    + koef*(mu*exp(-g00)-mu_C)*(t4h2_3*(v20R-2*v200+v20L) + th1_2*(v2R0-2*v200+v2L0))
                    + koef*tau*mu*exp(-g00)*(v1RR-v1RL-v1LR+v1LL)/(12*h1*h2)
                    ;*/
        av1_V100=1;
        av2_V200=1;
        B[3*m00+1] = v100;
        B[3*m00+2] = v200;
        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*mL0,ag_GL0);
        triplets.emplace_back(3*m00,3*mR0,ag_GR0);
        triplets.emplace_back(3*m00,3*m0L,ag_G0L);
        triplets.emplace_back(3*m00,3*m0R,ag_G0R);
        triplets.emplace_back(3*m00,3*mL0+1,ag_V1L0);
        triplets.emplace_back(3*m00,3*mR0+1,ag_V1R0);
        triplets.emplace_back(3*m00,3*m0L+2,ag_V20L);
        triplets.emplace_back(3*m00,3*m0R+2,ag_V20R);

        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        /*triplets.emplace_back(3*m00+1,3*mR0+1,av1_V1R0);
        triplets.emplace_back(3*m00+1,3*mL0+1,av1_V1L0);
        triplets.emplace_back(3*m00+1,3*m0R+1,av1_V10R);
        triplets.emplace_back(3*m00+1,3*m0L+1,av1_V10L);
        triplets.emplace_back(3*m00+1,3*mL0,av1_GL0);
        triplets.emplace_back(3*m00+1,3*mR0,av1_GR0);*/

        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);
        /*triplets.emplace_back(3*m00+2,3*mR0+2,av2_V2R0);
        triplets.emplace_back(3*m00+2,3*mL0+2,av2_V2L0);
        triplets.emplace_back(3*m00+2,3*m0R+2,av2_V20R);
        triplets.emplace_back(3*m00+2,3*m0L+2,av2_V20L);
        triplets.emplace_back(3*m00+2,3*m0L,av2_G0L);
        triplets.emplace_back(3*m00+2,3*m0R,av2_G0R);*/
    }
  }
  //build matrix, solve with eigen
  matrix.setFromTriplets (triplets.begin(), triplets.end());
  eigen_solver_t solver(matrix);

  if (solver.info() != Eigen::Success)
  {
      fprintf(OUT, "Can not build preconditioner"); fclose(OUT);
      return -1;
  }
  eigen_vector_t ress = solver.solveWithGuess(B, VEC);
  if (solver.info() != Eigen::Success)
  {
      fprintf(OUT, "Failed to solve the system with Eigen, n=%d", n); fclose(OUT);
      n_c[4*ij]=0; n_c[4*ij+1]=0; n_c[4*ij+2]=0; n_c[4*ij+3]=0;
      n_l[4*ij]=0; n_l[4*ij+1]=0; n_l[4*ij+2]=0; n_l[4*ij+3]=0; t_len[ij]=0;
      return 0;
  }
  VEC=ress;

  R*=1000;
  for (int i=0; i<Dim; i++)
  {
    if (kappa*exp(-VEC[3*i])>kap_C) kap_C=kappa*exp(-VEC[3*i]);
  }
  // solve matrixt*T=Bt
  for (int i=0; i<Dim; i++)
  {
    // fill matrixt
    if (st[i]==0) // inner knots
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        //m0R=i+M1+1;//M0R[i];
        //m0L=i-M1-1;//M0L[i];
        m0R=M0R[i];
        m0L=M0L[i];

        g00=VEC[3*m00];
        t00=T[m00]; tR0=T[mR0]; tL0=T[mL0]; t0R=T[m0R]; t0L=T[m0L];
        v100=VEC[3*m00+1]; v1R0=VEC[3*mR0+1]; v10R=VEC[3*m0R+1]; v1L0=VEC[3*mL0+1]; v10L=VEC[3*m0L+1]; // try to take 3*i+3*r+1 instead of 3*i+3*r+4 & 3*i-3*l+1 instead of 3*i-3*l-2
        v200=VEC[3*m00+2]; v2R0=VEC[3*mR0+2]; v20R=VEC[3*m0R+2]; v2L0=VEC[3*mL0+2]; v20L=VEC[3*m0L+2];

        a_T00 = cv + 2*kap_C*th1_2 + 2*kap_C*th2_2;
        a_TR0 = cv*th14*(v100+v1R0) - kap_C*th1_2;
        a_TL0 = -cv*th14*(v100+v1L0) - kap_C*th1_2;
        a_T0R = cv*th24*(v200+v20R) - kap_C*th2_2;
        a_T0L = -cv*th24*(v200+v20L) - kap_C*th2_2;

        Bt[m00] = cv*t00 + cv*t00*(th14*(v1R0-v1L0)+th24*(v20R-v20L)) + (kappa*exp(-g00)-kap_C)*(th1_2*(tR0-2*t00+tL0)+th2_2*(t0R-2*t00+t0L))
                                                       - R*t00*(th12*(v1R0-v1L0)+th22*(v20R-v20L))
                                                       - (1./6)*mu*exp(-g00)*(th1_2*(v1R0-v1L0)*(v1R0-v1L0)+th2_2*(v20R-v20L)*(v20R-v20L))
                                                       - (tau/(3*h1*h2))*mu*exp(-g00)*(v1R0-v1L0)*(v20R-v20L)
                                                       + 0.5*mu*exp(-g00) * (th2_2*(v20R-v20L)*(v20R-v20L) + th1_2*(v1R0-v1L0)*(v1R0-v1L0))
                                                       + 0.25*tau*mu*exp(-g00) * ((v2R0-v2L0)/h1+(v10R-v10L)/h2) * ((v2R0-v2L0)/h1+(v10R-v10L)/h2)
                                                       ;//+ tau*cas_ft(t_hat,X[i],Y[i],R,mu,cv,kappa);
        triplets_t.emplace_back(m00,m00,a_T00);
        triplets_t.emplace_back(m00,mR0,a_TR0);
        triplets_t.emplace_back(m00,mL0,a_TL0);
        triplets_t.emplace_back(m00,m0R,a_T0R);
        triplets_t.emplace_back(m00,m0L,a_T0L);
    }
    else if (st[i]==1) // left border flow start
    {
        m00=i;
        a_T00=1;
        Bt[m00] = cas_theta(t_hat,X[i],Y[i],h1,tt);
        triplets_t.emplace_back(m00,m00,a_T00); // matrixt is Dim*Dim, T is Dim, i=0,..Dim-1; first coord is matrixt's string, second coord is matrixt's row
    }
    else if (st[i]==2) // left border no flow
    {
        m00=i; mR0=i+1; //-dtheta/dx1
        a_T00=-1;
        a_TR0=1;
        Bt[m00] = 0;//cas_theta(t_hat,X[i],Y[i],h1,tt);
        triplets_t.emplace_back(m00,m00,a_T00);
        triplets_t.emplace_back(m00,mR0,a_TR0);
    }
    else if (st[i]==3) // right border
    {
        m00=i; mL0=i-1;
        a_T00=1;
        a_TR0=-1;
        Bt[m00] = 0;//cas_theta(t_hat,X[i],Y[i],h1,tt);
        triplets_t.emplace_back(m00,m00,a_T00);
        triplets_t.emplace_back(m00,mL0,a_TR0);
    }
    else if (st[i]==4) // down border
    {
        m00=i; m0R=M0R[i];
        a_T00=-1;
        a_TR0=1;
        Bt[m00] = 0;//cas_theta(t_hat,X[i],Y[i],h1,tt);
        triplets_t.emplace_back(m00,m00,a_T00);
        triplets_t.emplace_back(m00,m0R,a_TR0);
    }
    else if (st[i]==5) // up border
    {
        m00=i; m0L=M0L[i];
        a_T00=1;
        a_TR0=-1;
        Bt[m00] = 0;//cas_theta(t_hat,X[i],Y[i],h1,tt);
        triplets_t.emplace_back(m00,m00,a_T00);
        triplets_t.emplace_back(m00,m0L,a_TR0);
    }
    else if (st[i]==6) // left down corner (test)
    {
        m00=i; m0R=M0R[i];
        a_T00=-1;
        a_TR0=1;
        Bt[m00] = 0;//cas_theta(t_hat,X[i],Y[i],h1,tt);
        triplets_t.emplace_back(m00,m00,a_T00);
        triplets_t.emplace_back(m00,m0R,a_TR0);
    }
    else if (st[i]==7) // right down corner (test)
    {
        m00=i; m0R=M0R[i]; t00=T[m00];
        a_T00=1;
        //a_TR0=1;
        Bt[m00] = t00;
        triplets_t.emplace_back(m00,m00,a_T00);
        //triplets_t.emplace_back(m00,m0R,a_TR0);
    }
    else if (st[i]==8) // left up corner (test) flow start
    {
        m00=i; m0L=M0L[i];
        a_T00=1;
        a_TR0=-1;
        Bt[m00] = 0;//cas_theta(t_hat,X[i],Y[i],h1,tt);
        triplets_t.emplace_back(m00,m00,a_T00);
        triplets_t.emplace_back(m00,m0L,a_TR0);
    }
    else if (st[i]==9) // left up corner (test) no flow
    {
        m00=i; m0L=M0L[i];
        a_T00=1;
        a_TR0=-1;
        Bt[m00] = 0;//cas_theta(t_hat,X[i],Y[i],h1,tt);
        triplets_t.emplace_back(m00,m00,a_T00);
        triplets_t.emplace_back(m00,m0L,a_TR0);
    }
    else if (st[i]==10) // right up corner (test)
    {
        m00=i; m0L=M0L[i]; t00=T[m00];
        a_T00=1;
        //a_TR0=-1;
        Bt[m00] = t00;
        triplets_t.emplace_back(m00,m00,a_T00);
        //triplets_t.emplace_back(m00,m0L,a_TR0);
    }
    else if (st[i]==11) // inner left up corner
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00];
        t00=T[m00]; tR0=T[mR0]; t0R=T[m0R]; tL0=T[mL0]; t0L=T[m0L];
        v100=VEC[3*m00+1]; v1R0=VEC[3*mR0+1]; v10R=VEC[3*m0R+1]; v1L0=VEC[3*mL0+1]; v10L=VEC[3*m0L+1];
        v200=VEC[3*m00+2]; v2R0=VEC[3*mR0+2]; v20R=VEC[3*m0R+2]; v2L0=VEC[3*mL0+2]; v20L=VEC[3*m0L+2];
        v100=0;
        v200=0;
        a_T00 = cv + 2*kap_C*th1_2 + 2*kap_C*th2_2;
        a_TR0 = cv*th14*(v100+v1R0) - kap_C*th1_2;
        a_TL0 = -cv*th14*(v100+v1L0) - kap_C*th1_2;
        a_T0R = cv*th24*(v200+v20R) - kap_C*th2_2;
        a_T0L = -cv*th24*(v200+v20L) - kap_C*th2_2;

        Bt[m00] = cv*t00 + cv*t00*(th14*(v1R0-v1L0)+th24*(v20R-v20L)) + (kappa*exp(-g00)-kap_C)*(th1_2*(tR0-2*t00+tL0)+th2_2*(t0R-2*t00+t0L))
                                                       - R*t00*(th12*(v1R0-v1L0)+th22*(v20R-v20L))
                                                       - (1./6)*mu*exp(-g00)*(th1_2*(v1R0-v1L0)*(v1R0-v1L0)+th2_2*(v20R-v20L)*(v20R-v20L))
                                                       - (tau/(3*h1*h2))*mu*exp(-g00)*(v1R0-v1L0)*(v20R-v20L)
                                                       + 0.5*mu*exp(-g00) * (th2_2*(v20R-v20L)*(v20R-v20L) + th1_2*(v1R0-v1L0)*(v1R0-v1L0))
                                                       + 0.25*tau*mu*exp(-g00) * ((v2R0-v2L0)/h1+(v10R-v10L)/h2) * ((v2R0-v2L0)/h1+(v10R-v10L)/h2);
        triplets_t.emplace_back(m00,m00,a_T00);
        triplets_t.emplace_back(m00,mR0,a_TR0);
        triplets_t.emplace_back(m00,mL0,a_TL0);
        triplets_t.emplace_back(m00,m0R,a_T0R);
        triplets_t.emplace_back(m00,m0L,a_T0L);
    }
  }
  //build matrixt, solve with eigen
  matrixt.setFromTriplets (triplets_t.begin(), triplets_t.end());
  eigen_solver_t solvert(matrixt);

  if (solvert.info() != Eigen::Success)
  {
      fprintf(OUT, "Can not build preconditioner"); fclose(OUT);
      return -1;
  }
  eigen_vector_t res (Dim);
  res = solvert.solveWithGuess(Bt, T);
  if (solvert.info() != Eigen::Success)
  {
      fprintf(OUT, "Failed to solve the system with Eigen, n=%d", n); fclose(OUT);
      n_c[4*ij]=0; n_c[4*ij+1]=0; n_c[4*ij+2]=0; n_c[4*ij+3]=0;
      n_l[4*ij]=0; n_l[4*ij+1]=0; n_l[4*ij+2]=0; n_l[4*ij+3]=0; t_len[ij]=0;
      return 0;
  }
  T=res;
  /*if (n==250 || n==500 || n==1000 || n==1500 || n==2000 || n==2500 || n==3000 || n==3500 || n==4000 || n==4500 || n==5000 || n==6000 || n==7000 || n==8000 || n==9000 || n==10000 || n==15000 || n==20000 || n==25000 || n==30000 || n==40000 || n==50000) {ress_st=ress; res_st=res; dif=-1;}
  if (n==260 || n==520 || n==1020 || n==1520 || n==2020 || n==2520 || n==3020 || n==3520 || n==4020 || n==4520 || n==5020 || n==6020 || n==7020 || n==8020 || n==9020 || n==10020 || n==15020 || n==20020 || n==25020 || n==30020 || n==40020 || n==50020)
  {
      for (int i=0; i<Dim; i++)
      {
          if (dif<fabs(VEC[3*i]-ress_st[3*i])) dif=fabs(VEC[3*i]-ress_st[3*i]);
          if (dif<fabs(VEC[3*i+1]-ress_st[3*i+1])) dif=fabs(VEC[3*i+1]-ress_st[3*i+1]);
          if (dif<fabs(VEC[3*i+2]-ress_st[3*i+2])) dif=fabs(VEC[3*i+2]-ress_st[3*i+2]);
          if (dif<fabs(T[i]-res_st[i])) dif=fabs(T[i]-res_st[i]);
      }
      //char file1[20];    sprintf(file1, "%d", n);    strcat(file1,txt); OUT1 = fopen(file1, "w");   fprintf(OUT1, "%lf\n",dif);
        //                                                                fclose(OUT1);
  }*/
  //maxg=-1; maxu1=-1; maxu2=-1; maxt=-1;
  t_hat+=tau; mu_C=0; kap_C=0; n++; t_len[ij]=n*tau;
  if ((n-1)%10==0)
  {
      char file1[20];    sprintf(file1, "%d", n);    strcat(file1,txt); OUT1 = fopen(file1, "w");   for (int i=0; i<Dim; i++) //fprintf(OUT1, "%lf\n",VEC[3*i]);
                                                                        fprintf(OUT1, "%lf\n",T[i]);
                                                                        fclose(OUT1);
  }
}
// solve end =========================================================================================================
fprintf(OUT, "1");

//normi
maxg=0; maxu1=0; maxu2=0; maxt=0;
sumu1=0; sumu2=0; sumg=0; sumt=0;
for (int i=0; i<Dim; i++)
{
    if (maxg<fabs(VEC[3*i])) maxg=fabs(VEC[3*i]);
    if (maxu1<fabs(VEC[3*i+1])) maxu1=fabs(VEC[3*i+1]);
    if (maxu2<fabs(VEC[3*i+2])) maxu2=fabs(VEC[3*i+2]);
    if (maxt<fabs(T[i])) {maxt=fabs(T[i]);}// fprintf(OUT, " %d/%d ", i,st[i]);}
    sumg+=(VEC[3*i])*(VEC[3*i]);
    sumu1+=(VEC[3*i+1])*(VEC[3*i+1]);
    sumu2+=(VEC[3*i+2])*(VEC[3*i+2]);
    sumt+=(T[i])*(T[i]);
}
n_c[4*ij]=maxu1; n_c[4*ij+1]=maxu2; n_c[4*ij+2]=maxg; n_c[4*ij+3]=maxt;
n_l[4*ij]=sqrt(h1*h2*sumu1); n_l[4*ij+1]=sqrt(h1*h2*sumu2); n_l[4*ij+2]=sqrt(h1*h2*sumg); n_l[4*ij+3]=sqrt(h1*h2*sumt);
fclose(OUT);
return 0;
}

int Sxema_hcond_vgas_2_flowp_test(double mu, int N,int M1, int M2, int Dim,
                      Str_param_gas *str_param_gas, Str_setka *str_setka, P_she *p_she,
                      int *st, double *X, double *Y, int *M0R, int *M0L,
                      double *n_c, double *n_l, double *n_w, int ij, int K, int cas, double *t_len, double tv, double tg, double tt)
{
printf("hcond_vgas_2_flowp_test \n");

double  h1, h2, tau, R, kappa, cv, eps;
tau = str_param_gas->T/N;//p_she->tau;
h1 = str_param_gas->X1/M1;//p_she->h1;
h2 = str_param_gas->X2/M2;//p_she->h2;
R = str_param_gas->R;
kappa = str_param_gas->kappa;
cv = str_param_gas->cv;
eps = str_param_gas->eps;

double EX2 = str_param_gas->EX2;
// assigning functions
double (*cas_u1)(double , double , double, double , double );
double (*cas_u2)(double , double , double );
double (*cas_g)(double , double , double, double , double );
double (*cas_theta)(double , double , double, double , double );

if (cas==4)// for DEBUGGING TEST 1-2 functions u1, u2, g and theta are used
    {
        cas_u1=u1_1;
        cas_u2=u2_1;
        cas_g=g_1;
        cas_theta=theta_1;
    }

//declare vectors
eigen_vector_t VEC (3*Dim);
eigen_vector_t T (Dim);
//VEC consists of g, v1, v2: VEC = (G_00, V1_00,V2_00, G_10, V1_10, V2_10, ..., G_M10, V1_M10, V2_M10, G_01, V1_01, V2_01, G_11, V1_11, V2_11,...)
//knots' coord
int m00, mR0, m0R, mL0, m0L;
//initial conditions---------------------------------------
for (int i=0; i<Dim; i++)
{
    VEC[3*i] = cas_g(0,X[i],Y[i],h1,tg);
    VEC[3*i+1] = cas_u1(0,X[i],Y[i],h1,tv);
    VEC[3*i+2] = cas_u2(0,X[i],Y[i]);
    T[i] = cas_theta(0,X[i],Y[i],h1,tt);
}

//const
double th16=tau/(6*h1);
double th26=tau/(6*h2);
double th14=tau/(4*h1);
double th24=tau/(4*h2);
double th12=tau/(2*h1);
double th22=tau/(2*h2);
double t4h1_3=4*tau/(3*h1*h1);
double t4h2_3=4*tau/(3*h2*h2);
double th1_2=tau/(h1*h1);
double th2_2=tau/(h2*h2);

char txt[] = ".txt";
char file[20];
sprintf(file, "%d", ij);
strcat(file,txt);
FILE *OUT;
OUT = fopen(file, "w");
fprintf(OUT, "1 ");
FILE *OUT1;
double koef=1.;
double mu_C=0, kap_C=0;
double t_hat=tau;
//function values in the grid konts
double g00=0, gR0=0, g0R=0, gL0=0, g0L=0;
double t00=0, tR0=0, t0R=0, tL0=0, t0L=0;
double v100=0, v1R0=0, v10R=0, v1L0=0, v10L=0, v1RR=0, v1LL=0, v1RL=0, v1LR=0;
double v200=0, v2R0=0, v20R=0, v2L0=0, v20L=0, v2RR=0, v2LL=0, v2RL=0, v2LR=0;
double g2R0=0, g3R0=0, v12R0=0, v13R0=0, v22R0=0, v23R0=0; // extra knots values for border equations
//matrix coefficients for corresponding components
double ag_G00=0, ag_GR0=0, ag_GL0=0, ag_G0R=0, ag_G0L=0, ag_V100=0, ag_V200=0, ag_V1R0=0, ag_V1L0=0, ag_V20R=0, ag_V20L=0;
double av1_GR0=0, av1_GL0=0, av1_V100=0, av1_V1R0=0, av1_V1L0=0, av1_V10R=0, av1_V10L=0, av1_V20R=0, av1_V20L=0;
double av2_G0R=0, av2_G0L=0, av2_V200=0, av2_V2R0=0, av2_V2L0=0, av2_V20R=0, av2_V20L=0, av2_V1R0=0, av2_V1L0=0;
double a_T00=0, a_TR0=0, a_TL0=0, a_T0R=0, a_T0L=0;

//solve=========================================================================================================================================
int n=0; double dif=1;
double maxg0=-1, maxu10=-1, maxu20=-1, maxt0=-1;
double maxg=-1, maxu1=-1, maxu2=-1, maxt=-1;
double sumu1=0, sumu2=0, sumg=0, sumt=0; double eta=0.1;
eigen_vector_t ress_st (3*Dim);
eigen_vector_t res_st (Dim);
while (fabs(dif)>eps && n<50100)
{
    dif=-1;
    //declare matrix and right part
eigen_matrix_t matrix (3*Dim, 3*Dim);
eigen_matrix_t matrixt (Dim, Dim);
eigen_vector_t B (3*Dim);
eigen_vector_t Bt (Dim);

std::vector<eigen_triplet_t> triplets;
std::vector<eigen_triplet_t> triplets_t;
// firstly solve: matrix*VEC = B
// then solve: matrixt*T = Bt
  for (int i=0; i<Dim; i++)
  {
    if (mu*exp(-VEC[3*i])>mu_C) mu_C=mu*exp(-VEC[3*i]);
  }

    R*=0.001;
    // solve matrix*VEC=B
  for (int i=0; i<Dim; i++)
  {
    //fill matrix
    if (st[i]==0) // inner knots
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*mR0]; g0R=VEC[3*m0R]; gL0=VEC[3*mL0]; g0L=VEC[3*m0L];
        t00=T[m00]; tR0=T[mR0]; tL0=T[mL0]; t0R=T[m0R]; t0L=T[m0L];
        v100=VEC[3*m00+1]; v1R0=VEC[3*mR0+1]; v10R=VEC[3*m0R+1]; v1L0=VEC[3*mL0+1]; v10L=VEC[3*m0L+1]; v1RR=VEC[3*m0R+4]; v1LL=VEC[3*m0L-2]; v1RL=VEC[3*m0L+4]; v1LR=VEC[3*m0R-2];
        v200=VEC[3*m00+2]; v2R0=VEC[3*mR0+2]; v20R=VEC[3*m0R+2]; v2L0=VEC[3*mL0+2]; v20L=VEC[3*m0L+2]; v2RR=VEC[3*m0R+5]; v2LL=VEC[3*m0L-1]; v2RL=VEC[3*m0L+5]; v2LR=VEC[3*m0R-1];

        ag_G00 = 1. + 2*eta*tau*th1_2 + 2*eta*tau*th2_2;
        ag_GL0 = -th14*(v100+v1L0) - eta*tau*th1_2;
        ag_GR0 = th14*(v100+v1R0) - eta*tau*th1_2;
        ag_G0L = -th24*(v200+v20L) - eta*tau*th2_2;
        ag_G0R = th24*(v200+v20R) - eta*tau*th2_2;
        ag_V1L0 = -th12;
        ag_V20L = -th22;
        ag_V1R0 = th12;
        ag_V20R = th22;

        av1_V100 = 1. + koef*2*mu_C*t4h1_3 + koef*2*mu_C*th2_2;
        av1_V1R0 = th16*(v100 + v1R0) - koef*mu_C*t4h1_3;
        av1_V1L0 = -th16*(v100 + v1L0) - koef*mu_C*t4h1_3;
        av1_V10R = th24*(v200 + v20R) - koef*mu_C*th2_2;
        av1_V10L = -th24*(v200 + v20L) - koef*mu_C*th2_2;
        av1_GL0 = -R*th12*t00;
        av1_GR0 = R*th12*t00;

        av2_V200 = 1. + koef*2*mu_C*th1_2 + koef*2*mu_C*t4h2_3;
        av2_V20R = th26*(v200 + v20R) - koef*mu_C*t4h2_3;
        av2_V20L = -th26*(v200 + v20L) - koef*mu_C*t4h2_3;
        av2_V2R0 = th14*(v100 + v1R0) - koef*mu_C*th1_2;
        av2_V2L0 = -th14*(v100 + v1L0) - koef*mu_C*th1_2;
        av2_G0L = -R*th22*t00;
        av2_G0R = R*th22*t00;

        B[3*m00] = g00 + g00*(th14*(v1R0-v1L0) + th24*(v20R-v20L));// + tau*cas_fg(t_hat, X[i],Y[i]);

        B[3*m00+1] = v100 - th12*R*(tR0-tL0) + th24*v100*(v20R-v20L)
                    + koef*(mu*exp(-g00)-mu_C)*(t4h1_3*(v1R0-2*v100+v1L0) + th2_2*(v10R-2*v100+v10L))
                    + koef*tau*mu*exp(-g00)*(v2RR-v2RL-v2LR+v2LL)/(12*h1*h2)
                    ;//+ tau*cas_fu1(t_hat, X[i],Y[i], R,mu, koef);

        B[3*m00+2] = v200 - th22*R*(t0R-t0L) + th14*v200*(v1R0-v1L0)
                    + koef*(mu*exp(-g00)-mu_C)*(t4h2_3*(v20R-2*v200+v20L) + th1_2*(v2R0-2*v200+v2L0))
                    + koef*tau*mu*exp(-g00)*(v1RR-v1RL-v1LR+v1LL)/(12*h1*h2)
                    ;//+ tau*cas_fu2(t_hat, X[i],Y[i], R,mu, koef);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*mL0,ag_GL0);
        triplets.emplace_back(3*m00,3*mR0,ag_GR0);
        triplets.emplace_back(3*m00,3*m0L,ag_G0L);
        triplets.emplace_back(3*m00,3*m0R,ag_G0R);
        triplets.emplace_back(3*m00,3*mL0+1,ag_V1L0);
        triplets.emplace_back(3*m00,3*mR0+1,ag_V1R0);
        triplets.emplace_back(3*m00,3*m0L+2,ag_V20L);
        triplets.emplace_back(3*m00,3*m0R+2,ag_V20R);

        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+1,3*mR0+1,av1_V1R0);
        triplets.emplace_back(3*m00+1,3*mL0+1,av1_V1L0);
        triplets.emplace_back(3*m00+1,3*m0R+1,av1_V10R);
        triplets.emplace_back(3*m00+1,3*m0L+1,av1_V10L);
        triplets.emplace_back(3*m00+1,3*mL0,av1_GL0);
        triplets.emplace_back(3*m00+1,3*mR0,av1_GR0);

        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);
        triplets.emplace_back(3*m00+2,3*mR0+2,av2_V2R0);
        triplets.emplace_back(3*m00+2,3*mL0+2,av2_V2L0);
        triplets.emplace_back(3*m00+2,3*m0R+2,av2_V20R);
        triplets.emplace_back(3*m00+2,3*m0L+2,av2_V20L);
        triplets.emplace_back(3*m00+2,3*m0L,av2_G0L);
        triplets.emplace_back(3*m00+2,3*m0R,av2_G0R);
    }
    else if (st[i]==1) // left border flow start
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*mR0];  g2R0=VEC[3*mR0+3]; g3R0=VEC[3*mR0+6];
        v100=VEC[3*m00+1]; v1R0=VEC[3*mR0+1]; v12R0=VEC[3*mR0+4]; v13R0=VEC[3*mR0+7];
        //v100=0;
        ag_G00 = 1.;// - th12*v100;
        ag_GR0 = 0;//th12*v1R0;
        ag_V1R0=0;//2*th12;

        av1_V100=1;
        av2_V200=1;

        B[3*m00] = cas_g(t_hat,X[i],Y[i],h1,tg);//g00 + th12*g00*(v1R0-v100) + th12*((g2R0*v12R0-2*gR0*v1R0+g00*v100) - 0.5*(g3R0*v13R0-2*g2R0*v12R0+gR0*v1R0))
                                            //+ th12*(2-g00)*(v12R0-2*v1R0+v100 - 0.5*(v13R0-2*v12R0+v1R0));
        B[3*m00+1] = cas_u1(t_hat,X[i],Y[i],h1,tv);
        B[3*m00+2] = 0;//cas_u2(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*mR0,ag_GR0);
        triplets.emplace_back(3*m00,3*mR0+1,ag_V1R0);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);

    }
    else if (st[i]==2) // right border
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*mL0]; g2R0=VEC[3*mL0-3]; g3R0=VEC[3*mL0-6];
        v100=VEC[3*m00+1]; v1R0=VEC[3*mL0+1]; v12R0=VEC[3*mL0-2]; v13R0=VEC[3*mL0-5]; //here R means L in order to save memory
        //v100=0;
        ag_G00 = 1. + 2*th12*v100 - eta*tau*th1_2;
        ag_GL0 = -2*th12*v100 + eta*tau*th1_2;
        //ag_V1L0=-2*th12;

        B[3*m00] = g00;// + th12*g00*(v100-v1R0) - th12*((g00*v100-2*gR0*v1R0+g2R0*v12R0) - 0.5*(gR0*v1R0-2*g2R0*v12R0+g3R0*v13R0))
                                            //- th12*(2-g00)*(v100-2*v1R0+v12R0 - 0.5*(v1R0-2*v12R0+v13R0));
        av1_V100=1;
        av1_V1R0=-1;
        av2_V200=1;
        B[3*m00+1] = 0;//cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = 0;//cas_u2(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*mL0,ag_GL0);
        //triplets.emplace_back(3*m00,3*mL0+1,ag_V1L0);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+1,3*mL0+1,av1_V1R0);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);

    }
    else if (st[i]==3) // down border
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*m0R];  g2R0=VEC[3*m0R+3*(m0R-m00)]; g3R0=VEC[3*m0R+6*(m0R-m00)]; //here R0 means 0R in order to save memory
        v200=VEC[3*m00+2]; v2R0=VEC[3*m0R+2]; v22R0=VEC[3*m0R+3*(m0R-m00)+2]; v23R0=VEC[3*m0R+6*(m0R-m00)+2];
        v200=0;
        ag_G00 = 1. - th22*v200 + eta*tau*th2_2;
        ag_G0R = th22*v2R0 - eta*tau*th2_2;
        ag_V200=-2*th22;
        ag_V20R=2*th22;

        B[3*m00] = g00 + th22*g00*(v2R0-v200) + th22*((g2R0*v22R0-2*gR0*v2R0+g00*v200) - 0.5*(g3R0*v23R0-2*g2R0*v22R0+gR0*v2R0))
                                            + th22*(2-g00)*(v22R0-2*v2R0+v200 - 0.5*(v23R0-2*v22R0+v2R0));
        av1_V100=1;
        //av1_V1R0=1;
        av2_V200=1;
        B[3*m00+1] = 0;//cas_u1(t_hat,X[i],Y[i],h1,tv);
        B[3*m00+2] = 0;//cas_u2(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*m0R,ag_G0R);
        triplets.emplace_back(3*m00,3*m00+2,ag_V200);
        triplets.emplace_back(3*m00,3*m0R+2,ag_V20R);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        //triplets.emplace_back(3*m00+1,3*m0R+1,av1_V1R0);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);
    }
    else if (st[i]==4) // up border
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*m0L];  g2R0=VEC[3*m0L-3*(m00-m0L)]; g3R0=VEC[3*m0L-6*(m00-m0L)]; //here R0 means 0R in order to save memory
        v200=VEC[3*m00+2]; v2R0=VEC[3*m0L+2]; v22R0=VEC[3*m0L-3*(m00-m0L)+2]; v23R0=VEC[3*m0L-6*(m00-m0L)+2]; // also here R means L in order to save memory
        v200=0;
        ag_G00 = 1. + th22*v200 - eta*tau*th2_2;
        ag_G0L = -th22*v2R0 + eta*tau*th2_2;
        ag_V200=2*th22;
        ag_V20L=-2*th22;

        B[3*m00] = g00 + th22*g00*(v200-v2R0) - th22*((g2R0*v22R0-2*gR0*v2R0+g00*v200) - 0.5*(g3R0*v23R0-2*g2R0*v22R0+gR0*v2R0))
                                            - th22*(2-g00)*(v22R0-2*v2R0+v200 - 0.5*(v23R0-2*v22R0+v2R0));
        av1_V100=1;
        //av1_V1R0=-1;
        av2_V200=1;
        B[3*m00+1] = 0;//cas_u1(t_hat,X[i],Y[i],h1,tv);
        B[3*m00+2] = 0;//cas_u2(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*m0L,ag_G0L);
        triplets.emplace_back(3*m00,3*m00+2,ag_V200);
        triplets.emplace_back(3*m00,3*m0L+2,ag_V20L);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        //triplets.emplace_back(3*m00+1,3*m0L+1,av1_V1R0);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);
    }
    else if (st[i]==5) // left down corner (test)
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00];

        ag_G00 = 1;
        av1_V100=1;
        av2_V200=1;

        B[3*m00] = g00;
        B[3*m00+1] = 0;//cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = 0;//cas_u2(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);
    }
    else if (st[i]==6) // right down corner (test)
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*m0R];  g2R0=VEC[3*m0R+3*(m0R-m00)]; g3R0=VEC[3*m0R+6*(m0R-m00)]; //here R0 means 0R in order to save memory
        v200=VEC[3*m00+2]; v2R0=VEC[3*m0R+2]; v22R0=VEC[3*m0R+3*(m0R-m00)+2]; v23R0=VEC[3*m0R+6*(m0R-m00)+2];
        //v200=0;
        ag_G00 = 1. - th22*v200;
        ag_G0R = th22*v200;
        //ag_V200=-2*th22;
        //ag_V20R=2*th22;

        B[3*m00] = gR0;// + th22*g00*(v2R0-v200) + th22*((g2R0*v22R0-2*gR0*v2R0+g00*v200) - 0.5*(g3R0*v23R0-2*g2R0*v22R0+gR0*v2R0))
                         //                   + th22*(2-g00)*(v22R0-2*v2R0+v200 - 0.5*(v23R0-2*v22R0+v2R0));
        av1_V100=1;
        //av1_V1R0=-1;
        av2_V200=1;
        B[3*m00+1] = 0;//cas_u1(t_hat,X[i],Y[i],h1,tv);
        B[3*m00+2] = 0;//cas_u2(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*m0R,ag_G0R);
        triplets.emplace_back(3*m00,3*m00+2,ag_V200);
        triplets.emplace_back(3*m00,3*m0R+2,ag_V20R);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        //triplets.emplace_back(3*m00+1,3*mL0+1,av1_V1R0);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);
    }
    else if (st[i]==7) // left up corner (test)
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00];

        ag_G00 = 1;
        av1_V100=1;
        av2_V200=1;

        B[3*m00] = g00;
        B[3*m00+1] = 0;//cas_u1(t_hat,X[i],Y[i]);
        B[3*m00+2] = 0;//cas_u2(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);
    }
    else if (st[i]==8) // right up corner (test)
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        m0R=M0R[i];
        m0L=M0L[i];
        g00=VEC[3*m00]; gR0=VEC[3*m0L];  g2R0=VEC[3*m0L-3*(m00-m0L)]; g3R0=VEC[3*m0L-6*(m00-m0L)]; //here R0 means 0R in order to save memory
        v200=VEC[3*m00+2]; v2R0=VEC[3*m0L+2]; v22R0=VEC[3*m0L-3*(m00-m0L)+2]; v23R0=VEC[3*m0L-6*(m00-m0L)+2]; // also here R means L in order to save memory
        //v200=0;
        ag_G00 = 1. + th22*v200;
        ag_G0L = -th22*v200;
        //ag_V200=2*th22;
        //ag_V20L=-2*th22;

        B[3*m00] = gR0;// + th22*g00*(v200-v2R0) - th22*((g2R0*v22R0-2*gR0*v2R0+g00*v200) - 0.5*(g3R0*v23R0-2*g2R0*v22R0+gR0*v2R0))
                         //                   - th22*(2-g00)*(v22R0-2*v2R0+v200 - 0.5*(v23R0-2*v22R0+v2R0));
        av1_V100=1;
        //av1_V1R0=-1;
        av2_V200=1;
        B[3*m00+1] = 0;//cas_u1(t_hat,X[i],Y[i],h1,tv);
        B[3*m00+2] = 0;//cas_u2(t_hat,X[i],Y[i]);

        triplets.emplace_back(3*m00,3*m00,ag_G00);
        triplets.emplace_back(3*m00,3*m0L,ag_G0L);
        triplets.emplace_back(3*m00,3*m00+2,ag_V200);
        triplets.emplace_back(3*m00,3*m0L+2,ag_V20L);
        triplets.emplace_back(3*m00+1,3*m00+1,av1_V100);
        //triplets.emplace_back(3*m00+1,3*mL0+1,av1_V1R0);
        triplets.emplace_back(3*m00+2,3*m00+2,av2_V200);
    }
  }
  //build matrix, solve with eigen
  matrix.setFromTriplets (triplets.begin(), triplets.end());
  eigen_solver_t solver(matrix);

  if (solver.info() != Eigen::Success)
  {
      fprintf(OUT, "Can not build preconditioner");
      return -1;
  }
  eigen_vector_t ress = solver.solveWithGuess(B, VEC);
  if (solver.info() != Eigen::Success)
  {
      fprintf(OUT, "Failed to solve the system with Eigen, n=%d", n); fclose(OUT);
      n_c[4*ij]=0; n_c[4*ij+1]=0; n_c[4*ij+2]=0; n_c[4*ij+3]=0;
      n_l[4*ij]=0; n_l[4*ij+1]=0; n_l[4*ij+2]=0; n_l[4*ij+3]=0; t_len[ij]=0;
      return 0;
  }
  VEC=ress;

  R*=1000;
  for (int i=0; i<Dim; i++)
  {
    if (kappa*exp(-VEC[3*i])>kap_C) kap_C=kappa*exp(-VEC[3*i]);
  }
  // solve matrixt*T=Bt
  for (int i=0; i<Dim; i++)
  {
    //T[i] = cas_theta(t_hat,X[i],Y[i]);
    // fill matrixt
    if (st[i]==0) // inner knots
    {
        m00=i;
        mR0=i+1;
        mL0=i-1;
        //m0R=i+M1+1;//M0R[i];
        //m0L=i-M1-1;//M0L[i];
        m0R=M0R[i];
        m0L=M0L[i];

        g00=VEC[3*m00];
        t00=T[m00]; tR0=T[mR0]; tL0=T[mL0]; t0R=T[m0R]; t0L=T[m0L];
        v100=VEC[3*m00+1]; v1R0=VEC[3*mR0+1]; v10R=VEC[3*m0R+1]; v1L0=VEC[3*mL0+1]; v10L=VEC[3*m0L+1]; // try to take 3*i+3*r+1 instead of 3*i+3*r+4 & 3*i-3*l+1 instead of 3*i-3*l-2
        v200=VEC[3*m00+2]; v2R0=VEC[3*mR0+2]; v20R=VEC[3*m0R+2]; v2L0=VEC[3*mL0+2]; v20L=VEC[3*m0L+2];

        a_T00 = cv + 2*kap_C*th1_2 + 2*kap_C*th2_2;
        a_TR0 = cv*th14*(v100+v1R0) - kap_C*th1_2;
        a_TL0 = -cv*th14*(v100+v1L0) - kap_C*th1_2;
        a_T0R = cv*th24*(v200+v20R) - kap_C*th2_2;
        a_T0L = -cv*th24*(v200+v20L) - kap_C*th2_2;

        Bt[m00] = cv*t00 + cv*t00*(th14*(v1R0-v1L0)+th24*(v20R-v20L)) + (kappa*exp(-g00)-kap_C)*(th1_2*(tR0-2*t00+tL0)+th2_2*(t0R-2*t00+t0L))
                                                       - R*t00*(th12*(v1R0-v1L0)+th22*(v20R-v20L))
                                                       - (1./6)*mu*exp(-g00)*(th1_2*(v1R0-v1L0)*(v1R0-v1L0)+th2_2*(v20R-v20L)*(v20R-v20L))
                                                       - (tau/(3*h1*h2))*mu*exp(-g00)*(v1R0-v1L0)*(v20R-v20L)
                                                       + 0.5*mu*exp(-g00) * (th2_2*(v20R-v20L)*(v20R-v20L) + th1_2*(v1R0-v1L0)*(v1R0-v1L0))
                                                       + 0.25*tau*mu*exp(-g00) * ((v2R0-v2L0)/h1+(v10R-v10L)/h2) * ((v2R0-v2L0)/h1+(v10R-v10L)/h2)
                                                       ;//+ tau*cas_ft(t_hat,X[i],Y[i],R,mu,cv,kappa);
        triplets_t.emplace_back(m00,m00,a_T00);
        triplets_t.emplace_back(m00,mR0,a_TR0);
        triplets_t.emplace_back(m00,mL0,a_TL0);
        triplets_t.emplace_back(m00,m0R,a_T0R);
        triplets_t.emplace_back(m00,m0L,a_T0L);
    }
    else if (st[i]==1) // left border flow start
    {
        m00=i;
        a_T00=1;
        Bt[m00] = cas_theta(t_hat,X[i],Y[i],h1,tt);
        triplets_t.emplace_back(m00,m00,a_T00); // matrixt is Dim*Dim, T is Dim, i=0,..Dim-1; first coord is matrixt's string, second coord is matrixt's row
    }
    else if (st[i]==2) // right border
    {
        m00=i; mL0=i-1;
        a_T00=1;
        a_TR0=-1;
        Bt[m00] = 0;//cas_theta(t_hat,X[i],Y[i],h1,tt);
        triplets_t.emplace_back(m00,m00,a_T00);
        triplets_t.emplace_back(m00,mL0,a_TR0);
    }
    else if (st[i]==3) // down border
    {
        m00=i; m0R=M0R[i];
        a_T00=-1;
        a_TR0=1;
        Bt[m00] = 0;//cas_theta(t_hat,X[i],Y[i],h1,tt);
        triplets_t.emplace_back(m00,m00,a_T00);
        triplets_t.emplace_back(m00,m0R,a_TR0);
    }
    else if (st[i]==4) // up border
    {
        m00=i; m0L=M0L[i];
        a_T00=1;
        a_TR0=-1;
        Bt[m00] = 0;//cas_theta(t_hat,X[i],Y[i],h1,tt);
        triplets_t.emplace_back(m00,m00,a_T00);
        triplets_t.emplace_back(m00,m0L,a_TR0);
    }
    else if (st[i]==5) // left down corner (test)
    {
        m00=i; m0R=M0R[i];
        a_T00=-1;
        a_TR0=1;
        Bt[m00] = 0;//cas_theta(t_hat,X[i],Y[i],h1,tt);
        triplets_t.emplace_back(m00,m00,a_T00);
        triplets_t.emplace_back(m00,m0R,a_TR0);
    }
    else if (st[i]==6) // right down corner (test)
    {
        m00=i; m0R=M0R[i];
        a_T00=-1;
        a_TR0=1;
        Bt[m00] = 0;//cas_theta(t_hat,X[i],Y[i],h1,tt);
        triplets_t.emplace_back(m00,m00,a_T00);
        triplets_t.emplace_back(m00,m0R,a_TR0);
    }
    else if (st[i]==7) // left up corner (test) flow start
    {
        m00=i; m0L=M0L[i];
        a_T00=1;
        a_TR0=-1;
        Bt[m00] = 0;//cas_theta(t_hat,X[i],Y[i],h1,tt);
        triplets_t.emplace_back(m00,m00,a_T00);
        triplets_t.emplace_back(m00,m0L,a_TR0);
    }
    else if (st[i]==8) // right up corner (test)
    {
        m00=i; m0L=M0L[i];
        a_T00=1;
        a_TR0=-1;
        Bt[m00] = 0;//cas_theta(t_hat,X[i],Y[i],h1,tt);
        triplets_t.emplace_back(m00,m00,a_T00);
        triplets_t.emplace_back(m00,m0L,a_TR0);
    }
  }
  //build matrixt, solve with eigen
  matrixt.setFromTriplets (triplets_t.begin(), triplets_t.end());
  eigen_solver_t solvert(matrixt);

  if (solvert.info() != Eigen::Success)
  {
      fprintf(OUT, "Can not build preconditioner");
      return -1;
  }
  eigen_vector_t res (Dim);
  res = solvert.solveWithGuess(Bt, T);
  if (solvert.info() != Eigen::Success)
  {
      fprintf(OUT, "Failed to solve the system with Eigen, n=%d", n); fclose(OUT);
      n_c[4*ij]=0; n_c[4*ij+1]=0; n_c[4*ij+2]=0; n_c[4*ij+3]=0;
      n_l[4*ij]=0; n_l[4*ij+1]=0; n_l[4*ij+2]=0; n_l[4*ij+3]=0; t_len[ij]=0;
      return 0;
  }
  T=res;
  if (n==500 || n==1000 || n==2000 || n==3000 || n==4000 || n==5000 || n==6000 || n==7000 || n==8000 || n==9000 || n==10000 || n==15000 || n==20000 || n==25000 || n==30000 || n==40000 || n==50000) {ress_st=ress; res_st=res; dif=-1;}
  if (n==520 || n==1020 || n==2020 || n==3020 || n==4020 || n==5020 || n==6020 || n==7020 || n==8020 || n==9020 || n==10020 || n==15020 || n==20020 || n==25020 || n==30020 || n==40020 || n==50020)
  {
      for (int i=0; i<Dim; i++)
      {
          if (dif<fabs(VEC[3*i]-ress_st[3*i])) dif=fabs(VEC[3*i]-ress_st[3*i]);
          if (dif<fabs(VEC[3*i+1]-ress_st[3*i+1])) dif=fabs(VEC[3*i+1]-ress_st[3*i+1]);
          if (dif<fabs(VEC[3*i+2]-ress_st[3*i+2])) dif=fabs(VEC[3*i+2]-ress_st[3*i+2]);
          if (dif<fabs(T[i]-res_st[i])) dif=fabs(T[i]-res_st[i]);
      }
      char file1[20];    sprintf(file1, "%d", n);    strcat(file1,txt); OUT1 = fopen(file1, "w");   fprintf(OUT1, "%lf\n",dif);
                                                                        fclose(OUT1);
  }
  //maxg=-1; maxu1=-1; maxu2=-1; maxt=-1;
  t_hat+=tau; mu_C=0; kap_C=0; n++; t_len[ij]=n*tau;
  /*if ((n-1)%5==0)
  {
      char file1[20];    sprintf(file1, "%d", n);    strcat(file1,txt); OUT1 = fopen(file1, "w");   for (int i=0; i<Dim; i++) fprintf(OUT1, "%lf\n",T[i]);
                                                                        fclose(OUT1);
  }*/

}
// solve end =========================================================================================================
fprintf(OUT, "1");

//normi
maxg=0; maxu1=0; maxu2=0; maxt=0;
sumu1=0; sumu2=0; sumg=0; sumt=0;
for (int i=0; i<Dim; i++)
{
    if (maxg<fabs(VEC[3*i])) maxg=fabs(VEC[3*i]);
    if (maxu1<fabs(VEC[3*i+1])) maxu1=fabs(VEC[3*i+1]);
    if (maxu2<fabs(VEC[3*i+2])) maxu2=fabs(VEC[3*i+2]);
    if (maxt<fabs(T[i])) {maxt=fabs(T[i]);}// fprintf(OUT, " %d/%d ", i,st[i]);}
    sumg+=(VEC[3*i])*(VEC[3*i]);
    sumu1+=(VEC[3*i+1])*(VEC[3*i+1]);
    sumu2+=(VEC[3*i+2])*(VEC[3*i+2]);
    sumt+=(T[i])*(T[i]);
}
n_c[4*ij]=maxu1; n_c[4*ij+1]=maxu2; n_c[4*ij+2]=maxg; n_c[4*ij+3]=maxt;
n_l[4*ij]=sqrt(h1*h2*sumu1); n_l[4*ij+1]=sqrt(h1*h2*sumu2); n_l[4*ij+2]=sqrt(h1*h2*sumg); n_l[4*ij+3]=sqrt(h1*h2*sumt);
fclose(OUT);
return 0;
}

int Sxema_hcond_vgas_2_mymethod(double mu, int N,int M1, int M2, int Dim,
                      Str_param_gas *str_param_gas, Str_setka *str_setka, P_she *p_she,
                      int *st, double *X, double *Y, int *M0R, int *M0L,
                      double *n_c, double *n_l, double *n_w, int ij, int K, int cas, double *delta_mas)
{

return 0;
}
