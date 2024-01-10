#include "input.h"
#include "tex_tab.h"
#include <stdio.h>

void tex_tab(double mu, Str_param_gas *str_param_gas, Str_setka *str_setka, double *norm_c,double *norm_l, double *t_nm, char *filename, int i)
{
int n,m,M10, M20, N0,n_max,m_max,k_dr;
M10=str_setka->M10;
M20=str_setka->M20;
N0=str_setka->N0;
n_max=str_setka->n_max;
m_max=str_setka->m_max;
k_dr=str_setka->k_dr;
double h1,h2,tau;

FILE *fi1 = fopen(filename,"a");
fprintf(fi1,"\\begin{center}\n");
  switch (i){
  case 0: fprintf(fi1,"table norms of the ERROR in C,$\\,$ L2 for u1 and times\n");break;
  case 1: fprintf(fi1,"table norms of the ERROR in C,$\\,$ L2 for u2 and times\n");break;
  case 2: fprintf(fi1,"table norms of the ERROR in C,$\\,$ L2 for g and times\n");break;
  case 3: fprintf(fi1,"table norms of the ERROR in C,$\\,$ L2 for $\\theta$ and times\n");break;
  }
fprintf(fi1,"\\\\[2.0ex]  \n");

fprintf(fi1,"\\begin{tabular}{|p{0.6in}|");
  for (m=0;m<m_max;m++) fprintf(fi1,"p{0.8in}|");
  fprintf(fi1,"} \\hline\n");
  fprintf(fi1,"$\\tau\\setminus h_1, h_2$ ");
  h1=str_param_gas->X1/M10;
  h2=str_param_gas->X2/M20;
  for (m=0;m<m_max;m++)
    {
        fprintf(fi1,"& $h_1=%.5lf \,$ $h_2=%.5lf$",h1, h2);
        h1/=k_dr; h2/=k_dr;
    }
  fprintf(fi1,"\\\\ \\hline\n");

  tau=str_param_gas->T/N0;
  for(n = 0; n < n_max; n++)
    {
      fprintf(fi1,"\n");
      fprintf(fi1,"$%.5lf$ & ",tau);
      for(m = 0; m < m_max; m++)
        {
          fprintf(fi1, "$%10.3e$ $%10.3e$ $%10.3e$ ",norm_c[n*m_max*4+m*4+i],norm_l[n*m_max*4+m*4+i],t_nm[n*m_max+m]);//%10.3e
          if(m < m_max-1)
             fprintf(fi1,"&");
        }
      fprintf(fi1," \\\\ \\hline");
      tau/=k_dr;
     }
  fprintf(fi1,"\n\\end{tabular}\\\\[20pt]\n");

  fprintf(fi1,"\\end{center}\n");
  fclose(fi1);
}
void tex_tab_nested_grid(double mu, Str_param_gas *str_param_gas, Str_setka *str_setka, double *norm_c,double *norm_l, char *filename, int i)
{
    int n,m,M10, M20, N0,n_max,m_max,k_dr;
M10=str_setka->M10;
M20=str_setka->M20;
N0=str_setka->N0;
n_max=str_setka->n_max;
m_max=str_setka->m_max;
k_dr=str_setka->k_dr;
double h1,h2,tau;
tau=str_param_gas->T/N0;
h1=str_param_gas->X1/M10;
h2=str_param_gas->X2/M20;

FILE *fi1 = fopen(filename,"a");
fprintf(fi1,"\\begin{center}\n");
  switch (i){
  case 0: fprintf(fi1,"table norms of the ERROR in C,$\\,$ L2 for u1 and times\n");break;
  case 1: fprintf(fi1,"table norms of the ERROR in C,$\\,$ L2 for u2 and times\n");break;
  case 2: fprintf(fi1,"table norms of the ERROR in C,$\\,$ L2 for g and times\n");break;
  case 3: fprintf(fi1,"table norms of the ERROR in C,$\\,$ L2 for $\\theta$ and times\n");break;
  }
fprintf(fi1,"\\\\[2.0ex]  \n");

fprintf(fi1,"\\begin{tabular}{|p{0.6in}|");
  for (n=0;n<n_max+1;n++) fprintf(fi1,"p{0.8in}|");
  fprintf(fi1,"} \\hline\n");
  fprintf(fi1,"$\k\\setminus \\tau, h $ ");
  for (n=0;n<n_max+1;n++)
    {
        fprintf(fi1,"& %d ",n);
    }
  fprintf(fi1,"\\\\ \\hline\n");

  for(m = 0; m < m_max; m++)
    {
      fprintf(fi1,"\n");
      fprintf(fi1,"$\\tau=%.5lf \,$ $h=%.5lf$ &",tau, h2);
        tau/=k_dr; h2/=k_dr;
      for(n = 0; n < n_max+1; n++)
        {
          fprintf(fi1, "$%10.3e$ $%10.3e$ ",norm_c[m*(n_max+1)*4+n*4+i],norm_l[m*(n_max+1)*4+n*4+i]);//%10.3e
          if(n < n_max)
             fprintf(fi1,"&");
        }
      fprintf(fi1," \\\\ \\hline");
     }
  fprintf(fi1,"\n\\end{tabular}\\\\[20pt]\n");

  fprintf(fi1,"\\end{center}\n");
  fclose(fi1);
}

void tex_tab_flowp(double mu, Str_param_gas *str_param_gas, Str_setka *str_setka, double *norm_c,double *norm_l, double *t_nm, char *filename, int i)
{
    int n,m,M10, M20, N0,n_max,m_max,k_dr, var;
M10=str_setka->M10;
M20=str_setka->M20;
N0=str_setka->N0;
var=str_setka->var_const;
k_dr=str_setka->k_dr;
double h1,h2,tau;
tau=str_param_gas->T/N0;
h1=str_param_gas->X1/M10;
h2=str_param_gas->X2/M20;

FILE *fi1 = fopen(filename,"a");
fprintf(fi1,"\\begin{center}\n");
  switch (i){
  case 0: fprintf(fi1,"table norms of the ERROR in C,$\\,$ L2 for u1 and times\n");break;
  case 1: fprintf(fi1,"table norms of the ERROR in C,$\\,$ L2 for u2 and times\n");break;
  case 2: fprintf(fi1,"table norms of the ERROR in C,$\\,$ L2 for g and times\n");break;
  case 3: fprintf(fi1,"table norms of the ERROR in C,$\\,$ L2 for $\\theta$ and times\n");break;
  }
fprintf(fi1,"\\\\[2.0ex]  \n");

fprintf(fi1,"\\begin{tabular}{|p{0.6in}|");
  for (n=0;n<var;n++) fprintf(fi1,"p{0.8in}|");
  fprintf(fi1,"} \\hline\n");
  fprintf(fi1,"$\\tilde{v}\\setminus \\tilde{g} $ ");
  for (n=0;n<var;n++)
    {
        fprintf(fi1,"& %d ",n+1);
    }
  fprintf(fi1,"\\\\ \\hline\n");

  for(m = 0; m < var; m++)
    {
      fprintf(fi1,"\n");
      fprintf(fi1,"$%d$& ", m+1);
      for(n = 0; n < var; n++)
        {
          fprintf(fi1, "$%10.3e$ ",t_nm[m*var+n]);//%10.3e
          if(n < var-1)
             fprintf(fi1,"&");
        }
      fprintf(fi1," \\\\ \\hline");
     }
  fprintf(fi1,"\n\\end{tabular}\\\\[20pt]\n");

  fprintf(fi1,"\\end{center}\n");
  fclose(fi1);
}
