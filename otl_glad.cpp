#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "input.h"
#include "tex_init.h"
#include "tex_end.h"
#include "hcond_vgas_2.h"
#include "setka.h"
#include "tex_tab.h"
int main(int argc, char **argv)
{
//----------------------------------------------------------
//Assigning filename the value specified on the command line
//----------------------------------------------------------
FILE *LOG;
    int exit_code = 0;
    char filename_dat[80];
    char filename_tex[80];
    char dat[] = ".dat";
    char tex[] = ".tex";
    char txt[] = ".txt";
    char log[] = "log";
    strcat(log,txt);
    LOG = fopen(log, "w");
    if((argc < 2) || (argc > 2)) {
        exit_code = argc;
    printf("wrong input arguements, %d", argc);
    fprintf(LOG, "wrong input arguements, %d", argc);
    return exit_code;
    }
    strcpy(filename_dat,argv[1]);
    strcpy(filename_tex,filename_dat);
    strcat(filename_dat,dat);
    strcat(filename_tex,tex);
//----------------------------------------------------------
//   присвоение значений элементам структур str_param_gas &
// str_setka, где хранятся параметры дифференциальной и
//                   разностной задач
//----------------------------------------------------------
    Str_param_gas str_param_gas;
    Str_setka str_setka;
    P_she p_she;
    exit_code = input(&str_param_gas,&str_setka,filename_dat);
    if(exit_code == -1) {
        exit_code = -1;
        fprintf(LOG, "%d", exit_code);
    return exit_code;
    }
    else if(exit_code == -2){
        exit_code = -2;
        fprintf(LOG, "%d", exit_code);
    return exit_code;
    }
/*cases:
0: debugging test
1: debugging test 2
2: task 1.1 nested grids
3: task 1.1 nested grids 2
4: task 1.2 the flow problem
5: task 1.2 the flow problem 2
6: ----
*/
int cas=str_setka.cas;
int len=0;
int n_max=str_setka.n_max, m_max=str_setka.m_max;
int var_const=str_setka.var_const;
switch (cas){
case 0:
    len=n_max*m_max;
    break;
case 1:
    len=n_max*m_max;
    break;
case 2:
    len=(n_max+1)*m_max;
    break;
case 3:
    len=(n_max+1)*m_max;
    break;
case 4:
    len=var_const*var_const;
    break;
case 5:
    len=var_const*var_const;
    break;
case 6:
    len=1;
    break;
}

double *norm_c, *norm_l, *norm_w, *t_len, *mas_len;
double delta_mas=0;
norm_c=(double*)malloc(4*len*sizeof(double));
    if (norm_c == 0) {
        printf("Malloc error for norm_c \n");
        return 2;
    }
norm_l=(double*)malloc(4*len*sizeof(double));
    if (norm_c == 0) {
        printf("Malloc error for norm_l \n");
        return 2;
    }
t_len=(double*)malloc(len*sizeof(double));
    if (norm_c == 0) {
        printf("Malloc error for t_len \n");
        return 2;
    }
mas_len=(double*)malloc(len*sizeof(double));
    if (norm_c == 0) {
        printf("Malloc error for mas_len \n");
        return 2;
    }
for (int i=0; i<len; i++)
  {
      mas_len[i]=0; t_len[i]=0;
  }
for (int i=0; i<4*len; i++)
  {
      norm_l[i]=0; norm_c[i]=0;
  }

int N=str_setka.N0, M1=str_setka.M10, M2=str_setka.M20, k_dr=str_setka.k_dr;
int ij=0, K=0;
double mu=str_param_gas.mu;
int Dim=0, *st; //number of grids
double *X, *Y; int *M0R, *M0L; //arrays with size=Dim: st - knot status; X,Y - knot coordinates; M0R, M0L - knot neighbours

tex_init(filename_tex); // printing header of filename.tex
printf("M10=%d M20=%d N0=%d k_dr=%d n_max=%d m_max=%d\n",M1,M2, N,k_dr,n_max,m_max);
//for case 4
double v_start=5, g_start=0, t_start=303; //start flow characteristics
double tv=v_start, tg=g_start, tt=t_start; //tilde_v, tilde_g
switch (cas){
case 0:
  for (int n=0; n<n_max; n++)
  {
    for (int m=0; m<m_max; m++)
    {
        Dim = (M1+1)*(M2+1) - M1*M2/6;
        //we have rectangle 2x3, it has M1*M2 knots + M1 knots on the top border + M2 knots on the right border - M1*M2/6 knots in the square Omega_{01}
        //Dim=(M1+1)*(M2+1);
        printf("%d \n", Dim);
        st=(int*)malloc(Dim*sizeof(int));
        if (st == 0) {
            printf("Malloc error for st \n");
            return 2;
        }
        X=(double*)malloc(Dim*sizeof(double));
        if (X == 0) {
            printf("Malloc error for X \n");
            return 2;
        }
        Y=(double*)malloc(Dim*sizeof(double));
        if (Y == 0) {
            printf("Malloc error for Y \n");
            return 2;
        }
        M0R=(int*)malloc(Dim*sizeof(int));
        if (M0R == 0) {
            printf("Malloc error for M0R \n");
            return 2;
        }
        M0L=(int*)malloc(Dim*sizeof(int));
        if (M0L == 0) {
            printf("Malloc error for M0L \n");
            return 2;
        }

        p_she.h1=str_param_gas.X1/M1;
        p_she.h2=str_param_gas.X2/M2;
        p_she.tau=str_param_gas.T/N;

        fprintf(LOG, "start setka \n");
        if (Setka(Dim, M1, M2, &str_param_gas, &str_setka, &p_she, st, X, Y, M0R, M0L)!=0){ //fill arrays st, X, Y, M0R, M0L
                printf("Ошибка в процессе выполнения Setka \n");
                exit_code = 3;
            fprintf(LOG, "%d", exit_code);
            return exit_code;
            break;
            }
        printf("DEBUGGING TEST\n");
        fprintf(LOG, "start shema \n");
        clock_t BegClock, EndClock;
        BegClock = clock();
        if(Sxema_hcond_vgas_2(mu,N,M1,M2,Dim,&str_param_gas,&str_setka, &p_she, st, X, Y, M0R, M0L, norm_c,norm_l,norm_w, ij, K, cas, &delta_mas) != 0) {
                printf("Ошибка в процессе выполнения Sxema \n");
                exit_code = 2;
            fprintf(LOG, "%d", exit_code);
            return exit_code;
            break;
            }
        EndClock = clock();
        t_len[ij] = (double)(EndClock - BegClock)/CLOCKS_PER_SEC;
        free(st);
        free(X);
        free(Y);
        free(M0R);
        free(M0L);
        M1*=k_dr; M2*=k_dr; ij++;
    }
    N*=k_dr; M1=str_setka.M10; M2=str_setka.M20;
  }
  tex_tab(mu, &str_param_gas, &str_setka, norm_c, norm_l, t_len, filename_tex, 0);
  tex_tab(mu, &str_param_gas, &str_setka, norm_c, norm_l, t_len, filename_tex, 1);
  tex_tab(mu, &str_param_gas, &str_setka, norm_c, norm_l, t_len, filename_tex, 2);
  tex_tab(mu, &str_param_gas, &str_setka, norm_c, norm_l, t_len, filename_tex, 3);
  break;
case 1:
  for (int n=0; n<n_max; n++)
  {
    for (int m=0; m<m_max; m++)
    {
        Dim = (M1+1)*(M2+1) - M1*M2/6;
        //we have rectangle 2x3, it has M1*M2 knots + M1 knots on the top border + M2 knots on the right border - M1*M2/6 knots in the square Omega_{01}
        //Dim=(M1+1)*(M2+1);
        printf("%d \n", Dim);
        st=(int*)malloc(Dim*sizeof(int));
        if (st == 0) {
            printf("Malloc error for st \n");
            return 2;
        }
        X=(double*)malloc(Dim*sizeof(double));
        if (X == 0) {
            printf("Malloc error for X \n");
            return 2;
        }
        Y=(double*)malloc(Dim*sizeof(double));
        if (Y == 0) {
            printf("Malloc error for Y \n");
            return 2;
        }
        M0R=(int*)malloc(Dim*sizeof(int));
        if (M0R == 0) {
            printf("Malloc error for M0R \n");
            return 2;
        }
        M0L=(int*)malloc(Dim*sizeof(int));
        if (M0L == 0) {
            printf("Malloc error for M0L \n");
            return 2;
        }

        p_she.h1=str_param_gas.X1/M1;
        p_she.h2=str_param_gas.X2/M2;
        p_she.tau=str_param_gas.T/N;

        fprintf(LOG, "start setka \n");
        if (Setka(Dim, M1, M2, &str_param_gas, &str_setka, &p_she, st, X, Y, M0R, M0L)!=0){ //fill arrays st, X, Y, M0R, M0L
                printf("Ошибка в процессе выполнения Setka \n");
                exit_code = 3;
            fprintf(LOG, "%d", exit_code);
            return exit_code;
            break;
            }
        printf("DEBUGGING TEST 2\n");
        fprintf(LOG, "start shema \n");
        clock_t BegClock, EndClock;
        BegClock = clock();
        if(Sxema_hcond_vgas_2_mymethod(mu,N,M1,M2,Dim,&str_param_gas,&str_setka, &p_she, st, X, Y, M0R, M0L, norm_c,norm_l,norm_w, ij, K, cas, &delta_mas) != 0) {
                printf("Ошибка в процессе выполнения Sxema \n");
                exit_code = 2;
            fprintf(LOG, "%d", exit_code);
            return exit_code;
            break;
            }
        EndClock = clock();
        t_len[ij] = (double)(EndClock - BegClock)/CLOCKS_PER_SEC;
        free(st);
        free(X);
        free(Y);
        free(M0R);
        free(M0L);
        M1*=k_dr; M2*=k_dr; ij++;
    }
    N*=k_dr; M1=str_setka.M10; M2=str_setka.M20;
  }
  tex_tab(mu, &str_param_gas, &str_setka, norm_c, norm_l, t_len, filename_tex, 0);
  tex_tab(mu, &str_param_gas, &str_setka, norm_c, norm_l, t_len, filename_tex, 1);
  tex_tab(mu, &str_param_gas, &str_setka, norm_c, norm_l, t_len, filename_tex, 2);
  tex_tab(mu, &str_param_gas, &str_setka, norm_c, norm_l, t_len, filename_tex, 3);
    break;
case 2:
    double *VGth, *Tth, *Xth, *Yth;
  for (int m=0; m<m_max; m++) // original grid
  {
    Dim = (M1+1)*(M2+1) - M1*M2/6;
    VGth=(double*)malloc(3*Dim*sizeof(double)); //VG-results on tau-h grid
    if (VGth == 0) {
        printf("Malloc error for VGth \n");
        return 2;
    }
    Tth=(double*)malloc(Dim*sizeof(double)); //T-results on tau-h grid
    if (Tth == 0) {
        printf("Malloc error for Tth \n");
        return 2;
    }
    Xth=(double*)malloc(Dim*sizeof(double)); //X-coord on tau-h grid
    if (Xth == 0) {
        printf("Malloc error for Tth \n");
        return 2;
    }
    Yth=(double*)malloc(Dim*sizeof(double)); //Y-coord on tau-h grid
    if (Yth == 0) {
        printf("Malloc error for Tth \n");
        return 2;
    }
    for (int n=0; n<n_max+1; n++) // k=original, 1,2,3
    {
        Dim = (M1+1)*(M2+1) - M1*M2/6;
        //we have rectangle 2x3, it has M1*M2 knots + M1 knots on the top border + M2 knots on the right border - M1*M2/6 knots in the square Omega_{01}
        //Dim=(M1+1)*(M2+1);
        printf("%d \n", Dim);
        st=(int*)malloc(Dim*sizeof(int));
        if (st == 0) {
            printf("Malloc error for st \n");
            return 2;
        }
        X=(double*)malloc(Dim*sizeof(double));
        if (X == 0) {
            printf("Malloc error for X \n");
            return 2;
        }
        Y=(double*)malloc(Dim*sizeof(double));
        if (Y == 0) {
            printf("Malloc error for Y \n");
            return 2;
        }
        M0R=(int*)malloc(Dim*sizeof(int));
        if (M0R == 0) {
            printf("Malloc error for M0R \n");
            return 2;
        }
        M0L=(int*)malloc(Dim*sizeof(int));
        if (M0L == 0) {
            printf("Malloc error for M0L \n");
            return 2;
        }

        p_she.h1=str_param_gas.X1/M1;
        p_she.h2=str_param_gas.X2/M2;
        p_she.tau=str_param_gas.T/N;

        fprintf(LOG, "start setka \n");
        if (Setka(Dim, M1, M2, &str_param_gas, &str_setka, &p_she, st, X, Y, M0R, M0L)!=0){ //fill arrays st, X, Y, M0R, M0L
                printf("Ошибка в процессе выполнения Setka \n");
                exit_code = 3;
            fprintf(LOG, "%d", exit_code);
            return exit_code;
            break;
            }
        printf("TASK 1.1 NESTED GRIDS\n");
        fprintf(LOG, "start shema \n");

        if(Sxema_hcond_vgas_2_nested_grid(mu,N,M1,M2,Dim,&str_param_gas,&str_setka, &p_she, st, X, Y, M0R, M0L, norm_c,norm_l,norm_w, ij, K, cas, n, VGth, Tth, Xth, Yth) != 0) {
                printf("Ошибка в процессе выполнения Sxema \n");
                exit_code = 2;
            fprintf(LOG, "%d", exit_code);
            return exit_code;
            break;
            }
        free(st);
        free(X);
        free(Y);
        free(M0R);
        free(M0L);

        N*=2; M1*=2; M2*=2; ij++;
    }
    M1=str_setka.M10; M2=str_setka.M20; N=str_setka.N0;
    for (int j=0; j<m+1; j++) {M1*=k_dr; M2*=k_dr; N*=k_dr;}
    free(VGth); free(Tth); free(Xth); free(Yth);
  }
  fprintf(LOG, "start print \n");
  tex_tab_flowp(mu, &str_param_gas, &str_setka, norm_c, norm_l, t_len, filename_tex, 0);
  tex_tab_flowp(mu, &str_param_gas, &str_setka, norm_c, norm_l, t_len, filename_tex, 1);
  tex_tab_flowp(mu, &str_param_gas, &str_setka, norm_c, norm_l, t_len, filename_tex, 2);
  tex_tab_flowp(mu, &str_param_gas, &str_setka, norm_c, norm_l, t_len, filename_tex, 3);
  break;
case 3:
    printf("TASK 1.1 NESTED GRIDS 2\n");
    break;
case 4:
    //the problem will be solved with the only set of tau,h;
    Dim = (M1+1)*(M2+1) - M1*M2/6;
    //we have rectangle 2x3, it has M1*M2 knots + M1 knots on the top border + M2 knots on the right border - M1*M2/6 knots in the square Omega_{01}
    //Dim=(M1+1)*(M2+1);
    printf("%d \n", Dim);
    st=(int*)malloc(Dim*sizeof(int));
    if (st == 0) {
        printf("Malloc error for st \n");
        return 2;
    }
    X=(double*)malloc(Dim*sizeof(double));
    if (X == 0) {
        printf("Malloc error for X \n");
        return 2;
    }
    Y=(double*)malloc(Dim*sizeof(double));
    if (Y == 0) {
        printf("Malloc error for Y \n");
        return 2;
    }
    M0R=(int*)malloc(Dim*sizeof(int));
    if (M0R == 0) {
        printf("Malloc error for M0R \n");
        return 2;
    }
    M0L=(int*)malloc(Dim*sizeof(int));
    if (M0L == 0) {
        printf("Malloc error for M0L \n");
        return 2;
    }

    p_she.h1=str_param_gas.X1/M1;
    p_she.h2=str_param_gas.X2/M2;
    p_she.tau=str_param_gas.T/N;

    fprintf(LOG, "start setka \n");
        if (Setka_flowp(Dim, M1, M2, &str_param_gas, &str_setka, &p_she, st, X, Y, M0R, M0L)!=0){ //fill arrays st, X, Y, M0R, M0L
                printf("Ошибка в процессе выполнения Setka \n");
                exit_code = 3;
            fprintf(LOG, "%d", exit_code);
            return exit_code;
            break;
            }
    printf("TASK 1.2 THE FLOW PROBLEM\n");

    for (int v=0; v<var_const; v++)
    {
      for (int g=0; g<var_const; g++)
      {
        fprintf(LOG, "start shema \n");

        if(Sxema_hcond_vgas_2_flowp(mu,N,M1,M2,Dim,&str_param_gas,&str_setka, &p_she, st, X, Y, M0R, M0L, norm_c,norm_l,norm_w, ij, K, cas, t_len, tv, tg, tt) != 0) {
                printf("Ошибка в процессе выполнения Sxema \n");
                exit_code = 2;
            fprintf(LOG, "%d", exit_code);
            return exit_code;
            break;
            }

        ij++; tg+=0.5;
      }
      tv++; tg=g_start;
    }
    free(st); free(X); free(Y); free(M0R); free(M0L);
    fprintf(LOG, "start print \n");
    tex_tab_flowp(mu, &str_param_gas, &str_setka, norm_c, norm_l, t_len, filename_tex, 0);
    break;
case 5:
    printf("TASK 1.2, THE FLOW PROBLEM 2\n");
    break;
case 6:
    printf("TASK ??\n");
    break;
}

tex_end(filename_tex);
printf("press any key to continue\n");
fprintf(LOG, "end\n");
int end=0;
scanf("%d", &end);


free(norm_c);
free(norm_l);
free(t_len);
free(mas_len);
fclose(LOG);
return 0;
}

