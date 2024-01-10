#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "input.h"
int Setka(int Dim, int M1, int M2, Str_param_gas *str_param_gas, Str_setka *str_setka, P_she *p_she, int *st, double *X, double *Y, int *M0R, int *M0L)
{
double X1 = str_param_gas->X1;
double X2 = str_param_gas->X2;
double EX1 = str_param_gas->EX1;
double EX2 = str_param_gas->EX2;

double h1=X1/M1;//=p_she->h1;
double h2=X2/M2;//=p_she->h2;
EX1+=1; //use coordinates of block's right border
int i=0;
char txt[] = ".txt";
char file[20];
sprintf(file, "dim");
strcat(file,txt);
FILE *OUT;
OUT = fopen(file, "w");
for(int m2=0; m2<M2+1; m2++)
{
    for(int m1=0;m1<M1+1;m1++)
    {
        if (m2*h2<EX2+h2/3)
        {
            if (m1*h1>0 && m1*h1<X1 && m2*h2>0 && m2*h2<EX2) st[i]=0; //inner knots
            else if (m1*h1>EX1 && m1*h1<X1 && m2*h2>0 && m2*h2<EX2+h2/3) st[i]=0; //inner knots
            else if (m1*h1>EX1-h1/2 && m1*h1<EX1+h1/2 && m2*h2>EX2-h2/2) st[i]=11; //inner knot (EX1,EX2)
            else if (m1*h1<h1/2 && m2*h2>0 && m2*h2<EX2) st[i]=1; //left border
            else if (m1*h1>X1-h1/2 && m2*h2>0 && m2*h2<=EX2) st[i]=2; //right border
            else if (m1*h1>0 && m1*h1<X1 && m2*h2<h2/2) st[i]=3; //down border
            else if (m1*h1>0 && m1*h1<EX1 && m2*h2>EX2-h2/2) st[i]=4; //up border (y=EX2)
            else if (m1*h1<h1/2 && m2*h2<h2/2) st[i]=5; //left down corner
            else if (m1*h1>EX1-h1/2 && m2*h2<h2/2) st[i]=6; //right down corner
            else if (m1*h1<h1/2 && m2*h2>EX2-h2/2) st[i]=7; //left up corner (y=EX2)
            X[i]=m1*h1; Y[i]=m2*h2;
            if (m2*h2>0) M0L[i]=(m2-1)*(M1+1)+m1;
            else M0L[i]=-1;
            fprintf(OUT, "%d ", st[i]);
            if (m2*h2>EX2-h2/2 && m1*h1<EX1) M0R[i]=-1;
            else if (m2*h2<EX2-h2/2) M0R[i]=(m2+1)*(M1+1)+m1;
            else if (m2*h2>EX2-h2/2) M0R[i]=(M2/2)*(M1+1)+2*M1/3+1 + m1;//-M1/3;
            i++;
        }
        else if (m2*h2>EX2 && m1*h1>=EX1)
        {
            if (m1*h1>EX1 && m1*h1<X1 && m2*h2>EX2 && m2*h2<X2) st[i]=0; //inner knots
            else if (m1*h1<EX1+h1/2 && m2*h2>EX2 && m2*h2<X2) st[i]=1; //left border (x=EX1)
            else if (m1*h1>X1-h1/2 && m2*h2>EX2 && m2*h2<X2) st[i]=2; //right border
            else if (m1*h1>EX1 && m1*h1<X1 && m2*h2<h2/2) st[i]=3; //down border
            else if (m1*h1>EX1 && m1*h1<X1 && m2*h2>X2-h2/2) st[i]=4; //up border
            else if (m1*h1<EX1+h1/2 && m2*h2>X2-h2/2) st[i]=7; //left up corner
            else if (m1*h1>X1-h1/2 && m2*h2>X2-h2/2) st[i]=8; //right up corner
            X[i]=m1*h1; Y[i]=m2*h2;

            if (m2*h2<EX2+3*h2/2) M0L[i]=(M2/2)*(M1+1)+(m2-M2/2-1)*(2*M1/3+1)+m1;
            else M0L[i]=(M2/2)*(M1+1)+(m2-M2/2-1)*(2*M1/3+1)+m1;//-M1/3;
            fprintf(OUT, "%d ", st[i]);
            if /*(m2*h2<EX2+3*h2/2) M0R[i]=(M2/2-1)*(M1+1) +m1;
            else if*/ (m2*h2<X2) M0R[i]=(M2/2)*(M1+1)+(m2-M2/2+1)*(2*M1/3+1)+m1;//-M1/3;
            else M0R[i]=-1;
            i++;
            if (i>Dim) printf("%d, %m1=%d, m2=%d ", i, m1, m2);
        }
    }
    fprintf(OUT, "\n ");
}
printf("%d \n", i);
return 0;
}
int Setka_test(int Dim, int M1, int M2, Str_param_gas *str_param_gas, Str_setka *str_setka, P_she *p_she, int *st, double *X, double *Y, int *M0R, int *M0L)
{
double X1 = str_param_gas->X1;
double X2 = str_param_gas->X2;
double EX1 = str_param_gas->EX1;
double EX2 = str_param_gas->EX2;

double h1=X1/M1;//=p_she->h1;
double h2=X2/M2;//=p_she->h2;
EX1+=1; //use coordinates of block's right border
int i=0;
char txt[] = ".txt";
char file[20];
sprintf(file, "dim");
strcat(file,txt);
FILE *OUT;
OUT = fopen(file, "w");

for(int m2=0; m2<M2+1; m2++)
{
    for(int m1=0;m1<M1+1;m1++)
    {
        if (m1*h1>0 && m1*h1<X1 && m2*h2>0 && m2*h2<X2) st[i]=0; //inner knots
        else if (m1*h1<h1/3 && m2*h2>0 && m2*h2<X2) st[i]=1; //left border
        else if (m1*h1>X1-h1/3 && m2*h2>0 && m2*h2<X2) st[i]=2; //right border
        else if (m1*h1>0 && m1*h1<X1 && m2*h2==0) st[i]=3; //down border
        else if (m1*h1>0 && m1*h1<X1 && m2*h2>X2-h2/3) st[i]=4; //up border
        else if (m1*h1<h1/3 && m2*h2<h2/2) st[i]=5; //left down corner
        else if (m1*h1>X1-h1/3 && m2*h2<h2/3) st[i]=6; //right down corner
        else if (m1*h1<h1/3 && m2*h2>X2-h2/3) st[i]=7; //left up corner
        else if (m1*h1>X1-h1/3 && m2*h2>X2-h2/3) st[i]=8; //right up corner

        X[i]=m1*h1; Y[i]=m2*h2;
        if (m2*h2>0) M0L[i]=i-M1-1;//(m2-1)*(M1+1)+m1;
        else M0L[i]=-1;
        if (m2*h2<X2) M0R[i]=i+M1+1;//(m2+1)*(M1+1)+m1;
        else M0R[i]=-1;

        fprintf(OUT, "%d ", st[i]);
        //fprintf(OUT, "%lf ", Y[i]);
        i++;
    }
    fprintf(OUT, "\n ");
}
printf("%d \n", i);
fclose(OUT);
return 0;
}
//===============================================================================
int Setka_test2(int Dim, int M1, int M2, Str_param_gas *str_param_gas, Str_setka *str_setka, P_she *p_she, int *st, double *X, double *Y, int *M0R, int *M0L)
{
double X1 = str_param_gas->X1;
double X2 = str_param_gas->X2;
double EX1 = str_param_gas->EX1;
double EX2 = str_param_gas->EX2;

double h1=X1/M1;//=p_she->h1;
double h2=X2/M2;//=p_she->h2;
EX1+=1; //use coordinates of block's right border
int i=0;
char txt[] = ".txt";
char file[20];
sprintf(file, "dim");
strcat(file,txt);
FILE *OUT;
OUT = fopen(file, "w");

for(int m1=0;m1<M1+1;m1++)
{
    for(int m2=0; m2<M2+1; m2++)
    {
        if (m1*h1>0 && m1*h1<X1 && m2*h2>0 && m2*h2<X2) st[i]=0; //inner knots
        else if (m1*h1<h1/3 && m2*h2>0 && m2*h2<X2) st[i]=1; //left border
        else if (m1*h1>X1-h1/3 && m2*h2>0 && m2*h2<X2) st[i]=2; //right border
        else if (m1*h1>0 && m1*h1<X1 && m2*h2==0) st[i]=3; //down border
        else if (m1*h1>0 && m1*h1<X1 && m2*h2>X2-h2/3) st[i]=4; //up border
        else if (m1*h1<h1/3 && m2*h2<h2/2) st[i]=5; //left down corner
        else if (m1*h1>X1-h1/3 && m2*h2<h2/3) st[i]=6; //right down corner
        else if (m1*h1<h1/3 && m2*h2>X2-h2/3) st[i]=7; //left up corner
        else if (m1*h1>X1-h1/3 && m2*h2>X2-h2/3) st[i]=8; //right up corner

        X[i]=m1*h1; Y[i]=m2*h2;
        if (m2*h2>0) M0L[i]=(m1-1)*(M2+1)+m2;
        else M0L[i]=-1;
        if (m2*h2<X2) M0R[i]=(m1+1)*(M2+1)+m2;
        else M0R[i]=-1;

        fprintf(OUT, "%d ", st[i]);
        //fprintf(OUT, "%lf ", Y[i]);
        i++;
    }
    fprintf(OUT, "\n ");
}
printf("%d \n", i);
fclose(OUT);
return 0;
}

int Setka_flowp(int Dim, int M1, int M2, Str_param_gas *str_param_gas, Str_setka *str_setka, P_she *p_she, int *st, double *X, double *Y, int *M0R, int *M0L)
{
 double X1 = str_param_gas->X1;
double X2 = str_param_gas->X2;
double EX1 = str_param_gas->EX1;
double EX2 = str_param_gas->EX2;

double h1=X1/M1;//=p_she->h1;
double h2=X2/M2;//=p_she->h2;
EX1+=1; //use coordinates of block's right border
int i=0;
char txt[] = ".txt";
char file[20];
sprintf(file, "dim");
strcat(file,txt);
FILE *OUT;
OUT = fopen(file, "w");
for(int m2=0; m2<M2+1; m2++)
{
    for(int m1=0;m1<M1+1;m1++)
    {
        if (m2*h2<EX2+h2/3)
        {
            if (m1*h1>0 && m1*h1<X1 && m2*h2>0 && m2*h2<EX2) st[i]=0; //inner knots
            else if (m1*h1>EX1 && m1*h1<X1 && m2*h2>0 && m2*h2<EX2+h2/3) st[i]=0; //inner knots
            else if (m1*h1>EX1-h1/2 && m1*h1<EX1+h1/2 && m2*h2>EX2-h2/2) st[i]=11; //inner knot (EX1,EX2)
            else if (m1*h1<h1/2 && m2*h2>0 && m2*h2<EX2) st[i]=1; //left border flow start
            else if (m1*h1>X1-h1/2 && m2*h2>0 && m2*h2<=EX2) st[i]=3; //right border
            else if (m1*h1>0 && m1*h1<X1 && m2*h2<h2/2) st[i]=4; //down border
            else if (m1*h1>0 && m1*h1<EX1 && m2*h2>EX2-h2/2) st[i]=5; //up border (y=EX2)
            else if (m1*h1<h1/2 && m2*h2<h2/2) st[i]=6; //left down corner
            else if (m1*h1>EX1-h1/2 && m2*h2<h2/2) st[i]=7; //right down corner
            else if (m1*h1<h1/2 && m2*h2>EX2-h2/2) st[i]=8; //left up corner (y=EX2) flow start
            X[i]=m1*h1; Y[i]=m2*h2;
            if (m2*h2>0) M0L[i]=(m2-1)*(M1+1)+m1;
            else M0L[i]=-1;
            fprintf(OUT, "%d ", st[i]);
            if (m2*h2>EX2-h2/2 && m1*h1<EX1) M0R[i]=-1;
            else if (m2*h2<EX2-h2/2) M0R[i]=(m2+1)*(M1+1)+m1;
            else if (m2*h2>EX2-h2/2) M0R[i]=(M2/2)*(M1+1)+2*M1/3+1 + m1;//-M1/3;
            i++;
        }
        else if (m2*h2>EX2 && m1*h1>=EX1)
        {
            if (m1*h1>EX1 && m1*h1<X1 && m2*h2>EX2 && m2*h2<X2) st[i]=0; //inner knots
            else if (m1*h1<EX1+h1/2 && m2*h2>EX2 && m2*h2<X2) st[i]=2; //left border (x=EX1) no flow
            else if (m1*h1>X1-h1/2 && m2*h2>EX2 && m2*h2<X2) st[i]=3; //right border
            else if (m1*h1>EX1 && m1*h1<X1 && m2*h2<h2/2) st[i]=4; //down border
            else if (m1*h1>EX1 && m1*h1<X1 && m2*h2>X2-h2/2) st[i]=5; //up border
            else if (m1*h1<EX1+h1/2 && m2*h2>X2-h2/2) st[i]=9; //left up corner no flow
            else if (m1*h1>X1-h1/2 && m2*h2>X2-h2/2) st[i]=10; //right up corner
            X[i]=m1*h1; Y[i]=m2*h2;

            if (m2*h2<EX2+3*h2/2) M0L[i]=(M2/2)*(M1+1)+/*(m2-M2/2-1)*(2*M1/3+1)*/+m1;//m2-M2/2-1
            else M0L[i]=(M2/2)*(M1+1)+(m2-M2/2-1)*(2*M1/3+1)+m1;
            fprintf(OUT, "%d ", st[i]);
            if /*(m2*h2<EX2+3*h2/2) M0R[i]=(M2/2-1)*(M1+1) +m1;
            else if*/ (m2*h2<X2) M0R[i]=(M2/2)*(M1+1)+(m2-M2/2+1)*(2*M1/3+1)+m1;
            else M0R[i]=-1;
            i++;
            if (i>Dim) printf("%d, %m1=%d, m2=%d ", i, m1, m2);
        }
    }
    fprintf(OUT, "\n ");
}
printf("%d \n", i);
return 0;
}
