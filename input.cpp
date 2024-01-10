#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "input.h"
// input.h

int input(Str_param_gas *str_param_gas,Str_setka *str_setka,char *filename)
{
    FILE* input;
//    char filename[]="perenos.dat";
    char name[100];
    float X1, X2, EX1, EX2, T, mu, kappa, cv, R, eps; //EX = Empty X
    int mu_dr, mu_max, N0, M10, M20, k_dr, m_max, n_max, var_const, cas;
    if ((input = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "Can't open %s\n", filename);
        return -2;
    }
            if (fgets(name,100,input) == 0) {
                printf("Wrong input: the first line \n");
                fclose(input);
                return -1;
            }
            if (fgets(name,100,input) == 0) {
                printf("Wrong input: the second line \n");
                fclose(input);
                return -1;
            }
            if (fscanf(input, "%f", &X1) != 1) { //read X1
                printf("Wrong input X1  \n");
                fclose(input);
                return -1;
            }
            if (fgets(name,100,input ) == 0) {
                printf("Wrong input: text in line X1 \n");
                fclose(input);
                return -1;
            }
            if (fscanf(input, "%f", &X2) != 1) { //read X2
                printf("Wrong input X2  \n");
                fclose(input);
                return -1;
            }
            if (fgets(name,100,input ) == 0) {
                printf("Wrong input: text in line X2 \n");
                fclose(input);
                return -1;
            }
            if (fscanf(input, "%f", &EX1) != 1) { //read EX1
                printf("Wrong input EX1  \n");
                fclose(input);
                return -1;
            }
            if (fgets(name,100,input ) == 0) {
                printf("Wrong input: text in line EX1 \n");
                fclose(input);
                return -1;
            }
            if (fscanf(input, "%f", &EX2) != 1) { //read EX2
                printf("Wrong input EX2  \n");
                fclose(input);
                return -1;
            }
            if (fgets(name,100,input ) == 0) {
                printf("Wrong input: text in line EX2 \n");
                fclose(input);
                return -1;
            }
            if (fscanf(input, "%f", &T) != 1) { //read T
                printf("Wrong input T \n");
                fclose(input);
                return -1;
            }
            if (fgets(name,100,input ) == 0) {
                printf("Wrong input: text in line T \n");
                fclose(input);
                return -1;
            }
            if (fscanf(input, "%f", &mu) != 1) { //read mu
                printf("Wrong input mu \n");
                fclose(input);
                return -1;
            }
            if (fgets(name,100,input ) == 0) {
                printf("Wrong input: text in line mu \n");
                fclose(input);
                return -1;
            }
            if (fscanf(input, "%f", &kappa) != 1) { //read kappa
                printf("Wrong input kappa \n");
                fclose(input);
                return -1;
            }
            if (fgets(name,100,input ) == 0) {
                printf("Wrong input: text in line kappa \n");
                fclose(input);
                return -1;
            }
            if (fscanf(input, "%f", &cv) != 1) { //read cv
                printf("Wrong input cv \n");
                fclose(input);
                return -1;
            }
            if (fgets(name,100,input ) == 0) {
                printf("Wrong input: text in line cv \n");
                fclose(input);
                return -1;
            }
            if (fscanf(input, "%f", &R) != 1) { //read R
                printf("Wrong input R \n");
                fclose(input);
                return -1;
            }
            if (fgets(name,100,input ) == 0) {
                printf("Wrong input: text in line R \n");
                fclose(input);
                return -1;
            }
            if (fscanf(input, "%d", &mu_dr) != 1) { //read mu_dr
                printf("Wrong input mu_dr \n");
                fclose(input);
                return -1;
            }
            if (fgets(name,100,input ) == 0) {
                printf("Wrong input: text in line mu_dr \n");
                fclose(input);
                return -1;
            }
            if (fscanf(input, "%d", &mu_max) != 1) { //read mu_max
                printf("Wrong input mu_max \n");
                fclose(input);
                return -1;
            }
            if (fgets(name,100,input ) == 0) {
                printf("Wrong input: text in line mu_max \n");
                fclose(input);
                return -1;
            }
            if (fscanf(input, "%f", &eps) != 1) { //read eps
                printf("Wrong input eps \n");
                fclose(input);
                return -1;
            }
            if (fgets(name,100,input ) == 0) {
                printf("Wrong input: text in line eps \n");
                fclose(input);
                return -1;
            }
    printf("X1=%f X2=%f T=%f\n",X1,X2,T);
    printf("mu=%f kappa=%f cv=%f\n",mu,kappa,cv);
    printf("mu_dr=%d mu_max=%d eps=%f\n",mu_dr, mu_max, eps);
            if (fgets(name,100,input ) == 0) {
                printf("Wrong input: the first line after input gas parametrs \n");
                fclose(input);
                return -1;
            }
            if (fgets(name,100,input ) == 0) {
                printf("Wrong input: the second line after input gas parametrs \n");
                fclose(input);
                return -1;
            }
            if (fgets(name,100,input ) == 0) {
                printf("Wrong input: the third line after input gas parametrs \n");
                fclose(input);
                return -1;
            }
            if (fscanf(input, "%d", &N0) != 1) { //read N0
                printf("Wrong input: N0 \n");
                fclose(input);
                return -1;
            }
            if (fgets(name,100,input ) == 0) {
                printf("Wrong input: text in line N0 \n");
                fclose(input);
                return -1;
            }
            if (fscanf(input, "%d", &M10) != 1) { //read M10
                printf("Wrong input M10 \n");
                fclose(input);
                return -1;
            }
            if (fgets(name,100,input ) == 0) {
                printf("Wrong input: text in line M10 \n");
                fclose(input);
                return -1;
            }
            if (fscanf(input, "%d", &M20) != 1) { //read M20
                printf("Wrong input M20 \n");
                fclose(input);
                return -1;
            }
            if (fgets(name,100,input ) == 0) {
                printf("Wrong input: text in line M20 \n");
                fclose(input);
                return -1;
            }
            if (fscanf(input, "%d", &k_dr) != 1) { //read k_dr
                printf("Wrong input k_dr \n");
                fclose(input);
                return -1;
            }
            if (fgets(name,100,input ) == 0) {
                printf("Wrong input: text in line k_dr \n");
                fclose(input);
                return -1;
            }
            if (fscanf(input, "%d", &n_max) != 1) { //read n_max
                printf("Wrong input n_max \n");
                fclose(input);
                return -1;
            }
            if (fgets(name,100,input ) == 0) {
                printf("Wrong input: text in line n_max \n");
                fclose(input);
                return -1;
            }
            if (fscanf(input, "%d", &m_max) != 1) { //read m_max
                printf("Wrong input m_max \n");
                fclose(input);
                return -1;
            }
            if (fgets(name,100,input ) == 0) {
                printf("Wrong input: text in line m_max \n");
                fclose(input);
                return -1;
            }

            if (fscanf(input, "%d", &var_const) != 1) { //read var_const
                printf("Wrong input var_const \n");
                fclose(input);
                return -1;
            }
            if (fgets(name,100,input ) == 0) {
                printf("Wrong input: text in line var_const \n");
                fclose(input);
                return -1;
            }
            if (fscanf(input, "%d", &cas) != 1) { //read cas
                printf("Wrong input cas \n");
                fclose(input);
                return -1;
            }
            if (fgets(name,100,input ) == 0) {
                printf("Wrong input: text in line cas \n");
                fclose(input);
                return -1;
            }
    //printf("M10=%d M20=%d N0=%d k_dr=%d n_max=%d m_max=%d\n",M10,M20, N0,k_dr,n_max,m_max);
    str_param_gas->X1=(double) X1;
    str_param_gas->X2=(double) X2;
    str_param_gas->EX1=(double) EX1;
    str_param_gas->EX2=(double) EX2;

    str_param_gas->T=(double) T;
    str_param_gas->mu=(double) mu;
    str_param_gas->kappa=(double) kappa;
    str_param_gas->cv=(double) cv;
    str_param_gas->R=(double) R;
    str_param_gas->eps=(double) eps;
    str_param_gas->mu_dr= mu_dr;
    str_param_gas->mu_max= mu_max;
      str_setka->N0=N0;
      str_setka->M10=X1*M10;
      str_setka->M20=X2*M20;
      str_setka->k_dr=k_dr;
      str_setka->n_max=n_max;
      str_setka->m_max=m_max;
      str_setka->var_const=var_const;
      str_setka->cas=cas;
    fclose(input);
    return 0;
}
