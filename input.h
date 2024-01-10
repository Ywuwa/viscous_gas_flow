typedef struct
{
// zadanie parametrov gaza i oblasti:
// T - razmeri oblasti [0;T]
// X1, X2 - razmeri oblasti [0,X1]x[0,X2]
// mu - nachalnii parametr viazkosti
// mu_dr - koefficient droblenia viazkosti
// mu_max - chislo variantov viazkosti

double T;
double X1;
double X2;
double EX1;
double EX2;
double mu;
double kappa;
double cv;
double R;
int mu_dr;
int mu_max;
double eps;
} Str_param_gas;

typedef struct
{
// zadanie parametrov cetok

// k_dr - koefficient droblenia cetok

int N0;
int M10;
int M20;
int k_dr;
int n_max;
int m_max;
int var_const;
int cas;
} Str_setka;

typedef struct
{
double h1;
double h2;
double tau;
} P_she;

int input(Str_param_gas *str_param_gas, Str_setka *str_setka,char *filename);

