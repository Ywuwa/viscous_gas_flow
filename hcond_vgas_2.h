int Sxema_hcond_vgas_2(double mu, int N,int M1, int M2, int Dim,
                      Str_param_gas *str_param_gas, Str_setka *str_setka, P_she *p_she,
                      int *st, double *X, double *Y, int *M0R, int *M0L,
                      double *n_c, double *n_l, double *n_w, int ij, int K, int cas, double *delta_mas);

int Sxema_hcond_vgas_2_test(double mu, int N,int M1, int M2, int Dim,
                      Str_param_gas *str_param_gas, Str_setka *str_setka, P_she *p_she,
                      int *st, double *X, double *Y, int *M0R, int *M0L,
                      double *n_c, double *n_l, double *n_w, int ij, int K, int cas, double *delta_mas);

int Sxema_hcond_vgas_2_nested_grid(double mu, int N,int M1, int M2, int Dim,
                      Str_param_gas *str_param_gas, Str_setka *str_setka, P_she *p_she,
                      int *st, double *X, double *Y, int *M0R, int *M0L,
                      double *n_c, double *n_l, double *n_w, int ij, int K, int cas, int k, double *VGth, double *Tth, double *Xth, double *Yth);

int Sxema_hcond_vgas_2_flowp(double mu, int N,int M1, int M2, int Dim,
                      Str_param_gas *str_param_gas, Str_setka *str_setka, P_she *p_she,
                      int *st, double *X, double *Y, int *M0R, int *M0L,
                      double *n_c, double *n_l, double *n_w, int ij, int K, int cas, double *t_len, double tv, double tg, double tt);

int Sxema_hcond_vgas_2_flowp_test(double mu, int N,int M1, int M2, int Dim,
                      Str_param_gas *str_param_gas, Str_setka *str_setka, P_she *p_she,
                      int *st, double *X, double *Y, int *M0R, int *M0L,
                      double *n_c, double *n_l, double *n_w, int ij, int K, int cas, double *t_len, double tv, double tg, double tt);

int Sxema_hcond_vgas_2_mymethod(double mu, int N,int M1, int M2, int Dim,
                      Str_param_gas *str_param_gas, Str_setka *str_setka, P_she *p_she,
                      int *st, double *X, double *Y, int *M0R, int *M0L,
                      double *n_c, double *n_l, double *n_w, int ij, int K, int cas, double *delta_mas);
