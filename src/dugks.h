//extern real *rho, *ux, *uy, *uz;
extern real *h_rho, *h_ux, *h_uy, *h_uz;
extern real *h_f;

void init_dugks(void );
void dugks_step(int step);
void copy_vel(void );
void copy_ddf(void );
void upload_macro();
