void output(int step, real *r, real *u, real *v, real *w);
void output_2d(int step, real *r, real *u, real *v, real *w);
void dump(int step, const real *f);
void read_dump(int step, real *f);
void statistics_serial(const real rhg[], const real uxg[], const real uyg[], const real uzg[], const int nx, const int ny, const int nz, const int step);
