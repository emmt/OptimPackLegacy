

extern void lbfgs_(long *N, long *M, double *X, double *F, double *G,
		   int *DIAGCO, double *DIAG, long *IPRINT,
		   double *EPS, double *XTOL, double *W, long *IFLAG);

void lbfgs(long n, long m, double x[], double *f, double g[],
	   int diagco, double diag[],
	   double eps, double xtol, double w[], long *iflag)
{
  long iprint[2];
  iprint[0] = -1;
  iprint[1] = 0;
  lbfgs_(&n, &m, x, f, g, &diagco, diag, iprint, &eps, &xtol, w, iflag);
}
