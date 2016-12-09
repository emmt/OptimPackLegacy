extern __op_lbfgs;
/* PROTOTYPE
   void lbfgs(long n, long m, double array x, double array f, double array g,
              int diagco, double array diag, double eps, double xtol,
              double array w, long array iflag);
*/
/* __PROTOTYPE FORTRAN
   void LBFGS(long N, long M, double array X, double array F, double array G,
              int DIAGCO, double array DIAG, long array IPRINT,
              double EPS, double XTOL, double array W, long array IFLAG)
*/
func op_lbfgs_setup(n, m, eps=, xtol=)
{
  if (is_void(eps)) eps = 1e-8;
  if (is_void(xtol)) xtol = 1e-16;
  isave = array(long, 3);
  isave(1) = n;
  isave(2) = m;
  dsave = array(double, n*(2*m + 1) + 2*m);
  ws = [&isave, &dsave];
  return ws;
}

func op_lbfgs_task(ws) { return (*ws(1))(3); }

func op_lbfgs_next(x, &f, &g, ws, diag)
{
  local isave; eq_nocopy, isave, *ws(1);
  if (structof(isave) != long || numberof(isave) != 3)
    error, "corrupted workspace (ISAVE)";
  m = isave(1);
  n = isave(2);
  iflag = isave(3);
  
  local dsave; eq_nocopy, dsave, *ws(2);
  if (structof(dsave) != double || numberof(dsave) != n + 2*m*(n + 1))
    error, "corrupted workspace (DSAVE)";

  if (structof(x) != double || numberof(x) != n)
    error, "bad parameter array X";
  if (structof(f) != double || dimsof(f)(1))
    error, "bad function value F";
  if (structof(g) != double || numberof(g) != n)
    error, "bad gradient array G";

  //iprint = [-1, 0];
  diag = array(double, n);
  diagco = 0n;
  eps = 1e-8;
  xtol = 1e-16;
  
  __op_lbfgs, n, m, x, f, g, diagco, diag, /*iprint,*/ eps, xtol, dsave, iflag;
  isave(3) = iflag;
  if (iflag < 0) error, "line search error";
  return (iflag == 0 ? 3 : iflag);
  
//     IFLAG   is an INTEGER variable that must be set to 0 on initial entry
//             to the subroutine. A return with IFLAG<0 indicates an error,
//             and IFLAG=0 indicates that the routine has terminated without
//             detecting errors. On a return with IFLAG=1, the user must
//             evaluate the function F and gradient G. On a return with
//             IFLAG=2, the user must provide the diagonal matrix Hk0.
}
