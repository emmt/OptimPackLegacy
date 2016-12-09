extern __op_lbfgsb_wrapper;
/* PROTOTYPE
   long lbfgsb_wrapper(long n, long m, double array x, double array f,
                       double array g, double array l, double array u,
                       long array bnd, double factr, double pgtol,
                       char array task, char array csave, long array isave,
                       double array dsave, long iprint);
 */
func op_lbfgsb_setup(m, x, l, u, bnd, factr, pgtol)
{
  n = numberof(x);
  if (structof(x) != double) error, "X must be a double array";

  if (numberof(l) != n) error, "bad number of elements for L";
  if (structof(l) != double) error, "L must be a double array";

  if (numberof(u) != n) error, "bad number of elements for U";
  if (structof(u) != double) error, "U must be a double array";

  if (numberof(bnd) != n) error, "bad number of elements for BND";
  if (structof(bnd) != long) error, "BND must be a long array";

  ctask = csave = array(char, 61);
  ctask(1) = 'S';
  ctask(2) = 'T';
  ctask(3) = 'A';
  ctask(4) = 'R';
  ctask(5) = 'T';
  dsave = array(double, (2*m + 4)*n + (11*m + 8)*m + 29);
  isave = array(long, 3*n + 27);
  ws = [&m, &n, &l, &u, &bnd, &factr, &pgtol, &ctask, &csave, &isave, &dsave];
  job = __op_lbfgsb_wrapper(n, m, x, 0.0, x, *ws(3), *ws(4), *ws(5), *ws(6),
                            *ws(7), *ws(8), *ws(9), *ws(10), *ws(11), -1);
  if (job != 1) error, op_lbfgsb_task(ws);
  return ws;
}

func op_lbfgsb_msg(ws, csave) { return string(ws((csave ? 9 : 8))); }

func op_lbfgsb(&x, &f, &g, ws)
{
  m = *ws(1);
  n = *ws(2);
  if (numberof(x) != n || structof(x) != double)
    error, swrite(format="X must be a double array with %d elements", n);
  if (structof(f) != double || dimsof(f)(1))
    error, "F must be a double scalar";
  if (numberof(g) != n || structof(g) != double)
    error, swrite(format="G must be a double array with %d elements", n);
  return __op_lbfgsb_wrapper(n, m, x, f, g, *ws(3), *ws(4), *ws(5), *ws(6),
                             *ws(7), *ws(8), *ws(9), *ws(10), *ws(11), -1);
}


#if 0
#include "fits2.i"
#include "fft_utils.i"
img = fits_read("~/data/stsdas-testdata/data/sn1987a/sn1987a.fits");
psf = roll(fits_read("~/data/stsdas-testdata/data/sn1987a/sn1987a_psf.fits"));
psf *= 1.0/sum(psf);
x = deconv_lbfgsb(img,psf,maxiter=10,ndirs=5);
x2 =  deconv_vmlmb(img,psf,maxiter=10,verb=1)
#endif



