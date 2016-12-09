;;
;;local rosenbrock_nevals;
;;rosenbrock_nevals=0;
;;
;;func optim_test_rosenbrock(nil, start=, method=, ndirs=, frtol=)
;;/* DOCUMENT optim_test_rosenbrock, ...;
;;     Test optim_driver with Rosenbrock function.
;;     
;;   SEE ALSO: optim_test_rosenbrock_func. */
;;{
;;  x = is_void(start) ? [0.0, 0.0] : start;
;;  return optim_driver(optim_test_rosenbrock_func,
;;                      (is_void(start) ? [0.0, 0.0] : start),
;;                      method=method, fmin=0.0, verb=1, ndirs=ndirs,
;;                      frtol=frtol);
;;}
;;
;;func optim_test_rosenbrock_func(x, &g)
;;{
;;  ++rosenbrock_nevals;
;;  return f;
;;}
;;
;;func optim_test_quad(x, &g)
;;{
;;  u = (x - [1.0, 1.0, 1.0]);
;;  w = [3.0, 7.0, 2.0];
;;  g = 2.0*w*u;
;;  return sum(w*u*u);
;;}
;;
pro op_test, start=start, m=m, frtol=frtol
  n = 2L
  if n_elements(start) eq n then begin
    x = start
  end else begin
    x = [0D0, 0D0]
  end
  ws = op_vmlmb_setup(n, m, frtol=frtol, fatol=0d0)
  job = 1L
  iter = 0L
  eval = 0L

  while 1b do begin
    if job eq 1L then begin
      ;; Compute penalty and gradient.

      ;; Rosenbrock: f  = 100*(x2 - x1^2)^2 + (1 - x1)^2
      x1 = x(0)
      x2 = x(1)
      u = x2 - x1*x1
      v = 1D0 - x1
      f = 100D0*u*u + v*v
      g = [-400D0*u*x1 - 2D0*v, 200D0*u] ;
      
      eval = eval + 1L
    end

    if job ne 1L or eval eq 1L then begin
      if job eq 2L or job eq 3L then iter = iter + 1L
      print, iter, eval, f, x(0), x(1)
      if job ge 3L then return
    end

    ;; Call optimizer.
    job = op_vmlmb(x, f, g, ws)
  end
end
