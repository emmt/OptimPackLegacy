# List of things to do

## For release 1.5.0

* Safeguard the step when there are bounds.

* Fix Python interface.

* Backtrack (using Armijo's rule or a qudratic interpolation) when at least one
  bound constraint becomes active along the line search.  This is related to
  safeguarding the step.


## For release 1.6.0

* Use floating-point weights instead of `isfree` logical array (faster and more
  portable).

* Update documentation (use Doxygen?).

* Implement gradient-based convergence test (`gatol`/`grtol`).

* Implement damped BFGS updating, see Nocédal & Wright p. 537.

* Implement BLMVM.

* Saving best variables so far can be very cheap: it is sufficient to
  remember the best step length, and corresponding function value and
  gradient norm.  Not applicable with constraints.
