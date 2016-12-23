* Update documentation (use Doxygen?).

* Implement gradient-based convergence test.

* Use `fmin`.

* `gamma` should be computed (and applied) when there is a preconditioner.

* Backtrack (using Armijo's rule) when at least one bound constraint becomes
  active along the line search.

* Implement damped BFGS updating, see Nocédal & Wright p. 537.

* Scale initial gradient.  Perhaps use a typical value for the step norm

* Safeguard the step.

* Implement BLMVM.

* Saving best variables so far can be very cheap: it is sufficient to
  remember the best step length, and corresponding function value and
  gradient norm.  Not applicable with constraints.
