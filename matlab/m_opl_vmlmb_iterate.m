% function m_opl_vmlmb_iterate
%  task = m_opl_vmlmb_iterate(ws,xopt,f,grad,active,h);
%
%   X is a double precision array of length N.  On entry, X is an
%       approximation to the solution.  On exit with TASK = OPL_TASK_CONV, X
%       is the current approximation.
%
%   F is the address of a double precision variable.  On entry, F is the value
%       of the function at X.  On final exit, F is the function value at X.
%
%   G is a double precision array of length N.  On entry, G is the value of
%       the gradient at X.  On final exit, G is the value of the gradient at
%       X.
%
%   ISFREE is an optional integer array with length N provided by the caller
%       if the values in X have bounds.  If the parameters have no bounds,
%       ISFREE should be NULL (unconstrained minimization).  Otherwise,
%       elements set to zero in ISFREE indicate that the corresponding values
%       in X has reached a bound and should not be changed during the next
%       step because the gradient has the wrong sign (i.e.  the steepest
%       descent direction would violate the bound constraints):
%
%           ISFREE[i] = 0 if i-th value has a lower bound XLO[i]
%                         and X[i] = XLO[i] and G[i] >= 0
%                       0 if i-th value has an upper bound XHI[i]
%                         and X[i] = XHI[i] and G[i] <= 0
%                       1 otherwise
%
%       ISFREE needs only to be computed when TASK = OPL_TASK_FREEVARS.  If X
%       has (some) bounds, the caller is responsible for applying the bounds
%       to X before evaluating the function value F and the gradient G (i.e.,
%       when TASK = OPL_TASK_FG), e.g.:
%
%           if (X[i] < XLO[i]) X[i] = XLO[i];
%           if (X[i] > XHI[i]) X[i] = XHI[i];
%
%       If H is not specified (i.e., H is NULL) or if H[i] > 0 for all i such
%       that ISFREE[i] is non-zero, then ISFREE is left unchanged.
%
%   H is an optional double precision array with length N provided by the
%       caller and such that diag(H) is an approximation of the inverse of the
%       Hessian matrix.  If H is NULL, then the inverse of the Hessian is
%       approximated by a simple rescaling using Shanno & Phua formula.
%       Otherwise, if ISFREE is NULL, all elements of H must be strictly
%       greater than zero; else ISFREE[i] is set to zero if H[i] <= 0 (this is
%       the only case where ISFREE is modified).  
%
%   TASK is the value returned opl_vmlmb_iterate.  It
%       can have one of the following values:
%
%           OPL_TASK_FG - caller must evaluate the function and gradient at X
%               and call opl_vmlm_next.
%
%           OPL_TASK_FREEVARS - if variables are bounded, caller must
%               determine the set of free variables for the current variables
%               X and update IFREE accordingly.
%
%           OPL_TASK_NEWX - a new iterate has been computed.  The
%               approximation X, function F, and gradient G are available for
%               examination.
%
%           OPL_TASK_CONV - the search is successful.  The solution, function
%               value and gradient are available in X, F and G.
%
%           OPL_TASK_WARN - VMLMB is not able to satisfy the convergence
%               conditions.  The exit value of X contains the best
%               approximation found so far.  Warning message is given by
%               opl_vmlmb_get_reason(ws).
%
%           OPL_TASK_ERROR then there is an error in the input arguments.
%               Error message is  given by opl_vmlmb_get_reason(ws).
%

%
% See optimpacklegacy for explaination about VMLMB algorithm and its
% parameters
%
%	Definitions for optimization routines implemented in OptimPack
%	library.
%
%-----------------------------------------------------------------------------
%
%      Copyright (c) 2018, Ferreol SOULEZ.
%
%	This file is part of OptimPack.
%
%	OptimPack is  free software; you can redistribute  it and/or modify
%	it under the  terms of the GNU General  Public License as published
%	by the Free  Software Foundation; either version 2  of the License,
%	or (at your option) any later version.
%
%	OptimPack is  distributed in the hope  that it will  be useful, but
%	WITHOUT  ANY  WARRANTY;  without   even  the  implied  warranty  of
%	MERCHANTABILITY or  FITNESS FOR A PARTICULAR PURPOSE.   See the GNU
%	General Public License for more details.
%
%	You should have  received a copy of the  GNU General Public License
%	along with OptimPack (file  "LICENSE" in the top source directory);
%	if  not, write  to the  Free Software  Foundation, Inc.,  59 Temple
%	Place, Suite 330, Boston, MA 02111-1307 USA

