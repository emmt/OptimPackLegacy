% function m_opl_vmlmb_create
%
%  ws = m_opl_vmlmb_create(nparam, m, fatol, frtol,sftol, sgtol, sxtol, epsilon, delta);
%  Create the workspace ws for VMLMB.
%  
%   The arguments are:
%   N is the number of parameters.
%
%   M is the number of correction pairs to remember in order to compute the
%       limited memory variable metric (BFGS) approximation of the inverse of
%       the Hessian.  For large problems, M = 3 to 5 gives good results.  For
%       small problems, M should be less or equal N.  The larger is M (and N)
%       the more computer memory will be needed to store the workspace WS.
%
%   FRTOL is the relative error desired in the function (e.g.  FRTOL=1e-8).
%       Convergence occurs if the estimate of the relative error between F(X)
%       and F(XSOL), where XSOL is a local minimizer, is less or equal FRTOL.
%       FRTOL must have a non-negative floating point value.
%
%   FATOL is the absolute error desired in the function (e.g. FATOL=0.0).
%       Convergence occurs if the estimate of the absolute error between F(X)
%       and F(XSOL), where XSOL is a local minimizer, is less or equal FATOL.
%       FATOL must have a non-negative floating point value.
%
%   SFTOL, SGTOL, and SXTOL are tolerances for the line search subroutine (see
%       opl_csrch).   Recommended  values: SFTOL=0.001,  SGTOL=0.9,  SXTOL=0.1
%       (other values  may be more  suitable for highly  non-quadratic penalty
%       function).
%
%   DELTA is a small nonegative value used to compute a small initial step.
%
%   EPSILON is a small value, in the range [0,1), equals to the cosine of the
%       maximum angle between the search direction and the anti-gradient.  The
%       BFGS recursion is restarted, whenever the search direction is not
%       sufficiently "descending".
%
%   WS is a  workspace array. The caller must not modify the workspace array WS between calls
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

