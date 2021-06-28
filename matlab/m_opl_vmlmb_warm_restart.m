% function m_opl_vmlmb_warm_restart
%
%     task = m_opl_vmlmb_warm_restart(ws)
%
% This function initiates a new VMLM-B iteration assuming the same curvature
% even though the variables and objective function (and its gradient) may have
% changed.  Compared to `m_opl_vmlmb_restart`, the memorized L-BFGS model of
% the Hessian is kept.  Calling this function is not considered as a restart
% (hence the number of objective function calls, of iterations and of restarts
% are left unchanged).
%
% This function normally returns `OPL_TASK_FG` meaning that the caller shall
% compute the objective function and its gradient before calling
% `m_opl_vmlmb_iterate`.
