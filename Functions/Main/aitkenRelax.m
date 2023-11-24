function [sol_mp1_relaxed, res_m_nextIter, delta_m_nextIter, sol_m_NextIter] = ...
    aitkenRelax(sol_mp1, sol_m, idxSol, res_m, delta_m, negON, range)
%AITKENRELAX Apply Aitken delta^2 relaxation method with limited range of relaxation parameter
%   sol_mp1_damped = 
%       AITKENRELAX(sol_mp1, sol_m, idxSol, res_m, delta_m, negON, range)
%   returns the relaxed solution at the current iteration having used the
%   previous two iterations to improve the guess at the solution.
%
%   [sol_mp1_damped, res_m_nextCall] =
%       AITKENRELAX(sol_mp1, sol_m, idxSol, res_m, delta_m, negON, range)
%   also returns Delta d to be used as res_m in the next iter.
% 
%   [sol_mp1_damped, res_m_nextCall, delta_m_nextCall] =
%       AITKENRELAX(sol_mp1, sol_m, idxSol, res_m, delta_m, negON, range)
%   also returns the relaxation paramter to be used in next iter.
%
%   [sol_mp1_damped, res_m_nextCall, delta_m_nextCall, sol_m_NextCall] =
%       AITKENRELAX(sol_mp1, sol_m, idxSol, res_m, delta_m, negON, range)
%   also returns the solution to be used as sol_m in next iter.
%
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   sol_mp1   :     Solution at the current iteration
%   sol_m     :     Relaxed solution from the previous iteration
%   idxSol    :     Index of DOFs used to calculate the relaxation parameter
%   res_m     :     Delta d from the previous iteration 
%   delta_m   :     Relaxation parameter from previous iteration
%   negON     :     Boolean allowing or restricting the relaxation
%                   parameter to be negative (w<0).
%   range     :     2x1 vector with the minimum and maximum allowed
%                   relaxation parameters.
% 
%   --------------------------------------------------------------------
%   Output
%   --------------------------------------------------------------------
%   sol_mp1_relaxed     :   Relaxed solution at the current iteration
%   res_m_nextIter      :   Delta d to be used as res_m in the next iter
%   delta_m_nextIter    :   Relaxation paramter to be used in next iter
%   sol_m_nextIter      :   Solution to be used as sol_m in next iter


% Handle input.
if nargin < 7
    minRelax = 0.005;       maxRelax = 1;
else
    minRelax = range(1); maxRelax = range(2);
end

% Residual of current iterate, defined as difference in iterates.
% res_mp1 = sol_mp1(idxSol) - sol_m(idxSol);
res_mp1 = sol_mp1 - sol_m;

% Compute Aitken relaxation parameters.
aitkenPara = delta_m * (- res_m(idxSol)' * (res_mp1(idxSol) - res_m(idxSol)) / norm(res_mp1(idxSol) - res_m(idxSol),2)^2); 


% Constrain relaxation parameter.
if negON
    % If negatives aitken parameters are allowed, contrain the parameter
    % such that it is limited to the ranges [-maxRelax to -minRelax] and
    % [minRelax to maxRelax]
    delsign = sign(aitkenPara);
    if delsign == 0
        delsign = 1;
    end
    delta = max(min(maxRelax, abs(aitkenPara)),minRelax)*delsign;
else
    delta = max(min(maxRelax,aitkenPara),minRelax);
end

sol_mp1_relaxed = (1-delta) * sol_m + delta * sol_mp1;


% Prepare values for next call. Aitken relaxation is recursive! 
res_m_nextIter   = res_mp1;
delta_m_nextIter = delta;
sol_m_NextIter   = sol_mp1_relaxed; 

end
