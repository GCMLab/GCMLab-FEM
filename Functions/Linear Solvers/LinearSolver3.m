function [d, fE] = LinearSolver3(K, F, dfixed, free, fixed, beta)
%SOLVEMETHOD3 Solve for unknown vector using Penalty method
%   d = LinearSolver3(K,F)
%
%   ----------------------------------------------------------------------
%   Created by Endrina Rivas
%       endrina.rivas@uwaterloo.ca
%       Department of Civil Engineering
%       University of Waterloo
%       October 2015
%   ----------------------------------------------------------------------

% Partition stiffness matrix
KEE = K(fixed, fixed);
KEF = K(fixed, free);

% Apply boundary conditions
K(fixed,fixed) = eye(length(fixed))*beta;
F(fixed) = beta*dfixed;

% Solve for displacements
d = K\F;

% Apply displacement boundary conditions
d(fixed) = dfixed;

% Solve for unknown reaction forces
dE = d(fixed);
dF = d(free);
fE = KEE*dE + KEF*dF;

end