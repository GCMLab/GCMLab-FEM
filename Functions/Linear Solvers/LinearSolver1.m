function [d,fE] = LinearSolver1(K, F, dfixed, nDOF, free, fixed)
%SOLVEMTETHOD1 Solves for unknown vector using partition method
%   d = SolveMethod1(nodes,fixtempnode,fixtempvalue,K,F)
%
%   ----------------------------------------------------------------------
%   Created by Endrina Rivas
%       endrina.rivas@uwaterloo.ca
%       Department of Civil Engineering
%       University of Waterloo
%       November 2015
%   ----------------------------------------------------------------------

% Identify extra free nodes for LM method
if size(K,1) > nDOF
    free = [free nDOF+1:size(K,1)];
end

% Partition matrices
KEE = K(fixed, fixed);
KFF = K(free, free);
KEF = K(fixed, free);
KFE = KEF';
fF = F(free);

% Apply boundary conditions
dE = dfixed;

% Solve for unknown displacements
% This if statement allows for problems that do not have any fixed boundary
% conditions
if isempty(dE)
    dF = KFF\(fF);
else
    dF = KFF\(fF-KFE*dE);
end

% Update displacement vector
d(free, 1) = dF;
d(fixed, 1) = dE;

% Solve for reaction forces
% This if statement allows for problems that do not have any fixed boundary
% conditions
if isempty(dE)
    fE = [];
else
    fE = KEE*dE + KEF*dF;
end

end