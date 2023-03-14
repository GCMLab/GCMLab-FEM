function [d, fE] = LinearSolver1(K, F, dfixed, free, fixed, parallel)
%LINEARSOLVER1 Solves for unknown vector using partition method
% 	d = LINEARSOLVER1(K, F, dfixed, free, fixed) is the displacement 
% 	vector for a given stiffness matrix, K, and residual force vector, F.
% 	The free DOFs and fixed DOFs are given in the vectors free and fixed, 
% 	respectively. The displacements on the fixed DOFs are stored in the 
% 	vector dfixed. 
% 
% 	[d, fE] = LINEARSOLVER1(K, F, dfixed, free, fixed) also returns the 
% 	reaction forces on the fixed DOFs, provided as a column vector 
% 	(size nfixed x 1 in which nfixed is the number of fixed degrees of 
% 	freedom).
% 
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
% 	K:			Stiffness matrix (size )
% 	F: 			Column vector of residual forces (size ndof x 1)
% 	dfixed: 	Column vector of displacements on fixed degrees of 
% 				freedom (size nfixed x 1)
% 	free: 		Row vector of free degrees of freedom (size 1 x nfree in
% 				which nfree is the number of free degrees of freedom)
% 	fixed:		Row vector of fixed degrees of freedom (size 1 x nfixed)

% Identify extra free nodes
if size(K,1) > length(F)
    free = [free, length(F)+1:size(K,1)];
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
    dF = MatrixInvert(KFF,fF,parallel);
else
    dF = MatrixInvert(KFF,fF-KFE*dE,parallel);
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