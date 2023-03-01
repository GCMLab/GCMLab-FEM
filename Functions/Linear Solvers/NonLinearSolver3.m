function [d, fE] = NonLinearSolver3(K, F, dfixed, free, fixed, beta, parallel)
%NONLINEARSOLVER3 Solve for unknown vector using Penalty method
% 	d = NONLINEARSOLVER3(K, F, dfixed, free, fixed) is the column vector of 
% 	displacements (size ndof x 1 where ndof is the number of degrees of 
% 	freedom) for a given stiffness matrix (K) and residual force vector 
% 	(F). The free DOFs and fixed DOFs are given in the vectors free and 
% 	fixed, respectively. The displacements on the fixed DOFs are stored 
% 	in the vector dfixed. 
% 
% 	[d, fE] = NONLINEARSOLVER3(K, F, dfixed, free, fixed) also returns the 
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
% 	beta: 		Penalty parameter (large value)

% Partition stiffness matrix
KEE = K(fixed, fixed);
KEF = K(fixed, free);

% Apply boundary conditions
K(fixed,fixed) = eye(length(fixed))*beta;
F(fixed) = beta*dfixed;

% Solve for displacements
d = MatrixInvert(K,F,parallel);

% Apply displacement boundary conditions
d(fixed) = dfixed;

% Solve for unknown reaction forces
dE = d(fixed);
dF = d(free);
fE = KEE*dE + KEF*dF;

end