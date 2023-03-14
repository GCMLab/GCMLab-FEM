function [d, fE] = LinearSolver2(K, F, dfixed, free, fixed, parallel)
%LINEARSOLVER2 Solve for unknown vector by zero-ing rows and columns in 
%stiffness matrix that correspond to essential boundaries
% 	d = LINEARSOLVER2(K, F, dfixed, free, fixed) is the displacement 
% 	vector for a given stiffness matrix, K, and residual force vector, F.
% 	The free DOFs and fixed DOFs are given in the vectors free and fixed, 
% 	respectively. The displacements on the fixed DOFs are stored in the 
% 	vector dfixed. 
% 
% 	[d, fE] = LINEARSOLVER2(K, F, dfixed, free, fixed) also returns the 
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

% Partition stiffness matrix
KEE = K(fixed,fixed);
KEF = K(fixed,free);

% Apply boundary conditions
F(fixed) = dfixed;
F(free) = F(free) - K(fixed,free)'*dfixed;
K(fixed,:) = 0;
K(:,fixed) = 0;
K(fixed,fixed) = eye(length(fixed));

% solve for unknown displacements
d = MatrixInvert(K,F,parallel);

% apply displacement boundary conditions
d(fixed) = dfixed;

% solve for reaction forces
dE = d(fixed);
dF = d(free);
fE = KEE*dE + KEF*dF;

end

