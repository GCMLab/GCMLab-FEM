function [d,fE] = LinearSolver2(K, F, dfixed, free, fixed)
%SOLVEMETHOD2 Solve for unknown vector by zero-ing rows and columns in 
%stiffness matrix that correspond to essential boundaries
% Acknowledgements: Endrina Rivas

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
d = K\F;

% apply displacement boundary conditions
d(fixed) = dfixed;

% solve for reaction forces
dE = d(fixed);
dF = d(free);
fE = KEE*dE + KEF*dF;

end

