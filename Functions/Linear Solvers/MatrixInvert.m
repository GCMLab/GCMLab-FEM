function [d] = MatrixInvert(K,f,parallel_number)
%MATRIXINVERT Solves the linear system d = inv(K)*f
%   d = MatrixInvert(K,f,Control) is an ndof x 1 column vector of nodal 
%   displacements found by solving the system of linear equations 
%   d = inv(K)*f = K\f whether on a single core or in parallel

%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
% 	K:                  Partitioned stiffness matrix 
% 	f:                  Partitioned force/residual vector
%   parallel_number:    Number of cores on which to perform the inversion.
%                       1   - perform on a single core
%                       > 1 - creates a local pool and inverts in parallel
%                             if inverting in parallel, use as many cores
%                             as you have access to.
 
if parallel_number == 1
    d = K\f;
else
    % Create parallel processing pool
        parpool(parallel_number);
    % distribute to pool
        K = distributed(K);
        f = distributed(f);
    % Invert matrix
        d = K\f;
    % Return distributed to double
        d = gather(d);
end

