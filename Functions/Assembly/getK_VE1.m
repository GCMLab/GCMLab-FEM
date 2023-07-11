<<<<<<<< HEAD:Functions/Assembly/getK_TR1.m
function [K, R, Fint] = getK_TR1(~, ~, ~, ~, Fext, Fextnm1, Klin, ~, d, dnm1, ~, dt, ~,C,alpha)
%GETK_TR1 Stiffness matrix for iterative elastic case
%   [K, R, Fint] = GETK_TR1(Mesh, Quad, Material) returns the stiffness
========
function [K, R, Fint] = getK_VE1(Mesh, Quad, Material, ~, Fext, Fextnm1, Klin, M, d_m, dt, dtnm1,C,alpha)
%GETK_VE1 Stiffness matrix for iterative linear viscoelastic case
%   [K, R, Fint] = GETK_ELASTIC(Mesh, Quad, Material) returns the stiffness
>>>>>>>> 6b7568f5368fae3752142f4d29341d0e7ba7cdc1:Functions/Assembly/getK_VE1.m
%   matrix K, the residual vector R, and the internal force vector for the 
%   iterative solver where the problem uses a Kelvin-Voigt Model.
%   
%   Template file for other tangent matrix files 
%   --------------------------------------------------------------------
%   Accepted Inputs (in order)
%   --------------------------------------------------------------------
%   getK_elastic(Mesh, Quad, Material, Klin, M, d, dnm1, dnm2, stress, strain, dt, dtnm1)
%   Mesh:       Structure array with the following fields, may be updated
%               with new fields
%               .ne:    Total number of elements in the mesh
%               .nne:   Vector of number of nodes per element (size nsd x 1)
%               .nsd:   Number of spatial dimensions
%               .conn:  Array of element connectivity (size ne x nne)
%               .x:     Array of nodal spatial locations for
%                       undeformed mesh (size nn x nsd)
%               .DOF:   Array of DOF indices (size nn x nsd)
%               .nDOFe: Number of DOFs per element
%               .nDOF:  Total number of DOFs
%  
%   Quad:       Structure array with the following fields, may be updated
%               with new fields
%               .W:      Vector of quadrature weights (size nq x 1)      
%               .nq:     Number of quadrature points 
%               .Nq:     Cell array (size nq x 1) with shape functions  
%                        evaluated at each quadrature point
%               .dNdxiq: Cell array (size nq x 1) with derivative of shape 
%                        functions w.r.t. parent coordinates evaluated at 
%                        each quadrature point
% 
%   Material:   Structure array with the following fields, may be updated
%               with new fields
%               .t:         Material thickness
%
%   Fext:       External force vector at timestep n
%   Fextnm1:    External force vector at timestep n-1
%   Klin:       Linear elastic stiffness matrix
%   M:          Mass matrix
%   d:          unconverged degree of freedom vector at current timestep n and iteration
%   dnm1:       converged degree of freedom vector at timestep n-1
%   dnm2:       converged degree of freedom vector at timestep n-2
%   dt:         timestep size between timesteps n-1 and n
%   dtnm1:      timestep size between timesteps n-2 and n-1
%   alpha:       intagration parameter


% stiffness matrix in KV Model
K = alpha*Klin + C./dt;

% internal forces
Fint = C*(d_m.d-d_m.dnm1)./dt + Klin*(1-alpha)*d_m.dnm1 + alpha*Klin*d_m.d;

% residual
R = alpha*Fext + (1-alpha)*Fextnm1 - Fint;

end
