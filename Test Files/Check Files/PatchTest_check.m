function [disp_er, stress_er, reaction_er] = PatchTest_check(d, stress, Fext, Mesh)
%PATCHTEST_CHECK Calculates the error
%   [disp_er, stress_er, reaction_er] = PatchTest_check(d, stress, Fext)
%   calculates the error related to displacements, reaction forces, and
%   stress for using a patch test
%
%   ----------------------------------------------------------
%   Input
%   ----------------------------------------------------------
%   d:                  Displacement vectors
%   stress:             Nodal stress data
%   Fext:               External forces as the reactions
%   Mesh:               Mesh data structure
%
%   ----------------------------------------------------------
%   Output
%   ----------------------------------------------------------
%   disp_er:             Error related to Displacement vector
%   stress_er:           Error related to Nodal stress
%   reaction_er:         Error related to reaction forces

% Patch Test
%   ux = 0.002x
%   uy = -0.0006y
%   produces zero body forces and zero stresses except for ?(x)=2.
%   stress = [2*ones(1,Mesh.nn); zeros(1,Mesh.nn); zeros(1,Mesh.nn)]
%   Fext = [-0.5 0 0.5 0 0.5 0 -0.5 0 0 0 1 0 0 0 -1 0 0 0 0 0 0 0]'



% Calculate exact solutions
sigma = 2;
x = Mesh.x(:,1);
y = Mesh.x(:,2);
d_exact = zeros(2*Mesh.nn,1);
d_exact(1:2:end) = 0.002.*x;
d_exact(2:2:end) = -0.0006.*y;
Fext_exact = [-0.5 0 0.5 0 0.5 0 -0.5 0 0 0 1 0 0 0 -1 0 0 0 0 0 0 0]';
stress_exact = [sigma*ones(1,Mesh.nn); zeros(1,Mesh.nn); zeros(1,Mesh.nn)];

% Calculate the error
disp_er = norm(d - d_exact);
stress_er = sqrt(sum((stress - stress_exact).^2,'all'));
reaction_er = norm(Fext - Fext_exact);

end

