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
% ux = (1-nu)*t/E*x
% uy = (1-nu)*t/E*y
%   applied traction of t in both x and y directions
%   stress = [t*ones(1,Mesh.nn); t*ones(1,Mesh.nn); zeros(1,Mesh.nn)]

global E nu t

% Calculate exact solutions
sigma = t;
x = Mesh.x(:,1);
y = Mesh.x(:,2);
d_exact = zeros(2*Mesh.nn,1);
d_exact(1:2:end) = (1-nu)*t/E.*x;
d_exact(2:2:end) = (1-nu)*t/E.*y;
stress_exact = [sigma*ones(1,Mesh.nn); sigma*ones(1,Mesh.nn); zeros(1,Mesh.nn)];

% Calculate the error
disp_er = norm(d - d_exact);
stress_er = sqrt(sum((stress - stress_exact).^2,'all'));
reaction_er = sum(Fext);

end

