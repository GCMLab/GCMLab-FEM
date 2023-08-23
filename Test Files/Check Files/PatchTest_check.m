function [disp_er, stress_er, reaction_er] = PatchTest_check(d, stress, Fext, Mesh, Material, BC)
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
% ux = (1-nu)*traction/E*x
% uy = (1-nu)*traction/E*y
%   applied traction of traction in both x and y directions
%   stress = [traction*ones(1,Mesh.nn); traction*ones(1,Mesh.nn); zeros(1,Mesh.nn)]

% Calculate exact solutions
sigma = BC.traction;
x = Mesh.x(:,1);
y = Mesh.x(:,2);
d_exact = zeros(2*Mesh.nn,1);
d_exact(1:2:end) = (1-Material.Prop(1).nu)*BC.traction/Material.Prop(1).E0.*x;
d_exact(2:2:end) = (1-Material.Prop(1).nu)*BC.traction/Material.Prop(1).E0.*y;
stress_exact = [sigma*ones(1,Mesh.nn); sigma*ones(1,Mesh.nn); zeros(1,Mesh.nn)];

% Calculate the error
disp_er = abs(norm(d - d_exact));
stress_er = abs(sqrt(sum((stress - stress_exact).^2,'all')));
reaction_er = abs(sum(Fext));

end

