function [d_er, reaction_er] = Thermoelasticity_check(d, Fext, Mesh)
%THERMOELASTICITY_CHECK Calculates the error
%   [temp_er, stress_er, reaction_er] = Thermoelasticity_check(d, Fext, Mesh)
%   calculates the error related to displacements/temperatures, and reaction forces
%
%   ----------------------------------------------------------
%   Input
%   ----------------------------------------------------------
%   d:                  Displacement/temperature vector
%   Fext:               External forces as the reactions
%   Mesh:               Mesh data structure
%
%   ----------------------------------------------------------
%   Output
%   ----------------------------------------------------------
%   d_er:               Error related to displacement/temperature vector
%   reaction_er:        Error related to reaction forces

% Thermoelasticity test
% u:= (x1,x2) → -sin(pi*x1/2)*sin(pi*x2/2)
% v:= (x1,x2) → cos(pi*x1/2)*cos(pi*x2/2)
% T:= (x1,x2) → sin(pi*x1)*sin(pi*x2)

% Calculate exact solutions
d_exact = zeros(Mesh.nDOF,1);
aux_x = @(x) - sin(pi.*x(:,1)./2).*sin(pi.*x(:,2)./2)./1000;
aux_y = @(x) cos(pi.*x(:,1)./2).*cos(pi.*x(:,2)./2)./1000;
aux_t = @(x) sin(pi.*x(:,1)).*sin(pi.*x(:,2))./1000;

d_exact(1:3:end) = aux_x(Mesh.x);
d_exact(2:3:end) = aux_y(Mesh.x);
d_exact(3:3:end) = aux_t(Mesh.x);

% Calculate the error
d_er = abs(norm(d - d_exact));
reaction_er = abs(sum(Fext));

end

