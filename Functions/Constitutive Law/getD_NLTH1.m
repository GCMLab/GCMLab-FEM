function [D, alpha] = getD_NLTH1(nMat, Material, Mesh, dG)
%getD_NLTH1 Constitutive law for nonlinear diffusion problem
%   D = getD_NLTH1(nMat, Material, Mesh, dG) is the nonlinear diffusion 
%   tensor where the constitutive law is a nonlinear relation given as
%   k(T) = k0+αT^n. 
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   nMat:       Material type
%   Mesh:       Structure array with the following fields, may be updated
%               with new fields
%               .nsd:   Number of spatial dimensions
%   Material:   Structure array with the following fields, may be updated
%               with new fields
%               .k1: 
%                   (a) Initial diffusivity coefficient for isotropic case 
%                   (b) Initial diffusivity coefficient in x-direction for
%                       orthotropic case
%               .k2: 
%                   (a) Initial diffusivity coefficient in y-direction for
%                   orthotropic case
%               .a1 
%                   (a) Linear diffusivity coefficient for isotropic case 
%                   (b) Linear diffusivity coefficient in x-direction for
%                       orthotropic case
%               .k2: 
%                   (a) Linear diffusivity coefficient in y-direction for
%                   orthotropic case
%               .n: Power diffusivity coefficient
%   dG:         temperature at Gaussian point of element e for current
%               iteration
%   --------------------------------------------------------------------
%
%   Acknowledgements : Jonathan Zingaro

% Define Diffusivity Parameters
k1 = Material.Prop(nMat).k1; % For Isotropic and 1D Cases
if isfield(Material.Prop(nMat),'k2')
    k2 = Material.Prop(nMat).k2; % Orthotropic 2D Case
end

% Define Diffusivity Constant (α and n)
a1 = Material.Prop(nMat).a1; % For Isotropic and 1D Cases
if isfield(Material.Prop(nMat),'alpha2')
    a2 = Material.Prop(nMat).a2; % Orthotropic 2D Case
end

n = Material.Prop(nMat).n; % n Constant 

% Define Constitutive Matrix based on Isotropic or Orthotropic and
% Dimensionality of System
switch Mesh.nsd
    case 1
        D = k1+a1*dG.^n;
        alpha = a1;
    case 2
        switch Material.Dtype
            case 'ISO'
               alpha = a1*eye(2);
               D = k1*eye(2)+alpha*dG.^n;
            case 'ORTHO'
                k0 = [k1,0;0,k2];
                alpha = [a1,0;0,a2];
                D = k0+alpha*dG.^n;
        end
end
end
