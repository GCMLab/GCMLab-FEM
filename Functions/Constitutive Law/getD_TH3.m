function [D, alpha] = getD_TH3(nMat, Material, Mesh, dG)
%GETD_ST3 Constitutive Law for Nonlinear Diffusion Problem
%   D = GETD_ST3(nMat, Material, Mesh) is the nonlinear diffusion tensor
%   for a problem in which the constitutive law is a nonlinear relation
%   given by the following:
%   k(T) = k0+αT^n
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   nMat:       Material type
%   Mesh:       Structure array with the following fields, may be updated
%               with new fields
%               .nsd:   Number of spatial dimensions
%   Material:   Structure array with the following fields, may be updated
%               with new fields
%               .Dtype: constitutive law for two-dimensional elasticity
%               .E: Modulus of elasticity
%               .nu: Poisson's ratio
%   dG:         temperature at Gaussian point of element e
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
