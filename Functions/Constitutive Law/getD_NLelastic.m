function [D, Material] = getD_NLelastic(nMat, Material, Mesh, strain_e)
%GETD Elasticity tensor
%   D = GETD_NLelastic(E, nu, nsd) is the elasticity tensor for a problem 
%   in which the constitutive law is a non linear elastic relation given by
%   E = E0 + E1*I1^2
% 
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
%   strain_e:   element strains
%   --------------------------------------------------------------------

E0 = Material.Prop(nMat).E0;
E1 = Material.Prop(nMat).E1;
nu = Material.Prop(nMat).nu;

% constitutive matrix
switch Mesh.nsd
    case 1
        % strain invariant
        I1 = strain_e(1,1)^2;
        % elasticity modulus
        E = E0 + E1*I1^2;
        % constitutive matrix
        D = E;
    case 2    
        % strain invariant
        I1 = (strain_e(1,1)+strain_e(2,1))^2;
        % elasticity modulus
        E = E0 + E1*I1^2;
        % constitutive matrix
        switch Material.Dtype
            case 'PlaneStrain'
                D  = E/((1+nu)*(1-2*nu))*[1-nu  nu   0     ;
                                          nu    1-nu 0     ;
                                          0     0    0.5-nu];
            case 'PlaneStress'
                D  = E/(1-nu^2)*[1  nu 0       ;
                                 nu 1  0       ;
                                 0  0  (1-nu)/2];
            otherwise
                error('Material.Dtype is not correctly defined.')
        end
    case 3
        % strain invariant
        I1 = (strain_e(1,1)+strain_e(2,1)+ strain_e(3,1))^2;
        % elasticity modulus
        E = E0 + E1*I1^2;
        % constitutive matrix
        D = E/(1+nu)/(1-2*nu)*[1-nu nu nu 0 0 0;
                                nu 1-nu nu 0 0 0;
                                nu nu 1-nu 0 0 0;
                                0 0 0 (1-2*nu)/2 0 0;
                                0 0 0 0 (1-2*nu)/2 0;
                                0 0 0 0 0 (1-2*nu)/2];
end

% Update Young's modulus
Material.Prop(nMat).E = E;

end