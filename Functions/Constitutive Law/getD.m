function D = getD(nMat, Material, Mesh)
%GETD Elasticity tensor
%   D = GETD(nMat, Material, Mesh) is the elasticity tensor for a problem of 
%   spatial dimension, nsd, Young's modulus, E, and Poisson's ratio, nu. 
%   For a 1D problem, a scalar is returned. For a 2D problem, a 3x3 
%   matrix is returned, and for a 3D problem, a 6x6 matrix is returned. 
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
%   --------------------------------------------------------------------

E = Material.Prop(nMat).E0;
nu = Material.Prop(nMat).nu;

switch Mesh.nsd
    case 1
        D = E;
    case 2                 
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
        D = E/(1+nu)/(1-2*nu)*[1-nu nu nu 0 0 0;
                                nu 1-nu nu 0 0 0;
                                nu nu 1-nu 0 0 0;
                                0 0 0 (1-2*nu)/2 0 0;
                                0 0 0 0 (1-2*nu)/2 0;
                                0 0 0 0 0 (1-2*nu)/2];
end

end