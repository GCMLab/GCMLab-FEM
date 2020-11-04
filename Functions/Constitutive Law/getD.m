function D = getD(E, nu, nsd, Dtype)
%GETD Elasticity tensor
%   D = GETD(E, nu, nsd) is the elasticity tensor for a problem of 
%   spatial dimension, nsd, Young's modulus, E, and Poisson's ratio, nu. 
%   For a 1D problem, a scalar is returned. For a 2D problem, a 3x3 
%   matrix is returned, and a for a 3D problem, a 6x6 matrix is returned. 
% 
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   E:      Modulus of elasticity
%   nu:     Poisson's ratio
%   nsd:    Number of spatial dimensions

switch nsd
    case 1
        D = E;
    case 2                 
        switch Dtype
            case 'PlaneStrain'
                D  = E/((1+nu)*(1-2*nu))*[1-nu  nu   0     ;
                                          nu    1-nu 0     ;
                                          0     0    0.5-nu];
            case 'PlaneStress'
                D  = E/(1-nu^2)*[1  nu 0       ;
                                 nu 1  0       ;
                                 0  0  (1-nu)/2];
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