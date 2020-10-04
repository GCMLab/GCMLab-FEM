function D = getD(x, nsd, Material)
%GETD 

%   ----------------------------------------------------------------------
%   Created by Endrina Rivas & Matin Parchei Esfahani
%       endrina.rivas@uwaterloo.ca
%       Department of Civil Engineering
%       University of Waterloo
%       October 2015
% 	Last Updated: June 2016
%   ----------------------------------------------------------------------

E = Material.E(x);    

switch nsd
    case 1
        D = E;
    case 2
        nu = Material.nu(x);
                 
        switch Material.Dtype
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
        nu = Material.nu(x);
        D = E/(1+nu)/(1-2*nu)*[1-nu nu nu 0 0 0;
                                nu 1-nu nu 0 0 0;
                                nu nu 1-nu 0 0 0;
                                0 0 0 (1-2*nu)/2 0 0;
                                0 0 0 0 (1-2*nu)/2 0;
                                0 0 0 0 0 (1-2*nu)/2];
        
end

end