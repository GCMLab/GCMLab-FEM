function Mesh = NodeDOFs(Mesh, problemtype)
%NODEDOFS Define degrees of freedom (DOFs) in a mesh
% 	Mesh = NODEDOFS(Mesh) updates the structure array containing mesh
%   information with an array of degree of freedom indices
% 
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
% 	Mesh: structure array with the following fields,
% 		.nsd: 	number of spatial dimensions
% 		.nne: 	number of nodes per element
% 		.nn: 	total number of nodes
%       .problemtype:    
%           The type of problem being solved - sets the number of DoFs per
%           node. 1 - Equilibrium: nDoF = nsd, 2 - Diffusion: nDoF = 1,
%                 3 - Mixed: nDoF = nsd + 1
% 	
%   --------------------------------------------------------------------
%   Output
%   --------------------------------------------------------------------
% 	The function returns the Mesh structure array with new fields,
%       .nDOFn: number of DOFs per node
% 		.nDOFe:	number of DOFs per element
% 		.nDOF: 	total number of DOFs
% 		.DOF: 	array of DOF indices (size nn x nsd)

% Acknowledgemnents: Chris Ladubec, Matin Parchei-Esfahani
    switch problemtype
        case 1  % Equilibrium
            Mesh.nDOFn = Mesh.nsd;
        case 2  % Diffusion
            Mesh.nDOFn = 1;
        case 3  % Mixed
            Mesh.nDOFn = Mesh.nsd+1;
    end


    Mesh.nDOFe = Mesh.nne*Mesh.nDOFn;         % number of DOF per element
    Mesh.nDOF = Mesh.nn*Mesh.nDOFn;           % total number of DOF
    Mesh.DOF = zeros(Mesh.nn, Mesh.nDOFn); 
    


    for sd = 1:Mesh.nDOFn
       Mesh.DOF(:,sd) = (sd : Mesh.nDOFn : (Mesh.nDOF-(Mesh.nDOFn-sd)))';
    end

end