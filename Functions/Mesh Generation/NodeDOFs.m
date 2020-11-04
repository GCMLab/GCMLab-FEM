function Mesh = NodeDOFs(Mesh)
%NODEDOFS Define degrees of freedom (DOFs) in a mesh
% 	Mesh = NODEDOFS(Mesh) updates the structure array containing mesh
%           information with an array of degree of freedom indices
% 
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
% 	Mesh: structure array with the following fields,
% 		.nsd: 	number of spatial dimensions
% 		.nne: 	number of nodes per element
% 		.nn: 	total number of nodes
% 	
%   --------------------------------------------------------------------
%   Output
%   --------------------------------------------------------------------
% 	The function returns the Mesh structure array with new fields,
% 		.nDOFe:	number of DOFs per element
% 		.nDOF: 	total number of DOFs
% 		.DOF: 	array of DOF indices (size nn x nsd)

% Acknowledgemnents: Chris Ladubec, Matin Parchei-Esfahani

    Mesh.nDOFe = Mesh.nne*Mesh.nsd;         % number of DOF per element
    Mesh.nDOF = Mesh.nn*Mesh.nsd;           % total number of DOF
    Mesh.DOF = zeros(Mesh.nn, Mesh.nsd); 

    for sd = 1:Mesh.nsd
       Mesh.DOF(:,sd) = (sd : Mesh.nsd : (Mesh.nDOF-(Mesh.nsd-sd)))';
    end

end