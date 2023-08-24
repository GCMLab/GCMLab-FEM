function Mesh = NodeSetsCoupled(Mesh)
%NODESETSCOUPLED define sets of nodes and DOFs in the domain for coupled
% problems
%   Mesh = NODESETSCOUPLED(Mesh) updates the structure array containing mesh
%   information with relevant node sets for a coupled problem
% 
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   Mesh:
%       .nsd        number of spatial dimensions
%       .x          array of nodal spatial locations for undeformed mesh 
%                   (size nn x nsd in which nn is the number of nodes)
% 
%   --------------------------------------------------------------------
%   Output
%   --------------------------------------------------------------------
%   The function returns the Mesh structure array with new fields,
%       .left_dof_u       DOFs on the left edge of the domain
%       .left_dof_t       DOFs on the left edge of the domain
%       .right_dof_u      DOFs on the right edge of the domain
%       .right_dof_t      DOFs on the right edge of the domain
%       .xdofs_u          DOFs in the x-direction
%       .ydofs_u          DOFs in the y-direction
%       .zdofs_u          DOFs in the z-direction
%       .xdofs_t          DOFs in the x-direction
%       .ydofs_t          DOFs in the y-direction
%       .zdofs_t          DOFs in the z-direction
%   Two-dimensional meshes also contain the new fields,
%       .top_dof_u        DOFs on the top edge of the domain
%       .top_dof_ux       DOFs on the top boundary in the x-direction
%       .top_dof_uy       DOFs on the top boundary in the y-direction
%       .bottom_dof_u     DOFs on the bottom edge of the domain
%       .bottom_dof_ux    DOFs on the bottom boundary in the x-direction
%       .bottom_dof_uy    DOFs on the bottom boundary in the y-direction
%       .left_dof_ux      DOFs on the left boundary in the x-direction
%       .left_dof_uy      DOFs on the left boundary in the y-direction
%       .right_dof_ux     DOFs on the right boundary in the x-direction
%       .right_dof_uy     DOFs on the right boundary in the y-direction
%       .top_dof_t        DOFs on the top edge of the domain
%       .bottom_dof_t     DOFs on the bottom edge of the domain
%       .left_dof_t      DOFs on the left boundary of the domain
%       .right_dof_t     DOFs on the right boundary of the domain
%   Three-dimensional meshes also contain the new fields, 
%       .near_dof_u       DOFs on the nearest face of the domain
%       .near_dof_ux      DOFs on the near face in the x-direction
%       .near_dof_uy      DOFs on the near face in the y-direction
%       .near_dof_uz      DOFs on the near face in the z-direction
%       .far_dof_u        DOFs on the farthest face of the domain
%       .far_dof_ux       DOFs on the far face in the x-direction
%       .far_dof_uy       DOFs on the far face in the y-direction
%       .far_dof_uz       DOFs on the far face in the z-direction
%       .left_dof_uz      DOFs on the left face in the z-direction
%       .right_dof_uz     DOFs on the right face in the z-direction
%       .top_dof_uz       DOFs on the top face in the z-direction
%       .bottom_dof_uz    DOFs on the bottom face in the z-direction

%       .near_dof_t       DOFs on the nearest face of the domain
%       .far_dof_t        DOFs on the farthest face of the domain
%       .left_dof_t      DOFs on the left face in the z-direction
%       .right_dof_t     DOFs on the right face in the z-direction
%       .top_dof_t       DOFs on the top face in the z-direction
%       .bottom_dof_t    DOFs on the bottom face in the z-direction

switch Mesh.nsd
    case 1
      Mesh.left_dof_u = Mesh.left_dof(1:2:end);      
      Mesh.left_dof_t = Mesh.left_dof(2:2:end);      
      Mesh.right_dof_u = Mesh.right_dof(1:2:end);
      Mesh.right_dof_t = Mesh.right_dof(2:2:end);
      Mesh.xdofs_u = Mesh.xdofs(1:2:end);
      Mesh.ydofs_u = Mesh.ydofs(1:2:end);
      Mesh.zdofs_u = Mesh.zdofs(1:2:end);
      Mesh.xdofs_t = Mesh.xdofs(2:2:end);
      Mesh.ydofs_t = Mesh.ydofs(2:2:end);
      Mesh.zdofs_t = Mesh.zdofs(2:2:end);
      
    case 2
      
      Mesh.xdofs_u = 1:3:Mesh.nDOF;
      Mesh.ydofs_u = 2:3:Mesh.nDOF;
      Mesh.dofs_t = 3:3:Mesh.nDOF;

      Mesh.top_dof_ux = Mesh.top_dof(1:3:end);
      Mesh.top_dof_uy = Mesh.top_dof(2:3:end);
      Mesh.top_dof_u = sort([Mesh.top_dof_ux; Mesh.top_dof_uy]);
      
      Mesh.bottom_dof_ux = Mesh.bottom_dof(1:3:end);
      Mesh.bottom_dof_uy = Mesh.bottom_dof(2:3:end);
      Mesh.bottom_dof_u  = sort([Mesh.bottom_dof_ux; Mesh.bottom_dof_uy]);
      
      Mesh.left_dof_ux = Mesh.left_dof(1:3:end);
      Mesh.left_dof_uy = Mesh.left_dof(2:3:end);
      Mesh.left_dof_u  = sort([Mesh.left_dof_ux; Mesh.left_dof_uy]);
      
      Mesh.right_dof_ux = Mesh.right_dof(1:3:end);
      Mesh.right_dof_uy = Mesh.right_dof(2:3:end);
      Mesh.right_dof_u  = sort([Mesh.right_dof_ux; Mesh.right_dof_uy]);
      
      Mesh.top_dof_t = Mesh.top_dof(3:3:end);
      Mesh.bottom_dof_t = Mesh.bottom_dof(3:3:end);
      Mesh.left_dof_t = Mesh.left_dof(3:3:end);
      Mesh.right_dof_t = Mesh.right_dof(3:3:end);
        
    case 3
        
      Mesh.near_dof_ux = Mesh.near_dof(1:4:end);
      Mesh.near_dof_uy = Mesh.near_dof(2:4:end);
      Mesh.near_dof_uz = Mesh.near_dof(3:4:end);
      Mesh.near_dof_u = sort([Mesh.near_dof_ux; Mesh.near_dof_uy; Mesh.near_dof_uz]);
      
      Mesh.far_dof_ux = Mesh.far_dof(1:4:end);
      Mesh.far_dof_uy = Mesh.far_dof(2:4:end);
      Mesh.far_dof_uz = Mesh.far_dof(3:4:end);
      Mesh.far_dof_u = sort([Mesh.far_dof_ux; Mesh.far_dof_uy; Mesh.far_dof_uz]);
      
      Mesh.left_dof_uz = Mesh.left_dof(3:4:end);
      Mesh.right_dof_uz = Mesh.right_dof(3:4:end);
      Mesh.top_dof_uz = Mesh.top_dof(3:4:end);
      Mesh.bottom_dof_uz = Mesh.bottom_dof(3:4:end);

      Mesh.near_dof_t = Mesh.near_dof(4:4:end);
      Mesh.far_dof_t = Mesh.far_dof(4:4:end);
      Mesh.left_dof_t = Mesh.left_dof(4:4:end);
      Mesh.right_dof_t = Mesh.right_dof(4:4:end);
      Mesh.top_dof_t = Mesh.top_dof(4:4:end);
      Mesh.bottom_dof_t = Mesh.bottom_dof(4:4:end);

end


end