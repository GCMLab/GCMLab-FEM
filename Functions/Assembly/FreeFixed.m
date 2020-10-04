function BC = FreeFixed(BC, meshDOFs)
%% Identify free and fixed dofs

%   ----------------------------------------------------------------------
%   Created by Endrina Rivas
%       endrina.rivas@uwaterloo.ca
%       Department of Civil Engineering
%       University of Waterloo
%       March 2015
%	Last updated June 2016
%   ----------------------------------------------------------------------

BC.fixed = BC.fix_disp_dof;
BC.free = setdiff(meshDOFs, BC.fixed)';