function BC = FreeFixed(BC, meshDOFs)
%% Identify free and fixed dofs

BC.fixed = BC.fix_disp_dof;
BC.free = setdiff(meshDOFs, BC.fixed)';

end