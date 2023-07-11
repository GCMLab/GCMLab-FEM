function [Material, stiffnessmatrixfile_name, stressstrainfile_name, Control] = setMaterialModel(Material, Control)
%SETMATERIALMODEL - creates pointers to material model files
%           Sets Material
%    .ConstitutiveLawFile:           Constitutive model file
%    .StiffnessMatrixFile:           Stiffness matrix integration file
%    .StressStrainFile:              File used for postprocessing results
%    .PostProcessor:                 File which outputs results to vtk
%    .ProblemType:                   Type of problem, specifies DoFs per node
%           1 - Equilibrium / Deformation Problem
%           2 - Diffusion Problem
%           3 - Thermo-Elastic Problem (to be implemented)
%           

    switch Material.Model
        case 'LE1' % Linear elasticity
            Material.ConstitutiveLawFile = 'getD';
            Material.StiffnessMatrixFile = 'getK_LE1';
            Material.StressStrainFile = 'getStrain';
            Material.PostProcessor    = 'write2vtk_eqbm';
            Material.ProblemType      = 1;
            Material.Transient        = 0;
            Control.TimeCase = 'static'; 
        case 'LED1' % Linear elasticity and dynamic
            Material.ConstitutiveLawFile = 'getD';
            Material.StiffnessMatrixFile = 'getK_LE1_dynamic';
            Material.StressStrainFile = 'getStrain';
            Material.DampingFile = 'getC'; 
            Control.TimeCase = 'dynamic'; 
        case 'LET1' % Linear Elastic Transient
            Material.ConstitutiveLawFile = 'getD';
            Material.StiffnessMatrixFile = 'getK_LET1';
            Material.StressStrainFile = 'getStrain';
            Material.DampingFile = 'getC'; 
            Control.TimeCase = 'transient'; 
        case 'ST1' % Stiffening Model via 1st Invariant of strain
            Material.ConstitutiveLawFile = 'getD_ST1';
            Material.StiffnessMatrixFile = 'getK_ST1';
            Material.StressStrainFile = 'getStrain_ST1';
            Material.PostProcessor    = 'write2vtk_eqbm';
            Material.ProblemType      = 1;
            Material.Transient        = 0;
        case 'TR1' % Transient Linear Elastic
            Material.ConstitutiveLawFile = 'getD';
            Material.StiffnessMatrixFile = 'getK_TR1';
            Material.StressStrainFile = 'getStrain';
            Material.PostProcessor    = 'write2vtk_eqbm';
            Material.ProblemType      = 1;
            Material.Transient        = 1;
        case 'TH1' % Thermal Diffusion (Steady-state)
            Material.ConstitutiveLawFile = 'getD_TH1';
            Material.StiffnessMatrixFile = 'getK_TH1';
            Material.StressStrainFile = 'getFlux_TH1';
            Material.PostProcessor    = 'write2vtk_dfsn';
            Material.ProblemType      = 2;
            Material.Transient        = 0;
        case 'TH2' % Thermal Diffusion (Transient)
            Material.ConstitutiveLawFile = 'getD_TH1';
            Material.StiffnessMatrixFile = 'getK_TH1';
            Material.StressStrainFile = 'getFlux_TH1';
            Material.PostProcessor    = 'write2vtk_dfsn';
            Material.ProblemType      = 2;
            Material.Transient        = 1;
            Control.TimeCase = 'static'; 
        case 'ST2' % Softening Model via 1st Invariant of strain
            Material.ConstitutiveLawFile = 'getD_ST2';
            Material.StiffnessMatrixFile = 'getK_ST1';
            Material.StressStrainFile = 'getStrain_ST1';
			Material.PostProcessor    = 'write2vtk_eqbm';
            Material.ProblemType      = 1;
            Material.Transient        = 1;
            Control.TimeCase = 'static'; 
        case 'TR2' % Transient model with stiffening model via 1st invariant of strain
            Material.ConstitutiveLawFile = 'getD_ST1';
            Material.StiffnessMatrixFile = 'getK_TR2';
            Material.StressStrainFile = 'getStrain_ST1';
			Material.PostProcessor    = 'write2vtk_eqbm';
            Material.ProblemType      = 1;
            Material.Transient        = 1;
            Control.TimeCase = 'transient'; 
        case 'VE1' % Linear Viscoelastic Kelvin-Voigt Model
            Material.ConstitutiveLawFile = 'getD';
            Material.StiffnessMatrixFile = 'getK_VE1';
            Material.StressStrainFile = 'getStrain';
            Material.DampingFile = 'getC_VE1'; 
            Control.TimeCase = 'transient'; 
    end

    % Create function pointers
        % Tangent Matrix
        [~,stiffnessmatrixfile_name] = fileparts(Material.StiffnessMatrixFile);

        % Stress/strain 
        [~,stressstrainfile_name] = fileparts(Material.StressStrainFile);
        
    % D matrix is called inside assembly functions
end

