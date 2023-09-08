function [Material, stiffnessmatrixfile_name, stressstrainfile_name] = setMaterialModel(Material)
%SETMATERIALMODEL - creates pointers to material model files
%           Sets Material structure with the following fields
%    .ConstitutiveLawFile:           Constitutive model file
%    .StiffnessMatrixFile:           Stiffness matrix integration file
%    .DampingFile:                   Damping matrix integration file
%    .StressStrainFile:              File used for postprocessing results
%    .PostProcessor:                 File which outputs results to vtk
%    .ProblemType:                   Type of problem, specifies DoFs per node
%           1 - Equilibrium / Deformation Problem
%           2 - Diffusion Problem
%           3 - Thermo-Elastic Problem (to be implemented)
%    .TimeType:                      Type of time derivative associated
%                                    with the problem
%           0 - Static, no time derivatives
%           1 - First-order time derivative, rate of change/velocity, transient problems
%           2 - Second-order time derivative, acceleration
%    

    switch Material.Model
     %%%%%%%%%% Linear Elastic Models
            % Quasi-Static linear elasticity
            case 'LE1' 
                Material.ConstitutiveLawFile    = 'getD';
                Material.StiffnessMatrixFile    = 'getK_LE1';
                Material.StressStrainFile       = 'getStrain';
                Material.PostProcessor          = 'write2vtk_eqbm';
                Material.ProblemType            = 1;
                Material.TimeType               = 0;
            
            % Transient linear elasticity 
            case 'LET1' 
                Material.ConstitutiveLawFile    = 'getD';
                Material.StiffnessMatrixFile    = 'getK_LET1';
                Material.DampingFile            = 'getC'; 
                Material.StressStrainFile       = 'getStrain';
                Material.PostProcessor          = 'write2vtk_eqbm';
                Material.ProblemType            = 1;
                Material.TimeType               = 1;
                
            % Linear elasticity and dynamic     
            case 'LED1' 
                Material.ConstitutiveLawFile    = 'getD';
                Material.StiffnessMatrixFile    = 'getK_LED1';
                Material.DampingFile            = 'getC'; % Mass-based damping assumed
                Material.StressStrainFile       = 'getStrain';
                Material.PostProcessor          = 'write2vtk_eqbm';                
                Material.ProblemType            = 1;
                Material.TimeType               = 2;
                
            % Linear Viscoelastic Kelvin-Voigt Model    
            case 'VE1' 
                Material.ConstitutiveLawFile    = 'getD';
                Material.StiffnessMatrixFile    = 'getK_LET1';
                Material.DampingFile            = 'getC_VE1'; 
                Material.StressStrainFile       = 'getStrain';
                Material.PostProcessor          = 'write2vtk_eqbm';
                Material.ProblemType            = 1;
                Material.TimeType               = 1;
                
                
     %%%%%%%%%% Non-linear Models 
            % Stiffening Model via 1st Invariant of strain
            case 'ST1' 
                Material.ConstitutiveLawFile    = 'getD_ST1';
                Material.StiffnessMatrixFile    = 'getK_ST1';
                Material.DampingFile            = 'getC'; 
                Material.StressStrainFile       = 'getStrain_ST1';
                Material.PostProcessor          = 'write2vtk_eqbm';
                Material.ProblemType            = 1;
                Material.TimeType               = 0;
                
            % Softening Model via 1st Invariant of strain
            case 'ST2' 
                Material.ConstitutiveLawFile    = 'getD_ST2';
                Material.StiffnessMatrixFile    = 'getK_ST1';
                Material.DampingFile            = 'getC'; 
                Material.StressStrainFile       = 'getStrain_ST1';
                Material.PostProcessor          = 'write2vtk_eqbm';
                Material.ProblemType            = 1;
                Material.TimeType               = 0;
                
            % Transient model with stiffening model via 1st invariant of strain
            case 'TR2' 
                Material.ConstitutiveLawFile    = 'getD_ST1';
                Material.StiffnessMatrixFile    = 'getK_TR2';
                Material.DampingFile            = 'getC'; 
                Material.StressStrainFile       = 'getStrain_ST1';
                Material.PostProcessor          = 'write2vtk_eqbm';
                Material.ProblemType            = 1;
                Material.TimeType               = 1;
               
                
     %%%%%%%%%% Diffusion Models
            % Thermal Diffusion (Steady-state)
            case 'TH1' 
                Material.ConstitutiveLawFile    = 'getD_TH1';
                Material.StiffnessMatrixFile    = 'getK_TH1';
                Material.DampingFile            = 'getC'; 
                Material.StressStrainFile       = 'getFlux_TH1';
                Material.PostProcessor          = 'write2vtk_dfsn';
                Material.ProblemType            = 2;
                Material.TimeType               = 0;
            % Thermal Diffusion (Transient)
            case 'TH2' 
                Material.ConstitutiveLawFile    = 'getD_TH1';
                Material.StiffnessMatrixFile    = 'getK_TH1';
                Material.DampingFile            = 'getC'; 
                Material.StressStrainFile       = 'getFlux_TH1';
                Material.PostProcessor          = 'write2vtk_dfsn';
                Material.ProblemType            = 2;
                Material.TimeType               = 1;
            % Thermal Nonlinear Diffusion (Transient)
            case 'NLTH1'
                Material.ConstitutiveLawFile    = 'getD_NLTH1';
                Material.StiffnessMatrixFile    = 'getK_NLTH1';
                Material.DampingFile            = 'getC';
                Material.StressStrainFile       = 'getFlux_TH1';
                Material.PostProcessor          = 'write2vtk_dfsn';
                Material.ProblemType            = 2;
                Material.TimeType               = 1;
    end

    % Create function pointers
        % Tangent Matrix
        [~,stiffnessmatrixfile_name] = fileparts(Material.StiffnessMatrixFile);

        % Stress/strain 
        [~,stressstrainfile_name] = fileparts(Material.StressStrainFile);
        
    % D matrix is called inside assembly functions
end

