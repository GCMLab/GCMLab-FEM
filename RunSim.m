%% Clear variables and initialize code
    clearvars -global
    clear, clc, close all
    format compact
    tic;
    
    % Current directory
    curDir = pwd;
 
%% Input (these variables must be modified by the user)
    % Folder name where config files are stored
    DirFolder = 'Config Files';
    % Config files to run. Choose either 'all' or give the file name.
    FileList = 'NLElastic_Transient_2DPlate';
    
    
    % Directory for VTK Files (end with \)
    if ispc
        VTKFolder ='C:\Users\b3gee\Documents\Matlab Results\';
        %VTKFolder ='C:\Users\endri\Documents\Matlab Results\';
    else
        VTKFolder = '/Users/jonathanzingaro/Documents/MATLAB Results';
    end
     
    % output vtk files
    plot2vtk = 1;
    
    % output progress messages
    progress_on = 1;

%% Directories
    VTKFolder = fullfile(VTKFolder, DirFolder);
    FuncDir = fullfile(curDir, 'Functions');
    ConfigDir = fullfile(curDir, DirFolder);

%% Add validation folder and function folder to the search path
    % genpath adds all subfolders as well
    addpath(genpath(FuncDir));
    addpath(genpath(ConfigDir));  
    
%% Get list of config files
    % Cell array containing the names of the configuration files to 
    % run. This option selects all of the files in the folder
    if strcmp(FileList, 'all')
        ConfigFiles = dir([ConfigDir,'\*.m']);
        ConfigFiles = fullfile(ConfigDir,{ConfigFiles.name});
    else
        ConfigFiles = {FileList};
    end
    VTKDirs = fullfile(VTKFolder,ConfigFiles);
    numfiles = length(ConfigFiles);
    
%% Run
%try                     
    %% Run config files
    % Loop through every configuration file 
    for file = 1:numfiles
        % clear variables from previous simulation
        clearvars -except VTKDirs ConfigFiles...
                  curDir FuncDir  ConfigDir ...
                  file codeSubmitTime ...
                  exit_when_done print_log ...
                  plot2vtk progress_on AnalysisType

        clearvars -global

        % filename
        config_name_full = ConfigFiles{file};
        [~,config_name] = fileparts(config_name_full);

        vtk_dir = VTKDirs{file};
        
        if ~isfolder(vtk_dir) 
            mkdir(vtk_dir)
        end
        
        % run and time the simulation
        start_time = toc;
        run('Functions/Main/main_nonlinear')
        end_time = toc;

        disp(['run time: ' num2str(end_time - start_time)])
        close all
    end

% catch err
%     disp(err.message);
% 
%     errStack = struct2cell(err.stack);
%     errStackName = errStack(2,:);
%     errStackLine = errStack(3,:);
% 
%     disp([errStackName' errStackLine']);
%     disp(err.identifier);
% end
