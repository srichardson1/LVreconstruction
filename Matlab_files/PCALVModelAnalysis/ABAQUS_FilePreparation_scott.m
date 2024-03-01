%%%find the proper directory 
% LVWM_config;%%hao yu jue
LVWM_config_trial;
if ispc
    path(path, '.\HexMeshProcessing_not_using');
else
    path(path, './HexMeshProcessing_not_using')
end



abaqusDir_simulation = 'abaqusSimulation';
cd(resultDir);
if ~exist(abaqusDir_simulation, 'dir')
    mkdir(abaqusDir_simulation); 
    cd(abaqusDir_simulation);
    abaqusDir_simulation = pwd();
    cd(workingDir);
else
    cd(abaqusDir_simulation);
    abaqusDir_simulation = pwd();
    cd(workingDir);
end
    
abaqusDir_reconstruction = resultDir;


LV_WholeMesh_abaqusFilePreparation(abaqusDir_reconstruction,abaqusDir_simulation);



FiberGenerationDir = 'fibreGeneration'; 
cd(abaqusDir_reconstruction);
cd(FiberGenerationDir);
FiberGenerationDir = pwd(); %% this is the root folder for fibre generation
FiberResults = 'Results_fiber_60_45';
if ~exist(FiberResults, 'dir')
    mkdir(FiberResults);
end
cd(FiberResults);
FiberResults = pwd();
cd(workingDir);

a_fibre_endo = -60*pi/180; %% negative in the endocardium
a_fibre_epi = 60*pi/180; %% positive in the epicardium
a_sheet = 45.0*pi/180; %% same angle at both endocardium and epicardium
wall_thickness_calculated = 0;
LV_WholeMesh_FibreGeneration(a_fibre_endo, a_fibre_epi, a_sheet, ...
                                          wall_thickness_calculated, ...
                                          FiberGenerationDir,FiberResults,abaqusDir_simulation);
                                      
                                      



