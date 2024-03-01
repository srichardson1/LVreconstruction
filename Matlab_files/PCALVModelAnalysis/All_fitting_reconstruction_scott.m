global matlab_programme_folder_str
%% All_fitting_reconstruction
%% it will combine the following steps
%% LVWM_SALASeg_LongAxisAlignment
%% LVWM_SA_GuidPointsGeneration_LongAxisAlignment
%% LV_EndoFitting
%% LV_EpiFitting
%% LV_WholeMesh

cd(matlab_programme_folder_str);

%%%whole lv manual segmentation configuration
if ispc
    path(path, '.\segmentation');
    path(path, '.\BSplineFitting');
    path(path, '.\meshRelated');
end

if ismac
    path(path, './segmentation');
    path(path, './BSplineFitting');
    path(path, './meshRelated');
end

workingDir = pwd();

Pathname=char(matlab_programme_folder_str);
Filename=char('LVWM_config_scott.m');

projectConfig_dir = Pathname;
projectConfig_name = Filename;

cd(projectConfig_dir);
run(projectConfig_name);
cd(workingDir);



%%calling step by steps
disp('run LVWM_SALASeg_LongAxisAlignment');
LVWM_SALASeg_LongAxisAlignment(projectConfig_dir,projectConfig_name);

disp('run LVWM_SA_GuidPointsGeneration_LongAxisAlignment');
LVWM_SA_GuidPointsGeneration_LongAxisAlignment(projectConfig_dir,projectConfig_name);

disp('run LV_EndoFitting');
LV_EndoFitting(projectConfig_dir,projectConfig_name);

disp('run LV_EpiFitting');
LV_EpiFitting(projectConfig_dir,projectConfig_name);

disp('run LV_WholeMesh');
LV_WholeMesh(projectConfig_dir,projectConfig_name);

