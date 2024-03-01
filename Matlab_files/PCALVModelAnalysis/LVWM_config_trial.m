global total_number_of_SA_images matlab_programme_folder_str
%added on 20th Dec 2017
global MRIimageFlipB;
MRIimageFlipB = 0;
%%%%%%%%%%%%%%

%% restore default path
%restoredefaultpath;

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



%%load the patient config file


%if ispc 
%    [FileName,PathName,~] = uigetfile('..\Results\*.m');
%elseif ismac || isunix
%    [FileName,PathName,~] = uigetfile('./Results/*.m');
%end

% [FileName, PathName] = uigetfile( ...
%        {'*.m'}, ...
%         'Pick a file');

Pathname=char(matlab_programme_folder_str);
Filename=char('LVWM_config_scott.m');

projectConfig_dir = Pathname;
projectConfig_name = Filename;

cd(projectConfig_dir);
run(projectConfig_name);
cd(workingDir);


