%added on 20th Dec 2017
global MRIimageFlipB;
MRIimageFlipB = 0;
%%%%%%%%%%%%%%

%% restore default path
%restoredefaultpath;

cd('\Users\sr248u\OneDrive - University of Glasgow\Scott_working_directory\Core\Modelling\PCALVModelAnalysis');

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
if ispc 
    [FileName,PathName,~] = uigetfile('..\Results\*.m');
elseif ismac || isunix
    [FileName,PathName,~] = uigetfile('./Results/*.m');
end
% [FileName, PathName] = uigetfile( ...
%        {'*.m'}, ...
%         'Pick a file');
projectConfig_dir = PathName;
projectConfig_name = FileName(1:end-2);

cd(projectConfig_dir);
run(projectConfig_name);
cd(workingDir);


