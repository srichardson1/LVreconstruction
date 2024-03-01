global matlab_programme_folder_str

cd(matlab_programme_folder_str);

%%%whole lv manual segmentation configuration
if ispc
    path(path, '.\segmentation');
    path(path, '.\BSplineFitting');
    path(path, '.\meshRelated');
end

if ismac || isunix
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




seg_index = 1;
while 1
    list = {'1', '2', '3'};
    [imSelected, ~] = listdlg('ListString', list, 'SelectionMode', 'single',...
    'ListSize', [80,50], 'InitialValue',seg_index); %%current does not support change location
    %%cancer button will quit the loop
   
    if isempty(imSelected)
        break;
    end

   try 
    LVWM_LASegManual(imSelected,projectConfig_dir,projectConfig_name);
    seg_index = mod(seg_index, 3) + 1;
   catch 
       disp('need to re-segment');
       continue;
   end
   
end