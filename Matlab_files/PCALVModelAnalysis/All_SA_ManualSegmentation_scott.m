global matlab_programme_folder_str

cd(matlab_programme_folder_str);

if ispc
    path(path, '.\segmentation');
    path(path, '.\BSplineFitting');
    path(path, '.\meshRelated');
end

if ismac || ismac
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







%%figure out how many short axis images 
totalSASlices = 0;
for i = 1 : size(patientConfigs.SliceSpec, 1)
    if strcmp('SAcine', patientConfigs.SliceSpec(i,1).spec)
        totalSASlices = totalSASlices + 1;
    end
    
end


seg_index = 1;
while 1
    
    for i = 1 : totalSASlices
       list{i} = sprintf('%d', i); 
    end
    
    [imSelected, ~] = listdlg('ListString', list, 'SelectionMode', 'single',...
    'ListSize', [100,200], 'InitialValue',seg_index); %%current does not support change location
    %%cancer button will quit the loop
   
    if isempty(imSelected)
        break;
    end

   try 
    Function_LVWM_SASegManualUsingImpoint(imSelected,projectConfig_dir,projectConfig_name);
    seg_index = mod(seg_index, totalSASlices) + 1;
   catch 
       disp('need to re-segment');
       continue;
   end
   
end