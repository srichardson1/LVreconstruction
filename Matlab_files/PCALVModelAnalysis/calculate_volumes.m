global yes_or_no_checkbox period_time_str working_letter where_to_save_excel
% clear all; close all; clc;

% workingDir = pwd();
% [file,file_path] = uigetfile('.mat'); 


%cd(resultDir);
%load(C:\Users\staff\OneDrive - University of Glasgow\HV01\endDiastole\abaqusInputData);
%cd(workingDir);

endoNodes = abaqusInputData.endoNodes; 
epiNodes = abaqusInputData.epiNodes;

node = abaqusInputData.node(:,1:3);

nodes_x = node(endoNodes,1);
nodes_y = node(endoNodes,2);
nodes_z = node(endoNodes,3);

DT =  DelaunayTri(nodes_x,nodes_y,nodes_z);
[~, vol_endo] = convexHull(DT)

nodes_x = node(epiNodes,1);
nodes_y = node(epiNodes,2);
nodes_z = node(epiNodes,3);
DT =  DelaunayTri(nodes_x,nodes_y,nodes_z);
[~, vol_epi] = convexHull(DT)

vol_LV_wall = vol_epi - vol_endo



if yes_or_no_checkbox == true

folder_and_file=fullfile(where_to_save_excel,'volume_data.xls');
    

if ~exist(folder_and_file,'file')
initial_vals_to_save={'Patient ID' 'Early or End Diastole' 'Wall Volume' 'Epi Volume' 'Endo Volume'};
writecell(initial_vals_to_save,folder_and_file,'WriteMode','append')
end


all_vals_to_save={working_letter period_time_str vol_LV_wall vol_epi vol_endo};
writecell(all_vals_to_save,folder_and_file,'WriteMode','append')

end







