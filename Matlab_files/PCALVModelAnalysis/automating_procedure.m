%clear all
%clc
% this is the only part of the code that needs to be changed
%filename should not be changed
%insert patient id and the baseline/follow up information
%folder = 'C:\Users\sr248u\OneDrive - University of Glasgow\Desktop\growth_tensor_result';
%baseFileName = 'early_diastolic_growth_tensor.xls';

current_patient_number={'018'};
patient_number=018;
baseline={'c'};
follow_up={'b'};
diastole_time={'earlyDiastole'};
%--------------------------------------------------------------------------

current_patient=append('empa',current_patient_number);
fullFileName = fullfile(folder, baseFileName);


current_baseline=append(current_patient,baseline);
current_follow_up=append(current_patient,follow_up);

current_baseline_number=append(current_patient_number,baseline);
current_follow_up_number=append(current_patient_number,follow_up);

main_file='\Users\sr248u\OneDrive - University of Glasgow\Scott_working_directory\Core\Modelling\Results';
fibre_path='fibreGeneration\Results_fiber_60_45\fiberDir.txt';
sheet_path='fibreGeneration\Results_fiber_60_45\sheetDir.txt';
radDir_path='fibreGeneration\Results_fiber_60_45\radDir.txt';
cirDir_path='fibreGeneration\Results_fiber_60_45\cirDir.txt';

main_file_patient=fullfile(main_file,current_patient);



projectConfig_dir_first_cell=fullfile(main_file_patient,current_baseline_number);
projectConfig_dir_first=string(projectConfig_dir_first_cell);

ResultsDirAllScans_first_cell=fullfile(projectConfig_dir_first_cell,diastole_time);
ResultsDirAllScans_first=string(ResultsDirAllScans_first_cell);

fiberDir_first_cell=fullfile(ResultsDirAllScans_first_cell,fibre_path);
fiberDir_first_loader=string(fiberDir_first_cell);

sheetDir_first_cell=fullfile(ResultsDirAllScans_first_cell,sheet_path);
sheetDir_first_loader=string(sheetDir_first_cell);

radDir_first_cell=fullfile(ResultsDirAllScans_first_cell,radDir_path);
radDir_first_loader=string(radDir_first_cell);

cirDir_first_cell=fullfile(ResultsDirAllScans_first_cell,cirDir_path);
cirDir_first_loader=string(cirDir_first_cell);



%fiberDir_first = load(fiberDir_first_loader);
%sheetDir_first = load(sheetDir_first_loader);
%radDir_first = load(radDir_first_loader);
%cirDir_first = load(cirDir_first_loader);

%%=============================================================================================================
%%=============================================================================================================

projectConfig_dir_second_cell=fullfile(main_file_patient,current_follow_up_number);
projectConfig_dir_second=string(projectConfig_dir_second_cell);

ResultsDirAllScans_second_cell=fullfile(projectConfig_dir_second_cell,diastole_time);
ResultsDirAllScans_second=string(ResultsDirAllScans_second_cell);

fiberDir_second_cell=fullfile(ResultsDirAllScans_second_cell,fibre_path);
fiberDir_second_loader=string(fiberDir_second_cell);

sheetDir_second_cell=fullfile(ResultsDirAllScans_second_cell,sheet_path);
sheetDir_second_loader=string(sheetDir_second_cell);

radDir_second_cell=fullfile(ResultsDirAllScans_second_cell,radDir_path);
radDir_second_loader=string(radDir_second_cell);

cirDir_second_cell=fullfile(ResultsDirAllScans_second_cell,cirDir_path);
cirDir_second_loader=string(cirDir_second_cell);



%fiberDir_second = load(fiberDir_second_loader);
%sheetDir_second = load(sheetDir_second_loader);
%radDir_second = load(radDir_second_loader);
%cirDir_second = load(cirDir_second_loader);

%%=============================================================================================================
%%=============================================================================================================

%LVWM_GeoChangesWLmappedmesh





%%=============================================================================================================
%%=============================================================================================================

%all_vals=[patient_number meanFc mean_12 mean_13 mean_21 meanFr mean_23 mean_31 mean_32 meanFl i1 i2 i3];

%writematrix(all_vals,fullFileName,'WriteMode','append')




