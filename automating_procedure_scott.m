global folder_results_final_baseline folder_results_final_follow period_time_str current_patient where_to_save_excel yes_or_no_checkbox_growth


%clear all
%clc
% this is the only part of the code that needs to be changed
%filename should not be changed
%insert patient id and the baseline/follow up information
%folder = 'C:\Users\sr248u\OneDrive - University of Glasgow\Desktop\growth_tensor_result';
%baseFileName = 'early_diastolic_growth_tensor.xls';

%--------------------------------------------------------------------------


fibre_path='fibreGeneration\Results_fiber_60_45\fiberDir.txt';
sheet_path='fibreGeneration\Results_fiber_60_45\sheetDir.txt';
radDir_path='fibreGeneration\Results_fiber_60_45\radDir.txt';
cirDir_path='fibreGeneration\Results_fiber_60_45\cirDir.txt';



fiberDir_first_loader=fullfile(folder_results_final_baseline,fibre_path);
sheetDir_first_loader=fullfile(folder_results_final_baseline,sheet_path);
radDir_first_loader=fullfile(folder_results_final_baseline,radDir_path);
cirDir_first_loader=fullfile(folder_results_final_baseline,cirDir_path);



fiberDir_first = load(fiberDir_first_loader);
sheetDir_first = load(sheetDir_first_loader);
radDir_first = load(radDir_first_loader);
cirDir_first = load(cirDir_first_loader);

%%=============================================================================================================
%%=============================================================================================================

fiberDir_second_loader=fullfile(folder_results_final_follow,fibre_path);
sheetDir_second_loader=fullfile(folder_results_final_follow,sheet_path);
radDir_second_loader=fullfile(folder_results_final_follow,radDir_path);
cirDir_second_loader=fullfile(folder_results_final_follow,cirDir_path);




fiberDir_second = load(fiberDir_second_loader);
sheetDir_second = load(sheetDir_second_loader);
radDir_second = load(radDir_second_loader);
cirDir_second = load(cirDir_second_loader);

%%=============================================================================================================
%%=============================================================================================================

LVWM_GeoChangesWLmappedmesh_scott





%%=============================================================================================================
%%=============================================================================================================

%all_vals=[patient_number meanFc mean_12 mean_13 mean_21 meanFr mean_23 mean_31 mean_32 meanFl i1 i2 i3];

%writematrix(all_vals,fullFileName,'WriteMode','append')

if yes_or_no_checkbox_growth == true

folder_and_file_growth=fullfile(where_to_save_excel,'growth_tensor.xls');
    

if ~exist(folder_and_file_growth,'file')
initial_vals_to_save={'Patient Number' 'Early or End Diastole' '1' '2' '3' '4' '5' '6' '7' '8' '9' 'invariant 1' 'invariant 2' 'invariant 3'};
writecell(initial_vals_to_save,folder_and_file_growth,'WriteMode','append')
end


all_vals_to_save={current_patient period_time_str meanFc mean_12 mean_13 mean_21 meanFr mean_23 mean_31 mean_32 meanFl i1 i2 i3};
writecell(all_vals_to_save,folder_and_file_growth,'WriteMode','append')

end



