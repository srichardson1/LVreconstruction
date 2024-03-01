
patient_name = first_directory_name;

dicomDir = full_dicom_directory_computed;
resultDir = folder_results_final;

workingDir= '\Users\sr248u\OneDrive - University of Glasgow\Scott_working_directory\Core\Modelling\PCALVModelAnalysis';



abaqusDir_simulation = 'abaqusSimulation';
abaqusDir_reconstruction = resultDir;


cd(resultDir);
resultDir = pwd();
cd(workingDir);


if ~exist(resultDir, 'dir')
    mkdir(resultDir);
end

cd(resultDir);
if ~exist('fibreGeneration', 'dir')
    mkdir('fibreGeneration');
    cd(workingDir);
end

%%let's figure out whether we need to move basal plane or not
BasalMovingB = 0; %%if exist then it will try to move the first basal plane, otherwise it will not do anything

sampleN = 50;

%%%patient image data
patientIndex = 1;seriesIndex = 0;
patientConfigs(patientIndex,1).name = patient_name;
patientConfigs(patientIndex,1).TimeEndOfSystole = end_systole_time_for_config_file;
patientConfigs(patientIndex,1).TimeEndOfDiastole = end_diastole_time_for_config_file;
patientConfigs(patientIndex,1).TimeEarlyOfDiastole =early_diastole_time_for_config_file;
patientConfigs(patientIndex,1).sampleN = sampleN;
patientConfigs(patientIndex,1).SASliceDistance = 10; % mm this is not used, no need to change

patientConfigs(patientIndex,1).SASlicePositionApex = 7; %%%the starting slice position of apical region, which means slice 6 and 7 will considered to be at apical region
patientConfigs(patientIndex,1).totalLVOTSliceLocation = 3;
patientConfigs(patientIndex,1).totalTimeInstance = 25;
patientConfigs(patientIndex,1).timeInstanceSelected = patientConfigs(patientIndex,1).TimeEarlyOfDiastole;

% sampleN = patientConfigs(patientIndex,1).sampleN;

%this is for writing out the guide points when fitiing
sliceToBeSkipped = [];
sliceToBeSkippedLA = [];
basalSlices = [1 2 ];
middlSlices = [3 4 5 ];
apicaSlices = [6 7];

%%the following is for MI intensity window segmentation 
WW = 190; WC=120; %this is checked with Syngo viewer
intensityWindow = [ max(WC-WW/2, 0), WC+WW/2]; 
kcluster_centers = 3;
MIApex = 0;


%%not sure what to do with this, but leave it now
% scanIndex = 1;
basalMoving = 0;
basalMovingDistance = 0.0;
% if basalMoving == 1
%     cd(resultDir);
%     cd('..'); %%assume the adjustment for different scans is saved in the parent directory
%     load ResultsDirAllScans;
%     basalMovingDistance = ResultsDirAllScans(scanIndex).moveTowardsApex;
% end

%%SA basal 1: %SP A44.3
%%this is for end-diastole, in other time, such as in early-diastole and
%%mid-diastole, this slice does not cover the LV anymore, so needs to take
%%care when integrating other pathological information.

current_slice=first_sa_slice_for_config_file;

for j=1:total_SA_slices

first_sa_slice_for_config_file_str = string(current_slice);
first_sa_slice_composed=append('SA',first_sa_slice_for_config_file_str);
first_sa_slice_studyName=append('SA-plane-',first_sa_slice_for_config_file_str);


seriesIndex = seriesIndex + 1;
patientConfigs(patientIndex,1).dir(seriesIndex,1).studyDir = dicomDir;
patientConfigs(patientIndex,1).dirMidSA(seriesIndex,1).ImgDir = first_sa_slice_composed; %%%most base at end-diastole
patientConfigs(patientIndex,1).studyName(seriesIndex,1).studyName = first_sa_slice_studyName;
patientConfigs(patientIndex,1).SliceSpec(seriesIndex,1).spec = 'SAcine';

current_slice=current_slice+1;

end

 
%%%LVot 1
seriesIndex = seriesIndex + 1;
patientConfigs(patientIndex,1).dir(seriesIndex,1).studyDir = dicomDir;
patientConfigs(patientIndex,1).dirMidSA(seriesIndex,1).ImgDir = 'LVOT';
patientConfigs(patientIndex,1).studyName(seriesIndex,1).studyName = 'LAcine-LVOT-1';
patientConfigs(patientIndex,1).SliceSpec(seriesIndex,1).spec = 'LAcine_LVOT';

%%4CH
seriesIndex = seriesIndex + 1;
patientConfigs(patientIndex,1).dir(seriesIndex,1).studyDir = dicomDir;
patientConfigs(patientIndex,1).dirMidSA(seriesIndex,1).ImgDir = 'HLA';
patientConfigs(patientIndex,1).studyName(seriesIndex,1).studyName = 'LAcine-4CH-1';
patientConfigs(patientIndex,1).SliceSpec(seriesIndex,1).spec = 'LAcine_4CH';

%%1CH
seriesIndex = seriesIndex + 1;
patientConfigs(patientIndex,1).dir(seriesIndex,1).studyDir = dicomDir;
patientConfigs(patientIndex,1).dirMidSA(seriesIndex,1).ImgDir = 'VLA';
patientConfigs(patientIndex,1).studyName(seriesIndex,1).studyName = 'LAcine-1CH-1';
patientConfigs(patientIndex,1).SliceSpec(seriesIndex,1).spec = 'LAcine_1CH';

% SA distance will be constant. do not assume it, it will directly
% calculated from the dicom header
%LGE basal 1
%seriesIndex = seriesIndex + 1;
%patientConfigs(patientIndex,1).dir(seriesIndex,1).studyDir = dicomDir;
%patientConfigs(patientIndex,1).dirMidSA(seriesIndex,1).ImgDir = 'MI_SA';
%patientConfigs(patientIndex,1).studyName(seriesIndex,1).studyName = 'LGSeries';
%patientConfigs(patientIndex,1).SliceSpec(seriesIndex,1).spec = 'LGE_SA';
% 
% 
% %%the following is for MI intensity window segmentation 
% intensityWindow = [2342-551/2 2342+551/2]; 
% kcluster_centers = 3;
% MIApex = 0;
