function [sendo, sepi, abaqusInputData] = load_LVGeoData_1(result_dir, workingDir)

cd(result_dir);
load abaqusInputData;
%cd(workingDir);
[sendo, sepi] = extract_endo_epi_1(abaqusInputData);
%cd(workingDir);
