function [sendo, sepi, abaqusInputData] = load_LVGeoData_3(result_dir, workingDir)

cd(result_dir);
load abaqusInputData;
cd(workingDir);
[sendo, sepi] = extract_endo_epi_3(abaqusInputData);
cd(workingDir);
