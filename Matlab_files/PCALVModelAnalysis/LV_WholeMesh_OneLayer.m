clear all; close all; clc;
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

disp('run LV_WholeMesh for one layer mesh generation');


% LVWM_config;
if ~exist('projectConfig_name', 'var')
    LVWM_config; %% that is for standalone run
else
    cd(projectConfig_dir);
    run(projectConfig_name);
    cd(workingDir);
end

%%%% load parameter for one study
%%%% load parameter for one study
cd(resultDir);
if exist('SetParameters.m', 'file')
    SetParameters;
    cd(workingDir);
else
    disp('using the default mesh parameters');
    cd(workingDir);
    SetParametersDefault;
end
oneLayerDir = 'oneLayer';
cd(resultDir);
if ~exist(oneLayerDir, 'dir')
    mkdir(oneLayerDir);
    cd(oneLayerDir);
    oneLayerDir = pwd();
else
    cd(oneLayerDir);
    oneLayerDir = pwd();
end
cd(workingDir);

NwMesh = 1; % reset the mesh to be onelayer

%%%only for cubic spline
if(kS ~=4)
    disp('The program is ONLY for Cubic B-Spline with kS=4');
    disp('Please change kS in setparameters back to kS=4');
    stop;
end

 % c=======read geometry of the reference prolate spheriodal surface
disp('read geometry of the reference prolate spheriod');
[alpha0, w1, w2, umax, umin] = readGeometry_wholeLV(resultDir, prolateParametersFileName);

[uknot, vknot, AIJ_Endo, AIJ_Epi] =Fit_EndoEpi_Parameters(resultDir, Endo_FitparameterFileName, Epi_FitparameterFileName, ...
                               nptsu, nptsv, kS);

[uvnode, NNode,Nelement,ElementNode,neRegular, nnRegular, naReg] = PartitionUV(umax, w1, alpha0, NuMesh, NvMesh);                           
%%%%write out the mesh file
cd(oneLayerDir);
fid = fopen(FittedLVMeshFileName, 'w');
cd(workingDir);
fprintf(fid, 'TITLE = " 3D Mesh of a LV reconstructed from MR slices"\n');
fprintf(fid, 'VARIABLES = "x","y","z","u","v","w"\n');

%%%write out the reference mesh
ReferenceMesh(NwMesh, NuMesh, NvMesh, NNode,uvnode, Nelement, ElementNode, ...
                      w1, w2,alpha0,umax, fid);

%FittedLVMesh(w1, w2, alpha0, umax, NwMesh, NuMesh, NvMesh, NNode, Nelement, ElementNode,...
%                      uvnode,  nptsu, nptsv, uknot, vknot, kS, AIJ_Epi, AIJ_Endo, fid);
NeRegularApex = neRegular*NwMesh + naReg*NwMesh*NvMesh+1;
abaqusInputData = FittedLVMeshWithOutPut(w1, w2, alpha0, umax, NwMesh, NuMesh, NvMesh, NNode, Nelement,ElementNode, ...
                      uvnode,  nptsu, nptsv, uknot, vknot, kS, AIJ_Epi, AIJ_Endo,NeRegularApex, fid);
scaleTomm = 10; %%%since all data here after segmenation output it has been changed to be cm, while in Abaqus, better to be mm.
abaqusInputData.scaleTomm = 10;

% patient_name = empa024b
% projectConfig_name ='LVWM_config_empa024b_EarlyDia'
abaqus_input_name = split(projectConfig_name, '_');
abaqus_input_name = char(abaqus_input_name{end});
abaqusInputDataFileName = sprintf('%s_%s_abaqusInputData.mat', patient_name, abaqus_input_name);
cd(oneLayerDir);
save(abaqusInputDataFileName, 'abaqusInputData');
cd(workingDir);
                  
cd(oneLayerDir);
fidAba = fopen(AbaLVHexMeshFileName, 'w');
cd(workingDir);
 
abaqusInpGenerationForHexMesh(abaqusInputData, scaleTomm, fidAba); %%fidAba is closed in side the function       
              
                  
%%%plot guided spline
% PlotGuideSplines(resultDir, fid,outterGuidePointsFileName, innerGuidePointsFileName);
                  
fclose(fid);



%%need to output LV mesh in MRI coordinate system
cd(oneLayerDir);
fid = fopen('LVFittedMesh_MRICoor.dat', 'w');
cd(workingDir);

NodeMat = abaqusInputData.node;
ElemMat = abaqusInputData.elem;

nodeMat(:,2) = NodeMat(:,1)*abaqusInputData.scaleTomm;
nodeMat(:,3) = NodeMat(:,2)*abaqusInputData.scaleTomm;
nodeMat(:,4) = NodeMat(:,3)*abaqusInputData.scaleTomm;
nodeMatMRI = rotationBackToMRICoordinateSystemt(nodeMat,resultDir);

% TecplotHexMesh(nodeMatMRI, ElemMat,[],fid);
uvw = NodeMat(:,4:6);
TecplotHexMeshVec(nodeMatMRI, ElemMat,uvw,fid)
fclose(fid);