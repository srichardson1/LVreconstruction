%%%
clear all; 
close all; 
clc;

path(path, '.\segmentation');
path(path, '.\BSplineFitting');
path(path, '.\meshRelated');
workingDir = pwd();

%%load the patient config file
[FileName, PathName] = uigetfile( ...
       {'*.m'}, ...
        'Pick a file');%[FileName,PathName,~] = uigetfile('*.m');
projectConfig_dir = PathName;
projectConfig_name = FileName(1:end-2);

cd(projectConfig_dir);
run(projectConfig_name);
cd(workingDir)

%baselineScanIndex = 1;
%followUpScanIndex = 3;

baselineScanIndex = 3;
followUpScanIndex = 4;

%%then will only summarize the different values accoring to LGE intensity
% mean_remote_intensity = 0.642592;
% mean_MI_intensity = 0.84;
% mean_MVO = 0.95;

%Group dilation
%MI087
%mean_remote_intensity = 0.7526; %scan1
%mean_remote_intensity = 0.7601; %scan3

%MI102
%mean_remote_intensity = 0.7526; %scan1
%mean_remote_intensity = 0.6688; %scan3

% MI106
%mean_remote_intensity = 0.6671;%scan1
%mean_remote_intensity = 0.7079;%scan3

% MI108
%mean_remote_intensity = 0.6671;%scan1
%mean_remote_intensity = 0.7079;%scan3

% MI133
%mean_remote_intensity = 0.745422;%scan1
%mean_remote_intensity = 0.732116;%scan3

% MI167
%mean_remote_intensity = 0.382573;%scan1
%mean_remote_intensity = 0.6660;%scan3
%mean_remote_intensity = 0.6473;%scan1
%mean_remote_intensity = 0.6660;%scan3

%MI171
%mean_remote_intensity = 0.719467; %scan1
%mean_remote_intensity = 0.72915; %scan3

% MI190
%mean_remote_intensity = 0.757738;%scan1
%mean_remote_intensity = 0.686386;%scan3

%MI201
%mean_remote_intensity = 0.682421; %scan1
%mean_remote_intensity = 0.706659; %scan3

% No change group

% MI172
%mean_remote_intensity = 0.7431;%scan1
%mean_remote_intensity = 0.7466;%scan3

%MI179
%mean_remote_intensity = 0.7606;%scan1
%mean_remote_intensity = 0.7057;%scan3

% Shringake group

% MI132
%mean_remote_intensity = 0.724902;%scan1
%mean_remote_intensity = 0.771417;%scan3

%MI258
%mean_remote_intensity = 0.7101; %scan1
%mean_remote_intensity = 0.7029; %scan3

%MI265
%mean_remote_intensity = 0.74866;%scan1
%mean_remote_intensity = 0.732020;%scan3

% MI280
%mean_remote_intensity = 0.769496;%scan1
%mean_remote_intensity = 0.783905;%scan3

% MI288
%mean_remote_intensity = 0.8133;%scan1
mean_remote_intensity = 0.8190;%scan3

cd(projectConfig_dir);
fileNameToLoad = sprintf('DG_%s_to_%s.mat', ResultsDirAllScans(1,followUpScanIndex).folderName, ...
            ResultsDirAllScans(1,baselineScanIndex).folderName);
save(fileNameToLoad);
cd(workingDir)

%%now need to figure out the distance from the basal plane to the valvular
%%ring

II = [1 0 0;... 
      0 1 0; ...
      0 0 1];


for scanIndexRef = baselineScanIndex
    cd(projectConfig_dir);
    cd(ResultsDirAllScans(1,scanIndexRef).folderName);
    load abaqusInputData;
%
%   LV mesh with MI intensity
%
    load LGEMappedLVMeshforstrain;
    cd(workingDir);
    abaqusInputDataRef = abaqusInputData;
    clear abaqusInputData;
    
    cd(projectConfig_dir);
    cd(ResultsDirAllScans(1,scanIndexRef).folderName);
    cd(ResultsDirAllScans(1,scanIndexRef).folderNameFibre);
    fiberDir = load('fiberDir.txt');
    sheetDir = load('sheetDir.txt');
    radDir = load('radDir.txt');
    cirDir = load('cirDir.txt');
    cd(workingDir);

    for scanIndex = followUpScanIndex %scanIndexRef+1 :  totalScanNo %%only calculate the first and the last

        cd(projectConfig_dir);
        tecFileName = sprintf('tecplotDeform_%s_to_%s_mappedmesh.dat',ResultsDirAllScans(1,scanIndex).folderName, ...
            ResultsDirAllScans(1,scanIndexRef).folderName );
        tecfid = fopen(tecFileName,'w');
        
        tecFileName_strain = sprintf('tecplotStrain_%s_to_%s_mappedmesh.dat',ResultsDirAllScans(1,scanIndex).folderName, ...
            ResultsDirAllScans(1,scanIndexRef).folderName );
        tecStrainfid = fopen(tecFileName_strain,'w');
        
        cd(ResultsDirAllScans(1,scanIndex).folderName);
        load abaqusInputData;
        cd(workingDir);

        FF_r = []; EE_r = [];
        FF_c = []; EE_c = [];
        FF_l = []; EE_l = [];
        
        FF_f = []; EE_f = [];
        FF_s = []; EE_s = [];
        FF_n = []; EE_n = [];

        for i = 1 : size(abaqusInputData.elem,1)
            plist = abaqusInputData.elem(i,:);
            nodeRef_x =  abaqusInputDataRef.node(plist,1);
            nodeRef_y =  abaqusInputDataRef.node(plist,2);
            nodeRef_z =  abaqusInputDataRef.node(plist,3);

            nodeCur_x =  abaqusInputData.node(plist,1);
            nodeCur_y =  abaqusInputData.node(plist,2);
            nodeCur_z =  abaqusInputData.node(plist,3);
            FF = finiteStrainHex(nodeCur_x, nodeCur_y, nodeCur_z, ...
                               nodeRef_x, nodeRef_y, nodeRef_z);
 %
 %         get deformation gradient tensor node by node
 %
            F(1,1,i)=FF(1,1);
            F(1,2,i)=FF(1,2);
            F(1,3,i)=FF(1,3);
            
            F(2,1,i)=FF(2,1);
            F(2,2,i)=FF(2,2);
            F(2,3,i)=FF(2,3);
            
            F(3,1,i)=FF(3,1);
            F(3,2,i)=FF(3,2);
            F(3,3,i)=FF(3,3);
 %
 %          calculate det of deformation tensor 
 %
            J(i)=det(FF);
            
            EE = 1/2*(FF'*FF-II);

%             %%a simple way to figure out local direction 
%             cirDir = [ abaqusInputDataRef.node(plist(3),1) - abaqusInputDataRef.node(plist(2),1), ...
%                        abaqusInputDataRef.node(plist(3),2) - abaqusInputDataRef.node(plist(2),2), ...
%                        abaqusInputDataRef.node(plist(3),3) - abaqusInputDataRef.node(plist(2),3)];
%             cirDir=NormalizationVec(cirDir);
% 
%             logiDir = [ abaqusInputDataRef.node(plist(1),1) - abaqusInputDataRef.node(plist(2),1), ...
%                         abaqusInputDataRef.node(plist(1),2) - abaqusInputDataRef.node(plist(2),2), ...
%                         abaqusInputDataRef.node(plist(1),3) - abaqusInputDataRef.node(plist(2),3)];
%             logiDir=NormalizationVec(logiDir);
% 
%             raDir = cross(cirDir,logiDir);
%             logiDir = cross(raDir, cirDir);

            rdir = radDir(i,2:4);
            cdir = cirDir(i,2:4);
            ldir = cross(cdir, rdir);
            Rot = [cdir;rdir;ldir]; 
            FF_RCL = Rot*FF*(Rot');   EE_RCL = Rot*EE*(Rot');
            FF_c(i,1) = FF_RCL(1,1);  EE_c(i,1) = EE_RCL(1,1);
            FF_r(i,1) = FF_RCL(2,2);  EE_r(i,1) = EE_RCL(2,2);
            FF_l(i,1) = FF_RCL(3,3);  EE_l(i,1) = EE_RCL(3,3);
            
            F_cc(i,1)=FF_RCL(1,1);
            F_rc(i,1)=FF_RCL(2,1);
            F_lc(i,1)=FF_RCL(3,1);
            
            F_cr(i,1)=FF_RCL(1,2);
            F_rr(i,1)=FF_RCL(2,2);
            F_lr(i,1)=FF_RCL(3,2);            
            
            F_cl(i,1)=FF_RCL(1,3);
            F_rl(i,1)=FF_RCL(2,3);
            F_ll(i,1)=FF_RCL(3,3);
 %
 %          Green-Larange strain
 %
            E_cc(i,1)=EE_RCL(1,1);
            E_rr(i,1)=EE_RCL(2,2);
            E_ll(i,1)=EE_RCL(3,3);
            
            E_cr(i,1)=EE_RCL(1,2);
            E_cl(i,1)=EE_RCL(1,3);
            E_rl(i,1)=EE_RCL(2,3);
       
            Rot = [];
            fdir = fiberDir(i, 2:4);
            sdir = sheetDir(i, 2:4);
            ndir = cross(fdir,sdir);
            Rot = [fdir; sdir;ndir];
            FF_fsn = Rot*FF*(Rot');  EE_fsn = Rot*EE*(Rot');
            FF_f(i,1) = FF_fsn(1,1); EE_f(i,1) = EE_fsn(1,1);
            FF_s(i,1) = FF_fsn(2,2); EE_s(i,1) = EE_fsn(2,2);
            FF_n(i,1) = FF_fsn(3,3); EE_n(i,1) = EE_fsn(3,3);
            
            Lambda_f(i,1)=(2*EE_f(i,1)+1)^0.5; 
        
        end

        TecplotLVMeshWriteWithCentreData(tecfid,abaqusInputData,FF_f,FF_s,FF_n,FF_c, FF_r, FF_l);
        TecplotLVMeshWriteWithCentreDataStrain(tecStrainfid,abaqusInputData,EE_f,EE_s,EE_n,EE_c, EE_r, EE_l);
        

        fclose(tecfid);
        fclose(tecStrainfid);


       cd(projectConfig_dir);
       cd(ResultsDirAllScans(1,scanIndexRef).folderName);
       resultDir = pwd();
       cd(workingDir);
    
  %      NodeMat = abaqusInputData.node;
  %      ElemMat = abaqusInputData.elem;

  %      nodeMat(:,2) = NodeMat(:,1)*abaqusInputData.scaleTomm;
  %      nodeMat(:,3) = NodeMat(:,2)*abaqusInputData.scaleTomm;
  %      nodeMat(:,4) = NodeMat(:,3)*abaqusInputData.scaleTomm;
  %      nodeMatMRI = rotationBackToMRICoordinateSystemt(nodeMat,resultDir);
  %      nodeMat(:,1) = nodeMatMRI(:,2);
  %      nodeMat(:,2) = nodeMatMRI(:,3);
  %      nodeMat(:,3) = nodeMatMRI(:,4);
  %      abaqusInputDataInMRICoor = abaqusInputData;
  %      abaqusInputDataInMRICoor.node = nodeMat;
        
  %      cd(projectConfig_dir);
  %      tecFileName = sprintf('tecplotDeform_%s_to_%s_MRICoor_mappedmesh.dat', ResultsDirAllScans(1,scanIndex).folderName, ...
  %          ResultsDirAllScans(1,scanIndexRef).folderName);
  %      tecfid = fopen(tecFileName,'w');
  %      tecStrainFileName = sprintf('tecplotStrain_%s_to_%s_MRICoor_mappedmesh.dat', ResultsDirAllScans(1,scanIndex).folderName, ...
  %          ResultsDirAllScans(1,scanIndexRef).folderName);
  %      tecStrainfid = fopen(tecStrainFileName,'w');
  %      cd(workingDir);
  %      TecplotLVMeshWriteWithCentreData(tecStrainfid,abaqusInputDataInMRICoor,EE_f,EE_s,EE_n,EE_c, EE_r, EE_l);
  %  fclose(tecStrainfid);
  %
  %     Write Tecplot files in the reference configuration 
  %
        NodeMat = abaqusInputDataRef.node;
        ElemMat = abaqusInputDataRef.elem;

        nodeMat(:,2) = NodeMat(:,1)*abaqusInputData.scaleTomm;
        nodeMat(:,3) = NodeMat(:,2)*abaqusInputData.scaleTomm;
        nodeMat(:,4) = NodeMat(:,3)*abaqusInputData.scaleTomm;
        nodeMatMRI = rotationBackToMRICoordinateSystemt(nodeMat,resultDir,workingDir);
        nodeMat(:,1) = nodeMatMRI(:,2);
        nodeMat(:,2) = nodeMatMRI(:,3);
        nodeMat(:,3) = nodeMatMRI(:,4);
        
        abaqusInputDataInMRICoor = abaqusInputDataRef;
        abaqusInputDataInMRICoor.node = nodeMat;
        
        NodeMat_cur = abaqusInputData.node;
        nodeMat_cur(:,2) = NodeMat_cur(:,1)*abaqusInputData.scaleTomm;
        nodeMat_cur(:,3) = NodeMat_cur(:,2)*abaqusInputData.scaleTomm;
        nodeMat_cur(:,4) = NodeMat_cur(:,3)*abaqusInputData.scaleTomm;
        nodeMatMRI_cur = rotationBackToMRICoordinateSystemt(nodeMat_cur,resultDir,workingDir);
        nodeMat_cur(:,1) = nodeMatMRI_cur(:,2);
        nodeMat_cur(:,2) = nodeMatMRI_cur(:,3);
        nodeMat_cur(:,3) = nodeMatMRI_cur(:,4);
        abaqusInputDataInMRICoor_cur = abaqusInputData;
        abaqusInputDataInMRICoor_cur.node = nodeMat_cur;
        
        cd(projectConfig_dir);
        
        tecFileName = sprintf('tecplotDeform_%s_to_%s_MRICoor_mappedmesh.dat', ResultsDirAllScans(1,scanIndex).folderName, ...
            ResultsDirAllScans(1,scanIndexRef).folderName);
        tecfid = fopen(tecFileName,'w');
        cd(workingDir);
        TecplotLVMeshWriteWithCentreData(tecfid,abaqusInputDataInMRICoor,FF_f,FF_s,FF_n,FF_c, FF_r, FF_l);
        fclose(tecfid);

        
        cd(projectConfig_dir);
        
        tecFileName = sprintf('tecplotDeform_%s_to_%s_MRICoor_mappedmesh_cur.dat', ResultsDirAllScans(1,scanIndex).folderName, ...
            ResultsDirAllScans(1,scanIndexRef).folderName);
        tecfid = fopen(tecFileName,'w');
        cd(workingDir);
        TecplotLVMeshWriteWithCentreData(tecfid,abaqusInputDataInMRICoor_cur,FF_f,FF_s,FF_n,FF_c, FF_r, FF_l);
        fclose(tecfid);
        
        cd(projectConfig_dir);
        tecStrainFileName = sprintf('tecplotStrain_%s_to_%s_MRICoor_mappedmesh.dat', ResultsDirAllScans(1,scanIndex).folderName, ...
            ResultsDirAllScans(1,scanIndexRef).folderName);
        tecStrainfid = fopen(tecStrainFileName,'w');
        cd(workingDir);
        TecplotLVMeshWriteWithCentreData(tecStrainfid,abaqusInputDataInMRICoor,EE_f,EE_s,EE_n,EE_c, EE_r, EE_l);
        fclose(tecStrainfid);  
        
       cd(projectConfig_dir);
       tecuvwFileName = sprintf('tecplotDeform_%s_to_%s_uvwCoor_mappedmesh.dat', ResultsDirAllScans(1,scanIndex).folderName, ...
            ResultsDirAllScans(1,scanIndexRef).folderName);
       tecuvwfid = fopen(tecuvwFileName,'w');
       cd(workingDir);
       TecplotLVMeshWriteWithCentreDatauvw(tecuvwfid,abaqusInputDataRef,FF_f,FF_s,FF_n,FF_c, FF_r, FF_l);
       fclose(tecuvwfid);              
       LGEMappedelement=LGEMappedLVMeshforstrain.elemMat;
       LGEintensity=LGEMappedLVMeshforstrain.intensity;
          
       cd(resultDir);
       fid = fopen('LVMesh_LGE_Ref_configuration.dat','w');
       cd(workingDir);
       TecplotHexMeshRegions(nodeMatMRI,LGEMappedelement,LGEintensity,fid);
       fclose(fid);
       

       cd(resultDir);
       fid = fopen('LVMesh_LGE_Cur_configuration.dat','w');
       cd(workingDir);
       TecplotHexMeshRegions(nodeMatMRI_cur,LGEMappedelement,LGEintensity,fid);
       fclose(fid); 
       
       nodeMatuvw(:,1)=abaqusInputDataRef.node(:,4);
       nodeMatuvw(:,2)=abaqusInputDataRef.node(:,5);
       nodeMatuvw(:,3)=abaqusInputDataRef.node(:,6);
       
       cd(resultDir);
       fid = fopen('LVMesh_LGE_uvw_Ref_configuration.dat','w');
       cd(workingDir);
       TecplotHexMeshRegionsuvw(nodeMatuvw,LGEMappedelement,LGEintensity,fid);
       fclose(fid);       
       
%
%      calculate mean strain/deformation gradient and standard deviations in MI and remote zones
% 

cd(projectConfig_dir);
fileNameToLoad = sprintf('DG_%s_to_%s.mat', ResultsDirAllScans(1,followUpScanIndex).folderName, ...
            ResultsDirAllScans(1,baselineScanIndex).folderName);
load(fileNameToLoad);
cd(workingDir)

%%need to entre into the baseline scan to access the LGE mapping
cd(projectConfig_dir);
cd(ResultsDirAllScans(1,baselineScanIndex).folderName);
load LGEMappedLVMesh; %%LGEEachLVNode, LGEEachLVElem
load abaqusInputData;
cd(workingDir);


%%here we will only need to count the results out of apical top
cd(projectConfig_dir);
cd(ResultsDirAllScans(1,baselineScanIndex).folderName);
load XYZEndoFitted;
SetParameters;
cd(workingDir);
%
% element ID>24600 for apex region
%
%element_above_apex = neReg*NwMesh;
    element_above_apex =24600;
    LGEEachLVElem=LGEMappedLVMesh.intensity;

for el = 1 : size(abaqusInputData.elem,1)
    elem_nodes = abaqusInputData.elem(el,:);
    elem_xyz = abaqusInputData.node(elem_nodes,1:3);
    Nnode=8;
    for j=1:Nnode
        xelem(j)=elem_xyz(j,1);
        yelem(j)=elem_xyz(j,2);
        zelem(j)=elem_xyz(j,3);
    end
    P = [xelem',yelem',zelem'];
    %P(:)=[elem_xyz(:,1)',elem_xyz(:,2)',elem_xyz(:,3)'];
    dt =  DelaunayTri(P);
    [ch, elem_vol] = convexHull(dt);
    LVWallVol(el,1) = elem_vol;
end
kkk=0;
if kkk==1
h_FF_crl = figure();hold on;
subplot(231);
plot(LGEEachLVElem(1:element_above_apex),FF_c(1:element_above_apex),'.');
xlabel('LGE');ylabel('F_c_c');
subplot(232);
plot(LGEEachLVElem(1:element_above_apex),FF_r(1:element_above_apex),'.');
xlabel('LGE');ylabel('F_r_r');
subplot(233);
plot(LGEEachLVElem(1:element_above_apex),FF_r(1:element_above_apex),'.');
xlabel('LGE');ylabel('F_l_l');

h_FF_fsn = figure();hold on;
subplot(231);
plot(LGEEachLVElem(1:element_above_apex),FF_c(1:element_above_apex),'.');
xlabel('LGE');ylabel('F_f_f');
subplot(232);
plot(LGEEachLVElem(1:element_above_apex),FF_r(1:element_above_apex),'.');
xlabel('LGE');ylabel('F_s_s');
subplot(233);
plot(LGEEachLVElem(1:element_above_apex),FF_r(1:element_above_apex),'.');
xlabel('LGE');ylabel('F_n_n');

h_F_crl1 = figure();hold on;
subplot(231);
plot(LGEEachLVElem(1:element_above_apex),F_cr(1:element_above_apex),'.');
xlabel('LGE');ylabel('F_c_r');

subplot(232);
plot(LGEEachLVElem(1:element_above_apex),F_cl(1:element_above_apex),'.');
xlabel('LGE');ylabel('F_r_c');

subplot(233);
plot(LGEEachLVElem(1:element_above_apex),F_lc(1:element_above_apex),'.');
xlabel('LGE');ylabel('F_c_l');

h_F_crl2 = figure();hold on;
subplot(231);
plot(LGEEachLVElem(1:element_above_apex),F_lr(1:element_above_apex),'.');
xlabel('LGE');ylabel('F_l_c');

subplot(232);
plot(LGEEachLVElem(1:element_above_apex),F_rc(1:element_above_apex),'.');
xlabel('LGE');ylabel('F_r_l');

subplot(233);
plot(LGEEachLVElem(1:element_above_apex),F_rl(1:element_above_apex),'.');
xlabel('LGE');ylabel('F_l_r');

h_J_Lambda_f = figure();hold on;
subplot(2,2,1);
plot(LGEEachLVElem(1:element_above_apex),J(1:element_above_apex),'.');
xlabel('LGE');ylabel('J');
subplot(2,2,2);
plot(LGEEachLVElem(1:element_above_apex),Lambda_f(1:element_above_apex),'.');
xlabel('LGE');ylabel('\lambda_f');
else
end


mean_MVO = 1;

FF_c_remote=[];FF_c_MI=[];FF_c_MVO = [];
FF_r_remote=[];FF_r_MI=[];FF_r_MVO=[];
FF_l_remote=[];FF_l_MI=[];FF_l_MVO=[];

FF_f_remote = [];FF_f_MI = [];FF_f_MVO = [];
FF_s_remote = [];FF_s_MI = [];FF_s_MVO = []; 
FF_n_remote = [];FF_n_MI = [];FF_n_MVO = [];

F_cr_remote=[];F_cr_MI=[];F_cr_MVO = [];
F_rc_remote=[];F_rc_MI=[];F_rc_MVO = [];
F_cl_remote=[];F_cl_MI=[];F_cl_MVO = [];
F_lc_remote=[];F_lc_MI=[];F_lc_MVO = [];
F_rl_remote=[];F_rl_MI=[];F_rl_MVO = [];
F_lr_remote=[];F_lr_MI=[];F_lr_MVO = [];

LGE_remote = [];    LGE_MI = [];    LGE_MVO=[];
J_remote = [];      J_MI = [];      J_MVO=[];
Lambda_f_remote= [];Lambda_f_MI= [];Lambda_f_MVO= [];


vol_MVO = 0; vol_MI = 0; vol_Remote = 0;
for i = 1: element_above_apex
    if LGEEachLVElem(i)<=mean_remote_intensity
        FF_c_remote = [FF_c_remote;FF_c(i)];
        FF_r_remote = [FF_r_remote;FF_r(i)];
        FF_l_remote = [FF_l_remote;FF_l(i)];
        
        FF_f_remote = [FF_f_remote;FF_f(i)];
        FF_s_remote = [FF_s_remote;FF_s(i)];
        FF_n_remote = [FF_n_remote;FF_n(i)];
        
        F_cr_remote=[F_cr_remote;F_cr(i)];
        F_rc_remote=[F_rc_remote;F_rc(i)];
        F_cl_remote=[F_cl_remote;F_cl(i)];
        F_lc_remote=[F_lc_remote;F_lc(i)];
        F_rl_remote=[F_rl_remote;F_rl(i)];
        F_lr_remote=[F_lr_remote;F_lr(i)];
        
        J_remote=[J_remote;J(i)];
        Lambda_f_remote= [Lambda_f_remote;Lambda_f(i)];
 
         
        LGE_remote = [LGE_remote; LGEEachLVElem(i)];
        
        vol_Remote = vol_Remote + LVWallVol(i);
        
    elseif LGEEachLVElem(i)>=mean_MVO
        FF_c_MVO = [FF_c_MVO;FF_c(i)];
        FF_r_MVO = [FF_r_MVO;FF_r(i)];
        FF_l_MVO = [FF_l_MVO;FF_l(i)];
        
        FF_f_MVO = [FF_f_MVO;FF_f(i)];
        FF_s_MVO = [FF_s_MVO;FF_s(i)];
        FF_n_MVO = [FF_n_MVO;FF_n(i)];
        
        F_cr_MVO=[F_cr_MVO;F_cr(i)];
        F_rc_MVO=[F_rc_MVO;F_rc(i)];
        F_cl_MVO=[F_cl_MVO;F_cl(i)];
        F_lc_MVO=[F_lc_MVO;F_lc(i)];
        F_rl_MVO=[F_rl_MVO;F_rl(i)];
        F_lr_MVO=[F_lr_MVO;F_lr(i)];
        
        J_MVO=[J_MVO;J(i)];
        Lambda_f_MVO= [Lambda_f_MVO;Lambda_f(i)];
         
        LGE_MVO = [LGE_MVO; LGEEachLVElem(i)];
         
        vol_MVO = vol_MVO + LVWallVol(i);
        
    else
        FF_c_MI = [FF_c_MI;FF_c(i)];
        FF_r_MI = [FF_r_MI;FF_r(i)];
        FF_l_MI = [FF_l_MI;FF_l(i)];
        
        FF_f_MI = [FF_f_MI;FF_f(i)];
        FF_s_MI = [FF_s_MI;FF_s(i)];
        FF_n_MI = [FF_n_MI;FF_n(i)];
        
        F_cr_MI=[F_cr_MI;F_cr(i)];
        F_rc_MI=[F_rc_MI;F_rc(i)];
        F_cl_MI=[F_cl_MI;F_cl(i)];
        F_lc_MI=[F_lc_MI;F_lc(i)];
        F_rl_MI=[F_rl_MI;F_rl(i)];
        F_lr_MI=[F_lr_MI;F_lr(i)];
        
        J_MI=[J_MI;J(i)];
        Lambda_f_MI= [Lambda_f_MI;Lambda_f(i)];

        LGE_MI = [LGE_MI; LGEEachLVElem(i)];
         
        vol_MI = vol_MI + LVWallVol(i);
    end
    
end

if isempty(FF_c_MVO)
    FF_c_MVO = 0;
end
if isempty(FF_r_MVO)
    FF_r_MVO = 0;
end
if isempty(FF_l_MVO)
    FF_l_MVO = 0;
end
if isempty(FF_f_MVO)
    FF_f_MVO = 0;
end
if isempty(FF_s_MVO)
    FF_s_MVO = 0;
end
if isempty(FF_n_MVO)
    FF_n_MVO = 0;
end
if isempty(F_cr_MVO)
    F_cr_MVO = 0;
end
if isempty(F_rc_MVO)
    F_rc_MVO = 0;
end
if isempty(F_cl_MVO)
    F_cl_MVO = 0;
end
if isempty(F_lc_MVO)
    F_lc_MVO = 0;
end
if isempty(F_rl_MVO)
    F_rl_MVO = 0;
end
if isempty(F_lr_MVO)
    F_lr_MVO = 0;
end
if isempty(J_MVO)
    J_MVO = 0;
end
if isempty(Lambda_f_MVO)
    Lambda_f_MVO = 0;
end

kkk=0;
if kkk==1
figure(h_FF_crl);
subplot(234);
plot(LGE_MI,FF_c_MI, '.');
subplot(235);
plot(LGE_MI,FF_r_MI, '.');
subplot(236);
plot(LGE_MI,FF_l_MI, '.');

figure(h_FF_fsn);
subplot(234);
plot(LGE_MI,FF_f_MI, '.');
subplot(235);
plot(LGE_MI,FF_s_MI, '.');
subplot(236);
plot(LGE_MI,FF_n_MI, '.');

figure(h_F_crl1);
subplot(234);
plot(LGE_MI,F_cr_MI, '.');
subplot(235);
plot(LGE_MI,F_rc_MI, '.');
subplot(236);
plot(LGE_MI,F_cl_MI, '.');

figure(h_F_crl2);
subplot(234);
plot(LGE_MI,F_lc_MI, '.');
subplot(235);
plot(LGE_MI,F_rl_MI, '.');
subplot(236);
plot(LGE_MI,F_lr_MI, '.');

figure(h_J_Lambda_f);
subplot(2,2,3);
plot(LGE_MI,J_MI,'.');
subplot(2,2,4);
plot(LGE_MI,Lambda_f_MI,'.');

else
end

%%summarize and plot as a table

    
%row_name = {'remote'; 'MI'; 'MVO'; 'MI+MVO'};
%aveFF_c = [mean(FF_c_remote);mean(FF_c_MI);mean(FF_c_MVO); mean([FF_c_MI;FF_c_MVO ])];
%aveFF_r = [mean(FF_r_remote);mean(FF_r_MI);mean(FF_r_MVO); mean([FF_r_MI;FF_r_MVO ])];
%aveFF_l = [mean(FF_l_remote);mean(FF_l_MI);mean(FF_l_MVO); mean([FF_l_MI;FF_l_MVO ])];

%aveFF_f = [mean(FF_f_remote);mean(FF_f_MI);mean(FF_f_MVO); mean([FF_f_MI;FF_f_MVO ])];
%aveFF_s = [mean(FF_s_remote);mean(FF_s_MI);mean(FF_s_MVO); mean([FF_s_MI;FF_s_MVO ])];
%aveFF_n = [mean(FF_n_remote);mean(FF_n_MI);mean(FF_n_MVO); mean([FF_n_MI;FF_n_MVO ])];

%h_table = figure();
%FF_table = table(aveFF_c,aveFF_r,aveFF_l,aveFF_f, aveFF_s, aveFF_n,'RowNames',row_name);
%uitable('Data',FF_table{:,:}','RowName',FF_table.Properties.VariableNames,...
%    'ColumnName',FF_table.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);

row_name = {'remote'; 'MI'; 'MVO'; 'MI+MVO'};
aveF_cc = [mean(FF_c_remote);mean(FF_c_MI);mean(FF_c_MVO); mean([FF_c_MI;FF_c_MVO ])];
aveF_cl = [mean(F_cl_remote);mean(F_cl_MI);mean(F_cl_MVO); mean([F_cl_MI;F_cl_MVO ])];
aveF_cr = [mean(F_cr_remote);mean(F_cr_MI);mean(F_cr_MVO); mean([F_cr_MI;F_cr_MVO ])];
aveF_lc = [mean(F_lc_remote);mean(F_lc_MI);mean(F_lc_MVO); mean([F_lc_MI;F_lc_MVO ])];
aveF_ll = [mean(FF_l_remote);mean(FF_l_MI);mean(FF_l_MVO); mean([FF_l_MI;FF_l_MVO ])];
aveF_lr = [mean(F_lr_remote);mean(F_lr_MI);mean(F_lr_MVO); mean([F_lr_MI;F_lr_MVO ])];
aveF_rc = [mean(F_rc_remote);mean(F_rc_MI);mean(F_rc_MVO); mean([F_rc_MI;F_rc_MVO ])];
aveF_rl = [mean(F_rl_remote);mean(F_rl_MI);mean(F_rl_MVO); mean([F_rl_MI;F_rl_MVO ])];
aveF_rr = [mean(FF_r_remote);mean(FF_r_MI);mean(FF_r_MVO); mean([FF_r_MI;FF_r_MVO ])];
aveF_f = [mean(FF_f_remote);mean(FF_f_MI);mean(FF_f_MVO); mean([FF_f_MI;FF_f_MVO ])];
aveF_s = [mean(FF_s_remote);mean(FF_s_MI);mean(FF_s_MVO); mean([FF_s_MI;FF_f_MVO ])];
aveF_n = [mean(FF_n_remote);mean(FF_n_MI);mean(FF_n_MVO); mean([FF_n_MI;FF_f_MVO ])];
aveJ = [mean(J_remote);mean(J_MI);mean(J_MVO); mean([J_MI;J_MVO ])];
aveLambda_f = [mean(Lambda_f_remote);mean(Lambda_f_MI);mean(Lambda_f_MVO); mean([Lambda_f_MI;Lambda_f_MVO ])];

h_table1 = figure();
F_table1 = table(aveF_cc,aveF_cl,aveF_cr,aveF_lc,aveF_ll,aveF_lr,aveF_rc,aveF_rl,aveF_rr,aveF_f,aveF_s,aveF_n,aveJ,aveLambda_f,'RowNames',row_name);
uitable('Data',F_table1{:,:}','RowName',F_table1.Properties.VariableNames,...
    'ColumnName',F_table1.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);

row_name = {'remote'; 'MI'; 'MVO'; 'MI+MVO'};
stdF_cc = [std(FF_c_remote);std(FF_c_MI);std(FF_c_MVO); std([FF_c_MI;FF_c_MVO ])];
stdF_cl = [std(F_cl_remote);std(F_cl_MI);std(F_cl_MVO); std([F_cl_MI;F_cl_MVO ])];
stdF_cr = [std(F_cr_remote);std(F_cr_MI);std(F_cr_MVO); std([F_cr_MI;F_cr_MVO ])];
stdF_lc = [std(F_lc_remote);std(F_lc_MI);std(F_lc_MVO); std([F_lc_MI;F_lc_MVO ])];
stdF_ll = [std(FF_l_remote);std(FF_l_MI);std(FF_l_MVO); std([FF_l_MI;FF_l_MVO ])];
stdF_lr = [std(F_lr_remote);std(F_lr_MI);std(F_lr_MVO); std([F_lr_MI;F_lr_MVO ])];
stdF_rc = [std(F_rc_remote);std(F_rc_MI);std(F_rc_MVO); std([F_rc_MI;F_rc_MVO ])];
stdF_rl = [std(F_rl_remote);std(F_rl_MI);std(F_rl_MVO); std([F_rl_MI;F_rl_MVO ])];
stdF_rr = [std(FF_r_remote);std(FF_r_MI);std(FF_r_MVO); std([FF_r_MI;FF_r_MVO ])];
stdF_f = [std(FF_f_remote);std(FF_f_MI);std(FF_f_MVO); std([FF_f_MI;FF_f_MVO ])];
stdF_s = [std(FF_s_remote);std(FF_s_MI);std(FF_s_MVO); std([FF_s_MI;FF_s_MVO ])];
stdF_n = [std(FF_n_remote);std(FF_n_MI);std(FF_n_MVO); std([FF_n_MI;FF_n_MVO ])];
stdJ = [std(J_remote);std(J_MI);std(J_MVO); std([J_MI;J_MVO ])];
stdLambda_f = [std(Lambda_f_remote);std(Lambda_f_MI);std(Lambda_f_MVO); std([Lambda_f_MI;Lambda_f_MVO ])];

h_table2 = figure();
F_table2 = table(stdF_cc,stdF_cl,stdF_cr,stdF_lc,stdF_ll,stdF_lr,stdF_rc,stdF_rl,stdF_rr,stdF_f,stdF_s,stdF_n,stdJ,stdLambda_f,'RowNames',row_name);
uitable('Data',F_table2{:,:}','RowName',F_table2.Properties.VariableNames,...
    'ColumnName',F_table2.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);

cd(projectConfig_dir);
fileNameF = sprintf('AveF_component_%s_to_%s.xls', ResultsDirAllScans(1,followUpScanIndex).folderName, ...
            ResultsDirAllScans(1,baselineScanIndex).folderName);

variable = {'aveF_cc';'aveF_cl';'aveF_cr';'aveF_lc';'aveF_ll';'aveF_lr';'aveF_rc';'aveF_rl';'aveF_rr';'aveF_f';'aveF_s';'aveF_n';'aveJ';'aveLambda_f'};
remote =[mean(FF_c_remote);mean(F_cl_remote);mean(F_cr_remote);mean(F_lc_remote);mean(FF_l_remote);mean(F_lr_remote);mean(F_rc_remote);mean(F_rl_remote);mean(FF_r_remote);mean(FF_f_remote);mean(FF_s_remote);mean(FF_n_remote);mean(J_remote);mean(Lambda_f_remote)];
MI=[mean(FF_c_MI);mean(F_cl_MI);mean(F_cr_MI);mean(F_lc_MI);mean(FF_l_MI);mean(F_lr_MI);mean(F_rc_MI);mean(F_rl_MI);mean(FF_r_MI);mean(FF_f_MI);mean(FF_s_MI);mean(FF_n_MI);mean(J_MI);mean(Lambda_f_MI)];
MVO=[mean(FF_c_MVO);mean(F_cl_MVO);mean(F_cr_MVO);mean(F_lc_MVO);mean(FF_l_MVO);mean(F_lr_MVO);mean(F_rc_MVO);mean(F_rl_MVO);mean(FF_r_MVO);mean(FF_f_MVO);mean(FF_s_MVO);mean(FF_n_MVO);mean(J_MVO);mean(Lambda_f_MVO)]; 
MI_MVO = [mean([FF_c_MI;FF_c_MVO]);mean([F_cl_MI;F_cl_MVO]);mean([F_cr_MI;F_cr_MVO]);mean([F_lc_MI;F_lc_MVO]);mean([FF_l_MI;FF_l_MVO]);mean([F_lr_MI;F_lr_MVO]);mean([F_rc_MI;F_rc_MVO]);...
    mean([F_rl_MI;F_rl_MVO]);mean([FF_r_MI;FF_r_MVO]);mean([FF_f_MI;FF_f_MVO]);mean([FF_s_MI;FF_s_MVO]);mean([FF_n_MI;FF_n_MVO]);mean([J_MI;J_MVO]);mean([Lambda_f_MI;Lambda_f_MVO])];
F_table11=table(variable,remote,MI,MVO,MI_MVO);
writetable(F_table11,fileNameF);

fileNameF = sprintf('stdF_component_%s_to_%s.xls', ResultsDirAllScans(1,followUpScanIndex).folderName, ...
            ResultsDirAllScans(1,baselineScanIndex).folderName);

variable = {'stdF_cc';'stdF_cl';'stdF_cr';'stdF_lc';'stdF_ll';'stdF_lr';'stdF_rc';'stdF_rl';'stdF_rr';'stdF_f';'stdF_s';'stdF_n';'stdJ';'stdLambda_f'};
remote =[std(FF_c_remote);std(F_cl_remote);std(F_cr_remote);std(F_lc_remote);std(FF_l_remote);std(F_lr_remote);std(F_rc_remote);std(F_rl_remote);std(FF_r_remote);std(FF_f_remote);std(FF_s_remote);std(FF_n_remote);std(J_remote);std(Lambda_f_remote)];
MI=[std(FF_c_MI);std(F_cl_MI);std(F_cr_MI);std(F_lc_MI);std(FF_l_MI);std(F_lr_MI);std(F_rc_MI);std(F_rl_MI);std(FF_r_MI);std(FF_f_MI);std(FF_s_MI);std(FF_n_MI);std(J_MI);std(Lambda_f_MI)];
MVO=[std(FF_c_MVO);std(F_cl_MVO);std(F_cr_MVO);std(F_lc_MVO);std(FF_l_MVO);std(F_lr_MVO);std(F_rc_MVO);std(F_rl_MVO);std(FF_r_MVO);std(FF_f_MVO);std(FF_s_MVO);std(FF_n_MVO);std(J_MVO);std(Lambda_f_MVO)]; 
MI_MVO = [std([FF_c_MI;FF_c_MVO]);std([F_cl_MI;F_cl_MVO]);std([F_cr_MI;F_cr_MVO]);std([F_lc_MI;F_lc_MVO]);std([FF_l_MI;FF_l_MVO]);std([F_lr_MI;F_lr_MVO]);std([F_rc_MI;F_rc_MVO]);...
    std([F_rl_MI;F_rl_MVO]);std([FF_r_MI;FF_r_MVO]);std([FF_f_MI;FF_f_MVO]);std([FF_s_MI;FF_s_MVO]);std([FF_n_MI;FF_n_MVO]);std([J_MI;J_MVO]);std([Lambda_f_MI;Lambda_f_MVO])];
F_table11=table(variable,remote,MI,MVO,MI_MVO);
writetable(F_table11,fileNameF);
cd(workingDir)

if vol_MVO > 0 
h_box_FF = figure();
for i=1:length(FF_c_remote)
FF_c_Gvar_remote(i) = "F_cc_remote";
end
for i=1:length(FF_c_MI)
FF_c_Gvar_MI(i)=      "F_cc_MI    ";
end
for i=1:length(FF_c_MVO)
FF_c_Gvar_MVO(i)=     "F_cc_MVO   ";
end
FF_c_Gvar = [FF_c_Gvar_remote';FF_c_Gvar_MI';FF_c_Gvar_MVO'];
FF_c_group = [FF_c_remote;FF_c_MI;FF_c_MVO];
subplot(231);
boxplot(FF_c_group,FF_c_Gvar);

for i=1:length(FF_r_remote)
FF_r_Gvar_remote(i) = "F_rr_remote";
end
for i=1:length(FF_r_MI)
FF_r_Gvar_MI(i) =     "F_rr_MI    ";
end
for i=1:length(FF_r_MVO)
FF_r_Gvar_MVO(i) =    "F_rr_MVO   ";
end
FF_r_Gvar = [FF_r_Gvar_remote';FF_r_Gvar_MI';FF_r_Gvar_MVO'];
FF_r_group = [FF_r_remote;FF_r_MI;FF_r_MVO];
subplot(232);
boxplot(FF_r_group,FF_r_Gvar);

for i=1:length(FF_l_remote)
FF_l_Gvar_remote(i) = "F_ll_remote";
end
for i=1:length(FF_l_MI)
FF_l_Gvar_MI(i) =    "F_ll_MI    ";
end
for i=1:length(FF_l_MVO)
FF_l_Gvar_MVO(i) = "F_ll_MVO   ";
end
FF_l_Gvar = [FF_l_Gvar_remote';FF_l_Gvar_MI';FF_l_Gvar_MVO'];
FF_l_group = [FF_l_remote;FF_l_MI;FF_l_MVO];
subplot(233);
boxplot(FF_l_group,FF_l_Gvar);

for i=1:length(FF_f_remote)
FF_f_Gvar_remote(i) = "F_ff_remote";
end
for i=1:length(FF_f_MI)
FF_f_Gvar_MI(i) = "F_ff_MI    ";
end
for i=1:length(FF_f_MVO)
FF_f_Gvar_MVO(i) = "F_ff_MVO   ";
end
FF_f_Gvar = [FF_f_Gvar_remote';FF_f_Gvar_MI';FF_f_Gvar_MVO'];
FF_f_group = [FF_f_remote;FF_f_MI;FF_f_MVO];
subplot(234);
boxplot(FF_f_group,FF_f_Gvar);

for i=1:length(FF_s_remote)
FF_s_Gvar_remote(i) = "F_ss_remote";
end
for i=1:length(FF_s_MI)
FF_s_Gvar_MI(i) = "F_ss_MI    ";
end
for i=1:length(FF_s_MVO)
FF_s_Gvar_MVO(i) = "F_ss_MVO   ";
end
FF_s_Gvar = [FF_s_Gvar_remote';FF_s_Gvar_MI';FF_s_Gvar_MVO'];
FF_s_group = [FF_s_remote;FF_s_MI;FF_s_MVO];
subplot(235);
boxplot(FF_s_group,FF_s_Gvar);

for i=1:length(FF_n_remote)
FF_n_Gvar_remote(i) = "F_nn_remote";
end
for i=1:length(FF_n_MI)
FF_n_Gvar_MI(i) = "F_nn_MI    ";
end
for i=1:length(FF_n_MVO)
FF_n_Gvar_MVO(i) = "F_nn_MVO   ";
end
FF_n_Gvar = [FF_n_Gvar_remote';FF_n_Gvar_MI';FF_n_Gvar_MVO'];
FF_n_group = [FF_n_remote;FF_n_MI;FF_n_MVO];
subplot(236);
boxplot(FF_n_group,FF_n_Gvar);

else
h_box_FF = figure();
for i=1:length(FF_c_remote)
FF_c_Gvar_remote(i) = "F_cc_remote";
end
for i=1:length(FF_c_MI)
FF_c_Gvar_MI(i)=      "F_cc_MI    ";
end
FF_c_Gvar = [FF_c_Gvar_remote';FF_c_Gvar_MI'];
FF_c_group = [FF_c_remote;FF_c_MI];
subplot(231);
boxplot(FF_c_group,FF_c_Gvar);

for i=1:length(FF_r_remote)
FF_r_Gvar_remote(i) = "F_rr_remote";
end
for i=1:length(FF_r_MI)
FF_r_Gvar_MI(i) =     "F_rr_MI    ";
end
FF_r_Gvar = [FF_r_Gvar_remote';FF_r_Gvar_MI'];
FF_r_group = [FF_r_remote;FF_r_MI];
subplot(232);
boxplot(FF_r_group,FF_r_Gvar);

for i=1:length(FF_l_remote)
FF_l_Gvar_remote(i) = "F_ll_remote";
end
for i=1:length(FF_l_MI)
FF_l_Gvar_MI(i) =    "F_ll_MI    ";
end
FF_l_Gvar = [FF_l_Gvar_remote';FF_l_Gvar_MI'];
FF_l_group = [FF_l_remote;FF_l_MI];
subplot(233);
boxplot(FF_l_group,FF_l_Gvar);

for i=1:length(FF_f_remote)
FF_f_Gvar_remote(i) = "F_ff_remote";
end
for i=1:length(FF_f_MI)
FF_f_Gvar_MI(i) = "F_ff_MI    ";
end
FF_f_Gvar = [FF_f_Gvar_remote';FF_f_Gvar_MI'];
FF_f_group = [FF_f_remote;FF_f_MI];
subplot(234);
boxplot(FF_f_group,FF_f_Gvar);

for i=1:length(FF_s_remote)
FF_s_Gvar_remote(i) = "F_ss_remote";
end
for i=1:length(FF_s_MI)
FF_s_Gvar_MI(i) = "F_ss_MI    ";
end
FF_s_Gvar = [FF_s_Gvar_remote';FF_s_Gvar_MI'];
FF_s_group = [FF_s_remote;FF_s_MI];
subplot(235);
boxplot(FF_s_group,FF_s_Gvar);

for i=1:length(FF_n_remote)
FF_n_Gvar_remote(i) = "F_nn_remote";
end
for i=1:length(FF_n_MI)
FF_n_Gvar_MI(i) = "F_nn_MI    ";
end
FF_n_Gvar = [FF_n_Gvar_remote';FF_n_Gvar_MI'];
FF_n_group = [FF_n_remote;FF_n_MI];
subplot(236);
boxplot(FF_n_group,FF_n_Gvar);
end


[vol_MVO/(vol_MVO+vol_MI+vol_Remote)  (vol_MI+vol_MVO)/(vol_MVO+vol_MI+vol_Remote)]

%
%   Green-Larrange strain tensor
%

    
EE_c_remote=[];   EE_r_remote=[];   EE_l_remote=[];
EE_c_MVO = [];    EE_r_MVO=[];      EE_l_MVO=[];
EE_c_MI=[];       EE_r_MI=[];       EE_l_MI=[];
EE_f_remote = []; EE_s_remote = []; EE_n_remote = [];
EE_f_MVO = [];    EE_s_MVO = [];    EE_n_MVO = [];
EE_f_MI = [];     EE_s_MI = [];     EE_n_MI = [];

%LGE_remote = []; LGE_MI = []; LGE_MVO=[];
E_cr_remote = [];   E_cr_MI=[];     E_cr_MVO=[];
E_cl_remote = [];   E_cl_MI=[];     E_cl_MVO=[];
E_rl_remote = [];   E_rl_MI=[];     E_rl_MVO=[];

%vol_MVO = 0; vol_MI = 0; vol_Remote = 0;
for i = 1: element_above_apex
    if LGEEachLVElem(i)<=mean_remote_intensity
        EE_c_remote = [EE_c_remote;EE_c(i)];
        EE_r_remote = [EE_r_remote;EE_r(i)];
        EE_l_remote = [EE_l_remote;EE_l(i)];
        
        EE_f_remote = [EE_f_remote;EE_f(i)];
        EE_s_remote = [EE_s_remote;EE_s(i)];
        EE_n_remote = [EE_n_remote;EE_n(i)];
        
        E_cr_remote = [E_cr_remote;E_cr(i)];
        E_cl_remote = [E_cl_remote;E_cl(i)];
        E_rl_remote = [E_rl_remote;E_rl(i)];

    elseif LGEEachLVElem(i)>=mean_MVO
        
        EE_c_MVO = [EE_c_MVO;EE_c(i)];
        EE_r_MVO = [EE_r_MVO;EE_r(i)];
        EE_l_MVO = [EE_l_MVO;EE_l(i)];
        
        EE_f_MVO = [EE_f_MVO;EE_f(i)];
        EE_s_MVO = [EE_s_MVO;EE_s(i)];
        EE_n_MVO = [EE_n_MVO;EE_n(i)];
        
        E_cr_MVO = [E_cr_MVO;E_cr(i)];
        E_cl_MVO = [E_cl_MVO;E_cl(i)];
        E_rl_MVO = [E_rl_MVO;E_rl(i)];
        
    else
        
        EE_c_MI = [EE_c_MI;EE_c(i)];
        EE_r_MI = [EE_r_MI;EE_r(i)];
        EE_l_MI = [EE_l_MI;EE_l(i)];
        
        EE_f_MI = [EE_f_MI;EE_f(i)];
        EE_s_MI = [EE_s_MI;EE_s(i)];
        EE_n_MI = [EE_n_MI;EE_n(i)];
        
        E_cr_MI = [E_cr_MI;E_cr(i)];
        E_cl_MI = [E_cl_MI;E_cl(i)];
        E_rl_MI = [E_rl_MI;E_rl(i)];
    end
    
end

if isempty(EE_c_MVO)
    EE_c_MVO = 0;
end
if isempty(EE_r_MVO)
    EE_r_MVO = 0;
end
if isempty(EE_l_MVO)
    EE_l_MVO = 0;
end
if isempty(EE_f_MVO)
    EE_f_MVO = 0;
end
if isempty(EE_s_MVO)
    EE_s_MVO = 0;
end
if isempty(EE_n_MVO)
    EE_n_MVO = 0;
end
if isempty(E_cr_MVO)
    E_cr_MVO = 0;
end
if isempty(E_cl_MVO)
    E_cl_MVO = 0;
end
if isempty(E_rl_MVO)
    E_rl_MVO = 0;
end

kkk=0;
if kkk==1
h_EE_crl = figure();hold on;
subplot(231);
plot(LGEEachLVElem(1:element_above_apex),EE_c(1:element_above_apex),'.');
xlabel('LGE');ylabel('E_c_c');
subplot(232);
plot(LGEEachLVElem(1:element_above_apex),EE_r(1:element_above_apex),'.');
xlabel('LGE');ylabel('E_r_r');
subplot(233);
plot(LGEEachLVElem(1:element_above_apex),EE_r(1:element_above_apex),'.');
xlabel('LGE');ylabel('E_l_l');

h_EE_crl1 = figure();hold on;
subplot(231);
plot(LGEEachLVElem(1:element_above_apex),E_cr(1:element_above_apex),'.');
xlabel('LGE');ylabel('E_c_r');
subplot(232);
plot(LGEEachLVElem(1:element_above_apex),E_cl(1:element_above_apex),'.');
xlabel('LGE');ylabel('E_c_l');
subplot(233);
plot(LGEEachLVElem(1:element_above_apex),E_rl(1:element_above_apex),'.');
xlabel('LGE');ylabel('E_r_l');

h_EE_fsn = figure();hold on;
subplot(231);
plot(LGEEachLVElem(1:element_above_apex),EE_c(1:element_above_apex),'.');
xlabel('LGE');ylabel('E_f_f');
subplot(232);
plot(LGEEachLVElem(1:element_above_apex),EE_r(1:element_above_apex),'.');
xlabel('LGE');ylabel('E_s_s');
subplot(233);
plot(LGEEachLVElem(1:element_above_apex),EE_r(1:element_above_apex),'.');
xlabel('LGE');ylabel('E_n_n');

figure(h_EE_crl);
subplot(234);
plot(LGE_MI,EE_c_MI, '.');
subplot(235);
plot(LGE_MI,EE_r_MI, '.');
subplot(236);
plot(LGE_MI,EE_l_MI, '.');

figure(h_EE_crl1);
subplot(234);
plot(LGE_MI,E_cr_MI, '.');
subplot(235);
plot(LGE_MI,E_cl_MI, '.');
subplot(236);
plot(LGE_MI,E_rl_MI, '.');

figure(h_EE_fsn);
subplot(234);
plot(LGE_MI,EE_f_MI, '.');
subplot(235);
plot(LGE_MI,EE_s_MI, '.');
subplot(236);
plot(LGE_MI,EE_n_MI, '.');

else
end

%%summarize and plot as a table

%row_name = {'remote'; 'MI'; 'MVO'; 'MI+MVO'};
%aveEE_c = [mean(EE_c_remote);mean(EE_c_MI);mean(EE_c_MVO); mean([EE_c_MI;EE_c_MVO ])];
%aveEE_r = [mean(EE_r_remote);mean(EE_r_MI);mean(EE_r_MVO); mean([EE_r_MI;EE_r_MVO ])];
%aveEE_l = [mean(EE_l_remote);mean(EE_l_MI);mean(EE_l_MVO); mean([EE_l_MI;EE_l_MVO ])];

%aveEE_f = [mean(EE_f_remote);mean(EE_f_MI);mean(EE_f_MVO); mean([EE_f_MI;EE_f_MVO ])];
%aveEE_s = [mean(EE_s_remote);mean(EE_s_MI);mean(EE_s_MVO); mean([EE_s_MI;EE_s_MVO ])];
%aveEE_n = [mean(EE_n_remote);mean(EE_n_MI);mean(EE_n_MVO); mean([EE_n_MI;EE_n_MVO ])];

%h_table = figure();
%EE_table = table(aveEE_c,aveEE_r,aveEE_l,aveEE_f, aveEE_s, aveEE_n,'RowNames',row_name);
%uitable('Data',EE_table{:,:}','RowName',EE_table.Properties.VariableNames,...
%    'ColumnName',EE_table.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);

row_name = {'remote'; 'MI'; 'MVO'; 'MI+MVO'};
aveE_cc = [mean(EE_c_remote);mean(EE_c_MI);mean(EE_c_MVO); mean([EE_c_MI;EE_c_MVO ])];
aveE_rr = [mean(EE_r_remote);mean(EE_r_MI);mean(EE_r_MVO); mean([EE_r_MI;EE_r_MVO ])];
aveE_ll = [mean(EE_l_remote);mean(EE_l_MI);mean(EE_l_MVO); mean([EE_l_MI;EE_l_MVO ])];
aveE_cr = [mean(E_cr_remote);mean(E_cr_MI);mean(E_cr_MVO); mean([E_cr_MI;E_cr_MVO ])];
aveE_cl = [mean(E_cl_remote);mean(E_cl_MI);mean(E_cl_MVO); mean([E_cl_MI;E_cl_MVO ])];
aveE_rl = [mean(E_rl_remote);mean(E_rl_MI);mean(E_rl_MVO); mean([E_rl_MI;E_rl_MVO ])];
aveE_ff = [mean(EE_f_remote);mean(EE_f_MI);mean(EE_f_MVO); mean([EE_f_MI;EE_f_MVO ])];

h_table1 = figure();
E_table1 = table(aveE_cc,aveE_rr,aveE_ll,aveE_cr, aveE_cl, aveE_rl, aveE_ff,'RowNames',row_name);
uitable('Data',E_table1{:,:}','RowName',E_table1.Properties.VariableNames,...
    'ColumnName',E_table1.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);

row_name = {'remote'; 'MI'; 'MVO'; 'MI+MVO'};
stdE_cc = [std(EE_c_remote);std(EE_c_MI);std(EE_c_MVO); std([EE_c_MI;EE_c_MVO ])];
stdE_rr = [std(EE_r_remote);std(EE_r_MI);std(EE_r_MVO); std([EE_r_MI;EE_r_MVO ])];
stdE_ll = [std(EE_l_remote);std(EE_l_MI);std(EE_l_MVO); mean([EE_l_MI;EE_l_MVO ])];
stdE_cr = [std(E_cr_remote);std(E_cr_MI);std(E_cr_MVO); std([E_cr_MI;E_cr_MVO ])];
stdE_cl = [std(E_cl_remote);std(E_cl_MI);std(E_cl_MVO); std([E_cl_MI;E_cl_MVO ])];
stdE_rl = [std(E_rl_remote);std(E_rl_MI);std(E_rl_MVO); std([E_rl_MI;E_rl_MVO ])];
stdE_ff = [std(EE_f_remote);std(EE_f_MI);std(EE_f_MVO); std([EE_f_MI;EE_f_MVO ])];

h_table2 = figure();
E_table2 = table(stdE_cc,stdE_rr,stdE_ll,stdE_cr,stdE_cl,stdE_rl,stdE_ff,'RowNames',row_name);
uitable('Data',E_table2{:,:}','RowName',E_table2.Properties.VariableNames,...
    'ColumnName',E_table2.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);


cd(projectConfig_dir);
fileNameF = sprintf('AveE_component_%s_to_%s.xls', ResultsDirAllScans(1,followUpScanIndex).folderName, ...
            ResultsDirAllScans(1,baselineScanIndex).folderName);

variable = {'aveE_cc';'aveE_rr';'aveE_ll';'aveE_cr';'aveE_cl';'aveE_rl';'aveE_ff'};
remote =[mean(EE_c_remote);mean(EE_r_remote);mean(EE_l_remote);mean(E_cr_remote);mean(E_cl_remote);mean(E_rl_remote);mean(EE_f_remote)];
MI=[mean(EE_c_MI);mean(EE_r_MI);mean(EE_l_MI);mean(E_cr_MI);mean(E_cl_MI);mean(E_rl_MI);mean(EE_f_MI)];
MVO=[mean(EE_c_MVO);mean(EE_r_MVO);mean(EE_l_MVO);mean(E_cr_MVO);mean(E_cl_MVO);mean(E_rl_MVO);mean(EE_f_MVO)]; 
MI_MVO = [mean([EE_c_MI;EE_c_MVO]);mean([EE_r_MI;EE_r_MVO]);mean([EE_l_MI;EE_l_MVO]);mean([E_cr_MI;E_cr_MVO]);mean([E_cl_MI;E_cl_MVO]);mean([E_rl_MI;E_rl_MVO]);mean([EE_f_MI;EE_f_MVO])];
F_table22=table(variable,remote,MI,MVO,MI_MVO);
writetable(F_table22,fileNameF);

fileNameF = sprintf('stdE_component_%s_to_%s.xls', ResultsDirAllScans(1,followUpScanIndex).folderName, ...
            ResultsDirAllScans(1,baselineScanIndex).folderName);

variable = {'stdE_cc';'stdE_rr';'stdE_ll';'stdE_cr';'stdE_cl';'stdE_rl';'stdE_ff'};
remote =[std(EE_c_remote);std(EE_r_remote);std(EE_l_remote);std(E_cr_remote);std(E_cl_remote);std(E_rl_remote);std(EE_f_remote)];
MI=[std(EE_c_MI);std(EE_r_MI);std(EE_l_MI);std(E_cr_MI);std(E_cl_MI);std(E_rl_MI);std(EE_f_MI)];
MVO=[std(EE_c_MVO);std(EE_r_MVO);std(EE_l_MVO);std(E_cr_MVO);std(E_cl_MVO);std(E_rl_MVO);std(EE_f_MVO)]; 
MI_MVO = [std([EE_c_MI;EE_c_MVO]);std([EE_r_MI;EE_r_MVO]);std([EE_l_MI;EE_l_MVO]);std([E_cr_MI;E_cr_MVO]);std([E_cl_MI;E_cl_MVO]);std([E_rl_MI;E_rl_MVO]);std([EE_f_MI;EE_f_MVO])];
F_table22=table(variable,remote,MI,MVO,MI_MVO);
writetable(F_table22,fileNameF);

cd(workingDir)

stop

if vol_MVO > 0 
h_box_EE = figure();
for i=1:length(EE_c_remote)
EE_c_Gvar_remote(i) = "EE_c_remote";
end
for i=1:length(EE_c_MI)
EE_c_Gvar_MI(i)=      "EE_c_MI    ";
end
for i=1:length(EE_c_MVO)
EE_c_Gvar_MVO(i)=     "EE_c_MVO   ";
end
EE_c_Gvar = [EE_c_Gvar_remote';EE_c_Gvar_MI';EE_c_Gvar_MVO'];
EE_c_group = [EE_c_remote;EE_c_MI;EE_c_MVO];
subplot(231);
boxplot(EE_c_group,EE_c_Gvar);

for i=1:length(EE_r_remote)
EE_r_Gvar_remote(i) = "EE_r_remote";
end
for i=1:length(EE_r_MI)
EE_r_Gvar_MI(i) =     "EE_r_MI    ";
end
for i=1:length(EE_r_MVO)
EE_r_Gvar_MVO(i) =    "EE_r_MVO   ";
end
EE_r_Gvar = [EE_r_Gvar_remote';EE_r_Gvar_MI';EE_r_Gvar_MVO'];
EE_r_group = [EE_r_remote;EE_r_MI;EE_r_MVO];
subplot(232);
boxplot(EE_r_group,EE_r_Gvar);

for i=1:length(EE_l_remote)
EE_l_Gvar_remote(i) = "EE_l_remote";
end
for i=1:length(EE_l_MI)
EE_l_Gvar_MI(i) =    "EE_l_MI    ";
end
for i=1:length(EE_l_MVO)
EE_l_Gvar_MVO(i) = "EE_l_MVO   ";
end
EE_l_Gvar = [EE_l_Gvar_remote';EE_l_Gvar_MI';EE_l_Gvar_MVO'];
EE_l_group = [EE_l_remote;EE_l_MI;EE_l_MVO];
subplot(233);
boxplot(EE_l_group,EE_l_Gvar);

for i=1:length(EE_f_remote)
EE_f_Gvar_remote(i) = "EE_f_remote";
end
for i=1:length(EE_f_MI)
EE_f_Gvar_MI(i) = "EE_f_MI    ";
end
for i=1:length(EE_f_MVO)
EE_f_Gvar_MVO(i) = "EE_f_MVO   ";
end
EE_f_Gvar = [EE_f_Gvar_remote';EE_f_Gvar_MI';EE_f_Gvar_MVO'];
EE_f_group = [EE_f_remote;EE_f_MI;EE_f_MVO];
subplot(234);
boxplot(EE_f_group,EE_f_Gvar);

for i=1:length(EE_s_remote)
EE_s_Gvar_remote(i) = "EE_s_remote";
end
for i=1:length(EE_s_MI)
EE_s_Gvar_MI(i) = "EE_s_MI    ";
end
for i=1:length(EE_s_MVO)
EE_s_Gvar_MVO(i) = "EE_s_MVO   ";
end
EE_s_Gvar = [EE_s_Gvar_remote';EE_s_Gvar_MI';EE_s_Gvar_MVO'];
EE_s_group = [EE_s_remote;EE_s_MI;EE_s_MVO];
subplot(235);
boxplot(EE_s_group,EE_s_Gvar);

for i=1:length(EE_n_remote)
EE_n_Gvar_remote(i) = "EE_n_remote";
end
for i=1:length(EE_n_MI)
EE_n_Gvar_MI(i) = "EE_n_MI    ";
end
for i=1:length(EE_n_MVO)
EE_n_Gvar_MVO(i) = "EE_n_MVO   ";
end
EE_n_Gvar = [EE_n_Gvar_remote';EE_n_Gvar_MI';EE_n_Gvar_MVO'];
EE_n_group = [EE_n_remote;EE_n_MI;EE_n_MVO];
subplot(236);
boxplot(EE_n_group,EE_n_Gvar);

else
h_box_EE = figure();
for i=1:length(EE_c_remote)
EE_c_Gvar_remote(i) = "EE_c_remote";
end
for i=1:length(EE_c_MI)
EE_c_Gvar_MI(i)=      "EE_c_MI    ";
end
EE_c_Gvar = [EE_c_Gvar_remote';EE_c_Gvar_MI'];
EE_c_group = [EE_c_remote;EE_c_MI];
subplot(231);
boxplot(EE_c_group,EE_c_Gvar);

for i=1:length(EE_r_remote)
EE_r_Gvar_remote(i) = "EE_r_remote";
end
for i=1:length(EE_r_MI)
EE_r_Gvar_MI(i) =     "EE_r_MI    ";
end
EE_r_Gvar = [EE_r_Gvar_remote';EE_r_Gvar_MI'];
EE_r_group = [EE_r_remote;EE_r_MI];
subplot(232);
boxplot(EE_r_group,EE_r_Gvar);

for i=1:length(EE_l_remote)
EE_l_Gvar_remote(i) = "EE_l_remote";
end
for i=1:length(EE_l_MI)
EE_l_Gvar_MI(i) =    "EE_l_MI    ";
end
EE_l_Gvar = [EE_l_Gvar_remote';EE_l_Gvar_MI'];
EE_l_group = [EE_l_remote;EE_l_MI];
subplot(233);
boxplot(EE_l_group,EE_l_Gvar);

for i=1:length(EE_f_remote)
EE_f_Gvar_remote(i) = "EE_f_remote";
end
for i=1:length(EE_f_MI)
EE_f_Gvar_MI(i) = "EE_f_MI    ";
end
EE_f_Gvar = [EE_f_Gvar_remote';EE_f_Gvar_MI'];
EE_f_group = [EE_f_remote;EE_f_MI];
subplot(234);
boxplot(EE_f_group,EE_f_Gvar);

for i=1:length(EE_s_remote)
EE_s_Gvar_remote(i) = "EE_s_remote";
end
for i=1:length(EE_s_MI)
EE_s_Gvar_MI(i) = "EE_s_MI    ";
end
EE_s_Gvar = [EE_s_Gvar_remote';EE_s_Gvar_MI'];
EE_s_group = [EE_s_remote;EE_s_MI];
subplot(235);
boxplot(EE_s_group,EE_s_Gvar);

for i=1:length(EE_n_remote)
EE_n_Gvar_remote(i) = "EE_n_remote";
end
for i=1:length(EE_n_MI)
EE_n_Gvar_MI(i) = "EE_n_MI    ";
end
EE_n_Gvar = [EE_n_Gvar_remote';EE_n_Gvar_MI'];
EE_n_group = [EE_n_remote;EE_n_MI];
subplot(236);
boxplot(EE_n_group,EE_n_Gvar);
end



%       thresholdMI=0.25;
%       Nelement=length(LGEMappedelement(:,1));
%       Nnode=length(LGEMappedelement(1,:));
      
%       for i=1:Nelement
%        MIintensityelement(i)=0;
%        for j=1:Nnode
%        inode(j)=LGEMappedelement(i,j);
%        end
             
%        for j=1:Nnode
%        MIintensityelement(i)=MIintensityelement(i)+LGEintensity(inode(j));
%        end
%        MIintensityelement(i)=MIintensityelement(i)/Nnode;
%       end
       
%       kMI=0;
%       kREM=0;     
%       for i=1:Nelement
%           if MIintensityelement(i)>=thresholdMI
%             kMI=kMI+1;
%             MIzoneelement(kMI)=i;
%             MIzoneintensity(kMI)=MIintensityelement(i);
%           else
%             kREM=kREM+1;
%             REMzoneelement(kREM)=i;
%             REMzoneintensity(kREM)=MIintensityelement(i);
%           end
%       end
       
%       for i=1:kMI   
%       vol_MI(i) = 0;
%       for j=1:Nnode
%       inode(j)=LGEMappedelement(MIzoneelement(i),j);
%       end
    
%       for j=1:Nnode
%       xelem(j)=nodeMatMRI(inode(j),2);
%       yelem(j)=nodeMatMRI(inode(j),3);
%       zelem(j)=nodeMatMRI(inode(j),4);
%       end
%       P = [xelem',yelem',zelem'];
%       dt =  DelaunayTri(P);
%       [ch, vol_elem] = convexHull(dt);
%       vol_MI(i) = vol_MI(i) + vol_elem;
%       end
%       vol_MI_Total=sum(vol_MI);
       
%       for i=1:kMI   
%       FF_f_MI(i)=vol_MI(i)*FF_f(MIzoneelement(i));
%       end       
%       FF_f_MI_mean=sum(FF_f_MI)/vol_MI_Total;
%       FF_f_MI_dev=std(FF_f_MI)*(kMI-1)/vol_MI_Total; 
       
%       for i=1:kMI   
%       FF_s_MI(i)=vol_MI(i)*FF_s(MIzoneelement(i));
%       end       
%       FF_s_MI_mean=sum(FF_s_MI)/vol_MI_Total;
%       FF_s_MI_dev=std(FF_s_MI)*(kMI-1)/vol_MI_Total;
       
%       for i=1:kMI   
%       FF_n_MI(i)=vol_MI(i)*FF_n(MIzoneelement(i));
%       end       
%       FF_n_MI_mean=sum(FF_n_MI)/vol_MI_Total;
%       FF_n_MI_dev=std(FF_n_MI)*(kMI-1)/vol_MI_Total;
 
%       for i=1:kMI   
%       FF_c_MI(i)=vol_MI(i)*FF_c(MIzoneelement(i));
%       end       
%       FF_c_MI_mean=sum(FF_c_MI)/vol_MI_Total;
%       FF_c_MI_dev=std(FF_c_MI)*(kMI-1)/vol_MI_Total;   

%       for i=1:kMI   
%       FF_l_MI(i)=vol_MI(i)*FF_l(MIzoneelement(i));
%       end       
%       FF_l_MI_mean=sum(FF_l_MI)/vol_MI_Total;
%       FF_l_MI_dev=std(FF_l_MI)*(kMI-1)/vol_MI_Total;  
       
%       for i=1:kMI   
%       FF_r_MI(i)=vol_MI(i)*FF_r(MIzoneelement(i));
%       end       
%       FF_r_MI_mean=sum(FF_r_MI)/vol_MI_Total;
%       FF_r_MI_dev=std(FF_r_MI)*(kMI-1)/vol_MI_Total;   
          
%       for i=1:kREM   
%       vol_REM(i) = 0;
%       for j=1:Nnode
%       inode(j)=LGEMappedelement(REMzoneelement(i),j);
%       end
%       for j=1:Nnode
%       xelem(j)=nodeMatMRI(inode(j),2);
%       yelem(j)=nodeMatMRI(inode(j),3);
%       zelem(j)=nodeMatMRI(inode(j),4);
%       end
%       P = [xelem',yelem',zelem'];
%       dt =  DelaunayTri(P);
%       [ch, vol_elem] = convexHull(dt);
%       vol_REM(i) = vol_REM(i) + vol_elem;
%       end
%       vol_REM_Total=sum(vol_REM);
       
%       for i=1:kREM   
%       FF_f_REM(i)=vol_REM(i)*FF_f(REMzoneelement(i));
%       end
%       FF_f_REM_mean=sum(FF_f_REM)/vol_REM_Total;
%       FF_f_REM_dev=std(FF_f_REM)*(kREM-1)/vol_REM_Total;
          
%       for i=1:kREM   
%       FF_s_REM(i)=vol_REM(i)*FF_s(REMzoneelement(i));
%       end
%       FF_s_REM_mean=sum(FF_s_REM)/vol_REM_Total;
%       FF_s_REM_dev=std(FF_s_REM)*(kREM-1)/vol_REM_Total; 
       
%       for i=1:kREM   
%       FF_n_REM(i)=vol_REM(i)*FF_n(REMzoneelement(i));
%       end
%       FF_n_REM_mean=sum(FF_n_REM)/vol_REM_Total;
%       FF_n_REM_dev=std(FF_n_REM)*(kREM-1)/vol_REM_Total;          
       
%       for i=1:kREM   
%       FF_c_REM(i)=vol_REM(i)*FF_c(REMzoneelement(i));
%       end
%       FF_c_REM_mean=sum(FF_c_REM)/vol_REM_Total;
%       FF_c_REM_dev=std(FF_c_REM)*(kREM-1)/vol_REM_Total;   
       
%       for i=1:kREM   
%       FF_l_REM(i)=vol_REM(i)*FF_l(REMzoneelement(i));
%       end
%       FF_l_REM_mean=sum(FF_l_REM)/vol_REM_Total;
%       FF_l_REM_dev=std(FF_l_REM)*(kREM-1)/vol_REM_Total;  
       
%       for i=1:kREM   
%       FF_r_REM(i)=vol_REM(i)*FF_r(REMzoneelement(i));
%       end
%       FF_r_REM_mean=sum(FF_r_REM)/vol_REM_Total;
%       FF_r_REM_dev=std(FF_r_REM)*(kREM-1)/vol_REM_Total;  
       
%       cd(projectConfig_dir);
       
%       A={'FF_f_MI_mean',FF_f_MI_mean;'FF_n_MI_mean',FF_n_MI_mean;'FF_s_MI_mean',FF_s_MI_mean;...
%           'FF_c_MI_mean',FF_c_MI_mean;'FF_l_MI_mean',FF_l_MI_mean;'FF_r_MI_mean',FF_r_MI_mean;...
%           'FF_f_REM_mean',FF_f_REM_mean;'FF_n_REM_mean',FF_n_REM_mean;'FF_s_REM_mean',FF_s_REM_mean;...   
%           'FF_c_REM_mean',FF_c_REM_mean;'FF_l_REM_mean',FF_l_REM_mean;'FF_r_REM_mean',FF_r_REM_mean;...
%           'FF_f_MI_dev',FF_f_MI_dev;'FF_n_MI_dev',FF_n_MI_dev;'FF_s_MI_dev',FF_s_MI_dev;...
%           'FF_c_MI_dev',FF_c_MI_dev;'FF_l_MI_dev',FF_l_MI_dev;'FF_r_MI_dev',FF_r_MI_dev;...
%           'FF_f_REM_dev',FF_f_REM_dev;'FF_n_REM_dev',FF_n_REM_dev;'FF_s_REM_dev',FF_s_REM_dev;...
%           'FF_c_REM_dev',FF_c_REM_dev;'FF_l_REM_dev',FF_l_REM_dev;'FF_r_REM_dev',FF_r_REM_dev};      
      
%       FileName = sprintf('FF_%s_to_%s_statistic_results.xlsx', ResultsDirAllScans(1,scanIndex).folderName, ...
%            ResultsDirAllScans(1,scanIndexRef).folderName);
%       xlswrite(FileName,A,'w');
%       cd(workingDir);
%   
end
    end






