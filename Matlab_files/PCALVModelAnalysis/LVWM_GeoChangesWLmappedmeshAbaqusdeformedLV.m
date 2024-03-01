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

baselineScanIndex = 1;
followUpScanIndex = 2;
kk=1; % for standrd mesh
%kk=2; % for dense mesh
%%now need to figure out the distance from the basal plane to the valvular
%%ring

II = [1 0 0;... 
      0 1 0; ...
      0 0 1];


for scanIndexRef = baselineScanIndex
    cd(projectConfig_dir);
    cd(ResultsDirAllScans(1,scanIndexRef).folderName);
    load abaqusInputData;
    cd(workingDir);
    abaqusInputDataRef = abaqusInputData;
    clear abaqusInputData;
    
    cd(projectConfig_dir);
    cd(ResultsDirAllScans(1,scanIndexRef).folderName);
    cd(ResultsDirAllScans(1,scanIndexRef).folderNameFibre);
    %projectConfig_dir
    %ResultsDirAllScans(1,scanIndexRef).folderName
    %ResultsDirAllScans(1,scanIndexRef).folderNameFibre
    fiberDir = load('fiberDir.txt');
    sheetDir = load('sheetDir.txt');
    radDir = load('radDir.txt');
    cirDir = load('cirDir.txt');
    cd(workingDir);

    for scanIndex = followUpScanIndex %scanIndexRef+1 :  totalScanNo %%only calculate the first and the last

        cd(projectConfig_dir);
        
    %
    %   load x,y,z coordinates of the LV calculated by Abaqus
    %
        load 'coordinatesofLV_251.dat'
        load 'coordinatesofLV_252.dat'
        coordinatesofLV_25(:,:)=0.5*(coordinatesofLV_252(:,:)+coordinatesofLV_252(:,:));
    %
      
        tecFileName = sprintf('tecplotDG_%s_to_%s_mappedmeshAbaqus2.dat',ResultsDirAllScans(1,scanIndex).folderName, ...
            ResultsDirAllScans(1,scanIndexRef).folderName );
        tecfid = fopen(tecFileName,'w');
        
        cd(ResultsDirAllScans(1,scanIndex).folderName)
        %ResultsDirAllScans(1,scanIndex).folderName
        load abaqusInputData;
        cd(workingDir);
 
        FF_r = []; EE_r = [];
        FF_c = []; EE_c = [];
        FF_l = []; EE_l = [];
        
        FF_f = []; EE_f = [];
        FF_s = []; EE_s = [];
        FF_n = []; EE_n = [];
%
%     Update the LV node coordinates of Scan2 with the LV node coordinates of the LV produced by Abaqus 
% 
        for i = 1 : size(coordinatesofLV_25,1)
            abaqusInputData.node(i,1)=coordinatesofLV_25(i,2)/abaqusInputData.scaleTomm;
            abaqusInputData.node(i,2)=coordinatesofLV_25(i,3)/abaqusInputData.scaleTomm;
            abaqusInputData.node(i,3)=coordinatesofLV_25(i,4)/abaqusInputData.scaleTomm;
        end
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
            EE = 1/2*(FF'*FF-II);
%            
%             %%a simple way to figure out local direction 
%             cirDir = [ abaqusInputDataRef.node(plist(3),1) - abaqusInputDataRef.node(plist(2),1), ...
%                        abaqusInputDataRef.node(plist(3),2) - abaqusInputDataRef.node(plist(2),2), ...
%                        abaqusInputDataRef.node(plist(3),3) - abaqusInputDataRef.node(plist(2),3)];
%             cirDir=NormalizationVec(cirDir);
 
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
            
            Rot = [];
            fdir = fiberDir(i, 2:4);
            sdir = sheetDir(i, 2:4);
            ndir = cross(fdir,sdir);
            Rot = [fdir; sdir;ndir];
            FF_fsn = Rot*FF*(Rot');  EE_fsn = Rot*EE*(Rot');
            FF_f(i,1) = FF_fsn(1,1); EE_f(i,1) = EE_fsn(1,1);
            FF_s(i,1) = FF_fsn(2,2); EE_s(i,1) = EE_fsn(2,2);
            FF_n(i,1) = FF_fsn(3,3); EE_n(i,1) = EE_fsn(3,3);
 %
 %     Green-Larange strain tensor in cylindrical coordinate system
 %
            E_cc(i,1)=EE_RCL(1,1);
            E_rr(i,1)=EE_RCL(2,2);
            E_ll(i,1)=EE_RCL(3,3);
            
            E_cr(i,1)=EE_RCL(1,2);
            E_cl(i,1)=EE_RCL(1,3);
            E_rl(i,1)=EE_RCL(2,3);
 %
 %     generate three principal strains
 %
 %          aa(1,1)=E_cc(i,1);aa(1,2)=E_cr(i,1);aa(1,3)=E_cl(i,1);
 %          aa(2,1)=E_cr(i,1);aa(2,2)=E_rr(i,1);aa(2,3)=E_rl(i,1);  
 %          aa(3,1)=E_cl(i,1);aa(3,2)=E_rl(i,1);aa(3,3)=E_ll(i,1);
 %          strain123(i,:)=eig(aa);
 %
 %     take Ecc, Err and Ell as three principal strains for inverse problem
 %
           strain123(i,1)=E_cc(i,1);  
           strain123(i,2)=E_rr(i,1); 
           strain123(i,3)=E_ll(i,1); 
        end

        TecplotLVMeshWriteWithCentreData(tecfid,abaqusInputData,FF_f,FF_s,FF_n,FF_c, FF_r, FF_l);
 
        fclose(tecfid);
 
 %
 %     make strain data files for Inverse problem, element-based 
 %
 %      cd(projectConfig_dir);
 %      StrainforInverseFileName = sprintf('Strain3ForInverseProblemDeform_%s_to_%s_mappedmesh.dat',ResultsDirAllScans(1,scanIndex).folderName, ...
 %           ResultsDirAllScans(1,scanIndexRef).folderName );
 %      save(StrainforInverseFileName,'strain123','-ascii');
 %      cd(workingDir);
     
 %      for i = 1 : size(abaqusInputData.node,1)
       
 %      Dispxyz(i,1) = abaqusInputData.node(i,1) - abaqusInputDataRef.node(i,1);             
 %      Dispxyz(i,2) = abaqusInputData.node(i,2) - abaqusInputDataRef.node(i,2);
 %      Dispxyz(i,3) = abaqusInputData.node(i,3) - abaqusInputDataRef.node(i,3);
       
 %      end
 %
 %      make strain data files for Inverse problem, node-based 
 %
 %      cd(projectConfig_dir);
 %      StrainforInverseFileName = sprintf('Displacement3ForInverseProblemDeform_%s_to_%s_mappedmesh.dat',ResultsDirAllScans(1,scanIndex).folderName, ...
 %           ResultsDirAllScans(1,scanIndexRef).folderName );
 %      save(StrainforInverseFileName,'Dispxyz','-ascii');
 %      cd(workingDir);
           
       
       cd(projectConfig_dir);
       cd(ResultsDirAllScans(1,scanIndexRef).folderName);
       resultDir = pwd();
       cd(workingDir);
    
        NodeMat = abaqusInputData.node;
        ElemMat = abaqusInputData.elem;

        nodeMat(:,2) = NodeMat(:,1)*abaqusInputData.scaleTomm;
        nodeMat(:,3) = NodeMat(:,2)*abaqusInputData.scaleTomm;
        nodeMat(:,4) = NodeMat(:,3)*abaqusInputData.scaleTomm;
        nodeMatMRI = rotationBackToMRICoordinateSystemt(nodeMat,resultDir);
        nodeMat(:,1) = nodeMatMRI(:,2);
        nodeMat(:,2) = nodeMatMRI(:,3);
        nodeMat(:,3) = nodeMatMRI(:,4);
        abaqusInputDataInMRICoor = abaqusInputData;
        abaqusInputDataInMRICoor.node = nodeMat;
        
        cd(projectConfig_dir);
        tecFileName = sprintf('tecplotDG_%s_to_%s_MRICoor_mappedmeshAbaqus2.dat', ResultsDirAllScans(1,scanIndex).folderName, ...
            ResultsDirAllScans(1,scanIndexRef).folderName);
        tecfid = fopen(tecFileName,'w');
        cd(workingDir);
        TecplotLVMeshWriteWithCentreData(tecfid,abaqusInputDataInMRICoor,FF_f,FF_s,FF_n,FF_c, FF_r, FF_l);
        fclose(tecfid);
 
    end

end
    
    meanFc=mean(FF_c)
    meanFr=mean(FF_r)
    meanFl=mean(FF_l)
    stdFc=std(FF_c)
    stdFr=std(FF_r)
    stdFl=std(FF_l)



