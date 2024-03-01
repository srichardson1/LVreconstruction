%%%
%clear all; 
%close all; 
%clc;





path(path, '\Users\sr248u\OneDrive - University of Glasgow\Scott_working_directory\Core\Modelling\PCALVModelAnalysis\segmentation');
path(path, '\Users\sr248u\OneDrive - University of Glasgow\Scott_working_directory\Core\Modelling\PCALVModelAnalysis\BSplineFitting');
path(path, '\Users\sr248u\OneDrive - University of Glasgow\Scott_working_directory\Core\Modelling\PCALVModelAnalysis\meshRelated');
workingDir = pwd();

%%load the patient config file
%[FileName, PathName] = uigetfile( ...
%       {'*.m'}, ...
%        'Pick a file');%[FileName,PathName,~] = uigetfile('*.m');
%projectConfig_dir = PathName;
%projectConfig_name = FileName(1:end-2);

%cd(projectConfig_dir);
%run(projectConfig_name);
%cd(workingDir)

%ResultsDirAllScans="\Users\sr248u\OneDrive - University of Glasgow\Scott_working_directory\Core\Modelling\Results\empa002\002c\earlyDiastole";


baselineScanIndex = 1;
followUpScanIndex = 2;
kk=1; % for standrd mesh
%kk=2; % for dense mesh
%%now need to figure out the distance from the basal plane to the valvular
%%ring

II = [1 0 0;... 
      0 1 0; ...
      0 0 1];

%for scanIndexRef = baselineScanIndex



    
    

        
 %   projectConfig_dir="\Users\sr248u\OneDrive - University of Glasgow\Scott_working_directory\Core\Modelling\Results\empa002\002c\";
 %   ResultsDirAllScans="\Users\sr248u\OneDrive - University of Glasgow\Scott_working_directory\Core\Modelling\Results\empa002\002c\earlyDiastole";   
 %   fiberDir = load('C:\Users\sr248u\OneDrive - University of Glasgow\Scott_working_directory\Core\Modelling\Results\empa002\002c\earlyDiastole\fibreGeneration\Results_fiber_60_45\fiberDir.txt');
 %   sheetDir = load('\Users\sr248u\OneDrive - University of Glasgow\Scott_working_directory\Core\Modelling\Results\empa002\002c\earlyDiastole\fibreGeneration\Results_fiber_60_45\sheetDir.txt');
 %   radDir = load('\Users\sr248u\OneDrive - University of Glasgow\Scott_working_directory\Core\Modelling\Results\empa002\002c\earlyDiastole\fibreGeneration\Results_fiber_60_45\radDir.txt');
 %   cirDir = load('\Users\sr248u\OneDrive - University of Glasgow\Scott_working_directory\Core\Modelling\Results\empa002\002c\earlyDiastole\fibreGeneration\Results_fiber_60_45\cirDir.txt');


    projectConfig_dir=projectConfig_dir_first;
    ResultsDirAllScans=ResultsDirAllScans_first;   
    fiberDir = fiberDir_first;
    sheetDir = sheetDir_first;
    radDir = radDir_first;
    cirDir = cirDir_first;

      
    cd(projectConfig_dir);
    cd(ResultsDirAllScans);
    load abaqusInputData;
    cd(workingDir);
    abaqusInputDataRef = abaqusInputData;
    clear abaqusInputData;
    
    cd(projectConfig_dir);
    cd(ResultsDirAllScans);
    %cd(ResultsDirAllScans(1,scanIndexRef).folderNameFibre);
        cd(workingDir);

%   for scanIndex = followUpScanIndex %scanIndexRef+1 :  totalScanNo %%only calculate the first and the last

        if kk==1
        cd(projectConfig_dir);
        tecFileName = sprintf('tecplotDG_to_mappedmeshstand.dat');
        tecfid = fopen(tecFileName,'w');
        
        tecFileName_strain = sprintf('tecplotStrain_to_mappedmeshstand.dat');
        tecStrainfid = fopen(tecFileName_strain,'w');
        else
        cd(projectConfig_dir);
        tecFileName = sprintf('tecplotDG_%s_to_%s_mappedmeshdense.dat',ResultsDirAllScans);
        tecfid = fopen(tecFileName,'w');
        
        tecFileName_strain = sprintf('tecplotStrain_%s_to_%s_mappedmeshdense.dat',ResultsDirAllScans);
        tecStrainfid = fopen(tecFileName_strain,'w');
        end
            
       
        
    projectConfig_dir=projectConfig_dir_second;
    ResultsDirAllScans=ResultsDirAllScans_second;   
    fiberDir = fiberDir_second;
    sheetDir = sheetDir_second;
    radDir = radDir_second;
    cirDir = cirDir_second;
        
        
        
        
        cd(ResultsDirAllScans);
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
            
            FF_cr(i,1) = FF_RCL(1,2);
            FF_cl(i,1) = FF_RCL(1,3);
            
            FF_rc(i,1) = FF_RCL(2,1);
            FF_rl(i,1) = FF_RCL(2,3);
            
            FF_lc(i,1) = FF_RCL(3,1);
            FF_lr(i,1) = FF_RCL(3,2);            
            
            
            
            
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
        TecplotLVMeshWriteWithCentreDataStrain(tecStrainfid,abaqusInputData,EE_f,EE_s,EE_n,EE_c, EE_r, EE_l);
   
        fclose(tecfid);
        fclose(tecStrainfid);

        
        
        
           
       
       cd(projectConfig_dir);
       cd(ResultsDirAllScans);
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
        
        if kk==1
        cd(projectConfig_dir);
        tecFileName = sprintf('tecplotDG_to_MRICoor_mappedmeshstand.dat');
        tecfid = fopen(tecFileName,'w');
        tecStrainFileName = sprintf('tecplotStrain_to_MRICoor_mappedmeshstand.dat');
        tecStrainfid = fopen(tecStrainFileName,'w');
        cd(workingDir);
        TecplotLVMeshWriteWithCentreData(tecfid,abaqusInputDataInMRICoor,FF_f,FF_s,FF_n,FF_c, FF_r, FF_l)
        TecplotLVMeshWriteWithCentreData(tecStrainfid,abaqusInputDataInMRICoor,EE_f,EE_s,EE_n,EE_c, EE_r, EE_l);
        fclose(tecfid);
        fclose(tecStrainfid);
        else
         cd(projectConfig_dir);
        tecFileName = sprintf('tecplotDG_%s_to_%s_MRICoor_mappedmeshsdense.dat', ResultsDirAllScans);
        tecfid = fopen(tecFileName,'w');
        tecStrainFileName = sprintf('tecplotStrain_%s_to_%s_MRICoor_mappedmeshdense.dat', ResultsDirAllScans);
        tecStrainfid = fopen(tecStrainFileName,'w');
        cd(workingDir);
        TecplotLVMeshWriteWithCentreData(tecfid,abaqusInputDataInMRICoor,FF_f,FF_s,FF_n,FF_c, FF_r, FF_l)
        TecplotLVMeshWriteWithCentreData(tecStrainfid,abaqusInputDataInMRICoor,EE_f,EE_s,EE_n,EE_c, EE_r, EE_l);
        fclose(tecfid);
        fclose(tecStrainfid);      
        end 
 %   end

 

 
    meanFc=mean(FF_c)
    meanFr=mean(FF_r)
    meanFl=mean(FF_l)
    
    mean_12=mean(FF_cr)
    mean_13=mean(FF_cl)
    
    mean_21=mean(FF_rc)
    mean_23=mean(FF_rl)
    
    mean_31=mean(FF_lc)
    mean_32=mean(FF_lr)
    
    
    growth_tensor=[meanFc mean_12 mean_13; mean_21 meanFr mean_23; mean_31 mean_32 meanFl];
    growth_tensor_transpose=growth_tensor';
    
    i1=trace(growth_tensor*growth_tensor_transpose)
    i2=(((trace(growth_tensor_transpose))^2)-(trace(growth_tensor*growth_tensor)))/(2)
    i3=det(growth_tensor)
    
    
    
    stdFc=std(FF_c);
    stdFr=std(FF_r);
    stdFl=std(FF_l);



