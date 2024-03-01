%%now we need to consider how can we align all the 4 series of images togehter
%%assumption: the distance from the first basal slice to the valvular ring
%%does not change
%%there are other possible solutions: the ratio is fixed

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
        'Pick a file');
projectConfig_dir = PathName;
projectConfig_name = FileName(1:end-2);

cd(projectConfig_dir);
run(projectConfig_name);
cd(workingDir)

%%now need to figure out the distance from the basal plane to the valvular
%%ring
totalScanNo =  size(ResultsDirAllScans,2);

for scanIndex = 1 : totalScanNo
    cd(projectConfig_dir);
    cd(ResultsDirAllScans(1,scanIndex).folderName);
    load imDesired;
    cd(workingDir);
    
    LXSlice = LVOTSliceSorted(1).LXSlice;
    SXSlice =  SXSliceSorted(1).SXSlice;
 
    
    LXSliceSelected = LXSlice(timeInstanceSelected);
    SXSliceSelected = SXSlice(timeInstanceSelected);
    
    ImSA = SXSliceSelected.imData;
    ImSAInfo = infoExtract(SXSliceSelected.imInfo);
    
    ImLA = LXSliceSelected.imData;
    ImLAInfo = infoExtract(LXSliceSelected.imInfo);
    
    %%now we project the SA image into LA image
    [px, py]=findCrossLines(ImLA, ImSA, ImLAInfo, ImSAInfo);
    
    h = figure();hold on;
    imshow(ImLA,[]);hold on;
    plot(px,py,'ro','MarkerSize', 10); hold on
    line(px,py,'LineWidth', 2);
    
    
    %%%find out the two points, which will be used to define the annulus
    %%%ring 
    [xa1, ya1] = ginput(1);
    plot(xa1, ya1,'r+');
    [xa2, ya2] = ginput(1);
    plot(xa2, ya2,'r+');
    
    
%     lineVec = [xa2-xa1, ya2-ya1];
%     lineNor = [-lineVec(2) lineVec(1)];
%     lineNor = lineNor./sqrt(lineNor(1)^2 + lineNor(2)^2);
%     DistanceToAnnulusRing = dot(lineNor, [px(2)-xa2, py(2) - ya2]);
%     DistanceToAnnulusRing = abs(DistanceToAnnulusRing);
    DistanceToAnnulusRing = TwoLineSectionDistance(px,py,[xa1 xa2],[ya1 ya2]);
    
    ResultsDirAllScans(scanIndex).DistanceAnnulusRing = DistanceToAnnulusRing;
    ResultsDirAllScans(scanIndex).DistanceAnnulusRingMM = DistanceToAnnulusRing*ImLAInfo.PixelSpacing(1);
    
    if ImLAInfo.PixelSpacing(1) ~= ImLAInfo.PixelSpacing(2)
        disp('warning: pixel spacing in the x and y directions are different');
    end
    
    cd(projectConfig_dir);
    figureName = sprintf('scan%d.png',scanIndex);
    print(figureName, '-dpng');
    cd(workingDir);
    
    close(h);
   
    
end
%
%  look for max z
%
%maxDis = -1000;
%for scanIndex = 1 : totalScanNo
%    if maxDis<ResultsDirAllScans(scanIndex).DistanceAnnulusRingMM
%        maxDis = ResultsDirAllScans(scanIndex).DistanceAnnulusRingMM;
%    end
%end
%
%   look for min mean z
%
minDis = 1000;
for scanIndex = 1 : totalScanNo
    if minDis>ResultsDirAllScans(scanIndex).DistanceAnnulusRingMM
        minDis = ResultsDirAllScans(scanIndex).DistanceAnnulusRingMM;
    end
end

%for scanIndex = 1 : totalScanNo
%    ResultsDirAllScans(scanIndex).moveTowardsApex = maxDis-ResultsDirAllScans(scanIndex).DistanceAnnulusRingMM;
%end

for scanIndex = 1 : totalScanNo
    ResultsDirAllScans(scanIndex).moveTowardsApex = -minDis+ResultsDirAllScans(scanIndex).DistanceAnnulusRingMM;
end

cd(projectConfig_dir)
save ResultsDirAllScans ResultsDirAllScans;
cd(workingDir);

