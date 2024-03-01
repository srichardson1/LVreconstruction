function AbaqusInputToHexMeshWL(abaqusInputDir)

DirConfig;
% workingDir = pwd();
cd(abaqusInputDir);

epiB = -2; endoB = 2;innerB = 0;
Ver=load('Node.txt');
TMesh=load('Element.txt');
domainID = load('domainID.txt');

% EndoS1=load('./AbaqusInput/Endos1.txt');
%  filenameEndoS2='EndoS2.txt';
endofaceS1=[]; endofaceS2=[]; endofaceS3=[];
endofaceS4=[]; endofaceS5=[]; endofaceS6=[];
filenameEndoS1 = 'EndoS1.txt';
% endofaceS1_t = load(filenameEndoS1);
% endofaceS1 = endofaceS1_t(1):endofaceS1_t(3):endofaceS1_t(2);
endofaceS1 = readingFaces(filenameEndoS1);
endofaceS1 = endofaceS1';

filenameEndoS3 = 'EndoS3.txt';
% endofaceS3_t = load(filenameEndoS3);
% endofaceS3 = endofaceS3_t(1):endofaceS3_t(3):endofaceS3_t(2);
endofaceS3 = readingFaces(filenameEndoS3);
endofaceS3 = endofaceS3';

% endofaceS2_t=load(filenameEndoS2);
% endofaceS2 = endofaceS2_t(1):endofaceS2_t(3):endofaceS2_t(2);
% endofaceS2 = endofaceS2';
% endofaceS2 = readingFaces(filenameEndoS2);

%%%%base
basefaceS1=[]; basefaceS2=[]; basefaceS3=[];
basefaceS4=[]; basefaceS5=[]; basefaceS6=[];
filenameBaseS1='BaseS1.txt';
% basefaceS1_t = load(filenameBaseS1);
% basefaceS1 = basefaceS1_t(1):basefaceS1_t(3):basefaceS1_t(2);
basefaceS1 = readingFaces(filenameBaseS1);
basefaceS1 = basefaceS1';

%%%epi
epifaceS1=[]; epifaceS2=[]; epifaceS3=[];
epifaceS4=[]; epifaceS5=[]; epifaceS6=[];
filenameEpiS2='EpiS2.txt';
% epifaceS2_t = load(filenameEpiS2);
% epifaceS2 = epifaceS2_t(1):epifaceS2_t(3):epifaceS2_t(2);
epifaceS2 = readingFaces(filenameEpiS2);
epifaceS2 = epifaceS2';

filenameEpiS5='EpiS5.txt';
% epifaceS5_t = load(filenameEpiS5);
% epifaceS5 = epifaceS5_t(1):epifaceS5_t(3):epifaceS5_t(2);
epifaceS5 = readingFaces(filenameEpiS5);
epifaceS5 = epifaceS5';


%%%output to tecplot for check 
%fid = fopen('TecPlot_endoepiSurface.dat','w');
%fprintf(fid, 'TITLE = "endo surface plot" \n');
%fprintf(fid, 'VARIABLES = "x", "y", "z", "u", "v", "w" \n');
%totalEleQ = length(endofaceS1)+length(endofaceS3)+length(epifaceS2)+length(epifaceS5);
%fprintf(fid, 'ZONE T="Endo Mesh", N = %d, E=%d, F=FEPOINT, ET=QUADRILATERAL\n', size(Ver, 1), totalEleQ);
%for i = 1 : size(Ver,1)
%    fprintf(fid,'%f\t%f\t%f\t%f\t%f\t%f\n',Ver(i,2),Ver(i,3),Ver(i,4),0,0,0);
%end
%length(endofaceS1)
%for ii = 1 : length(endofaceS1)
%    fprintf(fid, '%d\t%d\t%d\t%d\n', TMesh(endofaceS1(i),2), TMesh(endofaceS1(i),3),TMesh(endofaceS1(i),4),TMesh(endofaceS1(i),5));
%TMesh(endofaceS1(ii),2)
%TMesh(endofaceS1(ii),3)
%TMesh(endofaceS1(ii),4)
%TMesh(endofaceS1(ii),5)
%end
%for i = 1 : length(endofaceS3)
%    fprintf(fid, '%d\t%d\t%d\t%d\n', TMesh(endofaceS3(i),2), TMesh(endofaceS3(i),6),TMesh(endofaceS3(i),7),TMesh(endofaceS3(i),3));
%end
%for i = 1 : length(epifaceS2)
%    fprintf(fid, '%d\t%d\t%d\t%d\n', TMesh(epifaceS2(i),6), TMesh(epifaceS2(i),9),TMesh(epifaceS2(i),8),TMesh(epifaceS2(i),7));
%end
%for i = 1 : length(epifaceS5)
%    fprintf(fid, '%d\t%d\t%d\t%d\n', TMesh(epifaceS5(i),4), TMesh(epifaceS5(i),8),TMesh(epifaceS5(i),9),TMesh(epifaceS5(i),5));
%end
%fclose(fid);




%%%%generage tetgen file 
fid=fopen('heart_real.1.node','w');
for ii = 1 : size(Ver,1)
    fprintf(fid, '%d\t%f\t%f\t%f\n', Ver(ii,1), Ver(ii,2), Ver(ii,3), Ver(ii,4));
end
fclose(fid);

fid=fopen('heart_real.1.ele','w');

for ii = 1 : size(TMesh,1)
    fprintf(fid,'%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n', TMesh(ii,1),TMesh(ii,2),TMesh(ii,3),TMesh(ii,4),TMesh(ii,5), ...
                                                                   TMesh(ii,6),TMesh(ii,7),TMesh(ii,8),TMesh(ii,9));
  
end
fclose(fid);

index = 0;
fid = fopen('heart_real.1.face','w');
for ii = 1 : length(endofaceS1)
    index = index+1;
    fprintf(fid,'%d\t%d\t%d\t%d\t%d\t%d\n', index, TMesh(endofaceS1(ii),2), TMesh(endofaceS1(ii),3),TMesh(endofaceS1(ii),4),TMesh(endofaceS1(ii),5),endoB);
end
for ii = 1 : length(endofaceS3)
    index = index+1;
    fprintf(fid,'%d\t%d\t%d\t%d\t%d\t%d\n', index, TMesh(endofaceS3(ii),2), TMesh(endofaceS3(ii),6),TMesh(endofaceS3(ii),7),TMesh(endofaceS3(ii),3),endoB);
end

%%%base
for ii = 1 : length(basefaceS1)
    index = index+1;
    fprintf(fid,'%d\t%d\t%d\t%d\t%d\t%d\n', index, TMesh(basefaceS1(ii),2),TMesh(basefaceS1(ii),3), TMesh(basefaceS1(ii),4),TMesh(basefaceS1(ii),5),innerB);
end

%%%%epi
for ii = 1 : length(epifaceS2)
    index = index+1;
    fprintf(fid,'%d\t%d\t%d\t%d\t%d\t%d\n', index, TMesh(epifaceS2(ii),6),TMesh(epifaceS2(ii),9), TMesh(epifaceS2(ii),8),TMesh(epifaceS2(ii),7),epiB);
end

for ii = 1 : length(epifaceS5)
    index = index+1;
    fprintf(fid,'%d\t%d\t%d\t%d\t%d\t%d\n', index, TMesh(epifaceS5(ii),4),TMesh(epifaceS5(ii),8), TMesh(epifaceS5(ii),9),TMesh(epifaceS5(ii),5),epiB);
end

fclose(fid);


cd(workingDir);


function facelist=readingFaces(filename)
fid = fopen(filename);
EndoS1=[];
while ~feof(fid)
    tline=fgetl(fid);
    faces=sscanf(tline,'%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f');
    if ~isempty(faces)
       EndoS1=[EndoS1 faces'];
    end
end
facelist=EndoS1;
fclose(fid);





