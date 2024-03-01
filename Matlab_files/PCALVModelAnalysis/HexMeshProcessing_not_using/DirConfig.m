global abaqusInputDir_input dataResultinput abaqusDirinput


abaqusInputDir = abaqusInputDir_input;
dataResult = dataResultinput; 
abaqusDir = abaqusDirinput;



workingDir = pwd();
cd(abaqusInputDir);
if ~exist(dataResult, 'dir')
    mkdir(dataResult);
end
cd(dataResult);
dataResult = pwd();
cd(workingDir);



