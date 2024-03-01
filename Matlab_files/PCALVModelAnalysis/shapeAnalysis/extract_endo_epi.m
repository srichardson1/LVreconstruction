function [sendo, sepi] = extract_endo_epi(abaqusInputData)

node = abaqusInputData.node;
%%%%%%%% node has entries
sendo = node(abaqusInputData.endoNodes,6);
sepi =  node(abaqusInputData.epiNodes,6);  