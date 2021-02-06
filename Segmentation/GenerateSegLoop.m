%% GenerateSegLoop for all cases in a folder
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2018
%% --------------------------------------------------------------
%% Will create the segmentation maps for mergetbls.m output for each case in a folder
%%% This function removes cells not in a segmentation map and relabels the CellNums
%%% from inform output and replaces them with the corresponding unique CellIDs
%% --------------------------------------------------------------
%% input: 
%%% Must input three variables to the code
%%%      1. wd -> a Clinical_Specimen directory
%% Usage: 
%%% wd2 =  'W:\Clinical_Specimen';
%% output: This function outputs the following for each image
%%% \inform_data\Component_Tiffs\M41_1_[33870,10031]_component_data_w_seg.tif
%%% This is a 13 layer image stack; the first 8 layers are from the orginal
%%% component_data; the last 5 correspond to -> 1. Tissue Seg; 2. NuclearSeg_1; 
%%% 3. NuclearSeg_2; 4. MembraneSeg_1; 5. MembraneSeg_2
%%--------------------------------------------------------------
%%
function[] = GenerateSegLoop(wd2)
% wd2 = 'E:\Clinical_Specimen\';
%
% get only clinical specimen sample names
%
fn = dir(wd2);
fd = fn(3:end);
ii = [fd.isdir];
fd = fd(ii);
samplenames = {fd(:).name};
%
ii = (contains(samplenames, 'Batch')...
    |contains(samplenames, 'tmp_inform_data')|...
    contains(samplenames, 'reject')|...
    contains(samplenames, 'Control')|...
    contains(samplenames, 'Flatfield'));
u1 = samplenames(~ii);
%
% loop through each sample and generate segmenation maps
%
for u2 = 1:length(u1)
    casenum = u1{u2};
    disp(casenum);
    GenerateOneSampleSeg(casenum, wd2);
end
end