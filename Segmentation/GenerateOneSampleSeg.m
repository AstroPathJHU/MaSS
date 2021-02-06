%% GenerateOneSeg for Vectra data given inform & merge tables output
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2018
%% --------------------------------------------------------------
%% Description
%%% This function will loop through each image in a case and call the
%%% GenerateOneSeg function on all of them
%% --------------------------------------------------------------
%% input: 
%%% Must input two variables to the code
%%%      1. casenum -->> SpecimenID
%%%      2. wd -->> Clinical_Specimen folder
%% Usage: 
%%% wd =  'W:\Clinical_Specimen\'; 
%%% casenum = 'M41_1';
%%%
%% output: For each image in the case it will generate the following output
%%% \inform_data\Component_Tiffs\M41_1_[33870,10031]_component_data_w_seg.tif
%%% This is a 13 layer image stack; the first 8 layers are from the orginal
%%% component_data; the last 5 correspond to -> 1. Tissue Seg; 2. NuclearSeg_1; 
%%% 3. NuclearSeg_2; 4. MembraneSeg_1; 5. MembraneSeg_2
%% --------------------------------------------------------------
%%
function [] = GenerateOneSampleSeg(casenum,wd2)
%
% get inform directory
%
wd = [wd2,casenum,'\inform_data'];
if exist(wd,'dir')
    %
    % check to make sure it exist before moving forward
    %
    pathtin = [wd,'\Phenotyped\Results\Tables\*table.csv'];
    fnamestin = dir(pathtin);
    nmtin = {fnamestin(:).name};
    %
    % get X and Y coords of images
    %
    regdata = extractBetween(nmtin,'[',']');
    regdata = cellfun(@(x)strsplit(x,','),regdata,'Uni',0);
    nl= vertcat(regdata{:});
    %
    % create coord table
    %
    CoordTable = table(nmtin',nl(:,1),nl(:,2),...
        'VariableNames',{'filename','ImageX','ImageY'});
    %
    % loop through each image to 'new' generate segmenation maps
    %
    parfor i1 = 1:height(CoordTable)
        fname = nmtin{i1};
        GenerateOneSeg(wd,fname, CoordTable(i1,:));
    end
end
end