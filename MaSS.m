%% function: MaSS (Merge a Single Sample)
%% --------------------------------------------------------------
%% Created by: Benjamin Green and Alex Szalay - Johns Hopkins - 01/09/2019
%% --------------------------------------------------------------
%% Description
% run cleanup protocol for multipass phenotyping of a single sample. The
% output .csv files for each image is placed into 
% *\inform output\Phenotyped\Results. Additional
% details located at the github repo: https://github.com/AstropathJHU/MaSS
%% --------------------------------------------------------------
%% input: 
%%% wd: working directory for inform data in the astropath workflow this is 
%%%  defined as: *\inform_data
%%%  the directory should contain the following subfolders
%%%     + Component_Tiffs (component_data.tiff for each image)
%%%     + Phenotyped (The *\Phenotyped folder should contain: all ABs 
%%%         analyzed in seperate folders, designated by their Target names)
%%%     | - + ABX1 (Every ABXN folder should have cell_seg_data.txt for
%%%                     each image)
%%%     | - + ABX2
%%%     | - + ABXN... 
%%%                EX. \inform_data\Phenotyped\CD8
%%%                EX. \inform_data\Phenotyped\Tumor
%%%                EX. \inform_data\Phenotyped\PDL1
%%%                EX. \inform_data\Phenotyped\PD1
%%%                EX. \inform_data\Phenotyped\FoxP3
%%%                EX. \inform_data\Phenotyped\CD163
%%% 
%%%  For segmentations; ABXN with the lowest Opal for each segmentation 
%%%     type should have the *_binary_seg_maps.tif images exported
%%% 
%%%  The binary segmentation maps should contain 4 layers;
%%%     Note: when saving the inform algorithm turn on all mask 
%%%     layers to be sure all segmentations are exported
%%% 
%%% sname: the sample name of the slide
%%% MergeConfig: full path to the merge configuration file
%%% logstring: the intial portion of the logstring (project;cohort;)
%% Usage:
%%% wd = 'W:\Clincial_Specimen\M1_1\inform_data'
%%% sname = 'M1_1'
%%% MergeConfig = '\\bki03\Clinical_Specimen_4\Batch\MergeConfig_01.xlsx'
%%% logstring = '1;2;'
%% output:
% *\inform_data\Phenotyped\Results\Tables\*_cleaned_phenotyped_table.csv &
% *\inform_data\Phenotyped\Results\tmp_ForFiguresTables\* _cleaned_phenotyped_table.mat 
% files for each image. Additional details located at the github repo: 
% https://github.com/AstropathJHU/MaSS
%% --------------------------------------------------------------
%% License: 
% Copyright (c) 2019 Benjamin Green, Alex Szalay.
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:

% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
% CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
% TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
% SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%% --------------------------------------------------------------
%%
function [err_val] = MaSS(wd, sname, MergeConfig, logstring)
%
filepath = fileparts(mfilename('fullpath'));
addpath(genpath(filepath))
%
if nargin < 4
    logstring = '';
    version = '0.01.001';
else
    l = strsplit(logstring, '-');
    if length(l) == 1
        version = '0.01.001';
        logstring = l{1};
    else 
        version = l{2};
        logstring = l{1};
    end
end
%
if strcmp(version,'0.0.2')
    imall = 1;
else 
    imall = 0;
end 
%
err_val = mywritetolog(wd, sname, logstring, '', 1, 'Tables', version);
e_code = err_handl(wd, sname, logstring, [], err_val, 'Tables', '');
if e_code == 1
    return
end
%
err_str = ['-wd: ', replace(wd, '\', '\\'), ' -sname: ', sname];
mywritetolog(wd, sname, logstring, err_str, 2, 'Tables');
%
% get Markers structure
%
err_str = ['parsing MergeConfig file: ', replace(MergeConfig, '\', '\\')];
mywritetolog(wd, sname, logstring, err_str, 2, 'Tables');
%
try
    %    
    [Markers, err_val] = createmarks(MergeConfig);
    %
catch
    err_val = 8;
    err_handl(wd, sname, logstring, [], err_val, 'Tables', '');
    return
end
%
e_code = err_handl(wd, sname, logstring, Markers, err_val, 'Tables', '');
if e_code == 1
    return
end
%
% get all filenames for a directory
%
err_str = 'get filenames and check image output';
mywritetolog(wd, sname, logstring, err_str, 2, 'Tables');
%
try
    %
    [filenms, sum_filenames, err_val] =  getfilenames(wd, Markers);
    %
catch
    err_val = 13;
    err_handl(wd, sname, logstring, Markers, err_val, 'Tables', '');
    return
end
%
e_code = err_handl(wd, sname, logstring, Markers, err_val, 'Tables', '');
if e_code == 1
    return
end
%
% run the file loop for the sample
%
err_str = ['merging ',num2str(length(filenms)), ' file(s)'];
mywritetolog(wd, sname, logstring, err_str, 2, 'Tables');
%
try
    %
    [e, nfiles, Markers] = ...
        fileloop(wd, sname, filenms, sum_filenames, Markers, logstring, imall);
    %
catch EM
    disp(EM.message);
    err_val = 15;
    err_handl(wd, sname, logstring, Markers, err_val, 'Tables', '');
    return
end
%
err_str = ['successfully merged ',num2str(nfiles),...
    ' file(s)'];
mywritetolog(wd, sname, logstring, err_str, 2, 'Tables');
if any(cell2mat(e))
    err_val = 15;
    err_handl(wd, sname, logstring, Markers, err_val, 'Tables', '');
    return
end
%
if length(filenms) == nfiles
    err_str = 'MaSS finished';
    mywritetolog(wd, sname, logstring, err_str, 2, 'Tables');
else
    err_str = 'MaSS failed';
    mywritetolog(wd, sname, logstring, err_str, 2, 'Tables');
end
%
end
%%