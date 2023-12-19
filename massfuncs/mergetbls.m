%% function:mergetbls; merge tables function for inform output
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description
% function takes in a filename, a data structure with marker information,
% and a directory location
% it generates a *\Tables\ directory and merged inform files
%% --------------------------------------------------------------
%%
function [fData, output, e_code, err_msg] = mergetbls(fname, sum_fname,... 
    Markers, wd, imall, i1, seg_markers, dep_markers)
%
fData = [];
err_msg = '';
%
% read in data
%
[C, units, e_code, output, err_msg] = readalltxt(fname, sum_fname,... 
    Markers, wd, i1, seg_markers, dep_markers);
%
if e_code ~= 0
    return
end
%
% units check
%
if any(~strcmp(units, 'pixels'))
    e_code = 17;
end
%
% select proper phenotypes
%
f = getphenofield(C, Markers, units);
%
% remove cells within X number of pixels in cells
%
d = getdistinct(f, Markers);
%
% removes double calls in hierarchical style
%
q = getcoex(d, Markers);
%
% get polygons from inform and remove other cells inside secondary 
% segmentation polygons
%
a = getseg(q, Markers);
if isa(a, 'double')
    e_code = 18;
    return
end
%
% determine expression markers by cells that are in X radius
%
fData = getexprmark(a, Markers);
%
% save the data
%
ii = parsave(fData, fname, Markers, wd, imall);
%
end
