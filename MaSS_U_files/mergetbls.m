function [fData] = mergetbls(fname,Markers,wd)
%
% read in data
%
[C, units] = readalltxt(fname,Markers,wd);
%
% select proper phenotypes
%
f = getphenofield(C, Markers, units);
%
% remove cells within X number of pixels in cells
%
d = getdistinct(f,Markers);
%
% removes double calls in hierarchical style
%
q = getcoex(d,Markers);
%
% get polygons from inform and remove other cells inside secondary 
% segmentation polygons
%
a = getseg(q,Markers);
%
% determine expression markers by cells that are in X radius
%
fData = getexprmark(a,Markers);
%
% save the data
%
ii = parsave(fData,fname,Markers,wd);
%
end