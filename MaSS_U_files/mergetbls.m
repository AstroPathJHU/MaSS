function [fData,logf] = mergetbls(fname,Markers,wd)
%
% read in data
%
logf ='read in data';
[C, p2, units] = readalltxt(fname,Markers,wd);
%
if p2 ~= length(Markers.all)
    fData = [];
    logf = 'path error in file';
    return
end
%
% select proper phenotypes
%
logf = [logf,';get phenotype fields'];
f = getphenofield(C, Markers, units);
%
% remove cells within X number of pixels in cells
%
logf = [logf,';find other cells'];
d = getdistinct(f,Markers);
%
% removes double calls in hierarchical style
%
logf = [logf,';resolve coexpressions of lineage markers'];
q = getcoex(d,Markers);
%
% get polygons from inform and remove other cells inside secondary 
% segmentation polygons
%
logf = [logf,';find segmentation, remove others in tumor'];
a = getseg(q,Markers);
%
% determine expression markers by cells that are in X radius
%
logf = [logf,';find expression marker expression'];
fData = getexprmark(a,Markers);
%
% save the data
%
logf = [logf,';save tables'];
ii = parsave(fData,fname,Markers,wd);
logf = [logf,';make table for this figure: ',ii];
%
end