function [fnames] = getallfils(wd)
%
% QA files only
%{
pp = [wd,'\Phenotyped\Results\QA_QC\Tables_QA_QC\*.csv'];
%}
% All files
%
pp = [wd,'\Phenotyped\PDL1\*cell_seg_data.txt'];
fils = dir(pp);
nm = {fils(:).name};
%nm = cellfun(@(x) extractBefore(x, '_cleaned'),nm,'Uni',0);
nm = cellfun(@(x) extractBefore(x, '_cell_seg'),nm,'Uni',0);
%
% inform files 
%
inffp = [wd,'\Phenotyped\PDL1'];
fnames.inffn = cellfun(@(x) [inffp, '\',x,'_cell_seg_data.txt'],nm,'Uni',0);
%
% component tiff names
%
imfp = [wd,'\Component_Tiffs'];
fnames.imfn = cellfun(@(x) [imfp, '\',x,'_component_data.tif'],nm,'Uni',0);
%
% segmentation tiff names
%
simfp = [wd,'\Phenotyped\PDL1'];
fnames.simfn = cellfun(@(x) [simfp, '\',x,'_binary_seg_maps.tif'],nm,'Uni',0);

end