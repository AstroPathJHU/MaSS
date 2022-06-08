%% function: getseg;
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description
% remove other cells with cell centers within non primary segmentation 
% boundaries; If a cell is on the boundary of the tissue the algorithm
% will close the cell to complete the calculation
%% --------------------------------------------------------------
%%
function [p] = getseg(q,Markers)
p = q;
%
% find the rows which do not have a membrane
%
mrows = isnan(p.fig.MeanMembraneDAPI);
nmr = p.fig(mrows,:);
p.fig = p.fig(~mrows,:);
%
% make global cellids for the image before segmentation correction
%
CellID = (1:1:height(p.fig))';
p.fig.Properties.VariableNames('CellID') = {'CellNum'};
p.fig = [table(CellID),p.fig];
%
% find the other cells before the tumor segmentation correction
%
orows = find(strcmp(p.fig.Phenotype, 'Other'));
%
% get segmentation outlines for all alternative segmentations
%
im3 = cell(length(Markers.altseg));
tcellids = cell(length(Markers.altseg));
trows = false(max(CellID),length(Markers.altseg));
%
for i1 = 1:length(Markers.altseg)
    %
    markalt = Markers.altseg{i1};
    idx = ismember(Markers.all, markalt);
    %
    SS = Markers.SegStatus(idx);
    s_markers_idx = Markers.SegStatus == SS & ismember(Markers.all,...
        Markers.lin)';
    s_markers = Markers.all(s_markers_idx);
    if ~isempty(Markers.add)
        iis = startsWith(Markers.add, s_markers);
        s_markers = [s_markers, Markers.add(iis)];
    end
    %
    % get folder and image names for altseg
    %
    fdname = [extractBefore(p.fname.folder,Markers.all{1}),...
        markalt,'\'];
    iname = [fdname,extractBefore(p.fname.name,...
        "]_cell_seg"),']_binary_seg_maps.tif'];
    if isempty(iname)
        iname = [fdname,extractBefore(p.fname.name,...
            "]_CELL_SEG"),']_binary_seg_maps.tif'];
    end
    %
    % get rows of altseg cells
    %
    trows(:,i1) = ismember(p.fig.Phenotype, s_markers);
    %
    % get inForm cellids of altseg cells
    %
    tcellids{i1} = double(p.fig.CellNum(trows(:,i1),:));
    %
    % read in the image for segmentation 
    %
    im = imread(iname,4);
    %
    % convert it to a linear index of labels
    %
    im = label2idx(im);
    im3{i1} = im(1,tcellids{i1});
end
%
%get filenames for 1ry seg images
%
iname = fullfile(p.fname.folder,p.fname.name);
iname = replace(iname, Markers.all{1}, Markers.seg{1});
iname = [extractBefore(iname,"]_cell_seg"),']_binary_seg_maps.tif'];
if isempty(iname)
    iname = [extractBefore(iname,"]_CELL_SEG"),']_binary_seg_maps.tif'];
end
%
% get cellids of 1ry seg cells
%
trowsall = sum(trows,2) > 0;
cellids = double(p.fig.CellNum(~trowsall,:));
%
% read in image convert to linear index 
%
im2 = imread(iname,4);
im4 = label2idx(im2);
%
% get image dimensions
%
s = {size(im2)};
p.size.T = 2; p.size.B = s{1}(1)-1; p.size.L = 2; p.size.R = s{1}(2)-1;
%
% Remove non 1ry cells
%
im4 = im4(1,cellids);
%
% get expected size of cell vector
%
s2 = max(CellID);
%
% compile the segmentation images
%
im5 = cell(1,s2);
%
% first input all altsegs
%
for i1 = 1:length(Markers.altseg)
    im5(CellID(trows(:,i1))) = im3{i1};
end
%
% next input the 1ry segmentations
%
im5(CellID(~trowsall)) = im4;
%
% fill the cells and get them in the proper format
%
[obj,s] = OrganizeCells(p,s,s2,im5);
%
% remove others in altsegs
%
objt = obj(trowsall);
objt = cat(1,objt{:});
%
X = p.fig.CellXPos(orows) + 1;
Y = p.fig.CellYPos(orows) + 1;
Oth = sub2ind(s{1},Y,X);
rows = ismember(Oth,objt);
%
% count number of deleted others
%
p.flags.segclean = sum(rows == 1);
%
% remove those others from the output p.fig
%
p.fig(orows(rows),:) = [];
%
% remove those others from the CellID vector
%
CellID(orows(rows)) = [];
%
% remove those others from the cell objects variable and save to p.obj
%
p.obj = obj(CellID);
%
% add on those cells without membranes to the bottom of p.fig
%
nmr.CellNum = nmr.CellID; 
p.fig = vertcat(p.fig,nmr);
%
% get new and final CellIDs
%
CellID = (1:1:height(p.fig))';
p.fig.CellID = CellID;
%
end
