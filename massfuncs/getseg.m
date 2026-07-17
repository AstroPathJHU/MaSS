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
% Load every segmentation owner's binary_seg_maps (layer 4). Polygon
% lookup uses the map matching the CellNum's parent SegStatus:
%   coexpression -> endsWith marker (row addcoex keeps)
%   plain lineage -> that marker
%   Other -> primary
%
[seg_maps, im_size, ok] = load_seg_maps(p, Markers);
if ~ok
    p = 18;
    return
end
s = {im_size};
p.size.T = 2; p.size.B = s{1}(1)-1; p.size.L = 2; p.size.R = s{1}(2)-1;
%
% get segmentation outlines for all alternative segmentations
%
im3 = cell(length(Markers.altseg));
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
    % get rows of altseg cells
    %
    trows(:,i1) = ismember(p.fig.Phenotype, s_markers);
    ids = double(p.fig.CellNum(trows(:,i1),:));
    phs = p.fig.Phenotype(trows(:,i1));
    %
    [polys, ok] = lookup_polys_by_kept_marker(ids, phs, Markers, seg_maps);
    if ~ok
        p = 18;
        return
    end
    im3{i1} = polys;
end
%
% primary-seg rows (not claimed by an altseg group)
%
trowsall = sum(trows,2) > 0;
cellids = double(p.fig.CellNum(~trowsall,:));
ph_pri = p.fig.Phenotype(~trowsall);
[im4, ok] = lookup_polys_by_kept_marker(cellids, ph_pri, Markers, seg_maps);
if ~ok
    p = 18;
    return
end
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

function [seg_maps, im_size, ok] = load_seg_maps(p, Markers)
%
ok = true;
im_size = [];
seg_maps = containers.Map('KeyType', 'double', 'ValueType', 'any');
%
owners = [Markers.seg(:); Markers.altseg(:)];
for i1 = 1:numel(owners)
    owner = owners{i1};
    SS = double(Markers.SegStatus(strcmp(Markers.all, owner)));
    iname = fullfile(p.fname.folder, p.fname.name);
    iname = replace(iname, Markers.all{1}, owner);
    iname2 = [extractBefore(iname, "]_cell_seg"), ']_binary_seg_maps.tif'];
    if isempty(extractBefore(iname, "]_cell_seg"))
        iname2 = [extractBefore(iname, "]_CELL_SEG"), ']_binary_seg_maps.tif'];
    end
    try
        im = imread(iname2, 4);
    catch
        ok = false;
        return
    end
    if isempty(im_size)
        im_size = size(im);
    end
    seg_maps(SS) = label2idx(im);
end
%
end

function parent_SS = cellnum_parent_segstatus(ph, Markers)
%
% SegStatus of the marker whose CellNum is on this row:
%   coexpression -> endsWith marker (addcoex keeps that row)
%   plain lineage -> that marker
%   Other -> primary
%
ph = char(string(ph));
if strcmp(ph, 'Other')
    parent_SS = double(Markers.SegStatus(strcmp(Markers.all, Markers.seg{1})));
    return
end
if ~isempty(Markers.add) && any(strcmp(string(Markers.add), string(ph)))
    ew = cellfun(@(x) endsWith(ph, x), Markers.all);
    matches = Markers.all(ew);
    if ~isempty(matches)
        [~, ix] = max(cellfun(@numel, matches));
        parent_SS = double(Markers.SegStatus(strcmp(Markers.all, matches{ix})));
        return
    end
end
ii = strcmp(Markers.all, ph);
if any(ii)
    parent_SS = double(Markers.SegStatus(ii));
else
    parent_SS = double(Markers.SegStatus(strcmp(Markers.all, Markers.seg{1})));
end
%
end

function [polys, ok] = lookup_polys_by_kept_marker(ids, phs, Markers, seg_maps)
%
ok = true;
polys = cell(1, numel(ids));
if isempty(ids)
    return
end
%
phs = cellstr(string(phs));
parent_SS = zeros(numel(ids), 1);
for i1 = 1:numel(ids)
    parent_SS(i1) = cellnum_parent_segstatus(phs{i1}, Markers);
end
%
uSS = unique(parent_SS)';
for SS1 = uSS
    if ~isKey(seg_maps, SS1)
        ok = false;
        return
    end
    im = seg_maps(SS1);
    sel = parent_SS == SS1;
    these = ids(sel);
    if any(these < 1 | these > numel(im))
        ok = false;
        return
    end
    polys(sel) = im(1, these);
end
%
end
