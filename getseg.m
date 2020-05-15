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
% get segmentation outlines for all alternative segmentations
%
im3 = cell(length(Markers.altseg), 1);
tcellids = cell(length(Markers.altseg),1);
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
    %
    % get folder and image names for altseg
    %
    fdname = [extractBefore(p.fname.folder,Markers.all{1}),...
        markalt,'\'];
    iname = [fdname,extractBefore(p.fname.name,'cell_seg_data.txt'),...
        'binary_seg_maps.tif'];
    %
    % get rows of altseg cells
    %
    trows(:,i1) = ismember(p.fig.Phenotype, s_markers);
    %
    xpos = p.fig.CellXPos(ismember(p.fig.Phenotype, s_markers));
    ypos = p.fig.CellYPos(ismember(p.fig.Phenotype, s_markers));
    %
    % get inForm cellids of altseg cells
    %
    tcellids{i1} = double(p.fig.CellNum(trows(:,i1),:));
    %
    % read in the image for segmentation
    %
    im = imread(iname,4);
    sz = size(im);
    ind = sub2ind(sz,ypos,xpos);
    %
    im = -1 .* (single(im) - 1);
    im = bwconncomp(im,4) ;
    cells = im.PixelIdxList;
    cells = cells(2:end);
    %
    im9 = cell(1,length(cells));
    for i2 = 1:length(cells)
        loc = ismember(ind, cells{i2});
        if find(loc)
            im9_ind = tcellids{i1}(loc);
            im9(im9_ind) = cells(i2); 
        end
    end
    %
    % convert it to a linear index of labels
    %
    im3{i1} = im9(1,tcellids{i1});
end
%
%get filenames for 1ry seg images
%
iname = fullfile(p.fname.folder,p.fname.name);
iname = replace(iname, Markers.all{1}, Markers.seg{1});
iname = replace(iname, 'cell_seg_data.txt','binary_seg_maps.tif');
%
% get cellids of 1ry seg cells
%
trowsall = sum(trows,2) > 0;
cellids = double(p.fig.CellNum(~trowsall,:));
xpos = p.fig.CellXPos(~trowsall,:);
ypos = p.fig.CellYPos(~trowsall,:);
%
% read in image convert to linear index 
%
im2 = imread(iname,4);
s = {size(im2)};
%
sz = size(im2);
ind = sub2ind(sz,ypos,xpos);
%
im2 = -1 .* (single(im2) - 1);
im2 = bwconncomp(im2,4) ;
cells = im2.PixelIdxList;
cells = cells(2:end);
%
im9 = cell(1,length(cells));
for i2 = 1:length(cells)
    loc = ismember(ind, cells{i2});
    if find(loc)
        im9_ind = cellids(loc);
        im9(im9_ind) = cells(i2);
    end
end
%
% get image dimensions
%
p.size.T = 2; p.size.B = s{1}(1)-1; p.size.L = 2; p.size.R = s{1}(2)-1;
%
% Remove non 1ry cells
%
im4 = im9(1,cellids);
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
on_mem = find(cellfun(@isempty, im5));
im5(on_mem) = [];
trowsall(on_mem) = []; 
p.fig(on_mem,:) = [];
%
% find the other cells before the tumor segmentation correction
%
orows = find(strcmp(p.fig.Phenotype, 'Other'));
CellID = (1:1:height(p.fig))';
s2 = max(CellID);
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
X = round(p.fig.CellXPos(orows));
Y = round(p.fig.CellYPos(orows));
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
%% function: OrganizeCells; 
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description
% close the cells on the edge of the tissue and organize the numbers
% inside the object cell to a circle so that the inside of the cell can
% be filled
%% --------------------------------------------------------------
%%
function [obj,s] = OrganizeCells(p,s,s2,im5)
%
% convert linear index to x,y
%
s = repmat(s,1,s2);
[y,x] = cellfun(@ind2sub,s,im5,'Uni',0);
%
% locate cells on edge of tissue & fill edges of cells on edge of image
%
X = cellfun(@min,x);
Y = cellfun(@min,y);
Y2 = cellfun(@max,y);
X2 = cellfun(@max,x);
rows = find(Y == p.size.T|Y2 == p.size.B|X == p.size.L|X2==p.size.R);
[x(rows),y(rows)] = cellfun(@(x,y)...
    filledges([x,y],p.size.T,p.size.B,p.size.L,p.size.R),...
    x(1,rows),y(1,rows),'Uni',0);
%
% find midpoints of the cell objects
%
mx = cellfun(@mean,x);
my = cellfun(@mean,y);
%
% get the angle of each point from the midpoint
%
ang = cellfun(@(x,y,mx,my)atan2(y-my,x-mx),...
    x,y,num2cell(mx),num2cell(my),'Uni',0);
%
% sorting the angles from the midpoint we can get indicies for a circle
%
[~,o] = cellfun(@sort,ang,'Uni',0);
%
% order the cells in the circle
%
x = cellfun(@(x,o)x(o),x,o,'Uni',0);
y = cellfun(@(y,o)y(o),y,o,'Uni',0);
%
% put the cells in reference coordinate rectangles
%
x1 = cellfun(@(x,y)x-y,x,num2cell(X),'Uni',0);
y1 = cellfun(@(x,y)x-y,y,num2cell(Y),'Uni',0);
%
% fill objects
%
obj = arrayfun(@(X,X2,Y,Y2)zeros(Y2-Y+1,X2-X+1),X,X2,Y,Y2,'Uni',0);
obj = cellfun(@fillobj,obj,y1,x1,'Uni',0);
[objy,objx] = cellfun(@(x)find(x==1),obj,'Uni',0);
%
% convert back to orginal image coordinates
%
objy = cellfun(@(x,y)x+y-1,objy,num2cell(Y),'Uni',0);
objx = cellfun(@(x,y)x+y-1,objx,num2cell(X),'Uni',0);
obj = cellfun(@sub2ind,s,objy,objx,'Uni',0);
%
% if the cell is very small then the obj is written as a column vector
% instead of a row vector
% to fix this find those cells
%
[~,r2] = cellfun(@size,obj,'Uni',0);
r2 = find(cat(1,r2{:})~=1);
%
% check if any of these cells exist
%
if ~isempty(r2)
    %
    % correct them to proper dimensions
    %
    for i4 = 1:length(r2)
        obj{r2(i4)} =  obj{r2(i4)}';
    end
end
end
%% function: filledges; 
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description
% this function takes in a cell object bon the edge of the image and 
% fills in the boundaries using T,B,L,R to define the edges as Top,
% Bottom, Left, Right
%% --------------------------------------------------------------
%%
function[x,y] = filledges(b,T,B,L,R)
w = 1;
tedge = [];
%
% get a vector that defines if the cell is on the top or bottom of the
% image
%
rows = b(:,2) == T |  b(:,2) == B;
edges = b(rows,:);
%
% find the edge where its pixel nieghborhood is missing a boundary
%
for i1 = 1:length(edges(:,1))
    x = edges(i1,:);
    neighbor=zeros(8,2);
    neighbor(1,:) = [x(1)-1, x(2)-1];
    neighbor(2,:) = [x(1),x(2)-1];
    neighbor(3,:) = [x(1)+1,x(2)-1];
    neighbor(4,:) = [x(1)-1,x(2)];
    neighbor(5,:) = [x(1)+1,x(2)];
    neighbor(6,:) = [x(1)-1,x(2)+1];
    neighbor(7,:) = [x(1),x(2)+1];
    neighbor(8,:) = [x(1)+1,x(2)+1];
    neighborhood(1) = sum(ismember(...
        neighbor([1,4,6],:),b,'rows'))>0;
    neighborhood(2) = sum(ismember(...
        neighbor([7,2],:),b,'rows'))>0;
    neighborhood(3) = sum(ismember(...
        neighbor([3,5,8],:),b,'rows'))>0;
    val = sum(neighborhood);
    if val == 1
        tedge(w,:) = x;
        w = w+1;
    end
end
%
% get a vector that defines if the cell is on the top or bottom of the
% image
%
rows = b(:,1) == L| b(:,1) ==  R;
edges = b(rows,:);
%
% find the edge where its pixel nieghborhood is missing a boundary
%
for i1 = 1:length(edges(:,1))
    x = edges(i1,:);
    neighbor=zeros(8,2);
    neighbor(1,:) = [x(1)-1, x(2)-1];
    neighbor(2,:) = [x(1),x(2)-1];
    neighbor(3,:) = [x(1)+1,x(2)-1];
    neighbor(4,:) = [x(1)-1,x(2)];
    neighbor(5,:) = [x(1)+1,x(2)];
    neighbor(6,:) = [x(1)-1,x(2)+1];
    neighbor(7,:) = [x(1),x(2)+1];
    neighbor(8,:) = [x(1)+1,x(2)+1];
    neighborhood(1) = sum(ismember(...
        neighbor([6,7,8],:),b,'rows'))>0;
    neighborhood(2) = sum(ismember(...
        neighbor([4,5],:),b,'rows'))>0;
    neighborhood(3) = sum(ismember(...
        neighbor([1,2,3],:),b,'rows'))>0;
    val = sum(neighborhood);
    if val == 1
        tedge(w,:) = x;
        w = w+1;
    end
end
%
% close the cell if it touches the edge
%
ns1 = [];ns2=[];
if ~isempty(tedge)&& length(tedge(:,1)) == 2 
    rows = tedge(:,1) == T | tedge(:,1) == B;
    %
    m1 = min(tedge(:,2))+1;
    m2 = max(tedge(:,2))-1;
    ns1(:,2) = m1:m2;
    %
    m1 = min(tedge(:,1))+1;
    m2 = max(tedge(:,1))-1;
    ns2(:,1) = m1:m2;
    %
    if sum(rows) == 0
        ns2(:,2) = tedge(1,2);
        ns = ns2;
    elseif sum(rows) == 2
        ns1(:,1) = tedge(1,1);
        ns = ns1;
    else
        ns1(:,1) = tedge(rows,1);
        ns2(:,2) = tedge(~rows,2);
        ns = vertcat(ns1,ns2);
    end
    %
    b = vertcat(b,ns);
    %
end
x = b(:,1);
y = b(:,2);
end
%% function: fillobj; 
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description
% fill the cell object 
%% --------------------------------------------------------------
%%
function[obj] = fillobj(obj,y,x)
s = size(obj);
x = x+1;
y = y+1;
obj(sub2ind(s,y,x)) = 1;
obj = imfill(obj);
end
%% function: getexprmark; 
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description
% determine the expression marker status of each cell; if the cell is in
% two cell membranes choose the cell that the expression marker call is
% closest to 
%% --------------------------------------------------------------
%%