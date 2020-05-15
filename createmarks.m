function[Markers] = createmarks(MergeConfig)
%
warning('off','MATLAB:table:ModifiedAndSavedVarnames')
%
try
    BIDtbl = readtable(MergeConfig);
    B = BIDtbl(:,{'Opal','Target',...
    'TargetType','CoexpressionStatus','SegmentationStatus',...
    'SegmentationHierarchy', 'ImageQA_QC', 'NumberofSegmentations'});
catch
    return
end
%
% remove DAPI row
%
dr = strcmp(B.Opal, 'DAPI');
B(dr,:) = [];
%
% start setting up Markers struct
%
Markers.Opals = B.Opal;
Markers.all = B.Target;
%
ii = strcmp('Tumor',B.ImageQA_QC);
Markers.all_original = Markers.all;
%
% change Tumor marker designation to 'Tumor'
%
if sum(ii) == 1
    Markers.all(ii) = {'Tumor'};
    Markers.Tumor{1} = 'Tumor';
elseif sum(ii) > 1
    Markers.Tumor =  ' MaSS must have only one "Tumor" designation';
    return
else
     Markers.Tumor{1} = '';
end
%
% get lineage and expression markers
%
LT = strcmp(B.TargetType,'Lineage');
Markers.lin = Markers.all(LT);
%
ET = strcmp(B.TargetType,'Expression');
Markers.expr = Markers.all(ET);
%
% get the markers with multiple segmentations, this will only be a
% capability on expression markers
%
nsegs = B.NumberofSegmentations;
if iscell(nsegs)
    nsegs = cellfun(@(x) str2double(x), nsegs, 'Uni',0);
    nsegs = cell2mat(nsegs);
end
if find(nsegs(~ET) > 1)
     Markers.nsegs =  [' MaSS can only handle expression markers',...
         ' with multiple segmentations'];
    return
end
Markers.nsegs = nsegs;
%
% Set up segmentation status to define number of segmentations and which is
% the primary segmentation
%
SS = B.SegmentationStatus;
if iscell(SS)
    SS = cell2mat(SS);
    SS = str2num(SS);
end
Markers.SegStatus = SS;
%
ii = nsegs == 1 & ~ismember(Markers.all,Markers.expr);
SS = SS(ii);
mn = Markers.all(ii);
%
% get number of different segmentations, remove markers with multiple
% segmentations from the contention
%
[~,y,~] = unique(SS);
ii = y(1);
%
Markers.seg = mn(ii);
%
Markers.altseg = cell(length(y)-1,1);
for i1 = 2:length(y)
    ii = y(i1);
    Markers.altseg(i1-1) = mn(ii);
end
%
% get coexpression status for lineage markers
%
CS = B.CoexpressionStatus(LT);
ii = ~strcmp(CS,'NA');
CS = CS(ii);
%
% track the corresponding target
%
TCS = Markers.lin(ii);
%
% get segmentation heirarchy
%
SH = B.SegmentationHierarchy;
if ~iscell(SH)
    SH = num2str(SH);
    SH = cellstr(SH);
end
Markers.SegHie = SH(LT);
%
% CS that is not NA in lineage markers; find which coexpressions are
% acceptable
%
Markers.add = [];
sego = [];
for i1 = 1:length(CS)
    %
    % get current target and opal
    %
    T = TCS{i1};
    ii = strcmp(T,Markers.all);
    o = Markers.Opals(ii);
    o = o{1};
    %
    % check them against rest of targets in coexpression
    %
    CStest = CS(~strcmp(TCS,T));
    TCStest = TCS(~strcmp(TCS,T));
    %
    for i2 = 1:length(CStest)
        o1 = CStest{i2};
        T1 = TCStest{i2};
        %
        % if the current target matches one of the targets in the rest 
        %
        if contains(o1,o)
            %
            % if the Markers.add is not empty; are both markers already
            % contained together
            %
            if ~isempty(Markers.add) && sum(contains(Markers.add,T)...
                & contains(Markers.add,T1))
                continue
            else
                track = length(Markers.add) + 1;
                Markers.add{track} = [T1,T];
                ii = strcmp(T1, Markers.lin);
                seg1 = Markers.SegHie(ii);
                ii = strcmp(T, Markers.lin);
                seg2 = Markers.SegHie(ii);
                %
                seg = max([str2num(seg1{1}),str2num(seg2{1})]);
                sego{track} = num2str(seg);
            end
        end
    end
end
%
Markers.SegHie = [Markers.SegHie;sego'];
%
% get coexpression status for expression markers
%
CS = B.CoexpressionStatus(ET);
for i1 = 1:length(CS)
    T = CS{i1};
    T = reshape(T,3,[])';
    [s,~] = size(T);
    x = arrayfun(@(x)contains(Markers.Opals,T(x,:)),1:s,'Uni',0);
    x = horzcat(x{:});
    Markers.Coex{i1} = sum(x,2);
end
%
%
%
Markers.Opals = cell2mat(Markers.Opals);
Markers.Opals = str2num(Markers.Opals);
%
Markers.all = Markers.all';
Markers.all_original = Markers.all_original';
Markers.lin = Markers.lin';
Markers.expr = Markers.expr';
Markers.nsegs = Markers.nsegs';
Markers.seg = Markers.seg';
Markers.altseg = Markers.altseg';
Markers.SegHie = Markers.SegHie';
end