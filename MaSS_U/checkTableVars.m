%% function: checkTableVars
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 07/07/2020
%% --------------------------------------------------------------
%% Description
% checks that the table input is correct 
%% --------------------------------------------------------------
%%
function [B, err_val] = checkTableVars(B)
%%
% check the table variables to be sure they are in the correct format for
% the code. If they are not convert them.
%%
%
err_val = 0;
%
% check the data type for Opal column
%
if isa(B.Opal,'double')
   %
   % if B.Opal is a 'double' convert to a string 
   %
   tmpopal = num2cell(B.Opal);
   tmpopal = cellfun(@(x) num2str(x), tmpopal, 'Uni', 0);
   ii = strcmp(tmpopal, 'NaN');
   %
   if sum(ii) > 1
      err_val = 2;
      ii = find(ii,1);
   end
   %
   tmpopal(ii) = {'DAPI'};
   ss = size(tmpopal);
   if ss(1) == 1
       B.Opal = tmpopal';
   else
       B.Opal = tmpopal;
   end
end
%
if ~isa(B.Opal, 'cell')
  err_val = 3;
  return
end
%
% check the data type for the coexpression status column
%
if isa(B.CoexpressionStatus,'double')
   %
   % if B.Opal is a 'double' convert to a string 
   %
   tmpCS = num2cell(B.CoexpressionStatus);
   tmpCS = cellfun(@(x) num2str(x), tmpCS, 'Uni', 0);
   %
   for i1 = 1:length(tmpCS)
       tmpCS_n = tmpCS{i1};
       if length(tmpCS_n) > 3
           ii = 3:3:length(tmpCS_n) - 1;
           t(1:length(tmpCS_n)) = char(0);
           t(ii) = ',';
           tmpCS_n = [tmpCS_n;t];
           tmpCS_n = reshape(tmpCS_n(tmpCS_n ~= 0),1,[]);
           tmpCS{i1} = tmpCS_n;
       end
   end
   %
   B.CoexpressionStatus = tmpCS;
   %
end
%
B.CoexpressionStatus = cellfun(@(x) replace(x, ',',''),...
      B.CoexpressionStatus, 'Uni',0);
%
if ~isa(B.Opal, 'cell')
    err_val = 4;
end
%
% remove the DAPI row
%
dr = strcmp(B.Opal, 'DAPI');
if sum(dr) ~= 1
    err_val = 5;
end
B(dr,:) = [];
%
% check the last 3 columns are all set as numeric
%
SS = B.SegmentationStatus;
if iscell(SS)
    %SS = cell2mat(SS);
    B.SegmentationStatus = str2double(SS);
end
%
SH = B.SegmentationHierarchy;
if ~iscell(SH)
    SH = num2str(SH);
    SH = cellstr(SH);
    B.SegmentationHierarchy = SH;
end
%
SS = B.NumberofSegmentations;
if iscell(SS)
    %SS = cell2mat(SS);
    B.NumberofSegmentations = str2double(SS);
end
%
if ~iscell(B.Compartment)
    err_val = 8;
    return
end
%
end