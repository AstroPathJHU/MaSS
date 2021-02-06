function [SW, EW] = resolveMultiLin(current_marker, Markers)
%
% get top color or highest numeric opal in the coexpression
%
SW = cellfun(@(x)startsWith(current_marker,x),Markers.all,'Uni',0);
SW = [SW{:}];
SW = find(SW) + 1;
%
% get bottom color or lowest numeric opal in the coexpression
%
EW = cellfun(@(x)endsWith(current_marker,x), Markers.all,'Uni',0);
EW = [EW{:}];
EW = find(EW) + 1;
%
end