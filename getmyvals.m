function [vals] = getmyvals(sim1,im)
%
sim2 = reshape(sim1, 1004,1344);
cells = label2idx(sim2);

vals = cellfun(@(x)sum(im(x,:),1),cells,'Uni',0);
ar = cellfun(@(x)length(x), cells,'Uni',0);
ar = [ar{:}];
%
cellids = 1:1:length(vals);
%
vals = vertcat(vals{:});
%
vals = [cellids',ar',vals];
%
end




