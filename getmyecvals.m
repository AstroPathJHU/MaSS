function [vals] = getmyecvals(sim1,im)
%
sim2 = reshape(sim1, 1004,1344);
cells = label2idx(sim2);
s = {[1004,1344]};
s2 = length(cells);
p.size.T = 2; p.size.B = s{1}(1)-1; p.size.L = 2; p.size.R = s{1}(2)-1;
%
[obj,~] = OrganizeCells(p,s,s2,cells);
%
vals = cellfun(@(x)sum(im(x,:),1),obj,'Uni',0);
ar = cellfun(@(x)length(x), obj,'Uni',0);
ar = [ar{:}];
%
cellids = 1:1:length(vals);
%
vals = vertcat(vals{:});
%
vals = [cellids',ar',vals];
%
end