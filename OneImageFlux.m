function [infm, my] = OneImageFlux(inff, imf, simf, imnum)
%
imnum = num2str(imnum);
fils = getinformfiles(inff);
% values will be in the same format as my
infm.memvals = getinformmemvals(fils, imnum);
infm.nucvals = getinformnucvals(fils, imnum);
infm.cytovals = getinformcytovals(fils, imnum);
infm.ecvals = getinformecvals(fils, imnum);
close all
% 
im = opencomps(imf);
sim = openseg(simf);
%
% each value will be in the format [cellid, area, total intensity, mean intensity]
% im = arrayfun(@(x)round(x,3),im(:,:),'Uni',1);
my.memvals = getmyvals(sim(:,3),im);
my.nucvals = getmyvals(sim(:,1),im);
my.ecvals = getmyecvals(sim(:,3),im);
my.cytovals = getmyvals(sim(:,2),im);
%
mkcytofigs(infm, my,imnum)
mkecfigs(infm, my,imnum,ii2)
mknucfigs(infm, my,imnum)
mkmemfigs(infm, my,imnum)
%
close all
my.ecvals(:,2:end) = (my.memvals(:,2:end) + my.nucvals(:,2:end) + my.cytovals(:,2:end));
mkecfigs2(infm,my,imnum)
mkecfigs3(infm, my, imnum,ii2)
%{
for i1 = 1:9
a1(:,i1) = int32(1000.*(my.ecvals(:,i1+1)+.0005));
%
a4(:,i1) = single(a1(:,i1))./1000;
end
my.ecvals(:,2:end) = (a4);
%
a1 = my.memvals(:,2:end);
a2 = my.nucvals(:,2:end);
a3 = my.cytovals(:,2:end);
%}
%my.ecvals(:,2:end) = (a1 + a2 + a3);
%mkecfigs4(infm,my,imnum)
end