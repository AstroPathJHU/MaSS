%% function: fillobj; 
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description
%%% fill the cell object 
%% --------------------------------------------------------------
%%
function[obj] = fillobj(obj,y,x)
s = size(obj);
x = x+1;
y = y+1;
obj(sub2ind(s,y,x)) = 1;
obj = imfill(obj);
end