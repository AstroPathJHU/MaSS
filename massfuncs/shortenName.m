%% function:shortenName; shorten the names of the inForm columns
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description
% If the names of the columns in inForm are longer than 63 char then
% matlab will only read the first 63 char of the header and the search
% strings need to be adjusted as such. This function takes in a string,
% checks if it is longer than 63 chars and, if it is, cuts the string at
% 63 chars
%% --------------------------------------------------------------
%%
function out = shortenName(n)
out = n;
if length(out) > 63
    out = extractBefore(out,64);
end
end