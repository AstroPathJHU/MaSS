%% extractcolorvec
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2018
%% --------------------------------------------------------------
%% Description
%%% extracts the color vector from the Colors column in the merge config 
%%% file
%% --------------------------------------------------------------
%%
function[err_val, mycol] = extractcolorvec(B)
%
err_val = 0;
%
% colors
%
if height(B) <= 7
    %%blue%green%yellow%red%orange%cyan%magenta%black%%
    mycolab = [0 1 0;
        1 1 0;
        1 0 0;
        0.91 0.41 0.17;
        0 1 1;
        1 0 1;];
    mycol.all = [0 0 1;
        mycolab(1:height(B)-1, :);
        0 0 0];
    %
elseif height(B) <= 10 && height(B) > 7
    %%blue%coral%green%yellow%red%orange%cyan%magenta%white%black%%
    mycolab = [1 .7529 .7961;
        0 1 0;
        1 1 0;
        1 0 0;
        0.91 0.41 0.17;
        0 1 1;
        1 0 1;
        1 1 1;];
    mycol.all = [0 0 1;
        mycolab(1:height(B)-1, :);
        0 0 0];
    %
else
    mycol.all = [];
    err_val = 1;
    %{
    error(['Error in ImageLoop > mkimageid: \n'...
        'Need to add color values for ',n,...
        ' color panel to mkimageid function'], class(n));
    %}
end
%
color_names = {'red','green','blue','cyan', ...
    'magenta','yellow','white','black','orange','coral'};
color_names2 = {'r','g','b','c','m','y','w','k','o','l'};
colors = [eye(3); 1-eye(3); 1 1 1; 0 0 0; 1 .7529 .7961; 0.91 0.41 0.17;];
%
if isa(B.Colors, 'cell')
    [ii,loc] = ismember(B.Colors, color_names);
    [ii1,loc1] = ismember(B.Colors, color_names2);
    ii = ii + ii1;
    loc = loc + loc1;
    %
    if sum(ii) ~= length(B.Colors)
        new_colors = B.Colors(~ii);
        new_colors = replace(new_colors, {'[',']',' '},'');
        new_colors = cellfun(@(x) strsplit(x, ','), new_colors, 'Uni', 0);
        new_colors = cellfun(@str2double, new_colors, 'Uni', 0);
        if any(cellfun(@length, new_colors)~=3)
            err_val = err_val + 2;
            return
        end
        new_colors = cell2mat(new_colors);
        if any(new_colors > 255)
            err_val = err_val + 2;
            return
        end
        loc(~ii) = (length(colors) + 1):...
            ((length(colors)) + (length(B.Colors) - sum(ii)));
        colors = [colors;new_colors];
    end
    %
    mycol.all = [colors(loc,:); 0, 0, 0];
    %
else
    err_val = err_val + 2;
    return
end
%
end

