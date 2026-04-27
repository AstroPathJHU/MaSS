%% validatecomponentstacks
%% --------------------------------------------------------------
%% Description
%%% preflight check that component_data.tif stacks are compatible with
%%% MergeConfig-derived marker opals before merge loop starts
%% --------------------------------------------------------------
%%
function [err_val, err_msg] = validatecomponentstacks(wd, filenms, Markers)
%
err_val = 0;
err_msg = '';
%
if isempty(filenms)
    return
end
%
STANDARD_OPALS_AFTER_DAPI = [480, 520, 540, 570, 620, 650, 690, 780];
n_standard_deck = 1 + numel(STANDARD_OPALS_AFTER_DAPI); % DAPI + 8 opals
expected_cols = 1 + length(Markers.Opals); % DAPI + active mergeconfig opals
%
bad = {};
for i1 = 1:length(filenms)
    %
    % convert *_cell_seg_data.txt -> *_component_data.tif name prefix
    %
    im_base = extractBefore(filenms(i1).name, ']_cell_seg_data');
    if isempty(im_base)
        im_base = extractBefore(filenms(i1).name, ']_CELL_SEG_DATA');
    end
    if isempty(im_base)
        bad{end+1} = sprintf('%s: could not parse image id from cell seg filename', filenms(i1).name); %#ok<AGROW>
        continue
    end
    %
    iname = [wd, '\Component_Tiffs\', im_base, ']_component_data.tif'];
    if ~exist(iname, 'file')
        bad{end+1} = sprintf('%s: component file missing', iname); %#ok<AGROW>
        continue
    end
    %
    try
        props = imfinfo(iname);
    catch
        bad{end+1} = sprintf('%s: component file unreadable/corrupt', iname); %#ok<AGROW>
        continue
    end
    %
    is_gray = strcmp({props.ColorType}, 'grayscale');
    gray_idx = find(is_gray);
    if isempty(gray_idx)
        bad{end+1} = sprintf('%s: no grayscale pages found', iname); %#ok<AGROW>
        continue
    end
    %
    % Ignore thumbnails/overview pages by only counting full-size pages
    h0 = props(gray_idx(1)).Height;
    w0 = props(gray_idx(1)).Width;
    full_size = arrayfun(@(k) props(k).Height == h0 && props(k).Width == w0, gray_idx);
    n_gray = sum(full_size);
    %
    ok_exact = (n_gray == expected_cols);
    ok_full_deck = (n_gray == n_standard_deck) && ...
        all(ismember(Markers.Opals, STANDARD_OPALS_AFTER_DAPI));
    %
    if ~(ok_exact || ok_full_deck)
        bad{end+1} = sprintf(['%s: grayscale full-size pages=%d, expected=%d ', ...
            '(DAPI+active opals), mergeconfig opals=%s'], ...
            iname, n_gray, expected_cols, mat2str(Markers.Opals(:)')); %#ok<AGROW>
    end
end
%
if ~isempty(bad)
    err_val = 20;
    n_show = min(5, numel(bad));
    shown = strjoin(bad(1:n_show), ' | ');
    if numel(bad) > n_show
        shown = [shown, ' | ... +', num2str(numel(bad)-n_show), ' more'];
    end
    err_msg = ['component stack preflight failed: ', shown];
end
%
end
