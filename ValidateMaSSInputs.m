%% ValidateMaSSInputs
%% --------------------------------------------------------------
%% Description
%%% Standalone preflight wrapper for MaSS/CreateImageQAQC inputs.
%%% Uses the same argument pattern:
%%%   ValidateMaSSInputs(wd, sname, MergeConfig, [logstring], [allimages])
%%% and validates MergeConfig parsing, inForm table discoverability, and
%%% component_data.tif channel compatibility before any heavy processing.
%% --------------------------------------------------------------
%%
function [err_val, report] = ValidateMaSSInputs(wd, sname, MergeConfig, logstring, allimages)
%
filepath = fileparts(mfilename('fullpath'));
addpath(genpath(filepath))
%
if nargin < 4
    logstring = '';
end
if nargin < 5
    allimages = 0; %#ok<NASGU> % accepted for CreateImageQAQC API parity
end
%
report = struct();
report.sample = sname;
report.mergeconfig = MergeConfig;
report.ok = false;
report.checks = {};
%
% parse merge config
[Markers, err_val] = createmarks(MergeConfig);
if err_val ~= 0
    report.checks{end+1} = sprintf('FAIL: createmarks err_val=%d', err_val); %#ok<AGROW>
    return
end
report.checks{end+1} = sprintf('PASS: parsed MergeConfig (%d active opal marker(s))', ...
    length(Markers.Opals)); %#ok<AGROW>
%
% locate fields/files exactly as MaSS does
[filenms, ~, err_val] = getfilenames(wd, Markers);
if err_val ~= 0
    report.checks{end+1} = sprintf('FAIL: getfilenames err_val=%d', err_val); %#ok<AGROW>
    return
end
if isempty(filenms)
    err_val = 13;
    report.checks{end+1} = 'FAIL: no matching inForm cell_seg_data files found'; %#ok<AGROW>
    return
end
report.checks{end+1} = sprintf('PASS: discovered %d image file(s) for merge', ...
    length(filenms)); %#ok<AGROW>
%
% component stack/channel preflight check
[err_val, err_msg] = validatecomponentstacks(wd, filenms, Markers);
if err_val ~= 0
    if isempty(err_msg)
        report.checks{end+1} = sprintf('FAIL: validatecomponentstacks err_val=%d', err_val); %#ok<AGROW>
    else
        report.checks{end+1} = sprintf('FAIL: %s', err_msg); %#ok<AGROW>
    end
    return
end
report.checks{end+1} = 'PASS: component_data.tif channels compatible with MergeConfig'; %#ok<AGROW>
%
% all checks passed
err_val = 0;
report.ok = true;
report.checks{end+1} = 'PASS: preflight complete (safe to run MaSS/CreateImageQAQC)'; %#ok<AGROW>
%
end
