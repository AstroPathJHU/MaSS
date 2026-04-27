function [use_parallel, pool, status_msg] = startrobustparpool(max_workers, task_label)
%
% Attempt to start/recover a local parallel pool.
% Falls back to serial execution if pool validation/start fails.
%
use_parallel = false;
pool = [];
status_msg = '';
%
if nargin < 1 || isempty(max_workers)
    max_workers = inf;
end
if nargin < 2
    task_label = 'task';
end
%
try
    pool = gcp('nocreate');
catch EM
    status_msg = ['Could not query existing parallel pool for ', task_label, ...
        '. Running serially. ', EM.message];
    return
end
%
if isempty(pool)
    try
        numcores = feature('numcores');
        requested_workers = min(numcores, max_workers);
        if isempty(requested_workers) || requested_workers < 1
            requested_workers = 1;
        end
        %
        pool = parpool('local', requested_workers);
    catch EM1
        %
        % Self-heal attempt for stale worker state, then retry once.
        %
        try
            pool = gcp('nocreate');
            if ~isempty(pool)
                delete(pool);
            end
            pool = parpool('local');
        catch EM2
            status_msg = ['Parallel pool unavailable for ', task_label, ...
                '. Running serially. First error: ', EM1.message, ...
                ' Retry error: ', EM2.message];
            pool = [];
            use_parallel = false;
            return
        end
    end
end
%
use_parallel = ~isempty(pool);
if use_parallel
    status_msg = ['Parallel pool ready for ', task_label, '.'];
else
    status_msg = ['Parallel pool unavailable for ', task_label, ...
        '. Running serially.'];
end
%
end
