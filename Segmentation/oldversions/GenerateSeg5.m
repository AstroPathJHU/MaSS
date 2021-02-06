wd2 = 'E:\Clinical_Specimen\';
%
fn = dir(wd2);
fd = fn(3:end);
ii = [fd.isdir];
fd = fd(ii);
samplenames = {fd(:).name};
%
ii = (contains(samplenames, 'Batch')...
    |contains(samplenames, 'tmp_inform_data')|...
    contains(samplenames, 'reject')|...
    contains(samplenames, 'Control'));
u1 = samplenames(~ii);
% 
for u2 = 1:length(u1)
    casenum = u1{u2};
    GenerateOneSampleSeg(casenum, wd2);
end
