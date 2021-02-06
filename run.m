function [] = run()
wd = '\\bki04\e$\Clinical_Specimen';
sample = 'M25_1';
path = [wd,'\',sample,'\inform_data'];
%
% fils -->> a data structure 
% fils.inffn -->> inform file names
% fils.imfn -->> component image file names
% fils.simfn -->> segmentation file names
%
fils = getallfils(path);
%
inffn = fils.inffn;
imfn = fils.imfn;
simfn = fils.simfn;
%
parfor i1 = 1:length(fils.inffn)
   inff = inffn{i1};
   imf = imfn{i1};
   simf = simfn{i1};
   vals = OneImageFlux(inff, imf, simf, i1); 
end
%
end