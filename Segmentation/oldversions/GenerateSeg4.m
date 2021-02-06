wd2 = 'E:\Clinical_Specimen\';
%
fn = dir(wd2);
fd = fn(3:end);
ii = [fd.isdir];
fd = fd(ii);
samplenames = {fd(:).name};
% ask for Control names seperate to keep as a control folder
%ii = contains(samplenames, 'Control');
%controlnames = samplenames(ii);
%samplenames = samplenames(ii);
ii = (contains(samplenames, 'Batch')...
    |contains(samplenames, 'tmp_inform_data')|...
    contains(samplenames, 'reject')|...
    contains(samplenames, 'Control'));
u1 = samplenames(~ii);
% 
formatspec = strcat({' %d16'},{' %s '},repmat('%d16 ',[1,3]),{ ' %s '},...
    repmat('%f32 ',[1,52]),{' '},repmat('%f32 ',[1,4]));
formatspec = formatspec{1};
%
%
ds.ImageLength = 1004;
ds.ImageWidth = 1344;
ds.Photometric = Tiff.Photometric.MinIsBlack;
ds.BitsPerSample   = 32;
ds.SamplesPerPixel = 1;
ds.SampleFormat = Tiff.SampleFormat.IEEEFP;
ds.RowsPerStrip    = 63;
ds.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
ds.Software = 'MATLAB';
%
for u2 = 47%47:length(u1)
    casenum = u1{u2};
    wd = [wd2,casenum,'\inform_data'];
    if exist(wd,'dir')
        wd1 = [wd,'\Phenotyped'];
        path.s1 = [wd1,'\PDL1\*maps.tif'];
        fnames.s1 = dir(path.s1);
        %
        path.s2 = [wd1,'\Tumor\*maps.tif'];
        fnames.s2 = dir(path.s2);
        %
        path.c = [wd,'\Component_Tiffs\*data.tif'];
        fnames.c = dir(path.c);
        %
        path.tin = [wd1,'\Results\Tables\*table.csv'];
        fnames.tin = dir(path.tin);
        %
        path.tout =  [wd,'\Component_Tiffs'];
        %
        fnames.tout = {fnames.tin(:).name};
        fnames.tout = replace(fnames.tout,'_cleaned_phenotype_table.csv','.csv');
        fnames.imout = replace(fnames.tout,'.csv','_seg.tif');
        fnames.imout2 = replace(fnames.tout,'.csv','_component_data_w_seg.tif');
        %
        regdata = extractBetween(fnames.tout,'[',']');
        regdata = cellfun(@(x)strsplit(x,','),regdata,'Uni',0);
        nl= vertcat(regdata{:});
        %
        CoordTable = table(fnames.tout',nl(:,1),nl(:,2),'VariableNames',{'filename','ImageX','ImageY'});
        c = zeros(1004,1344,8);
        %
        for i1 = 1:height(CoordTable)
            
            fname = strcat(fnames.tin(i1).folder,'\',fnames.tin(i1).name);
            f = readtable(fname,'Format',formatspec,'Delimiter',',');
            %
            T = [repmat(CoordTable(i1,:),height(f),1),f(:,{'CellID','CellNum',...
                'Phenotype','CellXPos','CellYPos'})];
            %
            s1 = imread(strcat(fnames.s1(i1).folder,'\',fnames.s1(i1).name),3);
            n1 = imread(strcat(fnames.s1(i1).folder,'\',fnames.s1(i1).name),2);
            n1 = reshape(n1,[],1);
            s1 = reshape(s1,[],1);
            %
            ii = table2array(T(~strcmp(T.Phenotype,'Tumor'),'CellNum'));
            ii2 = table2array(T(~strcmp(T.Phenotype,'Tumor'),'CellID'));
            %
            ii = double(ii);
            ii2 = double(ii2);
            %
            s1(~ismember(s1,ii)) = 0;
            n1(~ismember(n1,ii)) = 0;
            %ismember(ii,s1)
            %
            s3 = zeros(length(s1),1);
            n3 = zeros(length(n1),1);
            %            
            for i2 = 1:length(ii)
                v = ii(i2);
                s3(ismember(s1,v)) = ii2(i2);
                n3(ismember(n1,v)) = ii2(i2);
            end
            %
            %ismember(ii2,s3)
            s2 = imread(strcat(fnames.s2(i1).folder,'\',fnames.s2(i1).name),3);
            n2 = imread(strcat(fnames.s2(i1).folder,'\',fnames.s2(i1).name),2);
            t2 = imread(strcat(fnames.s2(i1).folder,'\',fnames.s2(i1).name),1);
            s2 = reshape(s2,[],1);
            n2 = reshape(n2,[],1);
            %
            ii = table2array(T(strcmp(T.Phenotype,'Tumor'),'CellNum'));
            ii2 = table2array(T(strcmp(T.Phenotype,'Tumor'),'CellID'));
            %
            ii = double(ii);
            ii2 = double(ii2);
            %
            s2(~ismember(s2,ii)) = 0;
            n2(~ismember(n2,ii)) = 0;
            %ismember(ii,s2)
            %
            s4 = zeros(length(s2),1);
            n4 = zeros(length(n2),1);
            %
            for i2 = 1:length(ii)
                v = ii(i2);
                s4(ismember(s2,v)) = ii2(i2);
                n4(ismember(n2,v)) = ii2(i2);
            end
            %
            %ismember(ii2,s4)
            %
            s5 = [single(s3),single(s4)];
            n5 = [single(n3),single(n4)];
            t2 = single(t2);
            %
            s5 = reshape(s5,1004,1344,2);
            n5 = reshape(n5,1004,1344,2);
            %imshow(s5(:,:,1))
            iname = strcat(path.tout,'\',fnames.imout{i1});
            %
            ii = Tiff(iname,'w');
            ii.setTag(ds);
            write(ii,t2(:,:,1));
            writeDirectory(ii)
            ii.setTag(ds);
            write(ii,n5(:,:,1));
            writeDirectory(ii)
            ii.setTag(ds);
            write(ii,n5(:,:,2));
            writeDirectory(ii)
            ii.setTag(ds);
            write(ii,s5(:,:,1));
            writeDirectory(ii)
            ii.setTag(ds);
            write(ii,s5(:,:,2));
            close (ii)
            %
            fname = strcat(path.tout,'\',fnames.tout{i1});
            T = T(:,{'filename','ImageX','ImageY','CellID','Phenotype',...
                'CellXPos','CellYPos'});
            %writetable(T, fname)
            %
            iname = strcat(path.tout,'\',fnames.imout2{i1});
            imin = strcat(fnames.c(i1).folder,'\',fnames.c(i1).name);
            c = readim(imin);
            %
            ii = Tiff(iname,'w');
            for i3 = 1:8
                d = c(:,:,i3);
                ii.setTag(ds);
                write(ii,d);
                writeDirectory(ii)
            end
            ii.setTag(ds);
            write(ii,t2(:,:,1));
            writeDirectory(ii)
            ii.setTag(ds)
            write(ii,n5(:,:,1));
            writeDirectory(ii)
            ii.setTag(ds)
            write(ii,n5(:,:,2));
            writeDirectory(ii)
            ii.setTag(ds)
            write(ii,s5(:,:,1));
            writeDirectory(ii)
            ii.setTag(ds)
            write(ii,s5(:,:,2));
            close(ii)
        end
    end
end
