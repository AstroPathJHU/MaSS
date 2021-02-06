wd2 = 'R:\Clinical_Specimen\';
%
fn = dir(wd2);
fd = fn(3:end);
ii = [fd.isdir];
fd = fd(ii);
samplenames = {fd(:).name};
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
for u2 = 47%1:length(u1)
    casenum = u1{u2};
    wd = [wd2,casenum,'\inform_data'];
    if exist(wd,'dir')
        wd1 = [wd,'\Phenotyped'];
        paths1 = [wd1,'\PDL1\*maps.tif'];
        fnamess1 = dir(paths1);
        fds1 = {fnamess1(:).folder};
        nms1 = {fnamess1(:).name};
        %
        paths2 = [wd1,'\Tumor\*maps.tif'];
        fnamess2 = dir(paths2);
        fds2 = {fnamess2(:).folder};
        nms2 = {fnamess2(:).name};
        %
        pathc = [wd1,'\PDL1\*data.tif'];
        fnamesc = dir(pathc);
        fdc = {fnamesc(:).folder};
        nmc = {fnamesc(:).name};
        %
        pathtin = [wd1,'\Results\Tables\*table.csv'];
        fnamestin = dir(pathtin);
        fdtin = {fnamestin(:).folder};
        nmtin = {fnamestin(:).name};
        %
        pathtout =  'W:\Segmentation\Comparison\M41_1\inForm_seg_v2';
        %
        inname = regData2.micronLoc;
        outname = regData2.regLoc;
        for x = 1:length(outname)
            inname2{x} = [num2str(inname{x}(1)),',',num2str(inname{x}(2))];
            outname2{x} = [num2str(outname{x}(1)),',',num2str(outname{x}(2))];
        end
        
        %
        fnamestout = replace(nmtin,inname2,outname2);
        fnamestout = replace(fnamestout,'_cleaned_phenotype_table.csv','.csv');
        fnamesimout = replace(fnamestout,'.csv','_seg.tif');
        fnamesimout2 = replace(fnamestout,'.csv','_component_data_w_seg.tif');
        %
        regdata = extractBetween(fnamestout,'[',']');
        regdata = cellfun(@(x)strsplit(x,','),regdata,'Uni',0);
        nl= vertcat(regdata{:});
        %
        CoordTable = table(fnamestout',nl(:,1),nl(:,2),...
            'VariableNames',{'filename','ImageX','ImageY'});
        c = zeros(1004,1344,8);
        %
        parfor i1 = 1:height(CoordTable)
            
            fname = [fdtin{i1},'\',nmtin{i1}];
            f = readtable(fname,'Format',formatspec,'Delimiter',',');
            %
            T = [repmat(CoordTable(i1,:),height(f),1),f(:,{'CellID','CellNum',...
                'Phenotype','CellXPos','CellYPos'})];
            %
            imnm1 = [fds1{i1},'\',nms1{i1}];
            %
            s1 = imread(imnm1,3);
            n1 = imread(imnm1,2);
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
            imnm2 = [fds2{i1},'\',nms2{i1}];
            %
            s2 = imread(imnm2,3);
            n2 = imread(imnm2,2);
            t2 = imread(imnm2,1);
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
            iname = [pathtout,'\',fnamesimout{i1}];
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
            fname = [pathtout,'\',fnamestout{i1}];
            T = T(:,{'filename','ImageX','ImageY','CellID','Phenotype',...
                'CellXPos','CellYPos'});
            writetable(T, fname)
            %{
            iname = [pathtout,'\',fnamesimout2{i1}];
            imin = [fdc{i1},'\',nmc{i1}];
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
            %}
        end
    end
end
