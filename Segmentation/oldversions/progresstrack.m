wd =  'S:\Clinical_Specimen_2';
fwpath = 'S:\flatw';
tpath =  'X:\Clinical_Specimen';

fn = dir(wd);
ii = [fn.isdir];
ii(1:2) = 0;
fd = fn(ii);
%
tmppath = [wd,'/tmp_inform_data'];
tmpfd = dir(tmppath);
tmpfd = tmpfd(3:end);
ii = [tmpfd.isdir];
tmpfd = tmpfd(ii);
ii = strcmp({tmpfd(:).name},'Project_Development');
tmpfd = tmpfd(~ii);
track = 1;
AB = cell(length(tmpfd)*2,1);
for i1 = 1:length(tmpfd)
    AB{track} = tmpfd(i1).name;
    AB{track+1} = [AB{track},'_inform_date'];
    track = track + 2;
end
%
ssname = [wd,'\Batch\samples_summary.xlsx'];
%
if ismember(ssname,{fn(:).name})
   ss = readtable(ssname);
else
    ss = array2table(zeros(0,27));
    ss.Properties.VariableNames = [{'Sample','Scan','Scan_date',...
        'im3','transfer_date','flatw','flatw_date'},AB,{'inform_all',...
        'inform_all_date','Merge_Tables','Merge_Tables_date','QC_Images',...
        'QC_ready_date','BKI04toBKI02','BKI04toBKI02_date'}];
end
%
samplenames = {fd(:).name};
ii = contains(samplenames, 'Control');
samplenames = samplenames(~ii);
ii = (contains(samplenames, 'Batch')...
    |contains(samplenames, 'tmp_inform_data')|...
    contains(samplenames, 'reject'));
samplenames = samplenames(~ii);
%
ss = [ss;[samplenames',repmat({''},length(samplenames),26)]];
%
s1 = height(ss);
ScanNum = cell(1);
transferdate = cell(1);
im3 = cell(1);
Scandate = cell(1);
flatwnum = cell(1);
flatwdate = cell(1);
infm = cell(s1,length(tmpfd));
infmd = cell(s1,length(tmpfd));
ifall =  cell(1);
ifalldate =  cell(1);
MergeTbls =  cell(1);
MergeTblsDate =  cell(1);
QCImagesdate = cell(1);
QCImages = cell(1);
%
%
for i1 = 1:height(ss)
   sname = samplenames{i1};
   Scan = dir(strcat(wd,'\',sname,'\im3\Scan*'));
   sid = {Scan(:).name};
   sid = vertcat(sid{:});
   sid = sort(sid,1,'ascend');
   ScanNum{i1} = sid(end);
   %
   Scanpath = strcat(wd,'\',sname,'\im3\Scan', num2str(ScanNum{i1}),'\');
   %
   R = getAnnotations(Scanpath,['Scan',num2str(ScanNum{i1})],sname);
   expectim3 = cellfun(@(x)erase(x,'_M2'),R,'Uni',0);
   expectim3num = length(unique(expectim3));
   %
   MSIpath = [Scanpath, 'MSI*'];
   MSIfolder = dir(MSIpath);
   transferdate{i1} = MSIfolder.date(1:11);
   %
   im3path = [Scanpath,'MSI/*.im3'];
   im3s = dir(im3path);
   actualim3num = length(im3s);
   im3{i1} = [num2str(actualim3num),'of',num2str(expectim3num)];
   %
   [~,idx] = sort(datenum({im3s(:).date},...
       'dd-mmm-yyyy hh:MM:ss'),1,'ascend');
   ii = idx == 1;
   Scandate{i1} = im3s(ii).date(1:11);
   %
   today = datetime();
   twodayago = today - 2;
   %
   if actualim3num == expectim3num||...
           twodayago >= datetime(transferdate{i1})
       flatwpath = [wd,'\',sname,'\im3\flatw'];
       if ~exist(flatwpath,'dir')
           command = ['\\bki02\c$\mcode\Im3Tools\doOneSample ',...
               wd,' ',fwpath,' ',sname];
           status = system(command);
       end
       flatw = dir([flatwpath,'\*.im3']);
       flatwnum{i1} = [num2str(length(flatw)),'of',num2str(actualim3num)];
       %
       [~,idx] = sort(datenum({flatw(:).date},...
           'dd-mmm-yyyy hh:MM:ss'),1,'ascend');
       ii = idx == 1;
       flatwdate{i1} = flatw(ii).date(1:11);
   end
   %
   informpath = [wd,'\',sname,'\inform_data'];
   %
   iffd = dir([informpath,'\Phenotyped\']);
   %
   insnum = zeros(length(tmpfd),1);
   expectedinform = zeros(length(tmpfd),1);
   iffdloc = zeros(length(tmpfd),1);
   dt = cell(length(tmpfd),1);
   trackinform = 0;
   for i2 = 1:length(tmpfd)
       tmpname = tmpfd(i2).name;
       [x,~]=ismember(tmpname,{iffd(:).name});
       %
       %if inform files do not yet exist in Clinical Specimen check the
       %tmp_data_folder if a new batch is ready to be moved for that
       %specimen
       %
       if x == 0
           %first find the inform files that are generated
           tmp2path = [tmpfd(i2).folder,'\',tmpname];
           tmp2 = dir(tmp2path);
           tmp2 = tmp2(3:end);
           ii = [tmp2.isdir];
           tmp2 = tmp2(ii);
           %loop through each file
           for i3 = 1:length(tmp2)
               numericfdspath = [tmp2(i3).folder,'\',tmp2(i3).name];
               %Find batch file to indicate that the inform is finished
               if isfile([numericfdspath,'\Batch.log'])
                   %now check the files and find out if any correspond to
                   %this case
                   cfiles = dir([numericfdspath,'\',sname,'*']);
                   if ~isempty(cfiles)
                       %trasnfer those files that are for this case
                       des1 = [informpath,'\Component_Tiffs'];
                       sor = [tmp2path,'\',sname,'*component_data.tif'];
                       [comps] = transferfls(sor,des1);
                       %
                       des = [informpath,'\Phenotyped\',tmpname];
                       sor = [tmp2path,'\',sname,'*binary_seg_maps.tif'];
                       [bin] = transferfls(sor,des);
                       %
                       sor = [tmp2path,'\',sname,'*.txt'];
                       [seg] = transferfls(sor,des);
                       %
                       sor = [tmp2path,'\Batch.log'];
                       copyfile(sor,des);
                       if comps
                           copyfile(sor,des1);
                       end
                       sor = [tmp2path,'\*.ifp'];
                       copyfile(sor,des);
                       if comps
                           copyfile(sor,des1);
                       end
                   end
               end
           end
       end 
      %
      %now check if that folder exists again and create trackers
      %
      [x,y] = ismember(tmpname,{iffd(:).name});
      iffdloc(i2) = y;
      %
      if x == 1
           trackinform = trackinform + 1;
           ins = iffd(y);
           inspath = [ins.folder,'\',ins.name];
           %
           %get number of files in Specified AB folder
           %
           insnames = dir([inspath,'\*seg_data.txt']);
           insnum(i2) = length(insnames);
           %
           %get number of files inForm had an error on to calculate
           %expected number of files
           %
           if isfile([inspath,'\Batch.log'])
               Batch = fileread([inspath,'\Batch.log']);
               InformErrors = length(strfind(Batch,[sname,'_[']))/2;
               expectedinform(i2) = actualim3num-InformErrors;
           else
               expectedinform(i2) = actualim3num;
           end
           %
           %make the number of files string 
           %
           infm{i1,i2} = [num2str(insnum(i2)),'of',num2str(expectedinform(i2))];
           %
           %find the most recent transfer date
           %
           [~,idx] = sort(datenum({insnames(:).date},...
               'dd-mmm-yyyy hh:MM:ss'),1,'ascend');
           ii = idx == 1;
           infmd{i1,i2} = insnames(ii).date(1:11);
           dt{i2} = insnames(ii).date(1:11);
      end
   end
   %
   %if trackinform is 6 all the inform data is done and move to Results
   %
   if trackinform == length(tmpfd)
       aifall = sum(insnum);
       exifall = sum(expectedinform);
       ifall{i1} = [num2str(aifall),'of',num2str(exifall)];
       [~,idx] = sort(datenum({iffd(iffdloc).date},...
           'dd-mmm-yyyy hh:MM:ss'),1,'ascend');
       ii = idx == 1;
       ifalldate{i1} = iffd(ii).date(1:11);
   end
   %
   %if Results exists:
   %
   Rfd = [informpath,'\Phenotyped\Results\Tables'];
   mergeroot = [informpath,'\Phenotyped'];
   if ~exist(Rfd,'dir') && trackinform == length(fd)
       %if directory does not exist but inform tracker ready then run merge
       %code
       addpath('\\BKI05\n$\bgcode\matlab\mergetbls\FileLoop')
       fileloop(mergeroot,Markers);
       rmpath('\\BKI05\n$\bgcode\matlab\mergetbls\FileLoop')
   end
   if exist(Rfd,'dir')
       %now if the results folder exists track
       Tablespath = [informpath,'\Phenotyped\Results\Tables*'];
       Tblsdate = dir(Tablespath);
       MergeTblsDate{i1} = Tblsdate.date(1:11);
       %
       Tablespath = [informpath,'\Phenotyped\Results\Tables\'];
       Tables = dir(Tablespath);
       Tablesnum = length(Tables);
       expectedTablesnum = min(insum);
       MergeTbls{i1} = [num2str(Tablesnum),'of',num2str(expectedTablesnum)];
   end
   QCfd = [informpath,'\Phenotyped\Results\ByImage'];
   if exist(Rfd,'dir') && exist(QCfd,'dir')
       %create QC images
        addpath('\\BKI05\n$\bgcode\matlab\mergetbls\ImageLoop')
        T = imageloop(Markers, mergeroot,sname);
        rmpath('\\BKI05\n$\bgcode\matlab\mergetbls\ImageLoop')
   end
   if exist(QCfd,'dir')
        %track QC images
        QCfd = [QCfd,'\QC_Lineage\All_Markers'];
        QCfiles = dir(QCfd);
        %
        %get date
        %
        [~,idx] = sort(datenum({QCfiles(:).date},...
           'dd-mmm-yyyy hh:MM:ss'),1,'ascend');
       ii = idx == 1;
       QCImagesdate{i1} = QCfiles(ii).date(1:11);
       %
       %get number of files
       %
       QCImages{i1} = length(QCfiles)/4;
       
   end
  % tfd = [tpath,'\',sname];
end

ss.Scan(1:length(ScanNum)) = ScanNum';
ss.transfer_date(1:length(transferdate)) = transferdate';
ss.im3(1:length(im3)) = im3';
ss.Scan_date(1:length(Scandate)) = Scandate';
ss.flatw(1:length(flatwnum)) = flatwnum';
ss.flatw_date(1:length(flatwnum)) = flatwdate';
for i1 = 1:length(tmpfd)
    dd = length(infm(:,i1));
    ss.(AB{i1})(1:dd) = infm(:,i1);
    id = [AB{i1},'_inform_date'];
    ss.(AB{i1})(1:dd) = infmd(:,i1);
end
ss.inform_all(1:length(ifall)) = ifall';
ss.inform_all_date(1:length(ifalldate)) = ifalldate';

ss.Merge_Tables(1:length(MergeTbls)) = MergeTbls';
ss.Merge_Tables_date(1:length(MergeTblsDate)) = MergeTblsDate';

ss.QC_Images(1:length(QCImages)) = QCImages';
ss.QC_ready_date(1:length(QCImagesdate)) = QCImagesdate';
%
%
writetable(ss,ssname)