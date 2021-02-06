function [infm] = getinformecvals(fils, imnum)
    %    
    Markers.all = {'PDL1','CD8','FoxP3','Tumor','PD1','CD163'};
    Markers.Opals = {'520','540','570','620','650','690'};
    %
    f1 = cellfun(@(x)['TotalEntireCell',x],Markers.Opals,'Uni',0);
    f2 = cellfun(@(x)['MeanEntireCell',x],Markers.Opals,'Uni',0);
    infm.ecvals = fils(:,[{'CellID','EntireCellArea'},f1,f2]);
    %
    for i1 = 1:length(Markers.Opals)
       a = infm.ecvals(:,f1{i1});
       a = table2array(a);
       b = infm.ecvals(:,'EntireCellArea');
       b = table2array(b);
       c = double(a)./double(b);
       %
       infm.ecvals.(['ExpectedMeanEntireCell',Markers.Opals{i1}]) = c(:);
    end
    %
    % create a figure
    %
    figname =  ['Mean EntireCell InForm Evaluation of Image ', imnum];
    Positions = {[0.10 0.65 0.25 0.25],[0.10 0.25 0.25 0.25],...
        [0.4 0.65 0.25 0.25], [0.4 0.25 0.25 0.25],...
        [0.7 0.65 0.25 0.25],[0.7 0.25 0.25 0.25]};
    %
    XX = figure('visible' , 'on', 'NumberTitle', 'off', 'Name',...
        figname);
    set(gcf, 'units','normalized','outerposition',[0 0 1 1],...
        'PaperUnits','inches');
    %
    for i1 = 1:length(Markers.Opals)
        op = Markers.Opals{i1};
        x = infm.ecvals.(['ExpectedMeanEntireCell', op]);
        y =  infm.ecvals.(['MeanEntireCell', op]);
         x(isnan(x)) = [];
        y(isnan(y)) = [];
        ii = abs(x - y) > 1;
         disp(num2str(find(ii)))
         disp(op)
         disp(num2str(x(ii)))
        disp(num2str(y(ii)))
         disp(' ')
        c = polyfit(x,y,1);
        y_est = polyval(c,x);
        axes('Position',Positions{i1});
        scatter(x,y);
        hold on
        plot(x,y_est,'r--','LineWidth',2)
        hold off
        xlabel({['Expected Mean EntireCell ', op],...
            'Calculated by: Inform ec Int / Inform Total Area'})
        ylabel(['Inform Mean EntireCell ', op])
        title(['EntireCell ', op,' Correlation Image ', imnum])
        ME = sum(y-y_est);
        MSE = sum((y-x).^2)/length(y);
        v1 = double(max(y) *(8/9));
        v2 = double(max(x) *(2/9));
        text(v2,v1,['MSE = ',num2str(MSE)])
        v1 = double(max(y)*(7/9));
        text(v2,v1,['slope = ',num2str(c(1))])
        box on
    end
     infm = infm.ecvals;
     print(XX,[figname,'.tif'],'-dtiff','-r0')
end