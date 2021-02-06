function mkcytofigs(infm, my, imnum)
%
Markers.Opals = {'520','540','570','620','650','690'};
%
    % create a figure
    %
    figname =  ['Mean Cytoplasm InForm vs. Mask Evaluation of Image ', imnum];
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
        x = my.cytovals(:,i1 + 3);
        ar = my.cytovals(:,2);
        x = x./ar;
        x = int32(x.*1000);
        x = single(x)./1000;
        y =  infm.cytovals.(['MeanCytoplasm', op]);
        x(isnan(y)) = [];
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
        xlabel({['Expected Mean Cytoplasm ', op],...
            'Calculated by: Intensity Under Mask'})
        ylabel(['Inform Mean Cytoplasm ', op])
        title(['Cytoplasm ', op,' Correlation Image ', imnum])
        ME = sum(y-y_est);
        y = int32(y.*1000);
        y = single(y)./1000;
        MSE = sum((y-x).^2)/length(y);
        v1 = double(max(y) *(8/9));
        v2 = double(max(x) *(2/9));
        text(v2,v1,['MSE = ',num2str(MSE)])
        v1 = double(max(y)*(7/9));
        text(v2,v1,['slope = ',num2str(c(1))])
        box on
    end
    print(XX,[figname,'.tif'],'-dtiff','-r0')
end
