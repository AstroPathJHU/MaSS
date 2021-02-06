function [F] = getAnnotations(root,scan,smp)
%%----------------------------------------------------
%% read XML file with grid annotations
%% Daphne Schlesinger, 2017
%% edited by Alex Szalay, 2018-0725
%% edited by Benjamin Green, 2018-1212
%%----------------------------------------------------
%% take the path for an xmlfile corresponding to a qptiff. 
%% The filename is written into a cell array, the output of a function: 
%%       F
%% 
%%------------------------------------------------------
    %
    % check what is the highest Scan number!!
    % 
    filepath = [root,smp,...
        '_', scan '_annotations.xml'];
    %
    XML = xmlread(filepath);
    annList = XML.getFirstChild;    
    ann = annList.item(5);
    %
    B = cell(1);
    track = 1;
    %
    % get rectangular annotations
    %
    for i1 = 1:2: ann.getLength - 1
        temp = ann.item(i1);
        s = temp.getAttribute('subtype');
        
       % try
       %     s = temp.getAttribute('subtype');
       % catch
       %     continue;
       % end
        if  s.equals("RectangleAnnotation")
            B{track} = temp;
            track = track + 1;
        end
    end
    %
    F = cell(1);
    track2 = 1;
    for i2 = 1 : length(B)
        %
        node = B{i2};
        history = node.item(7);
        histlastRef = history.getLength-2;
        histRef = history.item(histlastRef);
        %
        f =  histRef.item(3).getTextContent;
        t =  histRef.item(7).getTextContent;
        if strcmp(t, 'Acquired')
            f = char(f);
            F{track2} = f(5:end);
            track2 = track2 + 1;
        end
        %
    end
    %
end
