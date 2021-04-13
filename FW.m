function FW()
% =======================================================================
% Shows Food-Web's indexes distributions
%
% Alejandro Rozenfeld
% alejandro.rozenfeld@gmail.com
% =======================================================================

    % use uiimport to read FW's indexes file dbaseIndexes.csv into a Table.
    % and save it into dbaseIndexes.mat file.
    
    load('dbaseIndexes.mat')

        
    %Load food-webs fitting data...
    load('FWinfo_Real.mat') % ---> DatosFWs
    
    FWcodes   = cellfun(@(fw) fw{1}, DatosFWs', 'UniformOutput', false);                % Codes
    FWnNodes  = cell2mat(cellfun(@(fw) fw{2}(2), DatosFWs', 'UniformOutput', false));   % nNodes    
    FWnLinks  = cell2mat(cellfun(@(fw) fw{2}(3), DatosFWs', 'UniformOutput', false));   % nLinks    
    FWcats    = cell2mat(cellfun(@(fw) fw{3}, DatosFWs', 'UniformOutput', false));      % RandomR, ScaleFreeR, ExponR (Cathegories)    
    
    % erase FWs with fitting errors...
    FWsWithRgt1 = find(FWcats(:,1)>1.001 | FWcats(:,2)>1.001 | FWcats(:,3)>1.001);
    FWcats(FWsWithRgt1,:) = []; 
    FWcodes(FWsWithRgt1)  = [];
    FWnNodes(FWsWithRgt1)  = [];
    FWnLinks(FWsWithRgt1)  = [];
    %--------------------------------
    
    FW = table(FWcodes,FWcats,FWnNodes,FWnLinks);
    FW.Properties.VariableNames{1}='code';
    catsANDindexes = join(FW,dbaseIndexes);  %merge fitting and indexes tables...

    
    % ======== Filter food-webs to analyze: ==================== 
    % = 'coastal' or 'terrestrial' or 'freshwater' or 'marine' =
    ixCoastal     = find(strcmp(catsANDindexes.ecosystem,'coastal'));
    ixFreshWater  = find(strcmp(catsANDindexes.ecosystem,'freshwater'));
    ixMarine      = find(strcmp(catsANDindexes.ecosystem,'marine'));
    ixTerrestrial = find(strcmp(catsANDindexes.ecosystem,'terrestrial'));
    
    cantCoastal = numel(ixCoastal);
    cantFreshWt = numel(ixFreshWater);
    cantMarine  = numel(ixMarine);
    cantTerrest = numel(ixTerrestrial);
    
    totFWs = cantCoastal + cantFreshWt + cantMarine + cantTerrest;
    disp([    'coast: '   num2str(cantCoastal/totFWs*100) ...
          '% | terrest: ' num2str(cantTerrest/totFWs*100) ...
          '% | marine: '  num2str(cantMarine /totFWs*100) ...
          '% | freshwt: ' num2str(cantFreshWt/totFWs*100) ])
    
    
    ix2del = [1:size(catsANDindexes,1)];    %all indexes
   
    % ==== Select type of Food-Web to analyze (coastal, terrestrial, marine or freshwater) ===============
    ix2del(ixCoastal)     = [];            % keep coastal     indexes
%    ix2del(ixTerrestrial) = [];           % keep terrestrial indexes
%    ix2del(ixMarine)      = [];           % keep marine      indexes     
%    ix2del(ixFreshWater)  = [];           % keep freshwater  indexes
    % =======================================
    
    catsANDindexes(ix2del,:)=[];            %erase all but selected entries
   
  
    D2PG  = arrayfun(@(fw) norm(catsANDindexes{fw,2}(1:2)-[1,0]),[1:size(catsANDindexes,1)]'); %Distance 2 Pure Gaussian
    D2SF  = arrayfun(@(fw) norm(catsANDindexes{fw,2}(1:2)-[0,1]),[1:size(catsANDindexes,1)]'); %Distance 2 Pure ScaleFree
    D2PU  = arrayfun(@(fw) norm(catsANDindexes{fw,2}(1:2)-[0,0]),[1:size(catsANDindexes,1)]'); %Distance 2 Pure Uniform
    D2X   = [D2PG D2SF D2PU];
    Lats  = catsANDindexes{:,7};
    D2X = [D2X Lats];

    
    catsANDindexes.D2PG = D2PG;
    catsANDindexes.D2SF = D2SF;
    catsANDindexes.D2PU = D2PU;
    catsANDindexes.nFW = zeros(numel(D2PG),1);  %Aux. col. to calc FWs in bin
    
    colNames = catsANDindexes.Properties.VariableNames;
    colNames{7} = ['| ' colNames{7} ' |'];        
    
    
    
    
    hFigCorrel = figure('name','Correlations');
    hFigCorrel.WindowStyle = 'docked';
    
    
    pumX = uicontrol(hFigCorrel,'style','popupmenu','units','normalized',...
        'position',[0.75 0 0.2 0.05],...
        'string',{'D2PG', 'D2SF', 'D2PU', 'LAT'});
    

    ixYinfo = [3:4 6:7 9:21 24:26 27];
    pumY = uicontrol(hFigCorrel,'style','popupmenu','units','normalized',...
        'position',[0.05 0.95 0.2 0.05],...
        'string',colNames(ixYinfo));
        

    Bins = [5:100]';
    cellBins = arrayfun(@(num) num2str(num),Bins,'uniformoutput',false);
    pumBins = uicontrol(hFigCorrel,'style','popupmenu','units','normalized',...
        'position',[0.9 0.9 0.1 0.05],...
        'string',cellBins);
    
    pumX.Callback    = {@changePlot,hFigCorrel,pumX,pumY,pumBins,D2X,catsANDindexes(:,ixYinfo)};
    pumY.Callback    = {@changePlot,hFigCorrel,pumX,pumY,pumBins,D2X,catsANDindexes(:,ixYinfo)};
    pumBins.Callback = {@changePlot,hFigCorrel,pumX,pumY,pumBins,D2X,catsANDindexes(:,ixYinfo)};
    
    set(hFigCorrel,'WindowKeyPressFcn',{@KeyPressOnFig,pumX,pumY,pumBins,D2X,catsANDindexes(:,ixYinfo)});
    
    plot(D2X(:,pumX.Value),catsANDindexes{:,ixYinfo(pumY.Value)},'.r')
    
end

function changePlot(hObject, eventdata, hFigCorrel,pumX,pumY,pumBins,xPosibles,yPosibles)    
    fontSize = 25;

    X = xPosibles(:,pumX.Value);
    Y = yPosibles{:,pumY.Value};
    
    if pumY.Value == 4  %In case Y is Lat
        Y = abs(Y);    
    end
    
    xlims = [0 sqrt(2)];
    if pumX.Value == 4 % X = Lat
        xlims = [min(X) max(X)];
    end
    ylims = [min(Y) max(Y)];
    
    
    
    cantBins = str2num(pumBins.String{pumBins.Value});
    
    minX = min(X);
    maxX = max(X);    
    
    binStep = (maxX - minX) / cantBins;
    averYinBin=[];
     varYinBin=[];
    max_YinBin=[];
    min_YinBin=[];
    averXinBin=[];

    
    %Now Binning...
    for i = 1:cantBins,        
        minXinBin = minX + (i-1) * binStep;
        maxXinBin = minX + i     * binStep;
        ixInBin=find( minXinBin <= X & X < maxXinBin);
        if isempty(ixInBin),
            averXinBin(i) = (minXinBin + maxXinBin) / 2;
            averYinBin(i) = nan;
        else
            averXinBin(i) = nanmean(X(ixInBin));

            averYinBin(i) = nanmean(Y(ixInBin));
             varYinBin(i) =  nanvar(Y(ixInBin));
            max_YinBin(i) =     max(Y(ixInBin));
            min_YinBin(i) =     min(Y(ixInBin));
            
            if pumY.Value == 21  %If Y == nFWs selected
                averYinBin(i) = numel(ixInBin);  
            end
        end
    end
    
    
    ixNan = find(isnan(averYinBin));
    averXinBin(ixNan) = [];
    averYinBin(ixNan) = [];
     varYinBin(ixNan) = [];
    max_YinBin(ixNan) = [];
    min_YinBin(ixNan) = [];

    hold off;
    scatter(averXinBin, averYinBin,'sk','filled');
    errorbar(averXinBin, averYinBin,sqrt(varYinBin)/numel(varYinBin),'sk','color','r','markerfacecolor','k','markeredgecolor','k')
    xlim(xlims)

    rmfit;
    linearFit = ezfit(averXinBin,averYinBin,'a*x+b');        
    gaussFit  = ezfit(averXinBin,averYinBin,'gauss'); %'g+a*exp(-((x-m)^2)/(2*s^2))');  
    
%     %EXTRA FIT  ===========
    sqrFit = ezfit(averXinBin,averYinBin,'poly2');   %FIG D2PU
    showfit(sqrFit ,'fitcolor','r','fitlinewidth',1) %FIG D2PU
    linFit = ezfit(averXinBin(2:end),averYinBin(2:end),'a*x+b'); %FIG D2SF
    showfit(linFit ,'fitcolor','r','fitlinewidth',1)             %FIG D2SF
%     %======================
    hola = 1;
    

    rThr = 0.25;
    rThresh = 0.75 * max(linearFit.r,gaussFit.r);
    if rThresh >= rThr,
        cantFits = 0;
        if linearFit.r >= rThresh, 
           showfit(linearFit ,'fitcolor','r','fitlinewidth',1)
           cantFits = cantFits + 1;
        end
        if gaussFit.r  >= rThresh,
            showfit(gaussFit ,'fitcolor','b','fitlinewidth',1)
            cantFits = cantFits + 1;
        end

    end


    
    axH=gca;
    axH.YGrid='on';
    legend({'aver on bin', 'best fit'});
    axH.FontSize = fontSize;
    xlabel(pumX.String{pumX.Value},'fontsize',fontSize * 1.5);
    ylabel(pumY.String{pumY.Value},'interpreter','none','fontsize', fontSize * 1.5);
    hola = 1;
    
    
    
end


function [linearR gaussR extraR] = calcR(X,Y,cantBins)        
    
    minX = min(X);
    maxX = max(X);    
    
    binStep = (maxX - minX) / cantBins;
    averYinBin=[];
     varYinBin=[];
    max_YinBin=[];
    min_YinBin=[];
    averXinBin=[];
    
    %Ahora Binneo...
    for i = 1:cantBins,        
        minXinBin = minX + (i-1) * binStep;
        maxXinBin = minX + i     * binStep;
        ixInBin=find( minXinBin <= X & X < maxXinBin);
        if isempty(ixInBin),
            averXinBin(i) = (minXinBin + maxXinBin) / 2;
            averYinBin(i) = nan;
             varYinBin(i) = nan;
            max_YinBin(i) = nan;
            min_YinBin(i) = nan;                
        else
            averXinBin(i) = nanmean(X(ixInBin));

            averYinBin(i) = nanmean(Y(ixInBin));
             varYinBin(i) =  nanvar(Y(ixInBin));
            max_YinBin(i) =     max(Y(ixInBin));
            min_YinBin(i) =     min(Y(ixInBin));                
        end
    end

    linearFit = ezfit(averXinBin,averYinBin,'a*x+b');        
    gaussFit  = ezfit(averXinBin,averYinBin,'gauss'); %'g+a*exp(-((x-m)^2)/(2*s^2))');    
%     extraFit  = ezfit(averXinBin,averYinBin,'poly2');  %FIG RvsBin D2PU
    extraFit  = ezfit(averXinBin(2:end),averYinBin(2:end),'a*x+b');  %FIG RvsBin D2SF
    linearR = linearFit.r;
    gaussR  = gaussFit.r;
    extraR  = extraFit.r;
end



function KeyPressOnFig(obj,ev,pumX,pumY,pumBins,xPosibles,yPosibles)
    
    
    if ismember(ev.Character,{'b' 'B'}) % Showing linear.R and gauss.R vs cantBins 
        pumBinsValueAct = pumBins.Value;
        pumBins_cantItems = numel(pumBins.String);
        X = xPosibles(:,pumX.Value);
        Y = yPosibles{:,pumY.Value};
        if pumY.Value == 4  %if Y is Lat
            Y = abs(Y);
        end
        
        linearRs = [];
        gaussRs  = [];
        extraRs  = []; %FIG RvsBin D2PU
        cantsBins= [];
        for i = 1:pumBins_cantItems,   
            i
            cantBins = str2num(pumBins.String{i});
            [linearR gaussR extraR] = calcR(X,Y,cantBins);
            linearRs(end+1)  = linearR;
            gaussRs(end+1)   = gaussR;
            extraRs(end+1)   = extraR;
            cantsBins(end+1) = cantBins;
        end
        RsVScantBinsFIG = figure('name','RsVScantBins');
        RsVScantBinsFIG.WindowStyle = 'docked';
        if isempty(extraRs)
            plot(cantsBins,linearRs,'-sr',cantsBins,gaussRs,'-vb')
            legend({'linear' 'gauss'})
        else
            plot(cantsBins,linearRs,'-sr',cantsBins,gaussRs,'-vb',cantsBins,extraRs,'-vk') %FIG RvsBin D2PU
            legend({'linear' 'gauss' 'extraR'})
        end
        ca = gca;
        ca.XGrid='on';
        ca.XMinorGrid='on';
        ca.XMinorTick='on';
        
        
        xlabel('cant Bins')
        ylabel('R')
        
    elseif ismember(ev.Character,{'s' 'S'}),  %save .fig and .eps
        fileName = [pumY.String{pumY.Value} '_VS_' pumX.String{pumX.Value} '_' pumBins.String{pumBins.Value} 'bins'];
        savefig(obj,fileName)
        saveas (obj,fileName,'epsc')
        
    end
        
end