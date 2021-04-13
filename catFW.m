function catFWs()    
% =======================================================================
% Categorize food-webs
%
% Alejandro Rozenfeld
% alejandro.rozenfeld@gmail.com
% =======================================================================

    Thr=1;
    pathFiles = 'FoodWebs_CSV_Path';
    if exist('FWinfo_Real.mat','file')==0, %Create FWinfo_Real.mat if it doesn't exist
        files=dir([pathFiles 'fw*.csv']);
        cantFWs=numel(files);

        Rs=[]; %fwId gaussR powerR size
        DatosFWs={};
        for ixFW=1:cantFWs,   
            ixFW
            if ismember(ixFW,[234,249]), %Skip corrupted files
                continue
            end
            fname=files(ixFW).name;
            codIni = find(fname=='_',1) + 1;
            codFin = find(fname=='.',1) - 1;
            fwCode = fname(codIni:codFin);
            M=loadFW([pathFiles fname]);
            fw=logical(M+M'); %build symetric connectivity matrix (incoming and outgoing links)
            lambdas=sum(fw);
            lambdaMax=max(lambdas);
            nLinks = 0.5 * sum(lambdas);

            if lambdaMax < Thr, %skip food-webs with max connected sites lower than Thr links
                continue
            end

    
            hfig=figure('units','normalized','windowstyle','docked');

            hsp1=subplot(4,3,[1 9]);
            
            S=size(M,1); %num of Species in FW.
            histEdges=[0.5:1:S+1];
            h=histogram(lambdas,histEdges); 
            xlim([0 lambdaMax])

            nonZeroFreqs=find(h.Values>0);

            x=(h.BinEdges(nonZeroFreqs)+h.BinEdges(nonZeroFreqs+1))/2;
            y=h.Values(nonZeroFreqs);



            FileName=files(ixFW).name;

            fwCode=FileName(find(FileName=='_')+1 : find(FileName=='.')-1);

            posGuion=find(fwCode=='_');
            fwCode_ = [ fwCode(1:posGuion-1) '\_' fwCode(posGuion+1:end) ];

            title(['FW ' fwCode_ ' | S: ' num2str(S) ' | L: ' num2str(nLinks)  ' | \lambda_{max}: ' num2str(lambdaMax)])

            ylim([0 max(h.Values)+1])
            gaussFit=ezfit(x,y,...
                ['g+a*exp(-((x-m)^2)/(2*s^2));lin' ...
                 ';g=' num2str(min(h.Values)) ...
                 ';a=' num2str(max(h.Values))... 
                 ';m=' num2str(mean(h.Values))...
                 ';s=' num2str(sqrt(var(lambdas)))...
                 ]); 
                          
            gaussR = gaussFit.r;
            a=gaussFit.m(1);
            Sigma=abs(gaussFit.m(4));
            m=gaussFit.m(3);
            maxLambdas = lambdaMax;

            if (maxLambdas < m + Sigma) || m + Sigma/3 < 0 || h.Values(find(h.Values>0,1,'first')) == max(h.Values), %Looks like Uniform...
                gaussR=0;            
            end         
            
            

%             if gaussFit.r>0.15,            
                showfit(gaussFit,'fitcolor','b','fitlinewidth',1);
                hPane=findobj(allchild(hfig),'type','annotationpane'); 
                paneChilds=allchild(hPane);
                hGaussFitText=paneChilds{1}; 
                panePos=get(hGaussFitText,'position');
                
                hsp2=subplot(4,3,10);
                hsp2.Visible = 'off';
                sp2Pos=get(hsp2,'position');

                panePos(1)=sp2Pos(1);
                panePos(2)=sp2Pos(2)+sp2Pos(4)-panePos(4);
                set(hGaussFitText,'position',panePos); 
                set(hGaussFitText,'backgroundcolor','none', 'linestyle', 'none')
%             end
            hT=text(0,0,['gaussR: ' num2str(gaussR)]);

            subplot(hsp1)
            powerFit=ezfit(x,y,'power;log');
            
            powerR = powerFit.r;
            if powerFit.m(2) > 0  % if exponent > 0
                powerR = 0;
            end                        
            
%             if powerFit.r>0.15,            
                showfit(powerFit ,'fitcolor','r','fitlinewidth',1)
                hPane=findobj(allchild(hfig),'type','annotationpane'); 
                paneChilds=allchild(hPane);
                hPowerFitText=paneChilds{1}(1);
                panePos=get(hPowerFitText,'position');
                
                hsp3=subplot(4,3,11);
                hsp3.Visible = 'off';
                sp3Pos=get(hsp3,'position');

                panePos(1)=sp3Pos(1);
                panePos(2)=sp3Pos(2)+sp3Pos(4)-panePos(4);
                set(hPowerFitText,'position',panePos); 
                set(hPowerFitText,'backgroundcolor','none', 'linestyle', 'none')
%             end

            hT=text(0,0,['powerR: ' num2str(powerR)]);

            subplot(hsp1)
            expFit=ezfit(x,y,'exp;log');
            
            exponR = expFit.r;  
            if expFit.m(2) > 0 % Increasing Exponential
                exponR = 0;
            end
            
%             if expFit.r>0.15,            
                showfit(expFit,'fitcolor','k','fitlinewidth',1)
                hPane=findobj(allchild(hfig),'type','annotationpane'); 
                paneChilds=allchild(hPane);
                hExpFitText=paneChilds{1}(1);                  
                panePos=get(hExpFitText,'position');
                
                hsp4=subplot(4,3,12);
                hsp4.Visible='off';
                sp4Pos=get(hsp4,'position');

                panePos(1)=sp4Pos(1);
                panePos(2)=sp4Pos(2)+sp4Pos(4)-panePos(4);
                set(hExpFitText,'position',panePos);    
                set(hExpFitText,'backgroundcolor','none', 'linestyle', 'none')
%             end
            text(0,0,['expR: ' num2str(exponR)])

            nNodes=size(fw,1);
                                           
            fitFrame = getframe(hfig);
            close(hfig)
            
            %  Food Web Graph =======
            hGraphFig = figure('windowstyle','docked');
            G = graph(fw); 
            plot(G)
            graphFrame = getframe(hGraphFig);            
            close(hGraphFig);         
            
            DatosFWs{end+1}={fwCode, [ixFW nNodes nLinks], [gaussR powerR exponR], fitFrame, graphFrame };
            % =======================================

        end
        save('FWinfo_Real.mat','DatosFWs','-mat','-v7.3')
    else
        load('FWinfo_Real.mat','-mat','DatosFWs')
    end
end
    
