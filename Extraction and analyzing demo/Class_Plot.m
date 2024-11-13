classdef Class_Plot
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
        %% ---- initial ------
        function obj = Class_Plot()
        end
        
        %% ---- truss geometry----------
        function Figs = PlotBasicTruss(obj, Nodes, Elements)           
            Nodes4 = Nodes;
            Figs = figure;
            axis equal;
            set(gca,'Visible','off');
            hold on
            for tt = 1:size(Elements,1)
                plot( [Nodes4(Elements(tt,1),1),Nodes4(Elements(tt,2),1)],...
                    [Nodes4(Elements(tt,1),2),Nodes4(Elements(tt,2),2)],'r-','LineWidth',3)
            end            
            for tt = 1:size(Nodes4,1)
                plot( Nodes4(tt,1),Nodes4(tt,2),'b.', 'markersize', 40)
            end

            hold off
        end


        %% ----- Continuum plot ------

        function Plot_Element(obj, Nodes, Elements)
            figure;
            axis equal
            set(gca,'xtick',[],'xticklabel',[])
            set(gca,'ytick',[],'yticklabel',[])       
            Enum = size(Elements,1);
            set(gca,'Visible','off');
            hold on 
            colors = zeros(Enum,3);
            for i =1:Enum
                index = Elements(i,:);
                XX = Nodes(index,1);
                YY = Nodes(index,2);
                

                colors(i,:) = [0,1,1];
                patch('vertices',Nodes,'faces',Elements,'facecolor', color);
                
                set(h,'edgecolor',[0 0.5 0.5])
            end   
            hold off  
            rotate3d on
        end

        function Plot_Element_PrincipleStress(obj, Nodes, Elements, PS, option, scalarFactor, cb, mm)
            aa=-5;
            figure;
            axis equal
            set(gca,'xtick',[],'xticklabel',[])
            set(gca,'ytick',[],'yticklabel',[])   
            set(gca,'Visible','off'); 

            Nnum = size(Nodes,1);
            Enum = size(Elements,1);
            SS = scalarFactor;
            
            PSmax = SS*1;
            PSmin =  -1*SS*1;
            switch option
                case '1'
                    hold on
                    for i =1:Enum
                        index = Elements(i,:);
                        XX = Nodes(index,1);
                        YY = Nodes(index,2);
                        PStp = PS(:,i);
                        if PStp(1)>0 && ((PStp(1))>aa*(PStp(2)))
                            alfa = min(1,(PStp(1)/PSmax));
                            tp = alfa*[1,0,0] + (1-alfa)*[1,1,1];
                            fill(XX,YY,tp,'LineStyle','none');
                        else
                            alfa = min(1, (PStp(2)/PSmin));
                            tp = alfa*[0,0,1] + (1-alfa)*[1,1,1];
                            fill(XX,YY,tp,'LineStyle','none');
                        end
                    end
                    hold off
                    % color bar and color map settings
                    if cb>0
                        cmax = mm(2);
                        cmin = mm(1);
                        if cmin == cmax
                            cmin = cmax - 1;
                        end
                        nn = cb;
                        cc = [];
                        for i =1:nn/2
                            aa = (nn/2-i)/(nn/2-1);
                            cc(end+1,:)=[1,0,0]*aa + (1-aa)*[1,1,1];
                        end
                        for i =1:nn/2
                            aa = (1-(nn/2-i)/(nn/2-1));
                            cc(end+1,:)=[0,0,1]*aa + (1-aa)*[1,1,1];
                        end
                        colormap(cc);
                        TheBar = colorbar;
                        

                        c2=get(TheBar,'YTick');
                        cc2 = [];
                        for i=1:cb
                            cc2(i) = c2(1)+(c2(end)-c2(1))*(i-1)/(cb-1);
                        end
                        set(TheBar,'YTick',cc2);
                        c2 = get(TheBar, 'YTickLabel');
                        cc2={};                        
                        for i=1:ceil(cb/2)
                            tp=cmin+(0-cmin)*(i-1)/(ceil(cb/2)-1);
                            cc2(end+1)={num2str(tp,'%5.2e')};
                        end
                        for i=1+ceil(cb/2):cb
                            tp=(cmax-0)*(i-1-ceil(cb/2))/(cb-ceil(cb/2)-1);
                            cc2(end+1)={num2str(tp,'%5.2e')};
                        end
                        set(TheBar,'YTickLabel',cc2);
                        set(TheBar,'fontsize',14);
                                      
                        
                    end
                    
                case '2'
                    hold on
                    for i =1:Enum
                        index = Elements(i,:);
                        XX = Nodes(index,1);
                        YY = Nodes(index,2);
                        PStp = PS(:,i);
                        if PStp(1)<0
                            alfa=0;
                        else
                            alfa = min(1,(PStp(1)/PSmax));
                        end
                        tp = alfa*[1,0,0] + (1-alfa)*[1,1,1];
                        fill(XX,YY,tp,'LineStyle','none');
                    end
                    hold off        
                    figure;
                    axis equal
                    set(gca,'xtick',[],'xticklabel',[])
                    set(gca,'ytick',[],'yticklabel',[])                       
                    hold on
                    for i =1:Enum
                        index = Elements(i,:);
                        XX = Nodes(index,1);
                        YY = Nodes(index,2);
                        PStp = PS(:,i);
                        if PStp(1)>0
                            alfa=0;
                        else
                            alfa = min(1, (PStp(2)/PSmin));
                        end
                        tp = alfa*[0,0,1] + (1-alfa)*[1,1,1];
                        fill(XX,YY,tp,'LineStyle','none');                        
                    end
                    hold off                    
            end
            
        end
        

    end
end

