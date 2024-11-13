classdef Class_ThisCase

    properties
        
    end
    
    methods
        function obj = Class_ThisCase()
        end
        
        function Plot_Element_weights(obj, Nodes, Elements, ww,  cb)
            figure;
            axis equal
            set(gca,'xtick',[],'xticklabel',[])
            set(gca,'ytick',[],'yticklabel',[])   
            set(gca,'Visible','off'); 

            Nnum = size(Nodes,1);
            Enum = size(Elements,1);

            hold on
            for i =1:Enum
                index = Elements(i,:);
                XX = Nodes(index,1);
                YY = Nodes(index,2);
                alfa = ww(1,i);
                tp = alfa*[0,0,0] + (1-alfa)*[1,1,1];
                fill(XX,YY,tp,'LineStyle','none');

            end
            hold off
            % color bar settings
            if cb>0
                cmax = 1;
                cmin = 0;
                if cmin == cmax
                    cmin = cmax - 1;
                end
                nn = cb;
                cc = [];
                for i =1:nn
                    aa = (i-1)/(nn-1);
                    cc(end+1,:)=[0,0,0]*aa + (1-aa)*[1,1,1];
                end

                colormap(cc);
                TheBar = colorbar;

                % color bar
                c2=get(TheBar,'YTick');
                cc2 = [];
                for i=1:cb
                    cc2(i) = c2(1)+(c2(end)-c2(1))*(i-1)/(cb-1);
                end
                set(TheBar,'YTick',cc2);
                c2 = get(TheBar, 'YTickLabel');
                cc2={};                        
                for i=1:cb
                    tp=cmin+(cmax-cmin)*(i-1)/(cb-1);
                    cc2(end+1)={num2str(tp,'%5.2f')};
                end

                set(TheBar,'YTickLabel',cc2);
                set(TheBar,'fontsize',14);



            end
                    
            

        end
        

    end
end

