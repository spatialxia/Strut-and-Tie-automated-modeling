classdef Class_FEM_FrameTruss
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
        function obj = Class_FEM_FrameTruss()
        end
        
        %% 1
        function [U, SE, Length, Reaction] = Solve(obj,...
                Nodes, Elements, Young, Thick, Width, FF, GDof, allDof, fixedDof, freeDof)
            Enum = size(Elements,1);
            Fn = zeros(Enum,1);
            tp2 =zeros(6,Enum);
            
            [Stiff_G, Stiff_E, Length, TM] = obj.Stiff_TrussBar(GDof, allDof, Nodes, Elements, Young, Thick, Width);
            
            U = sparse(GDof,1);
            U(freeDof,:) = Stiff_G(freeDof,freeDof) \ FF(freeDof,:);
            U(fixedDof,:) = 0;        
            
            
            SEn = 0; % Axial force strain energy
            SEs = 0; % Shear Strain energy
            SEr = 0; % Bending strain energy
            for i = 1:Enum
                index = Elements(i,:);
                ElementDof = [(index(1)-1)*3+1,(index(1)-1)*3+2,(index(1)-1)*3+3, ...
                    (index(2)-1)*3+1,(index(2)-1)*3+2,(index(2)-1)*3+3];
                Ue = U(ElementDof);
                tp = Stiff_E(:,:,i)*Ue;
                Ft = TM(:,:,i)'*tp;
                tp2(:,i) = TM(:,:,i)'*tp;
                Ut = TM(:,:,i)'*Ue;
                Fn(i) = tp2(1,i);
                SEn = SEn + 0.5*Ut(1)*Ft(1) + 0.5*Ut(4)*Ft(4);
                SEs = SEs + 0.5*Ut(2)*Ft(2) + 0.5*Ut(5)*Ft(5);
                SEr = SEr + abs(0.5*Ut(3)*Ft(3)) + abs(0.5*Ut(6)*Ft(6));
            end
            
            SE =full( U' * Stiff_G * U);
            Reaction = tp2;
        end
        
        function [Stiff_G, Stiff_E, Length, TM] = Stiff_TrussBar_original(obj, GDof,allDof, Nodes, Elements, Young, Thick, Width)
            Enum = size(Elements,1);
            Stiff_G = sparse(GDof,GDof);
            Stiff_E = zeros(6,6,Enum);
            fcn = @plus;
            Length = zeros(Enum,1);
            TM = zeros(6,6,Enum);
            for i = 1:Enum
                index = Elements(i,:);
                NodesXY = [Nodes(index(1),1), Nodes(index(1),2), Nodes(index(2),1), Nodes(index(2),2)];
                xi = Nodes(index(1),1);
                yi = Nodes(index(1),2);
                xj = Nodes(index(2),1);
                yj = Nodes(index(2),2);
                ElementDof = [(index(1)-1)*3+1,(index(1)-1)*3+2,(index(1)-1)*3+3, ...
                    (index(2)-1)*3+1,(index(2)-1)*3+2,(index(2)-1)*3+3];          
                Length(i) = sqrt((xi-xj)^2 + (yi-yj)^2);
                cc = (xj-xi)/Length(i);
                ss = (yj-yi)/Length(i);
                TM(:,:,i) = [cc,ss,0, 0, 0, 0;
                              -1*ss, cc, 0, 0, 0, 0;
                              0, 0, 1, 0, 0, 0;
                              0, 0, 0, cc, ss, 0;
                              0, 0, 0, -1*ss, cc, 0;
                              0, 0, 0, 0, 0, 1];          
                b = Thick;
                h = Width;
                Area = b*h;
                II = b*(h^3)/12;
                ee= Young*(1e0);
                ee0=Young*(1e0);
                aa = Area;
                LL = Length(i);                
                Ke = [ee0*aa/LL, 0, 0, -1*ee0*aa/LL, 0, 0;
                      0, 12*ee*II/(LL^3), 6*ee*II/(LL^2), 0, -12*ee*II/(LL^3), 6*ee*II/(LL^2);
                      0, 6*ee*II/(LL^2), 4*ee*II/(LL), 0, -6*ee*II/(LL^2), 2*ee*II/(LL);
                      -1*ee0*aa/LL, 0, 0, ee0*aa/LL, 0, 0;
                      0, -12*ee*II/(LL^3), -6*ee*II/(LL^2), 0, 12*ee*II/(LL^3), -6*ee*II/(LL^2);
                      0, 6*ee*II/(LL^2), 2*ee*II/(LL), 0, -6*ee*II/(LL^2), 4*ee*II/(LL);];
                ttM = TM(:,:,i);
                Stiff_E(:,:,i) = ttM*Ke*(ttM');
                rows = zeros(6,6);
                cols = zeros(6,6);
                for iii= 1:6
                    rows(:,iii)=ElementDof';
                    cols(iii,:) = ElementDof;
                end  
                Stiff_G = fcn(Stiff_G,sparse(rows,cols,Stiff_E(:,:,i),GDof,GDof)); 
                
                LL = Length(i);
                AA =[1,0,0,0,0,0; 0,1,0,0,0,0; 0,0,1,0,0,0;
                     1,0,0,LL,0,0; 0,1,LL,0,LL^2,LL^3;
                     0,0,1,0,2*LL,3*LL^2];
                XX = 0; 
                YY = 0; 
                HH = [0,0,0,1,0,0;
                      0,0,0,0,-2*YY,-6*XX*YY] 
            end
        end
        
        function [StressTOP, StressBOT] = BeamStress_original(obj, Nodes, Elements, Young, Thick, Width, UU, nn, XX, YY ) 
            Enum = size(Elements,1);
            Length = zeros(Enum,1);
            TM = zeros(6,6,Enum);
            i = nn ;
            Em= Young;
            Eb= Young;

                index = Elements(i,:);
                NodesXY = [Nodes(index(1),1), Nodes(index(1),2), Nodes(index(2),1), Nodes(index(2),2)];
                xi = Nodes(index(1),1);
                yi = Nodes(index(1),2);
                xj = Nodes(index(2),1);
                yj = Nodes(index(2),2);
                ElementDof = [(index(1)-1)*3+1,(index(1)-1)*3+2,(index(1)-1)*3+3, ...
                    (index(2)-1)*3+1,(index(2)-1)*3+2,(index(2)-1)*3+3];          
                Length(i) = sqrt((xi-xj)^2 + (yi-yj)^2);
                cc = (xj-xi)/Length(i);
                ss = (yj-yi)/Length(i);
                TM(:,:,i) = [cc,ss,0, 0, 0, 0;
                              -1*ss, cc, 0, 0, 0, 0;
                              0, 0, 1, 0, 0, 0;
                              0, 0, 0, cc, ss, 0;
                              0, 0, 0, -1*ss, cc, 0;
                              0, 0, 0, 0, 0, 1];          
                b = Thick;
                h = Width;
                Area = b*h;


                LL = Length(i);
                AA =[1,0,0,0,0,0; 0,1,0,0,0,0; 0,0,1,0,0,0;
                     1,0,0,LL,0,0; 0,1,LL,0,LL^2,LL^3;
                     0,0,1,0,2*LL,3*LL^2];
                XX = XX*LL; 
                YY1 = YY*0.5*h; 
                HH1 = [0,0,0,1,0,0;
                      0,0,0,0,-2*YY1,-6*XX*YY1];
                BB1 = HH1*inv(AA);
                YY2 = YY*-0.5*h;
                HH2 = [0,0,0,1,0,0;
                      0,0,0,0,-2*YY2,-6*XX*YY2];
                BB2 = HH2*inv(AA);                
                Ut = TM(:,:,i)'*UU(ElementDof);                
                StressTOP = [Em;Eb].*(BB1*Ut);
                StressBOT = [Em;Eb].*(BB2*Ut );
            
        end
        
        
        
        %% 2
        function [U, SE, Length, Reaction, SE_member] = xySolve(obj,...
                Nodes, Elements, Young, Thick, Width, FF, GDof, allDof, fixedDof, freeDof)
            Enum = size(Elements,1);
            Fn = zeros(Enum,1);
            tp2 =zeros(6,Enum);
            SE_member = zeros(4,Enum);
            [Stiff_G, Stiff_E, Length, TM] = obj.Stiff_TrussBar(GDof, allDof, Nodes, Elements, Young, Thick, Width);
            
            U = sparse(GDof,1);
            U(freeDof,:) = Stiff_G(freeDof,freeDof) \ FF(freeDof,:);
            U(fixedDof,:) = 0;        
            
            
            SEn = 0; % Axial force strain energy
            SEs = 0; % Shear Strain energy
            SEr = 0; % Bending strain energy
            for i = 1:Enum
                index = Elements(i,:);
                ElementDof = [(index(1)-1)*3+1,(index(1)-1)*3+2,(index(1)-1)*3+3, ...
                    (index(2)-1)*3+1,(index(2)-1)*3+2,(index(2)-1)*3+3];
                Ue = U(ElementDof);
                tp = Stiff_E(:,:,i)*Ue;
                Ft = TM(:,:,i)'*tp;
                tp2(:,i) = TM(:,:,i)'*tp;
                Ut = TM(:,:,i)'*Ue;
                Fn(i) = tp2(1,i);
                
                SEn =  0.5*Ut(1)*Ft(1) + 0.5*Ut(4)*Ft(4);
                SEs =  0.5*Ut(2)*Ft(2) + 0.5*Ut(5)*Ft(5);
                SEr =  abs(0.5*Ut(3)*Ft(3)) + abs(0.5*Ut(6)*Ft(6));
                SEt = 0.5*Ue'*tp;
                SE_member(:,i) = [SEn;SEs;SEr;SEt];
            end
            
            SE =0.5*full( U' * Stiff_G * U);
            Reaction = tp2;
        end      
        
        function [Stiff_G, Stiff_E, Length, TM] = Stiff_TrussBar(obj, GDof,allDof, Nodes, Elements, Young, Thick, Width)
            Enum = size(Elements,1);
            Stiff_G = sparse(GDof,GDof);
            Stiff_E = zeros(6,6,Enum);
            fcn = @plus;
            Length = zeros(Enum,1);
            TM = zeros(6,6,Enum);
            for i = 1:Enum
                index = Elements(i,:);
                NodesXY = [Nodes(index(1),1), Nodes(index(1),2), Nodes(index(2),1), Nodes(index(2),2)];
                xi = Nodes(index(1),1);
                yi = Nodes(index(1),2);
                xj = Nodes(index(2),1);
                yj = Nodes(index(2),2);
                ElementDof = [(index(1)-1)*3+1,(index(1)-1)*3+2,(index(1)-1)*3+3, ...
                    (index(2)-1)*3+1,(index(2)-1)*3+2,(index(2)-1)*3+3];          
                Length(i) = sqrt((xi-xj)^2 + (yi-yj)^2);
                cc = (xj-xi)/Length(i);
                ss = (yj-yi)/Length(i);
                TM(:,:,i) = [cc,ss,0, 0, 0, 0;
                              -1*ss, cc, 0, 0, 0, 0;
                              0, 0, 1, 0, 0, 0;
                              0, 0, 0, cc, ss, 0;
                              0, 0, 0, -1*ss, cc, 0;
                              0, 0, 0, 0, 0, 1];          
                b = Thick;
                h = Width;
                Area = b*h;
                II = b*(h^3)/12;
                ee= Young*(1e0);
                ee0=Young*(1e0);
                aa = Area;
                LL = Length(i);                
                Ke = [ee0*aa/LL, 0, 0, -1*ee0*aa/LL, 0, 0;
                      0, 12*ee*II/(LL^3), 6*ee*II/(LL^2), 0, -12*ee*II/(LL^3), 6*ee*II/(LL^2);
                      0, 6*ee*II/(LL^2), 4*ee*II/(LL), 0, -6*ee*II/(LL^2), 2*ee*II/(LL);
                      -1*ee0*aa/LL, 0, 0, ee0*aa/LL, 0, 0;
                      0, -12*ee*II/(LL^3), -6*ee*II/(LL^2), 0, 12*ee*II/(LL^3), -6*ee*II/(LL^2);
                      0, 6*ee*II/(LL^2), 2*ee*II/(LL), 0, -6*ee*II/(LL^2), 4*ee*II/(LL);];
                ttM = TM(:,:,i);
                Stiff_E(:,:,i) = ttM*Ke*(ttM');
                rows = zeros(6,6);
                cols = zeros(6,6);
                for iii= 1:6
                    rows(:,iii)=ElementDof';
                    cols(iii,:) = ElementDof;
                end  
                Stiff_G = fcn(Stiff_G,sparse(rows,cols,Stiff_E(:,:,i),GDof,GDof)); 
            end
        end

        function [StressTOP, StressBOT] = BeamStress(obj, Nodes, Elements, Young, Thick, Width, UU, nn, XX, YY ) 
            Enum = size(Elements,1);
            Length = zeros(Enum,1);
            TM = zeros(6,6,Enum);
            i = nn ;
            Em= Young;
            Eb= Young;

                index = Elements(i,:);
                xi = Nodes(index(1),1);
                yi = Nodes(index(1),2);
                xj = Nodes(index(2),1);
                yj = Nodes(index(2),2);
                ElementDof = [(index(1)-1)*3+1,(index(1)-1)*3+2,(index(1)-1)*3+3, ...
                    (index(2)-1)*3+1,(index(2)-1)*3+2,(index(2)-1)*3+3];          
                Length(i) = sqrt((xi-xj)^2 + (yi-yj)^2);
                cc = (xj-xi)/Length(i);
                ss = (yj-yi)/Length(i);
                TM(:,:,i) = [cc,ss,0, 0, 0, 0;
                              -1*ss, cc, 0, 0, 0, 0;
                              0, 0, 1, 0, 0, 0;
                              0, 0, 0, cc, ss, 0;
                              0, 0, 0, -1*ss, cc, 0;
                              0, 0, 0, 0, 0, 1];          
                b = Thick;
                h = Width;
                Area = b*h;


                LL = Length(i);
                AA =[1,0,0,0,0,0; 0,1,0,0,0,0; 0,0,1,0,0,0;
                     1,0,0,LL,0,0; 0,1,LL,0,LL^2,LL^3;
                     0,0,1,0,2*LL,3*LL^2];
                XX = XX*LL; 
                YY1 = YY*0.5*h; 
                HH1 = [0,0,0,1,0,0;
                      0,0,0,0,-2*YY1,-6*XX*YY1];
                BB1 = HH1*inv(AA);
                YY2 = YY*-0.5*h;
                HH2 = [0,0,0,1,0,0;
                      0,0,0,0,-2*YY2,-6*XX*YY2];
                BB2 = HH2*inv(AA);                
                Ut = TM(:,:,i)'*UU(ElementDof);                
                StressTOP = ([Em;Eb].*(BB1*Ut))';
                StressBOT = ([Em;Eb].*(BB2*Ut ))';
            
        end
         
    end
end

