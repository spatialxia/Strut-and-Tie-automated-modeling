classdef Class_FEM_TrussBar
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
        function obj = Class_FEM_TrussBar()
        end   
        
        %% Common truss analysis - equilibrium      
        function [U, SE, Length, Reaction, Reaction_s, SE_n , SE_s, Stiff_G] = Solve(...
                obj, Nodes, Elements, Young, Area, FF,GDof, allDof, fixedDof, freeDof)
            Enum = size(Elements,1);
            Fn = zeros(Enum,1);
            tp2 = zeros(2,Enum);
            tp3 = zeros(1,GDof);
            
            
            [Stiff_G, Stiff_E, Length, TM] = obj.Stiff_TrussBar(GDof,allDof, Nodes, Elements, Young, Area);
            
            U = sparse(GDof,1);
            U(freeDof,:) = Stiff_G(freeDof,freeDof) \ FF(freeDof,:);
            U(fixedDof,:) = 0;        
            
            SE_n = 0;
            for i = 1:Enum
                index = Elements(i,:);
                ElementDof = [(index(1)-1)*2+1,(index(1))*2, (index(2)-1)*2+1,(index(2))*2 ];    
                Ue = U(ElementDof);
                tp = Stiff_E(:,:,i)*Ue;
                tp2(:,i) = TM(:,:,i)'*tp;     
                SE_n = SE_n + 0.5*Ue'*Stiff_E(:,:,i)*Ue;
            end
            
            SE_s = 0;
                tp = [1,0,0,0; 0,1,0,0; 
                      0,0,1,0; 0,0,0,1];     
                  kse= eye(4)*1e-6.*tp;
            for i = 1:Enum
                index = Elements(i,:);
                ElementDof = [(index(1)-1)*2+1,(index(1))*2, (index(2)-1)*2+1,(index(2))*2 ];    
                Ue = U(ElementDof);    
                SE_s = SE_s + 0.5*Ue'*kse*Ue;
            end            

            SE =full( U' * Stiff_G * U);
            Reaction = tp2;
            Reaction_s = tp3;
        end
        
        function [Stiff_G, Stiff_E, Length, TM] = Stiffness(obj, GDof,allDof, Nodes, Elements, Young, Area)
            Enum = size(Elements,1);
            Stiff_G = sparse(GDof,GDof);
            Stiff_E = zeros(4,4,Enum);
            fcn = @plus;
            Length = zeros(Enum,1);
            TM = zeros(4,2,Enum);
            for i = 1:Enum
                index = Elements(i,:);
                NodesXY = [Nodes(index(1),1), Nodes(index(1),2), Nodes(index(2),1), Nodes(index(2),2)];
                xi = Nodes(index(1),1);
                yi = Nodes(index(1),2);
                xj = Nodes(index(2),1);
                yj = Nodes(index(2),2);
                ElementDof = [(index(1)-1)*2+1,(index(1))*2, (index(2)-1)*2+1,(index(2))*2 ];           
                Length(i) = sqrt((xi-xj)^2 + (yi-yj)^2);
                cc = (xj-xi)/Length(i);
                ss = (yj-yi)/Length(i);
                TM(:,:,i) = [cc,0;ss,0;0,cc;0,ss];
                ke = Young*Area/Length(i)*[1,-1;-1,1];
                ttM = TM(:,:,i);               
                Stiff_E(:,:,i) = ttM*ke*(ttM');
                
                tp = [1,0,0,0; 0,1,0,0; 
                      0,0,1,0; 0,0,0,1];
                
                tke = Stiff_E(:,:,i)+ eye(4)*1e-6.*tp;
                
                rows = zeros(4,4);
                cols = zeros(4,4);
                for iii= 1:4
                    rows(:,iii)=ElementDof';
                    cols(iii,:) = ElementDof;
                end  
                Stiff_G = fcn(Stiff_G,sparse(rows,cols,tke,GDof,GDof)); 
            end
        end
        
        %% own use
        % type one: spring on global dofs
        function [U, SE, Length, Reaction, Reaction_s, SE_n , SE_s,...
                f0val, fval, df0dx, dfdx] = xySolve(...
                obj,Nodes, Elements, Young, Area, FF,GDof, allDof, fixedDof, freeDof, var)
            Enum = size(Elements,1);
            Fn = zeros(Enum,1);
            tp2 = zeros(2,Enum);
            tp3 = zeros(1,GDof);    
  
            df0dx = zeros(GDof,1);
            dfdx =  zeros(1,GDof);
            
            [Stiff_G, Stiff_E_n, Stiff_E_s, Length, TM] = obj.xyStiff_TrussBar(GDof,allDof, Nodes, Elements, Young, Area, var);
            
            U = sparse(GDof,1);
            U(freeDof,:) = Stiff_G(freeDof,freeDof) \ FF(freeDof,:);
            U(fixedDof,:) = 0;        
            
            SE_n = 0;
            for i = 1:Enum
                index = Elements(i,:);
                ElementDof = [(index(1)-1)*2+1,(index(1))*2, (index(2)-1)*2+1,(index(2))*2 ];    
                Ue = U(ElementDof);
                tp = Stiff_E_n(:,:,i)*Ue;
                tp2(:,i) = TM(:,:,i)'*tp;     
                SE_n = SE_n + 0.5*Ue'*Stiff_E_n(:,:,i)*Ue;
            end
                        
            SE_s = 0;
            for i = 1:GDof
                SE_s = SE_s + 0.5*U(i)*U(i)*var(i);
                tp3(1,i) = U(i)*var(i);
            end

            SE =full( U' * Stiff_G * U);
            Reaction = tp2;
            Reaction_s = tp3;
            
  
            for j = 1:GDof
                tpdf0dx = 0;
                tpdf0dx2 = 0;
                
                for i = 1:Enum
                    index = Elements(i,:);
                    ElementDof = [(index(1)-1)*2+1,(index(1))*2, (index(2)-1)*2+1,(index(2))*2 ];    
                    if sum(ismember(ElementDof,j))==1
                        Ue = U(ElementDof);
                        index2 = find(ElementDof==j);
                        tpm = zeros(4,4);
                        tpm(index2,index2) = 1;
                        tpm2 = Stiff_E_n(:,:,i);
                        tpm3 = inv(Stiff_E_n(:,:,i)+Stiff_E_s(:,:,i));
                        tpdf0dx = tpdf0dx + Ue'*tpm2*tpm3*tpm*Ue;
                        tpdf0dx2 = tpdf0dx2 + (0.5*tpm3*tpm*Ue)'*tpm2*Ue + 0.5*(Ue')*tpm2*(tpm3*tpm*Ue);
                    end
                end
                df0dx(j,1) = tpdf0dx; 
            end
            f0val = -1*SE_n;
            fval = -1;
        end
        
        function [Stiff_G, Stiff_E_n, Stiff_E_s, Length, TM] = xyStiff_TrussBar(obj, GDof,allDof, Nodes, Elements, Young, Area, var)
            Enum = size(Elements,1);
            Stiff_G = sparse(GDof,GDof);
            Stiff_E = zeros(4,4,Enum);
            Stiff_E2 = zeros(4,4,Enum);
            
            fcn = @plus;
            Length = zeros(Enum,1);
            TM = zeros(4,2,Enum);
            for i = 1:Enum
                index = Elements(i,:);
                NodesXY = [Nodes(index(1),1), Nodes(index(1),2), Nodes(index(2),1), Nodes(index(2),2)];
                xi = Nodes(index(1),1);
                yi = Nodes(index(1),2);
                xj = Nodes(index(2),1);
                yj = Nodes(index(2),2);
                ElementDof = [(index(1)-1)*2+1,(index(1))*2, (index(2)-1)*2+1,(index(2))*2 ];           
                Length(i) = sqrt((xi-xj)^2 + (yi-yj)^2);
                cc = (xj-xi)/Length(i);
                ss = (yj-yi)/Length(i);
                TM(:,:,i) = [cc,0;ss,0;0,cc;0,ss];
                ke = Young*Area/Length(i)*[1,-1;-1,1];
                ttM = TM(:,:,i);               
                Stiff_E(:,:,i) = ttM*ke*(ttM');
                
                tp = [var(ElementDof(1)),0,0,0; 0,var(ElementDof(2)),0,0; 
                      0,0,var(ElementDof(3)),0; 0,0,0,var(ElementDof(4))];
                
                Stiff_E2(:,:,i) = tp;  
                tke = Stiff_E(:,:,i) + tp;
                
                rows = zeros(4,4);
                cols = zeros(4,4);
                for iii= 1:4
                    rows(:,iii)=ElementDof';
                    cols(iii,:) = ElementDof;
                end  
                Stiff_G = fcn(Stiff_G,sparse(rows,cols,tke,GDof,GDof)); 
            end
            
            Stiff_E_n = Stiff_E;
            Stiff_E_s = Stiff_E2;
        end
                 
        % type two: equilibrium truss-spring element - shear equilibrium
        function [U, SE, Length, Reaction, SE_n , SE_s, factors,...
                UN, SEN] = xySolve2(...
                obj,Nodes, Elements, Young, Area, FF,GDof, allDof, fixedDof, freeDof)
            Enum = size(Elements,1);
            Fn = zeros(Enum,1);
            tpR = zeros(4,Enum);
            UN = zeros(4,Enum);
  
            df0dx = zeros(GDof,1);
            dfdx =  zeros(1,GDof);
            
            [Stiff_G, Stiff_E, Length, TM1, TM2] = obj.xyStiff_TrussBar2(GDof,allDof, Nodes, Elements, Young, Area);
            
            U = sparse(GDof,1);
            U(freeDof,:) = Stiff_G(freeDof,freeDof) \ FF(freeDof,:);
            U(fixedDof,:) = 0;        
            
            SE_n = 0;
            SE_s = 0;
            SEN = 0;
            for i = 1:Enum
                index = Elements(i,:);
                ElementDof = [(index(1)-1)*2+1,(index(1))*2, (index(2)-1)*2+1,(index(2))*2 ];    
                Ue = U(ElementDof);

                               
                tp = Stiff_E(:,:,i)*Ue;
                tpR(:,i) = TM1(:,:,i)'*tp;
                
                
                tp = TM2(:,:,i)'*Ue;
                tp2 = [tpR(1,i),tpR(3,i)]';
                SEN = SEN + 0.5*tp(1:2)'*tp2;

                tp(1:2) = 0;
                UN(:,i) = TM2(:,:,i)*tp;
                
                
                SE_n = SE_n + abs(tpR(1,i)*Length(i));
                SE_s = SE_s + abs(tpR(2,i)*Length(i));
            end
                        
            factors = SE_n / (SE_n + SE_s);
            SE =full( U' * Stiff_G * U);
            Reaction = tpR;
            
        end
        
        function [Stiff_G, Stiff_E, Length, TM1, TM2] = xyStiff_TrussBar2(obj, GDof,allDof, Nodes, Elements, Young, Area)
            Enum = size(Elements,1);
            Stiff_G = sparse(GDof,GDof);
            Stiff_E = zeros(4,4,Enum);
            Stiff_E_0 = zeros(4,4,Enum);
            
            fcn = @plus;
            Length = zeros(Enum,1);
            TM1 = zeros(4,4,Enum);
            TM2 = zeros(4,3,Enum);
            for i = 1:Enum
                index = Elements(i,:);
                NodesXY = [Nodes(index(1),1), Nodes(index(1),2), Nodes(index(2),1), Nodes(index(2),2)];
                xi = Nodes(index(1),1);
                yi = Nodes(index(1),2);
                xj = Nodes(index(2),1);
                yj = Nodes(index(2),2);
                ElementDof = [(index(1)-1)*2+1,(index(1))*2, (index(2)-1)*2+1,(index(2))*2 ];           
                Length(i) = sqrt((xi-xj)^2 + (yi-yj)^2);
                LL = Length(i);
                cc = (xj-xi)/Length(i);
                ss = (yj-yi)/Length(i);
                TM1(:,:,i) = [cc,-ss,0,0; ss, cc, 0,0; 0,0,cc,-ss; 0,0,ss,cc];  
                TM2(:,:,i) = [cc,0,-1*ss; ss, 0, cc; 0, cc, ss; 0,ss,-1*cc];         
                EA = Young*Area/Length(i);
                ES = 1e-9;
                ke = [EA,-1*EA, 0; 0,0,ES; -EA, EA, 0; 0, 0, -ES];
                ttM1 = TM1(:,:,i); 
                ttM2 = TM2(:,:,i);
                keg = ttM1*ke*(ttM2');
                Stiff_E(:,:,i) = keg;
%                 Stiff_E_0(:,:,i) = ke;
                
                rows = zeros(4,4);
                cols = zeros(4,4);
                for iii= 1:4
                    rows(:,iii)=ElementDof';
                    cols(iii,:) = ElementDof;
                end  
                Stiff_G = fcn(Stiff_G,sparse(rows,cols,keg,GDof,GDof)); 
            end
            
        end
   
        %% considering densities
        function [U, SE, Length, Reaction, Reaction_s, NodalForces, SE_n , SE_s] = SolveD(...
                obj,Nodes, Elements, Young, Area, FF,GDof, allDof, fixedDof, freeDof, tpElements)
            Enum = size(Elements,1);
            Fn = zeros(Enum,1);
            tp1 = zeros(4,Enum);
            tp2 = zeros(2,Enum);
            tp3 = zeros(1,GDof);
            
            
            [Stiff_G, Stiff_E, Length, TM] = obj.Stiff_TrussBarD(GDof,allDof, Nodes, Elements, Young, Area,tpElements);
            
            U = sparse(GDof,1);
            U(freeDof,:) = Stiff_G(freeDof,freeDof) \ FF(freeDof,:);
            U(fixedDof,:) = 0;        
            
            SE_n = 0;
            for i = 1:Enum
                index = Elements(i,:);
                ElementDof = [(index(1)-1)*2+1,(index(1))*2, (index(2)-1)*2+1,(index(2))*2 ];    
                Ue = U(ElementDof);
                tp = Stiff_E(:,:,i)*Ue;
                tp1(:,i) = tp;
                tp2(:,i) = TM(:,:,i)'*tp;     
                SE_n = SE_n + 0.5*Ue'*Stiff_E(:,:,i)*Ue;
            end
            
            SE_s = 0;
            for i = 1:GDof
                SE_s = SE_s + 0.5*U(i)*1e-6*U(i)*1;
                tp3(1,i) = 1e-6*U(i)*1;
            end

            SE =full( U' * Stiff_G * U);
            NodalForces = tp1;
            Reaction = tp2;
            Reaction_s = tp3;
        end
        
        function [Stiff_G, Stiff_E, Length, TM] = Stiff_TrussBarD(obj, GDof,allDof, Nodes, Elements, Young, Area,tpElements)
            Enum = size(Elements,1);
            Stiff_G = sparse(GDof,GDof);
            Stiff_E = zeros(4,4,Enum);
            fcn = @plus;
            Length = zeros(Enum,1);
            TM = zeros(4,2,Enum);
            for i = 1:Enum
                tp=sum(ismember(tpElements,i));
                if tp==0
                    dd=1e-0;
                else
                    dd=1e0;
                end
                index = Elements(i,:);
                NodesXY = [Nodes(index(1),1), Nodes(index(1),2), Nodes(index(2),1), Nodes(index(2),2)];
                xi = Nodes(index(1),1);
                yi = Nodes(index(1),2);
                xj = Nodes(index(2),1);
                yj = Nodes(index(2),2);
                ElementDof = [(index(1)-1)*2+1,(index(1))*2, (index(2)-1)*2+1,(index(2))*2 ];           
                Length(i) = sqrt((xi-xj)^2 + (yi-yj)^2);
                cc = (xj-xi)/Length(i);
                ss = (yj-yi)/Length(i);
                TM(:,:,i) = [cc,0;ss,0;0,cc;0,ss];
                ke = Young*Area/Length(i)*[1,-1;-1,1];
                ttM = TM(:,:,i);               
                Stiff_E(:,:,i) = ttM*ke*(ttM')*dd;
                
                tp = [1,0,0,0; 0,1,0,0; 
                      0,0,1,0; 0,0,0,1];
                
                tke = Stiff_E(:,:,i);
                
                rows = zeros(4,4);
                cols = zeros(4,4);
                for iii= 1:4
                    rows(:,iii)=ElementDof';
                    cols(iii,:) = ElementDof;
                end  
                Stiff_G = fcn(Stiff_G,sparse(rows,cols,tke,GDof,GDof)); 
            end
        end
    end
end

