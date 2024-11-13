classdef Class_FEM_PlaneElement
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
        %% initial
        function obj = Class_FEM_PlaneElement()
        end
        
        %% methods 
        function [ U, SE, E_Stiff, Stress, Pri_Stress] = FEM_Solve(obj, Nodes,Elements,thickness,D,GDof,F,freeDof, fixedDof,option)
            [ G_Stiff, E_Stiff,SM ] = obj.FEM_Stiffness_2D( GDof, Nodes, Elements, D, thickness,option );
            U = sparse(GDof,1);
            U(freeDof,:) = G_Stiff(freeDof,freeDof) \ F(freeDof,:);
            U(fixedDof,:) = 0;
            SE =full( U' * G_Stiff * U);    
            
            
            Stress = zeros(3,size(Elements,1));
            Pri_Stress = zeros(2,size(Elements,1));
            for i = 1:size(Elements,1)
                index = Elements(i,:);
                ElementDof = [index index+size(Nodes,1)];
                nn = length(index);  
                Ue = U(ElementDof);
                tp = SM(:,:,i);
                Stress(:,i) =  tp*Ue; 
                tpSP = [Stress(1,i),Stress(3,i);
                        Stress(3,i),Stress(2,i)];
                DD=eig(tpSP);
                Pri_Stress(:,i) = [max(DD);min(DD)];
            end
            
        end
        
        function [ G_Stiff, E_Stiff, SM ] = FEM_Stiffness_2D(obj, GDof, Nodes, Elements, D, thickness,option )
            G_Stiff = sparse(GDof, GDof);
            E_Stiff = zeros(8,8,size(Elements,1));
            fcn = @plus;
            [GaussWeights, GaussLocations] = obj.FEM_GaussQuadrature(option);

            parfor i=1:size(Elements,1)
                index = Elements(i,:);
                ElementDof = [index index+size(Nodes,1)];
                nn = length(index);

                for q =1:size(GaussWeights,1)
                    GaussPoint = GaussLocations(q,:);
                    xi = GaussPoint(1);
                    eta = GaussPoint(2);

                    [Shape_Natural, Shape_Natural_Drtv] =  obj.FEM_PlaneStress_ShapeFunctionQ4( xi, eta );

                    [JacobinaMatrix, JacobinaMatrix_inv, XYdrtv] = obj.FEM_Jacobian_Iso(Nodes(index,:), Shape_Natural_Drtv);

                    B = zeros(3, 2*nn);
                    B(1,1:nn) = XYdrtv(:,1)';
                    B(2,nn+1:2*nn) = XYdrtv(:,2)';
                    B(3,1:nn) = XYdrtv(:,2)';
                    B(3,nn+1:2*nn) = XYdrtv(:,1)';

                    tp = B'*D*thickness*B*GaussWeights(q)*det(JacobinaMatrix);
                    E_Stiff(:,:,i) = E_Stiff(:,:,i) + tp;
                    
                    rows = zeros(8,8);
                    cols = zeros(8,8);
                    for iii= 1:8
                        rows(:,iii)=ElementDof';
                        cols(iii,:) = ElementDof;
                    end  
                    G_Stiff = fcn(G_Stiff,sparse(rows,cols,tp,GDof,GDof));                       
%                     G_Stiff(ElementDof,ElementDof) = G_Stiff(ElementDof,ElementDof) + tp; 

                end
            end
            
            SM = zeros(3,8, size(Elements,1));
            for i = 1:size(Elements,1)
                index = Elements(i,:);
                ElementDof = [index index+size(Nodes,1)];
                nn = length(index);  
                xi = 0.5;
                eta = 0.5;     
                [Shape_Natural, Shape_Natural_Drtv] =  obj.FEM_PlaneStress_ShapeFunctionQ4( xi, eta );

                [JacobinaMatrix, JacobinaMatrix_inv, XYdrtv] = obj.FEM_Jacobian_Iso(Nodes(index,:), Shape_Natural_Drtv);

                B = zeros(3, 2*nn);
                B(1,1:nn) = XYdrtv(:,1)';
                B(2,nn+1:2*nn) = XYdrtv(:,2)';
                B(3,1:nn) = XYdrtv(:,2)';
                B(3,nn+1:2*nn) = XYdrtv(:,1)';

                tp = D*thickness*B*1;
                SM(:,:,i) = tp;
                    
            end

        end   
        
        function [ weights, locations ] = FEM_GaussQuadrature(obj, option )
        switch option
            case '2-Points'
                locations = ...
                    [ -0.577350269189626 -0.577350269189626;
                       0.577350269189626 -0.577350269189626;
                       0.577350269189626  0.577350269189626;
                      -0.577350269189626  0.577350269189626];
                weights = [1;1;1;1];

            case '1-Points'
                locations =[0 0];
                weights = [4];
        end

        end
        
        function [ Shape_Natural, Shape_Natural_Drtv] = FEM_PlaneStress_ShapeFunctionQ4(obj, xi, eta )
            Shape_Natural = 1/4*[ (1-xi)*(1-eta); (1+xi)*(1-eta);
                (1+xi)*(1+eta); (1-xi)*(1+eta)];

            Shape_Natural_Drtv = 1/4* [-(1-eta), -(1-xi); 1-eta, -(1+xi);
                1+eta, 1+xi; -(1+eta), (1-xi)];
        end
        
        function [ JacobinaMatrix, JacobinaMatrix_inv, XYdrtv ] = FEM_Jacobian_Iso(obj, NodeCordinates, Shape_Natural_Drtv)

            JacobinaMatrix = NodeCordinates' * Shape_Natural_Drtv;
            JacobinaMatrix_inv = inv(JacobinaMatrix);
            XYdrtv = Shape_Natural_Drtv*JacobinaMatrix_inv;

        end


        
        function [ U, SE, E_Stiff, Stress, Pri_Stress] = SolveForTO(obj, Nodes,Elements,thickness,D,GDof,F,freeDof, fixedDof,option,density)            
            [ G_Stiff, E_Stiff,SM ] = obj.FEM_Stiffness_2D_TO( GDof, Nodes, Elements, D, thickness,option, density);
            U = sparse(GDof,1);
            U(freeDof,:) = G_Stiff(freeDof,freeDof) \ F(freeDof,:);
            U(fixedDof,:) = 0;
            SE =full( U' * G_Stiff * U);    
            
            
            Stress = zeros(3,size(Elements,1));
            Pri_Stress = zeros(2,size(Elements,1));
            for i = 1:size(Elements,1)
                index = Elements(i,:);
                ElementDof = [index index+size(Nodes,1)];
                nn = length(index);  
                Ue = U(ElementDof);
                tp = SM(:,:,i);
                Stress(:,i) =  tp*Ue; 
                tpSP = [Stress(1,i),Stress(3,i);
                        Stress(3,i),Stress(2,i)];
                DD=eig(tpSP);
                Pri_Stress(:,i) = [max(DD);min(DD)];
            end
            
        end
        
        function [ G_Stiff, E_Stiff, SM ] = FEM_Stiffness_2D_TO(obj, GDof, Nodes, Elements, D, thickness,option, density)
            G_Stiff = sparse(GDof, GDof);
            E_Stiff = zeros(8,8,size(Elements,1));
            fcn = @plus;
            [GaussWeights, GaussLocations] = obj.FEM_GaussQuadrature(option);

            parfor i=1:size(Elements,1)
                index = Elements(i,:);
                ElementDof = [index index+size(Nodes,1)];
                nn = length(index);
                dd = density(i);
                for q =1:size(GaussWeights,1)
                    GaussPoint = GaussLocations(q,:);
                    xi = GaussPoint(1);
                    eta = GaussPoint(2);

                    [Shape_Natural, Shape_Natural_Drtv] =  obj.FEM_PlaneStress_ShapeFunctionQ4( xi, eta );

                    [JacobinaMatrix, JacobinaMatrix_inv, XYdrtv] = obj.FEM_Jacobian_Iso(Nodes(index,:), Shape_Natural_Drtv);

                    B = zeros(3, 2*nn);
                    B(1,1:nn) = XYdrtv(:,1)';
                    B(2,nn+1:2*nn) = XYdrtv(:,2)';
                    B(3,1:nn) = XYdrtv(:,2)';
                    B(3,nn+1:2*nn) = XYdrtv(:,1)';

                    tp = B'*D*thickness*B*GaussWeights(q)*det(JacobinaMatrix);
                    E_Stiff(:,:,i) = E_Stiff(:,:,i) + tp;
                    
                    rows = zeros(8,8);
                    cols = zeros(8,8);
                    for iii= 1:8
                        rows(:,iii)=ElementDof';
                        cols(iii,:) = ElementDof;
                    end  
                    tp = tp*dd^3;
                    G_Stiff = fcn(G_Stiff,sparse(rows,cols,tp,GDof,GDof));                       
%                     G_Stiff(ElementDof,ElementDof) = G_Stiff(ElementDof,ElementDof) + tp; 

                end
            end
            
            SM = zeros(3,8, size(Elements,1));
            for i = 1:size(Elements,1)
                index = Elements(i,:);
                ElementDof = [index index+size(Nodes,1)];
                nn = length(index);  
                xi = 0.5;
                eta = 0.5;     
                [Shape_Natural, Shape_Natural_Drtv] =  obj.FEM_PlaneStress_ShapeFunctionQ4( xi, eta );

                [JacobinaMatrix, JacobinaMatrix_inv, XYdrtv] = obj.FEM_Jacobian_Iso(Nodes(index,:), Shape_Natural_Drtv);

                B = zeros(3, 2*nn);
                B(1,1:nn) = XYdrtv(:,1)';
                B(2,nn+1:2*nn) = XYdrtv(:,2)';
                B(3,1:nn) = XYdrtv(:,2)';
                B(3,nn+1:2*nn) = XYdrtv(:,1)';

                tp = D*thickness*B*1;
                SM(:,:,i) = tp;
                    
            end

        end   
         
        function [ U, SE, E_Stiff, Stress, Pri_Stress, Pri_angle, SM, Area] = FEM_Solve_2(obj,...
                Nodes,Elements,thickness,D,GDof,F,freeDof, fixedDof,option,density, TOelements)
            
            [ G_Stiff, E_Stiff,SM, Area ] = obj.FEM_Stiffness_2D_2(...
                GDof, Nodes, Elements, D, thickness,option,density, TOelements );
            
            U = sparse(GDof,1);
            U(freeDof,:) = G_Stiff(freeDof,freeDof) \ F(freeDof,:);
            U(fixedDof,:) = 0;
            SE =full( U' * G_Stiff * U);    
            

            [GaussWeights, GaussLocations] = obj.FEM_GaussQuadrature('2-Points');
            Stress = zeros(3,size(Elements,1));
            Pri_Stress = zeros(2,size(Elements,1));
            Pri_angle = zeros(size(Elements,1),1);
            TE = TOelements;
            for i = 1:size(Elements,1)
                index = Elements(i,:);
                ElementDof = [index index+size(Nodes,1)];
                nn = length(index);  
                Ue = U(ElementDof);
                tp = ismember(TE,i);
                if sum(tp)==0
                    dd = 1e-9;
                else
                    dd = 1 ;
                end            
                
%                 tp = SM(:,:,i);
%                 Stress(:,i) =  tp*Ue; 
                
                for q =1:size(GaussWeights,1)
                    GaussPoint = GaussLocations(q,:);
                    xi = GaussPoint(1);
                    eta = GaussPoint(2);
                    [Shape_Natural, Shape_Natural_Drtv] =  obj.FEM_PlaneStress_ShapeFunctionQ4( xi, eta );
                    [JacobinaMatrix, JacobinaMatrix_inv, XYdrtv] = obj.FEM_Jacobian_Iso(Nodes(index,:), Shape_Natural_Drtv);
                    B = zeros(3, 2*nn);
                    B(1,1:nn) = XYdrtv(:,1)';
                    B(2,nn+1:2*nn) = XYdrtv(:,2)';
                    B(3,1:nn) = XYdrtv(:,2)';
                    B(3,nn+1:2*nn) = XYdrtv(:,1)';
                    Stress(:,i) = Stress(:,i) + D*1*B*Ue*GaussWeights(q)*dd*1;
                end   
                
                tpSP = [Stress(1,i),Stress(3,i);
                        Stress(3,i),Stress(2,i)];
                DD=eig(tpSP);
                Pri_Stress(:,i) = [max(DD);min(DD)];
                
                angle = 0.5*atan(2*Stress(3,i)/(Stress(1,i)-Stress(2,i)));
                tpp = Stress(1,i)*(cos(angle))^2 + Stress(2,i)*(sin(angle))^2 + 2*Stress(3,i)*sin(angle)*cos(angle);
                if abs(tpp-max(DD))<0.00001*abs(max(DD))
                    Pri_angle(i) = angle;
                else
                    Pri_angle(i) = angle+deg2rad(90);
                end                
            end
            %}
            

        end
        

        function [ G_Stiff, E_Stiff,SM,Area] = FEM_Stiffness_2D_2(obj,...
                GDof, Nodes, Elements, D, thickness,option, density, TOelements)
            G_Stiff = sparse(GDof, GDof);
            E_Stiff = zeros(8,8,size(Elements,1));
            TE = TOelements;
            fcn = @plus;
            [GaussWeights, GaussLocations] = obj.FEM_GaussQuadrature(option);
            Area = zeros(size(Elements,1),1);
            parfor i=1:size(Elements,1)
                index = Elements(i,:);
                ElementDof = [index index+size(Nodes,1)];
                nn = length(index);
                tp = ismember(TE,i);
                if sum(tp)==0
                    dd = 1e-9;
                else
                    dd = 1 ;
                end
                
                for q =1:size(GaussWeights,1)
                    GaussPoint = GaussLocations(q,:);
                    xi = GaussPoint(1);
                    eta = GaussPoint(2);

                    [Shape_Natural, Shape_Natural_Drtv] =  obj.FEM_PlaneStress_ShapeFunctionQ4( xi, eta );

                    [JacobinaMatrix, JacobinaMatrix_inv, XYdrtv] = obj.FEM_Jacobian_Iso(Nodes(index,:), Shape_Natural_Drtv);

                    B = zeros(3, 2*nn);
                    B(1,1:nn) = XYdrtv(:,1)';
                    B(2,nn+1:2*nn) = XYdrtv(:,2)';
                    B(3,1:nn) = XYdrtv(:,2)';
                    B(3,nn+1:2*nn) = XYdrtv(:,1)';

                    tp = B'*D*thickness*B*GaussWeights(q)*det(JacobinaMatrix)*dd;
                    E_Stiff(:,:,i) = E_Stiff(:,:,i) + tp;
%                     G_Stiff(ElementDof,ElementDof) = G_Stiff(ElementDof,ElementDof) + tp; 
                    Area(i) = Area(i)+det(JacobinaMatrix)*GaussWeights(q);
                    
                    rows = zeros(8,8);
                    cols = zeros(8,8);
                    for iii= 1:8
                        rows(:,iii)=ElementDof';
                        cols(iii,:) = ElementDof;
                    end  
                    G_Stiff = fcn(G_Stiff,sparse(rows,cols,tp,GDof,GDof));                     
                    
                end
            end
            
            SM = zeros(3,8, size(Elements,1));
            for i = 1:size(Elements,1)
                tp = ismember(TE,i);
                if sum(tp)==0
                    dd = 1e-9;
                else
                    dd = 1;
                end                
                index = Elements(i,:);
                ElementDof = [index index+size(Nodes,1)];
                nn = length(index);  
                xi = 0.0;
                eta = 0.0;     
                [Shape_Natural, Shape_Natural_Drtv] =  obj.FEM_PlaneStress_ShapeFunctionQ4( xi, eta );

                [JacobinaMatrix, JacobinaMatrix_inv, XYdrtv] = obj.FEM_Jacobian_Iso(Nodes(index,:), Shape_Natural_Drtv);

                B = zeros(3, 2*nn);
                B(1,1:nn) = XYdrtv(:,1)';
                B(2,nn+1:2*nn) = XYdrtv(:,2)';
                B(3,1:nn) = XYdrtv(:,2)';
                B(3,nn+1:2*nn) = XYdrtv(:,1)';

                tp = D*thickness*B*dd;
                SM(:,:,i) = tp;
                    
            end

        end   

        function [equFF] = FEM_equNodalForces(obj, Nodes, Elements, globalDofs, Pri_Stress, U, StiffE)
            UU = U;
            PS = Pri_Stress;
            GDofs = globalDofs;
            aa=0.2;
            
            equFF = zeros(GDofs,1);
            
            Enum = size(Elements,1);

            for i = 1: Enum
                index = Elements(i,:);
                ElementDof = [index index+size(Nodes,1)];     
                nn = length(index);  
                Ue = UU(ElementDof);
                Ke = StiffE(:,:,i);
                if PS(1,i)>0 && (abs(PS(1,i))>aa*-1*(PS(2,i)))
                    equFF(ElementDof) = equFF(ElementDof) + Ke*Ue;
                end
                
               
            end
        end        
          
        %% methods 
        function [ U, SE, E_Stiff, Stress, Pri_Stress, Pri_angle, SM] = FEM_Solve_3(obj,...
                Nodes,Elements,thickness,D,GDof,F,freeDof, fixedDof,option, TO_Stress,TO_PriStress,aa,bb, TOStiff, TOUU)
            
            FF1 = F;
            [equFF] = obj.FEM_equNodalForces( Nodes, Elements, GDof, TO_PriStress, TOUU, TOStiff);
            
            [ G_Stiff, E_Stiff, SM ] = obj.FEM_Stiffness_2D(...
                GDof, Nodes, Elements, D, thickness, option);
            
            FF= aa*FF1 - bb*equFF;
            
            U = sparse(GDof,1);
            U(freeDof,:) = G_Stiff(freeDof,freeDof) \ FF(freeDof,:);
            U(fixedDof,:) = 0;
            SE =full( U' * G_Stiff * U);    
            
            Stress = zeros(3,size(Elements,1));
            Pri_Stress = zeros(2,size(Elements,1));
            Pri_angle = zeros(size(Elements,1),1);
            
            [GaussWeights, GaussLocations] = obj.FEM_GaussQuadrature('2-Points');
            for i = 1:size(Elements,1)
                index = Elements(i,:);
                ElementDof = [index index+size(Nodes,1)];
                nn = length(index);  
                Ue = U(ElementDof);

                for q =1:size(GaussWeights,1)
                    GaussPoint = GaussLocations(q,:);
                    xi = GaussPoint(1);
                    eta = GaussPoint(2);
                    [Shape_Natural, Shape_Natural_Drtv] =  obj.FEM_PlaneStress_ShapeFunctionQ4( xi, eta );
                    [JacobinaMatrix, JacobinaMatrix_inv, XYdrtv] = obj.FEM_Jacobian_Iso(Nodes(index,:), Shape_Natural_Drtv);
                    B = zeros(3, 2*nn);
                    B(1,1:nn) = XYdrtv(:,1)';
                    B(2,nn+1:2*nn) = XYdrtv(:,2)';
                    B(3,1:nn) = XYdrtv(:,2)';
                    B(3,nn+1:2*nn) = XYdrtv(:,1)';
                    Stress(:,i) = Stress(:,i) + D*B*1*Ue*GaussWeights(q)*1;
                end                
                
                tpSP = [Stress(1,i),Stress(3,i);
                        Stress(3,i),Stress(2,i)];
                DD=eig(tpSP);
                Pri_Stress(:,i) = [max(DD);min(DD)];
                
                angle = 0.5*atan(2*Stress(3,i)/(Stress(1,i)-Stress(2,i)));
                tpp = Stress(1,i)*(cos(angle))^2 + Stress(2,i)*(sin(angle))^2 + 2*Stress(3,i)*sin(angle)*cos(angle);
                if abs(tpp-max(DD))<0.001*abs(max(DD))
                    Pri_angle(i) = angle;
                else
                    Pri_angle(i) = angle+deg2rad(90);
                end
                    
            end
            
        end
                
    end
end

