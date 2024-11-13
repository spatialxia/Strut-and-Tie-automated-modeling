classdef Class_FEM_Mesh
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
        %% initial
        function obj = Class_FEM_Mesh()
        end
        
        %% 1. Generating meshes

        function [NewElements, NewNodes] = FRefine(obj, Elements, Nodes, option, factor)
            op = option; % 1: division; 2: maxsize
            ff = factor;
            E0 = Elements;
            N0 = Nodes;     
            E1 = [];
            N1 = [];
            
            iterN = 0; 
            iterE = 0;   
            switch op
                case 1
                    NumE0 = size(E0,1);
                    for i = 1:NumE0
                        index1 = E0(i,1);
                        index2 = E0(i,2);
                        cor1 = N0(index1,:);
                        cor2 = N0(index2,:);
                            iterN = iterN + 1;
                            x0 = (1-1)*(cor2(1,1)-cor1(1,1))/(ff) + cor1(1,1);
                            y0 = (1-1)*(cor2(1,2)-cor1(1,2))/(ff) + cor1(1,2);
                            if iterN>1
                                a=ismember(N1,[x0,y0],'rows');
                                b = find(a==1);
                                if sum(a)<1
                                    N1(iterN,:) = [x0,y0];
                                    jp1 = iterN;
                                else
                                    iterN = iterN -1;
                                    jp1=b;
                                end   
                            else
                                N1(iterN,:) = [x0,y0];
                                jp1 = iterN;                                
                            end

                        for j = 1:ff
                            iterN = iterN + 1;
                            x2 = (j)*(cor2(1,1)-cor1(1,1))/(ff) + cor1(1,1);
                            y2 = (j)*(cor2(1,2)-cor1(1,2))/(ff) + cor1(1,2);
                            if iterN>1
                                a=ismember(N1,[x2,y2],'rows');
                                b = find(a==1);
                                if sum(a)<1
                                    N1(iterN,:) = [x2,y2];
                                    jp2 = iterN;
                                else
                                    iterN = iterN -1;
                                    jp2=b;
                                end   
                            else
                                N1(iterN,:) = [x2,y2];
                                jp2 = iterN;                                
                            end
                                                  
                            
                            iterE = iterE+1;
                            E1(iterE,:)=[jp1, jp2];
                            
                            jp1 = jp2;
                        end 
                    end
                    
                case 2 
                    NumE0 = size(E0,1);
                    for i = 1:NumE0
                        index1 = E0(i,1);
                        index2 = E0(i,2);
                        cor1 = N0(index1,:);
                        cor2 = N0(index2,:);
                        length = sqrt((cor2-cor1)*(cor2-cor1)');
                        ffnum = ceil(length/ff);
                        
                        iterN = iterN + 1;
                        x0 = (1-1)*(cor2(1,1)-cor1(1,1))/(ffnum) + cor1(1,1);
                        y0 = (1-1)*(cor2(1,2)-cor1(1,2))/(ffnum) + cor1(1,2);
                        if iterN>1
                            a=[];
                            for jj=1:size(N1,1)
                                if abs(N1(jj,1)-x0) < 1e-3*ff && abs(N1(jj,2)-y0) < 1e-3*ff
                                    a(jj)=1;
                                else
                                    a(jj)=0;
                                end
                            end
%                             a=ismember(N1,[x0,y0],'rows');
                            b = find(a==1);
                            if sum(a)<1
                                N1(iterN,:) = [x0,y0];
                                jp1 = iterN;
                            else
                                iterN = iterN -1;
                                jp1=b;
                            end   
                        else
                            N1(iterN,:) = [x0,y0];
                            jp1 = iterN;                                
                        end  

                        for j = 1:ffnum
                            iterN = iterN + 1;
                            x2 = (j)*(cor2(1,1)-cor1(1,1))/(ffnum) + cor1(1,1);
                            y2 = (j)*(cor2(1,2)-cor1(1,2))/(ffnum) + cor1(1,2);
                            if iterN>1
                                a=[];
                                for jj=1:size(N1,1)
                                    if abs(N1(jj,1)-x2) < 1e-3*ff && abs(N1(jj,2)-y2) < 1e-3*ff
                                        a(jj)=1;
                                    else
                                        a(jj)=0;
                                    end
                                end                                
%                                 a=ismember(N1,[x2,y2],'rows');
                                b = find(a==1);
                                if sum(a)<1
                                    N1(iterN,:) = [x2,y2];
                                    jp2 = iterN;
                                else
                                    iterN = iterN -1;
                                    jp2=b;
                                end   
                            else
                                N1(iterN,:) = [x2,y2];
                                jp2 = iterN;                                
                            end
                                                  
                            
                            iterE = iterE+1;
                            E1(iterE,:)=[jp1, jp2];
                            
                            jp1 = jp2;
                        end                         
                    end
            end
            
            NewElements = E1;
            NewNodes = N1;
        end

    end
end

