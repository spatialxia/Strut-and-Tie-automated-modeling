classdef Class_Opt_Related
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = Class_Opt_Related()
        end
        
        %% 1. Conventional filter
        function [xval_new] = Convention(obj, Nodes, Elements, rmin, xval)
            %UNTITLED8 此处显示有关此函数的摘要
            %   此处显示详细说明
            n = size(Elements,1);
            Enum = n;
            centroid = zeros(n,2);
            minim=1e-4;
            for i=1:size(Elements,1)
                index = Elements(i,:);
                ElementDof = [index index+size(Nodes,1)];
                nn = length(index);
                centroid(i,:) = [sum(Nodes(index,1))/4, sum(Nodes(index,2))/4]; 
            end


            for i=1:Enum
                conv = 0;
                conv_value1 =0;
                for j =1:Enum
                    dis = ( (centroid(i,1)-centroid(j,1))^2 + (centroid(i,2)-centroid(j,2))^2)^(0.5);
                    conv = conv + max(0,rmin-dis);
                    conv_value1 = conv_value1 + max(0,rmin-dis)*xval(j,1);
                end
                xval_new(i,1) = conv_value1 / (conv);


            end

        end
        
        %% Personal
        function [ xval_new] = StressAveraging(obj, Nodes, Elements, rmin, Stress, PreStress)
            n = size(Elements,1);
            Enum = n;
            xval = Stress;
            PS = PreStress;
            aa=-5;
            
            centroid = zeros(n,2);
            Pstp = zeros(2,size(Elements,1));
            ss_new=zeros(3,size(Elements,1));
            for i=1:size(Elements,1)
                index = Elements(i,:);
                ElementDof = [index index+size(Nodes,1)];
                nn = length(index);
                centroid(i,:) = [sum(Nodes(index,1))/4, sum(Nodes(index,2))/4]; 
            end

            for i=1:Enum
                conv = 0;
                conv_value1=0;
                conv_value2=0;
                conv_value3=0;
                for j =1:Enum
                    dis = ( (centroid(i,1)-centroid(j,1))^2 + (centroid(i,2)-centroid(j,2))^2)^(0.5);
                    if PS(1,j)>0 && ((PS(1,j))>aa*(PS(2,j))) 
%                     if 1
                        tp1=xval(1,j);
                        tp2=xval(2,j);
                        tp3=xval(3,j);
                    else
                        tp1=0; tp2=0; tp3=0;
                    end                    
                    conv = conv + max(0,rmin-dis);
                    conv_value1 = conv_value1 + max(0,rmin-dis)*tp1;    
                    conv_value2 = conv_value2 + max(0,rmin-dis)*tp2;   
                    conv_value3 = conv_value3 + max(0,rmin-dis)*tp3;   
                end
                St1 = conv_value1 / (conv);
                St2 = conv_value2 / (conv);
                St3 = conv_value3 / (conv);
                tpSP = [St1,St3;
                        St3,St2];
                DD=eig(tpSP);
                Pstp(:,i) = [max(DD);min(DD)];  
                ss_new(:,i) = [St1,St2,St3];
            end
            
            xval_new=Pstp;
        end


        function [ww] = Weights(obj, Nodes, Elements, rmin, PreStress)
            n = size(Elements,1);
            aa=-5;
            Enum = n;
            PS = PreStress;
            centroid = zeros(n,2);
            Pstp = zeros(1,size(Elements,1));
            
            tpmax = max(max(PS));
            for i=1:Enum
                if PS(1,i)>0 && (abs(PS(1,i))>aa*(PS(2,i)))
                    Pstp(i) = PS(1,i)/tpmax;
                else
                    Pstp(i)=0;
                end
                

            end
            
            ww=Pstp;            
        end
        
    end
end

