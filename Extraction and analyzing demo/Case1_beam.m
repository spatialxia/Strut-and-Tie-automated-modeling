clear;clc;
%--------------------------------------------------------------------------------------%
% Created on 24.02.2020
% author:Yi Xia
% Program applicatiao: A critical evaluation of topology optimization results for strut-and-tie modeling of reinforced concrete
% Content:1.A TO result extraction method for systematically transforming the optimized topology into a truss-like structure; 
%         2.a critical evaluation of topology optimization results for strut-and-tie modeling of reinforced concrete
% Input: topology optimization results
% Output: truss-like structure and three evaluation indexes
% TRS,SR and STS index are three indexes for evaluating truss-like models
% TRS index: Tensile region similarity index, tensile stress fields of the original structure and TO results are compared
% SR index: Steel reinforcement ratio, the SR ratio is calculated as the volume fraction of steel with respect to the concrete volume
% STS index: Suitable truss structure index, it measures the degree to which the obtained truss-like structure

%-----------------------------------------------------------------------------------------%
%% library
% creating instances,import the required functions
xyThin = Class_Thinning2D();
xy_Plot = Class_Plot();
xy_Mesh = Class_FEM_Mesh();
xy_FEM = Class_FEM_PlaneElement();
xy_BEAMFEM=Class_FEM_FrameTruss();
xy_SPE = Class_Opt_Related();
xy_Case = Class_ThisCase();
xyThin = Class_Thinning2D;
xy_PlaneFEM = Class_FEM_PlaneElement();

%Import topology optimization results
FFF = imread('Case1_beam.jpg'); 

%% truss-like structure extraction
%-----------------%
%According to the topology optimization results, the structure size, the grid size, the fixed loading point and the 
%    boundary condition position, the truss-like structure is extracted.
%----------------- %

% Generating binary designs,first clear solid/void (0/1) binary designs based on TO density fields are generated
BWW = im2bw(FFF,0.9); % Use 0.9 as the threshold
BWWd = double(BWW); 

xx=ceil(150); % The number of columns of the density matrix is determined according to the structure size and mesh size.
yy=ceil(94);% The number of rows of the density matrix is determined according to the structure size and mesh size.
BWW2=xyThin.SimplifyPoints(BWW,[xx,yy]);
BWW2=flipud(BWW2);
IM0 = 1-BWW2;
IM = IM0;
    figure;
    imshow(BWW2)
    
% Thinning
pd0 = sum(sum(IM0));
dd=1;
count=0;
IM2 = IM *0;
IM2(1,1)=1;%The density of loading point and boundary condition is fixed.
IM2(1,150)=1;%The density of loading point and boundary condition is fixed.
IM2(94,94)=1;% The density of loading point and boundary condition is fixed.
while dd==1
    count=count+1;
    [NM, jd] = xyThin.Thin_improve( IM, 2,IM2);
    IM = NM;
    [NM, jd] = xyThin.Thin_improve( IM, 1,IM2);
    IM = NM;
    pd1 = sum(sum(NM));
    if pd1 == pd0
        dd=0;
    else
        pd0=pd1;
    end
end
figure;
imshow(1-NM)

specialNodes = [1,1;150,1; 94,94];% The density of loading point and boundary condition is fixed.

% Node and bar extraction
Nodes = xyThin.DetermineNodes(NM,specialNodes);
[Nodes] = xyThin.SimplifyNodes(Nodes, 1.5*2, specialNodes);


[Elements , NM2, TopoRelation] = xyThin.MatchingElement_TopoRelation(Nodes, NM, 3, specialNodes); %holes, large factor, loose control.
[Nodes2, TopoRelation2] =xyThin.ImproveNodes(Nodes, TopoRelation, 3, specialNodes);
[Nodes3, TopoRelation3] =xyThin.ImproveNodes(Nodes2, TopoRelation2, 5, specialNodes);

NodesT = Nodes3;
TopoRelationT = TopoRelation3;
[Elements] = xyThin.MatchingElement_TopoRelation_withoutHole(NodesT, NM, NM2, TopoRelationT);
[Elements] = xyThin.MatchingElement_Fitting(NodesT, NM, 2, Elements, TopoRelationT); %last lines connect, large factor, loose control.

TrussElements = Elements;
TrussNodes = NodesT;

    % plot generated truss-like structures
    xy_Plot.PlotBasicTruss(TrussNodes, TrussElements);

clearvars Elements Nodes;

%% Continuum FEM
% Finite Element Analysis : Stiffness, Load and Boundary Conditions
TOmesh = BWW2;

Nodes = zeros((xx+1)*(yy+1),2);
Elements = zeros(1,4);
TOelements = [];

count = 0;
numx = xx;
numy =yy;
for i =1:(numy+1)
    for j = 1:(numx+1)
        count = count + 1;
        Nodes(count,:) = [j-1, i-1]; 
    end
end
count = 0;
for i = 1:numy
    for j = 1:numx
        
        if j>11 && j<40 && i>11 && i<40
        else
            count = count +1;
            Elements(count,:) = [(i-1)*(numx+1) + j, (i-1)*(numx+1)+ j+1, (i)*(numx+1) + j+1, (i)*(numx+1) + j];
            if TOmesh(i,j) == 0
                TOelements(end+1) = count;
            end
        end
        
        
    end
end

ElementsTO = zeros(1,4);
for i = 1:size(TOelements,2)
    ElementsTO(i,:) = Elements(TOelements(i),:);
end

E_young = 20e9;% Elastic modulus
Poisson = 0.2;% Poisson 's ratio
D = E_young/(1-Poisson^2) * [1 Poisson 0; Poisson 1 0; 0 0 (1-Poisson)/2];
thickness = 0.3;
option = '2-Points';
NodeNumber = size(Nodes,1);
GDof = 2 * NodeNumber;
tp = []; 
count=0;
for i =1:(numy+1)
    for j = 1:(numx+1)
        count = count + 1;
        tptp = ismember(Elements,count);
        if sum(sum(tptp))==0
            tp(end+1) = count;
            tp(end+1) = count + NodeNumber;
        end        
    end
end

fixedDof = [1; 1+NodeNumber; 151 + NodeNumber];% Loading point and boundary condition setting

allDof = [1:GDof];
allDof = setdiff(allDof, tp);
freeDof = setdiff(allDof, fixedDof);

F = sparse(GDof, 1);
FFNode1 = intersect(find(abs(Nodes(:,1)-94)<0.001),find(abs(Nodes(:,2)-94)<0.001));
ffDof = [FFNode1+NodeNumber];% Structural freedom setting of boundary conditions
F(ffDof,1)=-3000000;% applied load

Enum = size(Elements,1);
density = zeros(Enum,1);
density(:)=1;

Nodes(:,1) = Nodes(:,1)/150*7;
Nodes(:,2) = Nodes(:,2)/94*4.7;
TrussNodes(:,1) = TrussNodes(:,1)/150*7;
TrussNodes(:,2) = TrussNodes(:,2)/94*4.7;

%% FEA
% Finite element analysis : The plane four-node element is used 
%      to solve the structural displacement and stress distribution.
[ U0, SE, E_Stiff0, Stress0, Pri_Stress0, Pri_ang0, StrainMatrix, ElementArea] = xy_PlaneFEM.FEM_Solve_2(...
    Nodes,Elements,thickness,D,GDof,F,freeDof, fixedDof,option,density, TOelements);
 
TO_Stress = Stress0;
TO_Pri_Stress = Pri_Stress0;

[ U1, SE, E_Stiff1, Stress1, Pri_Stress1, Pri_ang1, SM] = xy_PlaneFEM.FEM_Solve_3(...
                Nodes,Elements,thickness,D,GDof,F,freeDof, fixedDof,option, TO_Stress,TO_Pri_Stress,1,0,E_Stiff0,U0);
  
%% TRS index    
% Calculate the TRS index
LL=0;% LL is the total length of the extracted truss-like structure calculated by summing lengths of all its members
for i =1:size(TrussElements,1)
    index = TrussElements(i,:);
    v1 = [TrussNodes(index(1),1),TrussNodes(index(1),2)];
    v2 = [TrussNodes(index(2),1),TrussNodes(index(2),2)];
    LL = LL+norm(v2-v1);
end
rmin=0.5*33/LL; % 33 is the ratio of structure volume to thickness.rmin is the averaging radius,
[avg_Stress_Ori] = xy_SPE.StressAveraging(Nodes, Elements, rmin, Stress1, Pri_Stress1);
[avg_Stress_TO] = xy_SPE.StressAveraging(Nodes, Elements, rmin, TO_Stress, Pri_Stress0);

[ww1] = xy_SPE.Weights(Nodes, Elements, rmin, avg_Stress_TO);
[ww2] = xy_SPE.Weights(Nodes, Elements, rmin, avg_Stress_Ori);
img1=zeros(xx,yy);
img2=img1;
count=0;
for i = 1:numy
    for j = 1:numx
        
        if j>11 && j<40 && i>11 && i<40
            img1(j,i)=0;
            img2(j,i)=0;
        else
            count = count +1;
            img1(j,i)=ww1(count);
            img2(j,i)=ww2(count);
        end     
    end
end

r1=ssim(img1,img2,'Exponents',[1,1,1]);

TRS_index = r1;

    % plot principal stress and weights
    xy_Plot.Plot_Element_PrincipleStress(Nodes, Elements, Pri_Stress0,'1',1.9315e+08*0.5, 8, [1.9315e+08,-1.9315e+08]*0.5)
    xy_Case.Plot_Element_weights( Nodes, Elements, ww1,  5)

%% SR_index
% Calculate the SR index
tf=0;
fy=450e6;% Yield strength of rebar
tpp=[];
for j=1:Enum
    x1= Pri_Stress0(1,j);
    x2= Pri_Stress0(2,j);
    a2= Pri_ang0(j);
   
    if x1>0 && (x1)>-5*(x2) 
        TPSE = (x1/fy)*33/(size(Elements,1));
        tf=tf+TPSE;
    end
end  
tf = tf/33;
SR_index = tf;

%% STS_index
% The STS index is calculated by beam element.
NodesT = TrussNodes;
ElementsT = TrussElements;
[ElementsT,NodesT] = xy_Mesh.FRefine(ElementsT,NodesT, 2, 0.2);
Enum = size(Elements,1);
Thick = thickness;

GDof = 3*size(NodesT,1);
allDof = [1:GDof];

% Boundary condition application : including position determination and structural freedom fixing
FixNode1 = intersect(find(abs(NodesT(:,1)-0.04667)<1e-4),find(abs(NodesT(:,2)-0.05)<1e-4));
FixNode2 = intersect(find(abs(NodesT(:,1)-7)<1e-4),find(abs(NodesT(:,2)-0.05)<1e-4));
fixedDof = [(FixNode1-1)*3+1, (FixNode1-1)*3+2, (FixNode2-1)*3+2];
freeDof = setdiff(allDof, fixedDof);


% Loading : determination of loading point and application of force
FFNode1 = intersect(find(abs(NodesT(:,1)-4.387)<1e-3),find(abs(NodesT(:,2)-4.7)<1e-3));
FF = sparse(GDof, 1);
Fdof = [(FFNode1-1)*3+2];
FF(Fdof) = -3000000;
Height = rmin;

% Internal force calculation of member:including rod axial force and shear force
[U1, SE, Length, Reactions, SEmember ] = xy_BEAMFEM.xySolve(...
    NodesT, ElementsT, E_young, thickness, thickness*0.01, FF, GDof, allDof, fixedDof, freeDof);
tttt1=0;
for i=1:size(Length,1)
    ll=Length(i);
    tttt1 = tttt1+ abs(Reactions(1,i))/(abs(Reactions(1,i))+abs(Reactions(2,i)));
end
APPI=tttt1/size(Length,1);
STS_index = APPI;

