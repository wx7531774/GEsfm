%%run fountain data
clear;
addpath(genpath('RtToolbox'));
load('castle-c1.mat');%%%Load input dataset
strfile = 'result.txt';%% save the calculated exterior parameters
originaltic=tic;
%
[finalTriplets,v,weight] = buildTripletsMinimal_c(pointMatchesInliers,EN,Rc,width,hight);%%build optimal minimum cover connected image triplet set (OMCTS)
nodesNum = length(unique(finalTriplets));
%classify the OMCTS into collinear and nonlinear ones, the therold is set
%to be 0.17 
indexNonLineal = find(weight(:,1)>rad2deg(0.17));
indexLineal = find(weight(:,1)<=rad2deg(0.17));
%%% admm for triplets in  groupNonLineal
Pss = cell(size(finalTriplets,1), 1);
[Xs,Y,IS] = optimizeTripletsIRLS( EN,finalTriplets(indexNonLineal,1:3),true,false );
for i=1:size(indexNonLineal,1)
    Ns = indexNonLineal(i,1);
    [~,~,~,Psc]=getCameraRtFromOneTriplet( projectE(Xs{i}));
    Pss{Ns} = Psc;
end
%%%% solving collinear triplet by using common ground points' depth values
m_sd=zeros(size(indexLineal,1),1);
for i=1:size(indexLineal,1)
    Ns = indexLineal(i,1);
    [Psc,m_std]=computeRTfromOneTriplet( Rt,M,K,finalTriplets(Ns,:));
    Pss{Ns} = Psc;
    m_sd(i,1)=m_std;
    %
end
%% similarity transform all the triplets to global frame
for i=1:size(v,1)
    curEdge=v(i,:);
    triplet1 = curEdge(1);
    triplet2 = curEdge(2);
    %
    [C,ia,ib]=intersect(finalTriplets(triplet1,1:3),finalTriplets(triplet2,1:3));
    [S,R,T]=  FindsimilarityTrans({Pss{triplet2}{ib(1)},Pss{triplet2}{ib(2)}},{Pss{curEdge(1)}{ia(1)},Pss{curEdge(1)}{ia(2)}});
    %% transfer triplet2
    for j=1:3
        Pss{triplet2}{j}(1:3,1:3)=Pss{triplet2}{j}(1:3,1:3)*inv(R);
        Pss{triplet2}{j}(1:3,4)=S*R*Pss{triplet2}{j}(1:3,4)+T;
    end
end
%% get camera R t from triplets global frame
RC=cell(nodesNum,1);
for i=1:nodesNum
    [is,js]=find(finalTriplets==i);
    curRt=Pss{is(1)}{js(1)};
    RC{i,1}=curRt;
end
%% compute big E from Ri ti 
for i =1:nodesNum-1
    for j=i+1:nodesNum
        Rti = RC{i,1};
        Rtj = RC{j,1};
        Eij = Rti(1:3,1:3) *(getCrossM(Rti(1:3,4))-getCrossM(Rtj(1:3,4)) )* Rtj(1:3,1:3)';%%%xiaot
        Eji = Rtj(1:3,1:3)*(getCrossM(Rtj(1:3,4))-getCrossM(Rti(1:3,4)))*Rti(1:3,1:3)';
        %%
        Enew(3*i-2:3*i,3*j-2:3*j)=Eij;
        Enew(3*j-2:3*j,3*i-2:3*i)=Eji;
    end
end
%%vertification by using the corresponding sufficient conditions in the
%%corresponding paper
 [EnewC,Scale] = projectE(Enew);
[RC2,Scale] = getCameraRtFromBigE(EnewC);
%%output intial pose
saveGesfm(RC2,strfile);

%%draw Rt to show our results
figure(2);
for i=1:nodesNum
    Rti = RC2{i,1};
    scale = 0.01;
    drawCamera1(Rti, width, hight, 1500, scale, 1); axis equal;
end

%Evaluation fountaion herzjesu_ castle
dataset = 'D:/Convex Optimization orientation/data/castle';
RT = load(fullfile(dataset,'/groundtruth.txt'));
nodesNumg=length(RT(:,1))/4;
RT_g = cell(nodesNumg,1);
figure(3);
for i=1:nodesNumg
    RTI=RT(i*4-3:i*4,1:3);
    R=RTI(1:3,1:3);
    T=RTI(4,1:3);
    Rtg=[R',T'];
    RT_g{i} = Rtg;
    scale = 0.001;
    drawCamera1(Rtg, 2000, 1500, 1500, scale, 1); axis equal;  
end
%Rotation error
Ro=RC2{1,1}(1:3,1:3);
Rg=RT_g{1,1}(1:3,1:3);
Ro2g= Ro'*Rg;
R_go = cell(nodesNumg,1);
R_err=zeros(nodesNumg,1);
for i=1:nodesNumg
   Roi=RC2{i,1}(1:3,1:3);
   Rgi=RT_g{i,1}(1:3,1:3);
   Rg_ti=Roi*Ro2g;
   R_d=Rg_ti'*Rgi;
   a=trace(R_d)/3;
   R_err(i)=rad2deg(acos(a));
   R_go{i}=Rg_ti;
end
m_rerr=mean(R_err);

%Translation error

To=cell(nodesNumg,1);
Tg=cell(nodesNumg,1);
T_go = cell(nodesNumg,1);
for i=1:nodesNumg
   To{i}=RC2{i,1}(1:3,4);
   Tg{i}=RT_g{i,1}(1:3,4);
end
T_er=999999.99;
S_=0.0;R_=[1 0 0;0 1 0;0 0 1];T_=[0 0 0];
for i =1:4096
rn=randperm(nodesNumg,3);
xl1=Tg{rn(1)};xl2=Tg{rn(2)};xl3=Tg{rn(3)};
xr1=To{rn(1)};xr2=To{rn(2)};xr3=To{rn(3)};
[s,R,T]=calculte3Dtransformation(xl1,xl2,xl3,xr1,xr2,xr3);
te=zeros(nodesNumg,1);
for j=1:nodesNumg
   e=Tg{j}-s*(R*To{j}+T);
   te(j)=norm(e);
end
T_er1=mean(te);
if T_er1<T_er
    S_=s;R_=R;T_=T;
    T_er=T_er1;
end
end
te1=zeros(nodesNumg,1);
for j=1:nodesNumg
   
   T_go{j}=S_*(R_*To{j}+T_);
   e=Tg{j}-T_go{j};
   te1(j)=norm(e);
end

for i=1:nodesNumg
    Rtg1=[R_go{i},T_go{i}];
    scale = 0.001;
    DrawCamera2(Rtg1, 2000, 1500, 1500, scale, 1); axis equal;  
end

function [s,R,T]=calculte3Dtransformation(xl1,xl2,xl3,xr1,xr2,xr3)
   xlc=(xl1(1)+xl2(1)+xl3(1))/3; ylc=(xl1(2)+xl2(2)+xl3(2))/3;zlc=(xl1(3)+xl2(3)+xl3(3))/3;
   xrc=(xr1(1)+xr2(1)+xr3(1))/3; yrc=(xr1(2)+xr2(2)+xr3(2))/3;zrc=(xr1(3)+xr2(3)+xr3(3))/3;
   tlc=[xlc ylc zlc]';trc=[xrc yrc zrc]';
   %centerlized
   xl1c=xl1-tlc;xl2c=xl2-tlc;xl3c=xl3-tlc;
   xr1c=xr1-trc;xr2c=xr2-trc;xr3c=xr3-trc;
   
   s=(norm(xl1c)+norm(xl2c)+norm(xl3c))/(norm(xr1c)+norm(xr2c)+norm(xr3c));
   
   xr1cs=s*xr1c;xr2cs=s*xr2c;xr3cs=s*xr3c;
   %H=xr1cs*xl1c'+xr2cs*xl2c'+xr3cs*xl3c';
   H=xl1c*xr1cs'+xl2c*xr2cs'+xl3c*xr3cs';
   [u,d_,v]=svd(H);
   R=u*v';
   T= (tlc-s*R*trc)/s;
end