%%run fountain data
%%run fountain data
clear;
addpath(genpath('RtToolbox'));
load('castle-c1.mat');
% t2 =load('castle-P30-e.mat');
% t1 =load('castle-P30-gt.mat');
strfile = 'castle-P30-c1.txt';
% e=t1.EN(1:3,:);
% e3 =EN(1:3,:);
% e33 =t2.EN(1:3,:);
% K=[2759.48 0 1520.69;0 2764.16 1006.81 ;0 0 1 ];
originaltic=tic;
%
[finalTriplets,v,weight] = buildTripletsMinimal_c(pointMatchesInliers,EN,Rc,width,hight);
nodesNum = length(unique(finalTriplets));
%
indexNonLineal = find(weight(:,1)>rad2deg(0.15));
indexLineal = find(weight(:,1)<=rad2deg(0.15));
%%% admm for triplets in  groupNonLineal
Pss = cell(size(finalTriplets,1), 1);
[Xs,Y,IS] = optimizeTripletsIRLS( EN,finalTriplets(indexNonLineal,1:3),false,false );
for i=1:size(indexNonLineal,1)
    Ns = indexNonLineal(i,1);
    [~,~,~,Psc]=getCameraRtFromOneTriplet( projectE(Xs{i}));
    Pss{Ns} = Psc;
    %
end
% figure(1);
% for j = 1:length(Pss)
%     Pst = Pss{j};
%     for i=1:length(Pst)
%         Rti = Pst{i,1};
%         scale = 0.001;
%         drawCamera1(Rti, 3072, 2048, 2759/2, scale, 1);    axis equal;
%     end
% end
%
% K =[2563.565,0,1542.212;0,2563.565,1180.39;0,0,1];
% K =[1 0 0;0 1 0;0 0 1];
for i=1:size(indexLineal,1)
    Ns = indexLineal(i,1);
    Psc=computeRTfromOneTriplet( Rt,M,K,finalTriplets(Ns,:));
    Pss{Ns} = Psc;
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
% figure(5);
% for j = 1:length(Pss)
%     Pst = Pss{j};
%     for i=1:length(Pst)
%         Rti = Pst{i,1};
%         scale = 0.001;
%         drawCamera1(Rti, 3072, 2048, 2759/2, scale, 1);    axis equal;
%     end
% end
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
%%vertification
 %[EnewC,Scale] = projectE(Enew);
[RC2,Scale] = getCameraRtFromBigE(Enew);
%%output intial pose
%%saveGesfm(RC2,strfile);

%%draw Rt
figure(3);
for i=1:nodesNum
    Rti = RC2{i,1};
    scale = 0.01;
    drawCamera1(Rti, width, hight, K(1,1), scale, 1); axis equal;
    
end




