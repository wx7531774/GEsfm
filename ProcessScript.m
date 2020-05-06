clear
dataset = 'D:/Convex Optimization orientation/data/test_unibw/gesfm';
matfile = 'Unibw.mat';
%
measurements1 = load(fullfile(dataset,'/measurement.txt'));
measurements1(isnan(measurements1))=0;
% measurements2 = load(fullfile(dataset,'/sfm.json/gesfm/structures.txt'));
% measurements2(isnan(measurements2))=0;
M=measurements1;%%%xiaot measurements2 from bin
%
camsNum=size(measurements1,1)/2;
pointMatchesInliers=zeros(camsNum,camsNum);
%
intrinsics = dlmread(fullfile(dataset,'/intrinsic.txt'));
K = intrinsics(1:3,1:3);
width = intrinsics(4,1);
hight = intrinsics(4,2);
%%
EN=zeros(camsNum*3,camsNum*3);
Rt = cell(camsNum,camsNum);
Rc =cell(camsNum,camsNum);
relativeInfo = load(fullfile(dataset,'/relativeInfo_essential.txt'));
relativeInfo(:,1:2) = relativeInfo(:,1:2)+1;
for i=1:camsNum-1
    for j=i+1:camsNum
        [flag,idx]= ismember([i,j],relativeInfo(:,1:2),'rows');
        if flag
            %% openMVG
           Rj = vec2mat(relativeInfo(idx,13:21),3);
           tj = relativeInfo(idx,22:24)';
           cj = -(Rj' *tj);
           Rt{i,j} = [Rj,tj];
           Rc{i,j}=[Rj,cj];
            Eij =-getCrossM(cj)*Rj';
            Eji = Rj*getCrossM(cj);
            EN(3*i-2:3*i,3*j-2:3*j) = Eij;
            EN(3*j-2:3*j,3*i-2:3*i)= Eji;
           %% 
           pointMatchesInliers(i,j) = relativeInfo(idx,3);           
        end
    end
end      

[ EN ] = normalizeForbineusNorm( EN );
EN(isnan(EN)) = 0;

M=M(:,(sum(abs(M)>10^-5,1)>=4));
save(matfile,'pointMatchesInliers','EN','M','Rt','Rc','K','width','hight');

