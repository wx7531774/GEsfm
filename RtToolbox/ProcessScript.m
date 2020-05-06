clear

width = 3072;
hight = 2048;
% image_size = repmat([width,hight]',1,length(P));
% centers =repmat([width/2,hight/2]',1,length(P));
% points = u_uncalib.points;
% index = u_uncalib.index;
% measurements = zeros(2*length(points),u_uncalib.pointnr);
% for i = 1:length(points)
%     for j = 1:size(points{i},2)
%         measurements(2*i-1:2*i,index{i}(j)) = points{i}(1:2,j);
%     end
% end
dataset = 'E:/xiao.teng/sfmIPI/3gpsfm/dataset/fountain-P11/sfm.json';
measurements = load(fullfile(dataset,'/tracks/measurement.txt'));
measurements(isnan(measurements))=0;
camsNum=size(measurements,1)/2;


pointMatchesInliers=zeros(camsNum,camsNum);
% pointMatchesGround=cell(3,3,2);
M=measurements;

%%
EN=zeros(camsNum*3,camsNum*3);
Rt = cell(camsNum,camsNum);
relativeInfo = load(fullfile(dataset,'/relativeInfo_essential.txt'));
relativeInfo(:,1:2) = relativeInfo(:,1:2)+1;
for i=1:camsNum-1
    for j=i+1:camsNum
        [flag,idx]= ismember([i,j],relativeInfo(:,1:2),'row');
        if flag
            %%
           rij = vec2mat(relativeInfo(idx,13:21),3);
           tij = relativeInfo(idx,22:24)';
           Rt{i,j} = [rij,tij];
           %%
           Eij = vec2mat(relativeInfo(idx,4:12),3)';
%            Eij = getCrossM(tij)*rij;
           EN(3*i-2:3*i,3*j-2:3*j)=Eij;
           EN(3*j-2:3*j,3*i-2:3*i)=Eij';
           pointMatchesInliers(i,j) = relativeInfo(idx,3);
           
        end
    end
end      

% [ EN ] = normalizeForbineusNorm( EN );
% FN(isnan(FN)) = 0;
% pointMatchesInliers = converPointInlier(pointMatchesInliers);
M=M(:,(sum(abs(M)>10^-5,1)>=4));
save(' fountain-P11.mat','pointMatchesInliers','EN','M','Rt','width','hight');

