data1 =load('castle-P30-e.mat');
data2 =load('castle-P30-e0.mat');
data3 = load('castle-P30-gt.mat');
e1 = data1.EN;
e2 =data2.EN;
e3 = data3.EN;
e1(isnan(e1)) = 0;
e2(isnan(e2)) = 0;
e3(isnan(e3)) = 0;
% dlmwrite('castle-P30en-e.txt',data1.EN(7:9,:)');
% dlmwrite('castle-P30en-bin.txt',data2.EN(7:9,:)');
% dlmwrite('castle-P30en-gt.txt',data3.EN(7:9,:)');

% if max(abs(data1.EN-data2.EN),[],'all')>0.7
%     nerr=1;
% end
