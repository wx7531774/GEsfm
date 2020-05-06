function [Psc, m_std]= computeRTfromOneTriplet(Rt,M,K,C)

r1 = eye(3);
t1 =[0,0,0]';
r12 = Rt{ C(1,1),C(1,2)}(:,1:3);
t12 = Rt{ C(1,1),C(1,2)}(:,4);
r13 = Rt{ C(1,1),C(1,3)}(:,1:3);
t13 = Rt{ C(1,1),C(1,3)}(:,4);
%%
index1 = intersect(find(M(2*C(1,1),:) ~=0),find(M(2*C(1,2),:) ~=0));%%% Xin W
index2 = intersect(find(M(2*C(1,1),:) ~=0),find(M(2*C(1,3),:) ~=0));%%% Xin W
index = intersect(index1,index2);
matches = [M(2*C(1,1)-1:2*C(1,1),index);M(2*C(1,2)-1:2*C(1,2),index);M(2*C(1,3)-1:2*C(1,3),index)];
%
P1 =K*[r1,t1];
P2 =K*[r12,t12];
P3 = K*[r13,t13];
% DLT
lam = zeros(size(matches,2),1);
for j=1:size(matches,2)
    A12 = [matches(1,j)*P1(3,:) - P1(1,:); matches(2,j)*P1(3,:) - P1(2,:);matches(3,j)*P2(3,:) - P2(1,:);matches(4,j)*P2(3,:) - P2(2,:)];
    A13 = [matches(1,j)*P1(3,:) - P1(1,:); matches(2,j)*P1(3,:) - P1(2,:);matches(5,j)*P3(3,:) - P3(1,:);matches(6,j)*P3(3,:) - P3(2,:)];
    [U,D,V]=svd(A12);
    X1 =V(:,4)/V(end,4);
    [U,D,V]=svd(A13);
    X2 =V(:,4)/V(end,4);
   lam(j,1)= X2(3,1)/X1(3,1);
   if lam(j,1)<0
       nerr=1;
   end
end
if length(find(lam<0)) >size(matches,2)/10
    nerr=1;
end

B=removeOutliers(lam);%% to robustlt calculated the ratios of depths
m_std=std(lam);
mlam = mean(B);
if isnan(r13'*(t13/mlam))
    nerror=-1;
end
if mlam<0
    nerror=-1;
end

Psc = cell(size(C,2),1);
Psc{1}=[r1,-r1'*t1];
Psc{2}=[r12,-r12'*t12];
Psc{3}=[r13,-r13'*(t13/mlam)];
end

function B = removeOutliers(A)

avg=mean(A);
a_std=std(A);
B=A;
sas=100000;
while sas~=0
    C=B(abs(B-avg)<2*a_std);
    avg=mean(C);
    a_std=std(C);
    sas=size(B,1)-size(C,1);
    B=C;
end

end