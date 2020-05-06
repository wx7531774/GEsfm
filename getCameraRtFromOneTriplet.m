function [s,r,t,Ps] = getCameraRtFromOneTriplet(Ek)
[V,D]=eig(Ek);
dd=diag(D);
[~,inds]=sort(abs(dd),'descend');
d=dd(inds(1:6));
pd=d(d>0);
nd=d(d<0);
if length(pd)~=3 || length(nd)~=3
    nerror = 1;
end
DDP=diag(pd);
DD=diag(d);
VV=V(:,inds(1:6));
X=VV(:,d>0);
Y=VV(:,d<0);
%
[Is,dIS] = obtainsign(X,Y);

Uh=(X-Y*Is)*sqrt(0.5);
Vh=(X+Y*Is)*sqrt(0.5);
%% normalize block rotation matrix
%
Ps = cell(size(Vh,1)/3,1);
R=cell(size(Vh,1)/3,1);
Scale=zeros(size(Vh,1)/3,1);
T=cell(size(Vh,1)/3,1);
for i=1:(size(Vh,1)/3)
      Vi=Vh(i*3-2:i*3,1:3);
       aa=det(Vi);
       Scale(i)=nthroot(aa,3);
       R{i}=Vi/Scale(i);
       %
    tt=inv(Vi)*Uh(i*3-2:i*3,1:3)*DDP;
    tt=0.5*(tt-tt');
%
    T{i}=[tt(3,2) tt(1,3) tt(2,1)]';
    Ps{i} = [R{i},T{i}];%% R{i} is the rotation matrix in cv
end
s =Scale;
r =R;
t =T;

end

function [Is2,dIS]=obtainsign(Xs,Ys)
%
IS = cell(8,1);
IS{1}=[1,0,0;0,1,0;0,0,1];
IS{2}=[-1,0,0;0,1,0;0,0,1];
IS{3}=[1,0,0;0,-1,0;0,0,1];
IS{4}=[1,0,0;0,-1,0;0,0,-1];
IS{5}=[-1,0,0;0,-1,0;0,0,-1];
IS{6}=[-1,0,0;0,-1,0;0,0,1];
IS{7}=[1,0,0;0,1,0;0,0,-1];
IS{8}=[-1,0,0;0,1,0;0,0,-1];
dIS = zeros(8,1);
for i=1:8
    mls = IS{i};
    for j=1:size(Xs,1)/3
      A= Xs(3*j-2 : 3*j , :)+(Ys(3*j-2 : 3*j , :)*mls); 
      dIS(i,1) = dIS(i,1) + norm(diag(A' * A))/ norm((A' *A),'fro');        
%      dIS(i,1) = dIS(i,1) + sum(diag(A' * A))/ norm((A' *A),'fro');       
    end     
end
[~,indx]=sort(dIS,'descend');
Is2 = IS{indx(1)};    
end