function [ Ec,Scale] = ProjectE( E )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
[V,D]=eig(E);

dd=diag(D);
[~,inds]=sort(abs(dd),'descend');
d=dd(inds(1:6));
pd=d(d>0);
nd=d(d<0);

DD=diag(d);
VV=V(:,inds(1:6));
VP=VV(:,d>0);
VN=VV(:,d<0);
DDN=diag(nd);
DDP=diag(pd);

m_d2=zeros(6,1);
pdd=(pd-nd)/2;
m_d2(1:3)=pdd;
m_d2(4:6)=-pdd;
XD2=VP;YD2=VN;
Is = obtainsign(XD2,YD2);
XYS=(XD2+YD2*Is)*sqrt(0.5);
Scale=zeros(size(XYS,1)/3,1);
%% Vi
for i=1:size(XYS,1)/3
    V_i=XYS(i*3-2:i*3,1:3);
    %       aa=det(V_i);
    %       if aa<0
    %           V_i = -V_i;
    %       end
    [u,d_,v]=svd(V_i);
    temp=diag(d_);
    scale=(abs(temp(1))+abs(temp(2))+abs(temp(3)))/3;
    temp(1:3)=1;
    d1=diag(temp);
    E_=u*d1*v';
    
    E_=E_*scale;
    Vi(i*3-2:i*3,1:3)=E_;
    Scale(i) = scale;
end
Ui=(XD2-YD2*Is)*sqrt(0.5);
XD22=(Ui+Vi)*sqrt(0.5);
YD22=(Vi-Ui)*sqrt(0.5);
EC=[XD22 YD22]*diag(m_d2)*[XD22 YD22]';

EC = (EC+EC')/2;
for i=1:size(XYS,1)/3
    EC(i*3-2:i*3,i*3-2:i*3)=0;
end
Ec=EC;
norm(EC-E);

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


