function [EK,Et,IS] = optimizeTripletsIRLS( EN,Cf,show,doIRLS)%xiaot
%OPTIMZEOLSSON Summary of this function goes here
%   Detailed explanation goes here

if nargin<3
    show=true;
end
% Cf=[Cf;18 19 20;19 20 21];
nn=size(EN,1);
n=size(Cf,1);
Eh=cell(n,1);%% EK observation
GAMMAS=cell(n,1);
SITA=cell(n,1);
BK =  cell(n,1);
DK = cell(n,1);
EK = cell(n,1);%%EK
IS = cell(n,1);
%%%%
%t=0
for i =1:n    
    GAMMAS{i}=zeros(9,9);
    SITA{i}=zeros(9,9);
    %
    Eh{i}=getTrippleInds(EN, Cf(i,:) );
    BK{i} = Eh{i};
    DK{i} = Eh{i};    
end

%% t=1:1000 for admm
sings=zeros(n,30000);
wss=ones(n,1);
count=1;
if doIRLS
    iterssIrlss=20;
else
    iterssIrlss=1;
end
%
for irls=1:iterssIrlss
    for t=1:100%1000
        %%compute Et
        Et = solveEK_Y4(GAMMAS,SITA,BK,DK,nn,Cf,wss,Eh);
        %%compute Bkt and Dkt, and update
        for k=1:n
            Ekt = getTrippleInds(Et, Cf(k,:) );
            GAMMASt_1 = GAMMAS{k};
            Bkt = solveBK(Ekt,GAMMASt_1);
            BK{k} = Bkt;
            %
            SITAt_1 = SITA{k};
            [Dkt,Is] = solveDK(Ekt,SITAt_1);
            DK{k}=Dkt;
            IS{k}=Is;
            %%update
            GAMMASt = GAMMASt_1 + Bkt -Ekt;
            SITAt = SITAt_1 + Dkt -Ekt;
            GAMMAS{k} = GAMMASt;
            SITA{k} = SITAt;
            %
            EK{k} = Ekt;
            ttt=svd(Ekt);
            sings(k,count)= ttt(7)/ttt(6);
        end        
        count=count+1;
    end
    %% weight the triplet
    
end
if show 
        a = log10(mean(sings,1));
    a = a(1:1000);
    noiseVec = 1:length(a);
    p  =regress(a',[1./noiseVec;noiseVec;ones(size(noiseVec))]');%polyfit(noiseVec,a,2);%regress(a',[1./noiseVec;noiseVec;ones(size(noiseVec))]'); %polyfit(noiseVec,a,2);regress(a',[1./noiseVec;noiseVec;ones(size(noiseVec))]'); %
    newa = [];
    
    for i = 1:length(noiseVec)
        newa(i) = p(1)/noiseVec(i)+p(2)*noiseVec(i)+p(3);%p(1)*noiseVec(i)^2+p(2)*noiseVec(i)+p(3);% p(1)/noiseVec(i)+p(2)*noiseVec(i)+p(3);%
    end
    newa(800:end) = newa(800);
    figure('DefaultAxesFontSize',19), plot((  newa),'LineWidth',3),xlabel('Iterations'),ylabel('$\displaystyle log_{10}(\frac{\sigma_7}{\sigma_6})$','interpreter','latex')
    
end

end
function Y= solveEK_Y4(GAMMAS,SITA,BK,DK,nn,Cf,wss,Eh)
%
Y=zeros(nn);
W=zeros(nn);
param1 = 100;%% parameters setting in the paper
param2 = 0.01;%% parameters setting in the paper
%
for k=1:length(GAMMAS)
    [~,tr]=getTrippleInds(zeros(nn), Cf(k,:) );
    BKt_1 = BK{k};  DKt_1 = DK{k};
    GAMMASt_1 = GAMMAS{k}; SITAt_1 = SITA{k};
    Ekh = Eh{k};
    % vertified
    col = 0.01;
    if (~all(abs(BKt_1(:)-Ekh(:))<col))  || (~all(abs(DKt_1(:)-Ekh(:))<col))
        berror = 1;
    end
    Y(tr,tr)=Y(tr,tr)+param1*(BKt_1+GAMMASt_1)+param2*(DKt_1+SITAt_1)+wss(k)*Ekh;
    W(tr,tr)=W(tr,tr)+ (param1+param2+wss(k))*ones(9);
end
Y=Y./W;
Y(W==0)=0;
Y=(Y+Y')/2;
for i=1:nn/3
    Y(3*i-2:3*i,3*i-2:3*i)=0;
end

end
function [Dkt,Is] = solveDK(Ekt,SITAt_1)
DK1 = Ekt-SITAt_1;
for iteration = 1:1
    [V,D]=eig(DK1);
    astopb = norm(DK1,'fro');
    dd=diag(D);
    [~,inds]=sort(abs(dd),'descend');
    d=dd(inds(1:6));
    pd=d(d>0);
    nd=d(d<0);
    m_d2=zeros(6,1);
    if length(pd) ~= 3 || length(nd) ~= 3
        nerror =1;
    end
    m_d2(1:3)=pd;
    m_d2(4:6)=nd;
    VV=V(:,inds(1:6));
    XD2=VV(:,d>0);
    YD2=VV(:,d<0);

%     [V,D]=eig(DK1);
%     astopb = norm(DK1,'fro');
%     dd=diag(D);
%     [~,inds]=sort(dd,'descend');
%     [~,inds_]=sort(dd,'ascend');
%     pd=dd(inds(1:3));
%     nd=dd(inds_(1:3));
%     m_d2=zeros(6,1);
%     m_d2(1:3)=pd;
%     m_d2(4:6)=nd;
%     VV=V(:,inds(1:9));
%     VV1=V(:,inds_(1:9));
%     XD2=VV(:,1:3);
%     YD2=VV1(:,1:3);
    Is = obtainsign(XD2,YD2);
    XYS=(XD2+YD2*Is)*sqrt(0.5);
    %% Vi
    Vi = zeros(9,3);
    for i=1:3
        V_i=XYS(i*3-2:i*3,1:3);
%         di = det(V_i)
        [u,d_,v]=svd(V_i);
        temp=diag(d_);
        scale=(abs(temp(1))+abs(temp(2))+abs(temp(3)))/3;
        temp(1:3)=1;
        d1=diag(temp);
        E=u*d1*v';
        E=E*scale;
        Vi(i*3-2:i*3,1:3)=E;
    end
    Ui=(XD2-YD2*Is)*sqrt(0.5);
    XD22=(Ui+Vi)*sqrt(0.5);
    YD22=(Vi-Ui)*sqrt(0.5);
    %

    Dkt=[XD22 YD22]*diag(m_d2)*[XD22 YD22]';
        Dkt = (Dkt+Dkt')/2;
    for i=1:3
        Dkt(3*i-2:3*i , 3*i-2:3*i)=0;
    end
    astope = norm(Dkt, 'fro');
    if abs(astope-astopb)<0.001
        break;
    end
    DK1 = Dkt;    
    
end
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

function  Bkt = solveBK(Ekt,GAMMASt_1)
     res = Ekt - GAMMASt_1;
     [V,D]=eig(res);
dd=diag(D);
[~,inds]=sort(dd,'descend');
d_d=dd(inds,:);
UU=V(:,inds);
%
temp=d_d;
temp1(1)=(temp(1)-temp(9))/2;temp1(2)=(temp(2)-temp(8))/2;temp1(3)=(temp(3)-temp(7))/2;
temp1(4)=0;temp1(5)=0;temp1(6)=0;
temp1(7)=(temp(7)-temp(3))/2;temp1(8)=(temp(8)-temp(2))/2;temp1(9)=(temp(9)-temp(1))/2;
D=diag(temp1);
%
Bkt=UU*D*UU';
end

function [ma,tr]=getTrippleInds(BigMatrix, tripleInds )
tr=[tripleInds(1)*3-2:tripleInds(1)*3 tripleInds(2)*3-2:tripleInds(2)*3 tripleInds(3)*3-2:tripleInds(3)*3];
ma=BigMatrix(tr,tr);
end
