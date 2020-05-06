function [S,R,T ] = FindsimilarityTrans( Ps2,Ps1 )
%R{i} in Ps is the rotation matrix in cv
%  
% a0=[0,0,0]'; b0=[1,0,0]';c0=[0,1,0]';d0=[0,0,1]';
% e0=[0,0,0]'; f0=[1,0,0]';g0=[0,1,0]';h0=[0,0,1]';
% p=[0,0,0;1,0,0;0,1,0;0,0,1;-1,0,0;0,-1,0;0,0,-1]';
% 
% %
% r11 = Ps1{1}(1:3,1:3);
% t11 = Ps1{1}(1:3,4);
% r12 = Ps1{2}(1:3,1:3);
% t12 = Ps1{2}(1:3,4);
% u1=r11' *p +t11;%%[a1,b1,c1,d1]
% u2 =r12'*p +t12;%%[e1,f1,g1,h1] 
% %
% r21 = Ps2{1}(1:3,1:3);
% t21 = Ps2{1}(1:3,4);
% r22 = Ps2{2}(1:3,1:3);
% t22 = Ps2{2}(1:3,4);
% v1=r21' *p +t21;%%[a2,b2,c2,d2]
% v2 =r22'*p +t22;%%[e2,f2,g2,h2]
% % [sR,t]=umeyama([v1,v2],[u1,u2],true);
% % tform = pcregistericp(pointCloud([v1,v2]'),pointCloud([u1,u2]'),'Extrapolate',true);
% [regParams,Bfit,ErrorStats]=absor([v1,v2],[u1,u2],'doScale',1);
% S = regParams.s;
% R = regParams.R;
% T = regParams.t;
% %
% s1 = norm(t12-t11);
% s2 = norm(t22-t21);
% s = s1/s2;
% t21s = t21*s;
% t22s = t22*s;
% %;
% Rt = (inv(r11)*r21+inv(r12)*r22)/2;
% [U,D,V]=svd(Rt);
% Ra = U*V';
% %
% t11s = t11-Ra*t21s;
% t12s =t12- Ra*t22s;
% if max(abs(t11s-t12s)) > 0.1
%     nerr = 1;
% end
% l = [t11s;t12s];
% B=[1,0,0;0,1,0;0,0,1;1,0,0;0,1,0;0,0,1];
% t = inv(B'*B)*(B'*l);
% R=Ra;
% S =s;
% T=(t11s+t12s)/2;

R1L=Ps1{1}(1:3,1:3);
t1L=Ps1{1}(1:3,4);
R2L=Ps1{2}(1:3,1:3);
t2L=Ps1{2}(1:3,4);
%
R1R=Ps2{1}(1:3,1:3);
t1R=Ps2{1}(1:3,4);
R2R=Ps2{2}(1:3,1:3);
t2R=Ps2{2}(1:3,4);

Rst=inv(R1L)* R1R;
Rst1=inv(R2L)*R2R;
if max(abs(Rst-Rst1))>0.01
    nerror =1;
end
Rt = (Rst+Rst1)/2;
[U,D,V]=svd(Rt);
Ra = U*V';
%%
t2R=Ra*t2R;t1R=Ra*t1R;

A=[t1R(1),1,0,0;
   t1R(2),0,1,0;
   t1R(3),0,0,1;
   t2R(1),1,0,0;
   t2R(2),0,1,0;
   t2R(3),0,0,1;];
L=[t1L(1);
   t1L(2);
   t1L(3);
   t2L(1);
   t2L(2);
   t2L(3)];
AA=A'*A;AL=A'*L;
st=inv(AA)*AL;
v=A*st-L;
if max(abs(v)>0.1)
  nerr =1;
end
S =st(1,:);
if S<0
  nerr=1;
end
T=st(2:4,:);
if isnan(T)
    error=-1;
end
R=Ra;
end

