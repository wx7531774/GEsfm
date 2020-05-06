function DrawCamera2(Rc, w, h, f, scale, lineWidth)
% R is the rotation matrix in cv
% i.e R' is the orientation matrix in a global frame 
% 

% Xcamera = Rt * Xworld
%  f = f;
V= [...
0 f/2 0 0 -w/2  w/2 w/2 -w/2
0 0 f/2 0 -h/2 -h/2 h/2  h/2
0 0 0 f/2  f    f    f   f];

V = V*scale;

% V = transformPtsByRt(V, Rt, false);
% V vector in camera frame to V vector in world frame
R = Rc(1:3, 1:3);
C = Rc(1:3,4);
V = R' *V +repmat(C,1,size(V,2));


hold on;
plot3(V(1,[1 4]),V(2,[1 4]),V(3,[1 4]),'-b','LineWidth',lineWidth);
plot3(V(1,[1 3]),V(2,[1 3]),V(3,[1 3]),'-g','LineWidth',lineWidth);
plot3(V(1,[1 2]),V(2,[1 2]),V(3,[1 2]),'-r','LineWidth',lineWidth);
%
plot3(V(1,[1 5]),V(2,[1 5]),V(3,[1 5]),'-r','LineWidth',lineWidth);
plot3(V(1,[1 6]),V(2,[1 6]),V(3,[1 6]),'-r','LineWidth',lineWidth);
plot3(V(1,[1 7]),V(2,[1 7]),V(3,[1 7]),'-r','LineWidth',lineWidth);
plot3(V(1,[1 8]),V(2,[1 8]),V(3,[1 8]),'-r','LineWidth',lineWidth);

plot3(V(1,[5 6 7 8 5]),V(2,[5 6 7 8 5]),V(3,[5 6 7 8 5]),'-r','LineWidth',lineWidth);