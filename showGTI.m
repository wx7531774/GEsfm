%%% show ground truth camera
clear


dataset = 'D:/Convex Optimization orientation/data/test_building';
RT = load(fullfile(dataset,'/building.txt'));
nodesNum=length(RT(:,1))/4;
figure(4);
for i=1:nodesNum
    RTI=RT(i*4-3:i*4,1:3);
    R=RTI(1:3,1:3);
    T=RTI(4,1:3);
    Rti=[R,T'];
    scale = 0.0002;
    drawCamera1(Rti, 20, 15, 15, scale, 1); axis equal;  
end