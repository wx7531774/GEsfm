function saveGesfm(rt,strfile)
 n = length(rt);
 pose = zeros(n,12);
 for i=1:n
     rti = rt{i,1};
     pose(i,:) = [rti(1,1:3),rti(2,1:3),rti(3,1:3),rti(1:3,4)'];     
 end
%%
save(strfile,'pose','-ascii', '-double')
end