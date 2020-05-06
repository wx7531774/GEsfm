function [finalTriplets,v,testPassed,origgTime,currentTime] = buildTripletsAviodCollinearFast(pointMatchesInliers,EN,Rt,width,hight)

numcam=size(EN,2)/3;
epipolsx=zeros(numcam,numcam);
epipolsy=zeros(numcam,numcam);
inliersNum= pointMatchesInliers;
adjec=zeros(numcam);
for i=1:numcam-1
    for j=i+1:numcam
        if  inliersNum(i,j)>0
            adjec(i,j)=1/inliersNum(i,j);
            adjec(j,i)=1/inliersNum(i,j);
        end
    end
end
G=graph(adjec);
center=[width,hight]'/2;

T=minspantree(G);
tripletGraph=adjacency(T);
graphsTrips=cell(4,1);
graphsTrips{1}=graph(tripletGraph);
try
    for i=1:5
        X=setdiff(G.Edges,T.Edges);
        G=graph(X.EndNodes(:,1) , X.EndNodes(:,2),X.Weight);
        T=minspantree(G);
        tripletGraph=tripletGraph+adjacency(T);
        graphsTrips{i+1}=graph(tripletGraph);
    end
catch
    
end
[c]=extractTripletsFromViewingGraph(graph(tripletGraph));

tts=zeros(size(c,1),1);% collinearity
tripletsErrors=zeros(size(c,1),2);%%det rotation / translation error
for i=1:size(c,1)
    rt12=Rt{c(i,1),c(i,2)};
    rt13=Rt{c(i,1),c(i,3)};
    rt21=Rtij2ji(rt12);
    rt23=Rt{c(i,2),c(i,3)};
    rt31=Rtij2ji(rt13);
    rt32=Rtij2ji(rt23);
    tts(i)=getCollinearityMeasurement( rt12,rt13,rt21,rt23,rt31,rt32 );
   %%
    tripletsErrors(i,1)=getRotationMeasurement( rt12,rt13,rt21,rt23,rt31,rt32 );
    tripletsErrors(i,2)=getTranslationMeasurement( rt12,rt13,rt21,rt23,rt31,rt32 );
end
%%% xiaot
% interId =intersect( intersect(find(tts>rad2deg(0.14)), find(tripletsErrors(:,1)<1.1)),  find(tripletsErrors(:,2)<rad2deg(1)) ); 
% c=c(interId ,:);%% 0.17 rad
% tripletsErrors=tripletsErrors(interId,:);
% tts=tts(interId,:);

adjec=zeros(size(c,1));
for i=1:size(c,1)-1    
    for j=i+1:size(c,1)        
        tuple1=c(i,:);
        tuple2=c(j,:);        
        numequals=0;
        indd1=1;
        indd2=1;        
        for k=1:4
            if tuple1(indd1)==tuple2(indd2)
                numequals=numequals+1;
                indd1=indd1+1;
                indd2=indd2+1;
            elseif  tuple1(indd1)>tuple2(indd2)
                indd2=indd2+1;
            else
                indd1=indd1+1;
            end
            if indd1>3 || indd2>3
                break;
            end            
        end
        if numequals>=2
            adjec(i,j)=1;
            adjec(j,i)=1;
        end        
    end
end
Gt=graph(adjec);

%  tts2=( tripletsErrors(:,2)).^(-1);%%xiaot
tts2= tts / tripletsErrors(:,2);%%xiaot

[ firstGroup,Gt]=makeMinimalGraph([c,repmat(tts2,1,3)],1:size(c,1),Gt,numcam);

testPassed=true;

finalTriplets=c(firstGroup,:);
v = bfsearch(Gt,1,{'edgetonew'});

end