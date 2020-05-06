function [finalTriplets,v,weight] = buildTripletsMinimal_c(pointMatchesInliers,EN,Rc,width,hight)
numcam=size(EN,2)/3;
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
if length(conncomp(G))>1
    nerr=1;
end
%
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
% [c]=extractTripletsFromViewingGraph(G);%%%FOR internet data
%
tts=zeros(size(c,1),1);% collinearity
tripletsErrors=zeros(size(c,1),3);%%det rotation / translation error
for i=1:size(c,1)
    rc12=Rc{c(i,1),c(i,2)};
    rc13=Rc{c(i,1),c(i,3)};
    rc21=Rcij2ji(rc12);
    rc23=Rc{c(i,2),c(i,3)};
    rc31=Rcij2ji(rc13);
    rc32=Rcij2ji(rc23);
    tts(i)=getCollinearityMeasurement(  rc12,rc13,rc21,rc23,rc31,rc32 );
   %%
    tripletsErrors(i,1)=getRotationMeasurement(  rc12,rc13,rc21,rc23,rc31,rc32 );
    tripletsErrors(i,2)=getTranslationMeasurement(  rc12,rc13,rc21,rc23,rc31,rc32 );
    tripletsErrors(i,3)=max(getRotationMeasurement(  rc12,rc13,rc21,rc23,rc31,rc32 ),getTranslationMeasurement(  rc12,rc13,rc21,rc23,rc31,rc32 ));
end

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

   %tts2=(tripletsErrors(:,1)).^(-1);%%xiaot
   %tts2=(tripletsErrors(:,2)).^(-1);%%xiaot
   tts2=(tripletsErrors(:,3)).^(-1);%%xiaot
   %tts2= tts ./(tripletsErrors(:,3));%%xiaot

[ firstGroup,Gt]=makeMinimalGraph([c,repmat(tts2,1,3)],1:size(c,1),Gt,numcam);
finalTriplets=c(firstGroup,:);
v = bfsearch(Gt,1,{'edgetonew'});
weight = [tts(firstGroup,:),tripletsErrors(firstGroup,:)];
end

function [ output_args ] = getTranslationMeasurement( rc12,rc13,rc21,rc23,rc31,rc32)
%GETCOLLINEARITYMEASUREMENT Summary of this function goes here
%   calculating the intersection angles of triplet's edges
c12 = rc12(:,4); 
c13 = rc13(:,4);
theta1 = acosd(dot(c12,c13)/norm(c12)/norm(c13));
%
c21 = rc21(:,4); 
c23 = rc23(:,4);
theta2 = acosd(dot(c21,c23)/norm(c21)/norm(c23));
%
c31 = rc31(:,4); 
c32 = rc32(:,4);
theta3= acosd(dot(c31,c32)/norm(c31)/norm(c32));
output_args=abs(theta1+theta2+theta3-180);
end

function [ output_args ] = getRotationMeasurement( rc12,rc13,rc21,rc23,rc31,rc32)
%GETCOLLINEARITYMEASUREMENT Summary of this function goes here
%   Detailed explanation goes here
r12 = rc12(:,1:3); 
r23 = rc23(:,1:3);
r31 = rc31(:,1:3);

errR = r12*r23*r31;
err_r=(errR(1,1)+errR(2,2)+errR(3,3)-1)/2.0;
if (errR(1,1)+errR(2,2)+errR(3,3)-1)/2.0<-1.00
    err_r=-1;
end

if(errR(1,1)+errR(2,2)+errR(3,3)-1)/2.0>1.00
    err_r=1;
end
errAngle_rot = rad2deg(acos(err_r));
if isreal(errAngle_rot)==0
    aqsa=0;
end
% errNorm = norm((errR-eye(3,3)),'fro');
output_args=errAngle_rot;
end

function [ output_args ] = getCollinearityMeasurement(rc12,rc13,rc21,rc23,rc31,rc32)
%GETCOLLINEARITYMEASUREMENT Summary of this function goes here
%   Detailed explanation goes here
c12 = rc12(:,4); 
c13 = rc13(:,4);
theta1 = acosd(dot(c12,c13)/norm(c12)/norm(c13));
%
c21 = rc21(:,4); 
c23 = rc23(:,4);
theta2 = acosd(dot(c21,c23)/norm(c21)/norm(c23));
%
c31 = rc31(:,4); 
c32 = rc32(:,4);
theta3= acosd(dot(c31,c32)/norm(c31)/norm(c32));
output_args=min([theta1,theta2,theta3]);
if output_args >60
    nerror=1;
end
end

function Rcji = Rcij2ji(Rcij)
rij = Rcij(:,1:3);
cij = Rcij(:,4);
rji = rij';
cji = -(rij * cij);%%xiaot
Rcji = [rji,cji];
end





