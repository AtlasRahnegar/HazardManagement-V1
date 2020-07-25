clear all
clc
close all
tic;
m=22;%areas
n=66;%facilities
t=3;%candidates
O=zeros(m,t);%optimization cost
maxCapacity = 25000;
maxVictims = 300000;
C=ones(n,1)*maxCapacity;%zeros(n,1);%capacity
V=(rand([m,1])+1)*maxVictims;%zeros(m,1);%number of victims
A=zeros(2*m+n,m*t);
B=zeros(2*m+n,1);
S=zeros(1,m*t);%index of selected candidates
S2=zeros(m,t);%index of selected candidates in another format
SC=zeros(m,n); %selection score
dirty=zeros(size(V));
epse=0;
idealVic=((t*maxCapacity)./(V))';
counter=1;
MaxIter=200;

l1x=516169.37+50000*((1:11)/11);
l1y=3960176.8+20000*((1:2)/2);
loca=[];%location of areas

for i=1:length(l1x)
    for j=1:length(l1y)
        temp=[l1x(i);l1y(j)]+(rand([2,1])-0.5)*1000;
        loca=[loca,temp];
    end
end

l2x=516169.37+50000*((1:11)/11);
l2y=3960176.8+20000*((1:6)/6);
locf=[];%location of facilities

for i=1:length(l2x)
    for j=1:length(l2y)
        temp=[l2x(i);l2y(j)]+(rand([2,1])-0.5)*1000;
        locf=[locf,temp];
    end
end

while eps<0.9 && counter<MaxIter
    for i = 1 : size(loca,2)
        for j = 1 : size(locf,2)
            SC(i,j) = sqrt((loca(1,i)-locf(1,j))^2+(loca(2,i)-locf(2,j))^2);
        end
        if counter>1 && dirty(i)~=1
            SC(i,dirtyF)=Inf;
        end
        [ss,ii]=sort(SC(i,:));
        S2(i,1:t)=ii(1:t);
        O(i,1:t)=ss(1:t);
    end
    
    S=reshape(S2',[1,m*t]);
    ideal = 3000/max(max(O(:,:)));
    O=O./max(max(O(:,:)));
    O=O./10;
    for i = 1 : m
        A(i,(i-1)*t+1:i*t)=O(i,:);
        B(i,1)=ideal*V(i,1);
        A(m+i,(i-1)*t+1:i*t)=ones(1,t);
        B(m+i,1)=idealVic(i)*V(i,1);
    end
    for i = 1 : n
        A(2*m+i,:) =(S==i);
        CC(i,:)=(S==i);
        B(2*m+i,:) = C(i,1);
        d(i,:) = C(i,1);
    end
    CC=double(CC);
    d=double(d);
    options = optimoptions('lsqlin','Display','off');

    X = lsqlin(A,B,[CC;A(m+1:2*m,:)],[d;V],[],[],zeros(m*t,1),[],[],options);
    XP=reshape(X,[t,m])';
    
    for i = 1 :n
        resC(i)=CC(i,:)*X-C(i);
    end
    for i = 1 :m
        resV(i)=sum(XP(i,:))/V(i);
    end
    sum(resV(resV<0)<-10);
    sum(resV>1);
        goodness = resV./idealVic;

    [dirtyval,dirtyareas]=sort(goodness,'descend');
    dirty(goodness>0.7)=1;
    dirtyF=unique(reshape(S2(dirty==1,:),[],1));
    eps=min(goodness);
    counter=counter+1;
    mean(goodness)
    size(find(resC<-1),2)
end
toc;
