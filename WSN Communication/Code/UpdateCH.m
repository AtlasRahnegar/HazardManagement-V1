function [CH_F,net,Num_itr,FirstDegreeF,SecondDegreeF] = UpdateCH(R,D,p,pMin,E,Emax,net,cost,USE_MCH)

% R : range of a cluster
% D : State of nodes (Dead or alive)
% p : Cprop
% pMin : Minimum value of CHprop
% E    : Residual energy in all nodes
% Emax : Maximum energy each node can hold
% net  : Matrix holding position of each node and their clustering
% cost : cost for each node to be a cluster head

% CH : vector showing the final Cluster Heads
N = size(net,2); % number of nodes

% initializing ��l��
CH_prop = max([p*(E./Emax);pMin*ones(1,N)]);
clust_idx = zeros(1,N); % clusters indices
CH_F = false(1,N);    % Final CH
CH_tent = false(1,N); % Tentative CH
Num_itr = zeros(1,N); % Overhead iterations
Snbr1 = zeros(1,N);
FirstDegreeF = false(1,N);
SecondDegreeF = false(1,N);
% Clustering
for i=1:N
    FirstDegree{i}=zeros(1,N);
    SecondDegree{i}=zeros(1,N);
    
    if ~D(i) % to make sure node is not dead
        % Clustering range of this node
        Dist = sqrt(((net(2,:)-net(2,i)).^2) + ((net(3,:)-net(3,i)).^2));
        Snbr = Dist <= R(i);
        if USE_MCH
            FirstDegree{i} = Snbr;
            for iter1=1:N
                if Snbr(iter1)==1 && i~=iter1
                    Dist2 = sqrt(((net(2,:)-net(2,iter1)).^2) + ((net(3,:)-net(3,iter1)).^2));
                    Snbr2 = Dist2 <= R(iter1);
                    Snbr=Snbr|Snbr2;
                    Snbr1=Snbr1|Snbr2;
                end
            end
            
            
            SecondDegree{i}=Snbr1;
        end
        CH_prev = 0;
        
        % Repeat
        while CH_prev~=1
            tmp = rand; Num_itr(i) = Num_itr(i)+1;
            if sum((CH_tent&Snbr&(~D)))>0
                tmp = cost; tmp(~(CH_tent&Snbr&(~D)))= inf;
                [~,my_CH] = min(tmp);
                clust_idx(i) = my_CH;
                if my_CH==i
                    if CH_prop(i)==1
                        CH_F(i)=true;
                    else
                        CH_tent(i)=true;
                    end
                end
            elseif CH_prop(i)==1
                CH_F(i)=true;
                clust_idx(i)=i;
            elseif tmp <= CH_prop(i)
                CH_tent(i)=true;
            end
            CH_prev = CH_prop(i);
            CH_prop(i) = min(2*CH_prop(i),1);
        end
        
        % finalize
        if ~CH_F(i)
            if sum((CH_tent&Snbr&(~D)))>0
                tmp = cost; tmp(~(CH_tent&Snbr&(~D)))= inf;
                [~,my_CH] = min(tmp);
                clust_idx(i) = my_CH;
                
            else
                CH_tent(i)=true;
                CH_F(i)=true;
                clust_idx(i)=i;
            end
        end
    end
end
% for i=1:N
%     if ~D(i)
%         for iter2nd=1:length(FirstDegree)
%             tmpp=SecondDegree{i};
%             tmpp(FirstDegree{iter2nd})=false;
%             SecondDegree{i}=tmpp;
%         end
%     end
% end
if USE_MCH
    for i=1:N
        FirstDegreeF=FirstDegreeF+FirstDegree{i};
        SecondDegreeF=SecondDegreeF+SecondDegree{i};
    end
end

net(1,:) = clust_idx;
CH_F(D) = false;
FirstDegreeF(CH_F)=0;
SecondDegreeF(CH_F)=0;