
% Preallocation for energy calculations
E = Ei*ones(1,N);          % Energy left in each node
EH = zeros(1,num_rounds);
% Preallocation for dead nodes calculations
Ecrit = 0;                 % Critical energy left in node to call it alive
Ncrit = fix((95/100)*N);   % Critical number for dead nodes to stop simulation
Dead = false(1,N);         % Status of all nodes 0:Alive 1:Dead
DeadH = zeros(1,num_rounds);
Cover = zeros(1,num_rounds);
% Preallocation for Bits sent calculations
BitsH = zeros(1,num_rounds);
BitsLost = zeros(1,num_rounds);
AED = zeros(1,num_rounds);
BitsNode = zeros(1,N);
TotSentCH = 0; TotRecCH = 0;
TotSentSink = 0; TotRecSink = 0;

figure(1); set(gcf,'Position',[68,25,1347,875]);

% Simulating for each round
for r=1:num_rounds % iterating on each round
    
    % Updating alive status to all nodes
    Dead(E<=Ecrit) = true; E(Dead) = 0; DeadH(r)=sum(double(Dead));
    
    %%%% Scheduling %%%%
    ON = Scheduling(N,net,Rsense,USE_schedule,Dead); % Apply scheduling algorithm
    
    %%%% Choosing Clusters heads %%%%
    [CH,net,Num_itr,FirstDegree,SecondDegree] = UpdateCH(ranges,Dead,p,pMin,E,Ei,net,cost,USE_MCH);
    
    %%%%%%%% Energy calculations %%%%%%%
    EH(r) = sum(E); %get all energy left in all nodes
    %%% first : Clustering stage %%%
    E = E - ((NonCHpl*(Etrans+Eagg)).*Num_itr);
    %%% second : Steady state stage %%%
    % Cluster heads energy calculations
    numClust = sum(double(CH));
    D = sqrt((net(2,CH) - SX).^2 + (net(3,CH) - SY).^2);
    if USE_MCH
        E(CH) = E(CH) - (((Etrans+Eagg)*CHpl)+(Efs*CHpl*(D.^ 2)));%+(mean(NonCHpl)*Erec*round(N/numClust)));
        E(CH) = E(CH) - FirstDegree(CH).*Erec*mean(NonCHpl);
    else
        E(CH) = E(CH) - (((Etrans+Eagg)*CHpl)+(Efs*CHpl*(D.^ 2))+(mean(NonCHpl)*Erec*round(N/numClust)));
    end
    % Rest of nodes energy calculations
    tmp = net(2:3,~CH&~Dead&ON);
    rest = size(tmp,2);
    mD = zeros(1,rest);
    for i=1:rest, mD(i) = fun(tmp(1,i),tmp(2,i),net,CH,SX,SY); end
    if USE_MCH
        E(FirstDegree>1) = E(FirstDegree>1) - (mean(NonCHpl)*(Eagg));
    end
    E(~CH&~Dead&ON) = E(~CH&~Dead&ON) - ((mean(NonCHpl)*Etrans) + (Efs*CHpl*(mD.^2)) + ((Erec+Eagg)*CHpl));
    
    %%% sent bits considering bucket loss %%%%
    for i=1:N
        if CH(i)==true
            Di = sqrt(((net(2,i)-SX).^2) + ((net(3,i)-SY).^2));
            delay = (3.3e-9)*Di - (1.7e-8);
            AED(r) = AED(r) + delay;
            tmp = P_loss(Di,Ta);
            BitsNode(i) = round((1-tmp)*CHpl);
            BitsLost(r) = BitsLost(r) + round(tmp*CHpl);
            TotRecSink = TotRecSink + BitsNode(i);
        else
            idx2 = find(((net(1,:)==net(1,i)) & CH) & ~Dead & ON);
            Di = sqrt(((net(2,i)-net(2,idx2)).^2) + ((net(3,i)-net(3,idx2)).^2));
            if isempty(Di)
                BitsNode(i) = 0;
            else
                tmp = P_loss(Di,Ta);
                BitsNode(i) = round((1-tmp)*NonCHpl(i));
                BitsLost(r) = BitsLost(r) + round(tmp*NonCHpl(i));
                TotRecCH = TotRecCH + BitsNode(i);
            end
        end
    end
    if r==1
        BitsH(r) = sum(BitsNode);
    else
        BitsH(r) = BitsH(r-1) + sum(BitsNode);
    end
    TotSentCH   = TotSentCH + sum(NonCHpl(~CH&~Dead&ON));
    TotSentSink = TotSentSink + numClust*CHpl;
    BitsLost(r) = BitsLost(r)/(TotSentCH+TotSentSink);
    AED(r) = AED(r)/sum(CH);
    
    % Coverability calculation
    netAlive = net(:,~Dead&ON);
    if isempty(netAlive)
        Cover(r) = 0;
    else
        [Xcov1,Ycov1] = Make_cir(netAlive(2,1),netAlive(3,1),Rsense);
        for i=2:size(netAlive,2)
            [X,Y] = Make_cir(netAlive(2,i),netAlive(3,i),Rsense);
            [Xcov1, Ycov1] = polybool('union', Xcov1, Ycov1, X, Y);
        end
        [Xuncov, Yuncov] = polybool('subtraction', pgonx, pgony, Xcov1, Ycov1);
        idxnan = find(isnan(Xuncov)); idxnan = [0,idxnan,length(Xuncov)+1];
        A = 0;
        for i=1:length(idxnan)-1
            tmpx = Xuncov(idxnan(i)+1:idxnan(i+1)-1);
            tmpy = Yuncov(idxnan(i)+1:idxnan(i+1)-1);
            A = A + polyarea(tmpx, tmpy);
        end
        Cover(r) = 1 - (A/(W*L));
    end
    
    %%%% Showing updated net %%%%
    %net = DrawNet(net,N,CH,Dead,SX,SY,Xcov1,Ycov1,ON);
    title(['CH : Red  ----  Dead : Empty  ---  round (',num2str(r),')']);
    drawnow
    if sum(Dead)>=Ncrit,break;end % Stop simulation when 5% or less is alive
end
close(1)
figure(1); set(gcf,'Position',[664,140,1170,775]);
%net = DrawNet(net,N,CH,Dead,SX,SY,Xcov1,Ycov1,ON);
title(['CH : Red  ----  Dead : Empty  ---  round (',num2str(r),')']);
T = (Tsetup+Tss)*(0:r-1);
EH = EH(1:r); EHdis = EH;
DeadH = DeadH(1:r); AliveH = N-DeadH;
BitsH = BitsH(1:r);
BitsLost = BitsLost(1:r);
AED = AED(1:r);
Cover = Cover(1:r);

PDR_CH = TotRecCH/TotSentCH;
PDR_Sink = TotRecSink/TotSentSink;
% disp('')
% disp('Analysis:')
% disp('=========================')
% disp(['Packet delivery ratio to the cluster heads = ',num2str(100*PDR_CH),'%'])
% disp(['Packet delivery ratio to sink = ',num2str(100*PDR_Sink),'%'])
% disp('')

% Plotting analysis of network preformance

figure(2); hold on
set(gcf,'Position',[131 59 792 613]);
plot(T,EHdis,'-x');
xlabel('Time (s)'); ylabel('Energy (j)')
title('Residual energy');hold off
if USE_schedule==1,legend('Nodes always ON','Using Scheduling'); end
if USE_schedule==2,legend('Nodes always ON','Using Scheduling','Using Muti hop'); end
if USE_MCH==1,legend('Nodes always ON','Using Scheduling','Using Muti hop','Using Muti hop MCH'); end

if USE_schedule==0, saveas(gcf,'Residual energy.png'); end
if USE_schedule==1, saveas(gcf,'Residual energy-Single.png'); end
if USE_schedule==2, saveas(gcf,'Residual energy-Multi.png'); end
if USE_MCH==1, saveas(gcf,'Residual energy-MultiMCH.png'); end
figure(3); hold on
set(gcf,'Position',[298 66 792 613]);
plot(N-AliveH,'-x');
xlabel('Round'); ylabel('No of dead Sensors nodes')
title('Life time of sensor nodes');hold off
if USE_schedule==1,legend('Nodes always ON','Using Scheduling'); end
if USE_schedule==2,legend('Nodes always ON','Using Scheduling','Using Muti hop'); end
if USE_MCH==1,legend('Nodes always ON','Using Scheduling','Using Muti hop','Using Muti hop MCH'); end

if USE_schedule==0, saveas(gcf,'Life time.png'); end
if USE_schedule==1, saveas(gcf,'Life time-Single.png'); end
if USE_schedule==2, saveas(gcf,'Life time-Multi.png'); end
if USE_MCH==1, saveas(gcf,'Life time-MultiMCH.png'); end

figure(4); hold on
set(gcf,'Position',[511 61 792 613]);
plot(T,BitsH,'-x');
xlabel('Time (s)'); ylabel('Overall Throughput (no of packets)')
title(['Throughput (' num2str(N) ' Sensor nodes)']);hold off
if USE_schedule==1,legend('Nodes always ON','Using Scheduling'); end
if USE_schedule==2,legend('Nodes always ON','Using Scheduling','Using Muti hop'); end
if USE_MCH==1,legend('Nodes always ON','Using Scheduling','Using Muti hop','Using Muti hop MCH'); end

if USE_schedule==0, saveas(gcf,'Throughput.png'); end
if USE_schedule==1, saveas(gcf,'Throughput-Single.png'); end
if USE_schedule==2, saveas(gcf,'Throughput-Multi.png'); end
if USE_MCH==1, saveas(gcf,'Throughput-MultiMCH.png'); end

figure(5); hold on
set(gcf,'Position',[711 61 792 613]);
plot(Cover,'-x');
xlabel('Round'); ylabel('Coverability')
title('Ratio of covered area to total area');hold off
if USE_schedule==1,legend('Nodes always ON','Using Scheduling'); end
if USE_schedule==2,legend('Nodes always ON','Using Scheduling','Using Muti hop'); end
if USE_MCH==1,legend('Nodes always ON','Using Scheduling','Using Muti hop','Using Muti hop MCH'); end

if USE_schedule==0, saveas(gcf,'Coverage.png'); end
if USE_schedule==1, saveas(gcf,'Coverage-Single.png'); end
if USE_schedule==2, saveas(gcf,'Coverage-Multi.png'); end
if USE_MCH==1, saveas(gcf,'Coverage-MultiMCH.png'); end
