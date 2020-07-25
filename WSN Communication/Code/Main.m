% Clean memory and command window
clear,clc,close all,warning('off')

%% Parameters
R = 5;              % Max range for wireless transmition for each node
N = 50;             % Number of nodes
W = 30;              % length of the network
L = 30;              % width of the network
Ei = 3;              % Initial energy of each node (joules)
Rsense = 10;          % radius for sensing area around each node

rates = 80;          % Packet size for normal node per round (bits) for each type
CHpl = 3000;         % Packet size for cluster head per round (bits)
p = 5/100;           % desired percentage of cluster heads
num_rounds = 2000;   % Max Number of simulated rounds
Tsetup = 4;          % average Time in seconds taken in setup phase
Tss = 10;            % average Time in seconds taken in steady state phase
Etrans = 1.0000e-09; % Energy for transmitting one bit
Erec = 1.0000e-09;   % Energy for receiving one bit
Eagg = 1.0000e-09;   % Data aggregation energy
Efs = 5.00e-8;       % Energy of free space model amplifier
Marg = 15;           % Distance from the sink to sensors area


Ta = [10 0 2;20 0 4.459;30 0 3.75;40 1.38 7.14;50 2.09 20.37
      60 2.54 10.18;70 1.86 10;80 2.77 30;90 5.9 30.95;100 23.33 60.74];
m = 1;
pMin = 10^-4;        % Lowest possible CH_prop
if length(rates)~=m || length(R)~=m
    error('The length of rates and ranges is not equal to m please check parameters')
end

%% Building the WSN

% Position of sink
SX = -Marg; SY = L/2;

% 1st row: Clustering indexing
% 2nd: x-position, 3rd: y-position
net = [zeros(1,N) ; rand([1,N])*W ; rand([1,N])*L ; randi(m,[1,N])];
tmp = sqrt((net(2,:)-SX).^2 + (net(3,:)-SY).^2); 
[~,I] = sort(tmp,'descend'); 
net = net(:,I);
NonCHpl = zeros(1,N); ranges = zeros(1,N);
for i=1:m
    idx = net(end,:)==i;
    NonCHpl(idx) = rates(i);
    ranges(idx) = R(i);
end

% area to be sensed 要感測的區域
pgonx = [0 0 W W];
pgony = [0 L L 0];

% calculating costs
cost = zeros(1,N);
for i=1:N
    Dist = sqrt(((net(2,:)-net(2,i)).^2) + ((net(3,:)-net(3,i)).^2));
    d = sqrt(((SX-net(2,i)).^2) + ((SY-net(3,i)).^2));
    Snbr = Dist <= ranges(i);
    for iter1=1:N
            if Snbr(iter1)==1 && i~=iter1
                Dist2 = sqrt(((net(2,:)-net(2,iter1)).^2) + ((net(3,:)-net(3,iter1)).^2));
                Snbr2 = Dist2 <= ranges(iter1);
                Snbr=Snbr|Snbr2;
            end
    end
        
    cost(i) = sum(Dist(Snbr))/(sum(Snbr)-1)+d;
    % special case of isolated node
    if isnan(cost(i)),cost(i)=max(Dist);end 
end

%% Start communication simulation
USE_MCH=0;
% first run no scheduling
USE_schedule = 0;
simulate_WSN

% second run using scheduling
USE_schedule = 1;
simulate_WSN

USE_schedule = 2;
simulate_WSN

USE_MCH = 1;
USE_schedule = 2;
simulate_WSN
warning('on')

