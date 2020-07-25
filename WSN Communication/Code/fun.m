function d = fun(x,y,net,CH,SX,SY)
% This is a function to calculate distance between two nodes
if isempty(net(2,CH)) % very rare
    d=sqrt((x - SX).^2 + (y - SY).^2);
else
    idx = (net(2,:)==x) & (net(3,:)==y);
    cluster = net(1,:) == net(1,idx);
    idx = cluster & CH;
    
    d=sqrt((x - net(2,idx)).^2 + (y - net(3,idx)).^2);
end

if isempty(d),d=0;end