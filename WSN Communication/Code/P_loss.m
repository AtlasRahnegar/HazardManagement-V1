function loss_ratio = P_loss(Dist,T)

minV = interp1(T(:,1),T(:,2),Dist,'linear','extrap');
maxV = interp1(T(:,1),T(:,3),Dist,'linear','extrap');

loss_ratio = (1/100)*(minV + (maxV-minV).*rand);
