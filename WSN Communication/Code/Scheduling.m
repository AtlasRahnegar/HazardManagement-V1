function ON = Scheduling(N,net,Rsense,USE,D)

ON = true(1,size(net,2)); % initialize all nodes to be ON-duty
if USE>=1
    % randomely order nodes for checking eligability rule
    % to simulate the random different delay for each node (like in paper)
    % i like rendom counter
    ss = 1:N; ss = ss(randperm(N,uint16(N*0.1)));%random order
    
    for i=ss
        [Xn, Yn] = Make_cir(net(2,i),net(3,i),Rsense); % sensing area of node in hand
        
        % Alive ON Neighbouring rang of this node
        Dist = sqrt(((net(2,:)-net(2,i)).^2) + ((net(3,:)-net(3,i)).^2));% cal dis R
        Snbr = find((Dist <= Rsense) & (Dist > 0) & ~D & ON);%the node must on,within the Radias and  not detect itself
        
        if ~isempty(Snbr)
            % The area neighbours cover 重疊 Xcov1, YCov1
            [Xcov1, Ycov1] = Make_cir(net(2,Snbr(1)),net(3,Snbr(1)),Rsense);
            poly1 = polyshape(Xcov1, Ycov1);
            
            for j=2:length(Snbr)
                %                 Coverage for each neighbor node
                [X,Y] = Make_cir(net(2,Snbr(j)),net(3,Snbr(j)),Rsense);
                poly2 = polyshape(X, Y);
                poly1 = union(poly1,poly2);
                
                %                 [Xcov1, Ycov1] = polybool('union', Xcov1, Ycov1, X, Y);
                
                %                 Coverage for 2nd degree neighbor nodes
                if USE>1
                    Dist2 = sqrt(((net(2,:)-net(2,Snbr(j))).^2) + ((net(3,:)-net(3,Snbr(j))).^2));
                    Snbr2 = find((Dist2 <= Rsense) & (Dist2 > 0) & ~D & ON);%the node must on,within the Radias and  not detect itself
                    if ~isempty(Snbr2)
                        for k=1:length(Snbr2)
                            [X,Y] = Make_cir(net(2,Snbr2(k)),net(3,Snbr2(k)),Rsense);
                            %                             try
                            poly2 = polyshape(X, Y);
                            poly1 = union(poly1,poly2);
                            
                            %                             [Xcov1, Ycov1] = polybool('union', Xcov1, Ycov1, X, Y);
                            %                             catch e
                            %                             end
                        end
                    end
                end
            end
            
            % The area lost if this node is off ?法 polyarea 做6角形
            poly2 = polyshape(Xn, Yn);
            poly1 = subtract(poly1,poly2);
            %             [Xuncov, Yuncov] = polybool('subtraction', Xn, Yn, Xcov1, Ycov1);
            %             A = polyarea(Xuncov, Yuncov);
            A = area(poly1) ;
            % Check eligability to be off 圓面積比較
            if A/(pi*Rsense^2)<0.4 %0.4 lower you have more nodes
                ON(i) = false;
            end
            
            
        end
    end
end



