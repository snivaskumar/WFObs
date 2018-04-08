function [x,d,p, F,D,E,G,H,Q,R, tur,n, x_est,x_unest, P_unest] = subsystem_turbine(strucObs,sol_in, p,Fk,Bk,C,QQ,RR, tur,state,turbLocArray, Subsys_length,RD, Sk1k1);
% [x,d, F,D,G,H,Q,R,y, tur,n, x_est,x_unest, P_unest] = subsystem_turbine(Fk,Bk,C,QQ,RR,yy, tur,stateLocArray,turbLocArray, Sk1k1);

% Model Decomposition
% tic
% aa  = double(abs(Fk)>1e-4);
% p   = symrcm(aa);
% A   = Fk(p,p);
% Bk  = Bk(p,:);
% Ck  = C(:,p);
% state = stateLocArray(p,:);
A       = Fk;
Bk      = Bk;
Ck      = C;
state   = state;
% toc
n = length(Fk);

if sol_in.k == 1
    turbine = turbLocArray;
    d = cell(tur,1);
    % tic
    parfor j = 1:tur
        for i = 1:n
            d{j}(i) = sqrt( (state(i,1)-turbine(j,1))^2 + (state(i,2)-turbine(j,2))^2 );
        end
    end
    % toc
    x = cell(tur,1);
    x_pos = cell(tur,1);

    RD;
    if Subsys_length <= 5
        sub_len = Subsys_length;
        Subsys_length = Subsys_length*RD;
    else
        Subsys_length = Subsys_length;
    end
    Subsys_length;

    % tic
    parfor j = 1:tur
        for i = 1:n
            if d{j}(i)<= (Subsys_length)
                x{j} = [x{j},i];
                x_pos{j} = [x_pos{j};state(i,:)];
            end
        end
        x{j} = x{j}';
    end
    % toc

    opu = [strucObs.obs_array_locu.x;strucObs.obs_array_locu.y]';
    opv = [strucObs.obs_array_locv.x;strucObs.obs_array_locv.y]';
    
%     figure, plot(state(:,2),state(:,1),'sb') 
%     for i = 1:tur
%         hold on, plot(turbine(i,2),turbine(i,1),'s','LineWidth',3)
%         hold on, plot(x_pos{i}(:,2),x_pos{i}(:,1),'*')
%     end
%     for i = 1:length(opu)
%         hold on, plot(opu(i,2),opu(i,1),'*k','LineWidth',2)
%     end
%     for i = 1:length(opv)
%         hold on, plot(opv(i,2),opv(i,1),'*r','LineWidth',2)
%     end
%     hold off
%     xlabel('y-direction'),ylabel('x-direction');
%     titl = sprintf('Decomposing the flow field into sub-systems (%dD)',sub_len);
%     title(titl),legend('Flow fields','Turbine 1','Sub-System 1','Turbine 2','Sub-System2');    
    
    x_est = x{1};
    for i = 2:tur
        x_est = union( x_est, x{i} );
    end
    x_unest = setdiff([1:n],x_est)';
    size(A);
    length(x_est);
    length(x_unest);

    d = cell(tur,1);
    % tic
    for i = 1:tur
        clear rtmp ctmp tmp
        tmp         = setdiff(x_est,x{i});
        [rtmp ctmp] = find( A(x{i},tmp) );
        ctmp        = unique(ctmp);
        d{i}        = [d{i}, tmp(ctmp)];
    end
    % toc
    for i = 1:tur
        d{i} = unique(d{i});
    end
    strucObs.subsystem.x = x;           strucObs.subsystem.d = d;
    strucObs.subsystem.x_est = x_est;   strucObs.subsystem.x_unest = x_unest;
end
x       = strucObs.subsystem.x;         d       = strucObs.subsystem.d;
x_est   = strucObs.subsystem.x_est;     x_unest = strucObs.subsystem.x_unest;
%% Sub-Systems

F   = cell(tur,1);                 % Local A
H   = cell(tur,1);                 % Local C
D   = cell(tur,1);                 % Local Internal Input 
E   = cell(tur,1); 
G   = cell(tur,1);                 % Local External Input
Q   = cell(tur,1);
R   = cell(tur,1);
y   = cell(tur,1);
% tic
for i = 1:tur
    F{i}    = A( x{i},x{i} );
    D{i}    = A( x{i},d{i} );
    E{i}    = A( x{i},x_unest );
    if strucObs.superOptimize == 1
        factor          = strucObs.superOptimizeFactor;
        indices         = find( abs(E{i})<factor );
        E{i}(indices)   = 0;
        E{i}            = sparse(E{i});
% % %         len_x       = length(x{i});
% % %         len_xunest  = length(x_unest);
% % %         E{i}        = sparse(zeros(len_x,len_xunest));
% % %         E{i}(1:min(len_x,len_xunest),1:min(len_x,len_xunest)) = sparse(diag(diag(A( x{i},x_unest ))));
    end
    G{i}    = Bk( x{i},: );
    H{i}    = Ck( :,x{i} );
    R{i}    = RR( :,: );
%     y{i}    = yy( : );
    Q{i}    = QQ( x{i},x{i} );
end
% toc
P_unest = [];
% P_unest = A( x_unest,x_unest )*Sk1k1( x_unest,x_unest )*A( x_unest,x_unest )' ...
%                              + A( x_unest,x_unest )*Sk1k1( x_unest,x_est )*A( x_unest,x_est )' ...
%                              + (A( x_unest,x_unest )*Sk1k1( x_unest,x_est )*A( x_unest,x_est )')' ...
%                              + A( x_unest,x_est )*Sk1k1( x_est,x_est )*A( x_unest,x_est )' ...
%                              + QQ( x_unest,x_unest );

end