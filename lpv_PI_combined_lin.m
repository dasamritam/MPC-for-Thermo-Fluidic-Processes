t_complete=tic;
%-------------------------------------------------------------------
% loads the state space matrices where the,
% control inputs := Q_N1, Q_N2, Q_N3, H1, H4
% disturbance := Q_in, Q_RD, Q_o1, Q_o2, Q_o3
%-------------------------------------------------------------------
A_th = Lin.A ;

%N_input = input('Enter 1 for 5 inputs and 2 for 2 input actuation ');
B_b=[Lin.B Lin.f];

%% Discretisation- 4th order Taylor series expansion
A_d = (1/24 * A_th^4 *Ts^4 ) + (1/6 * A_th^3 *Ts^3 )+ (1/2 * A_th^2 *Ts^2 ) + ( A_th*Ts ) + eye(size(A_th,1));
B_d = (1/24 * A_th^3*B_b *Ts^4 ) + (1/6 * A_th^2*B_b *Ts^3 )+ (1/2 * A_th^1*B_b *Ts^2 ) + ( B_b*Ts );


 
%% reference creation
nx = size(A_d,1);

nu=size(Lin.B,2);

total_Time = 0:Ts:T;

%%
dxh=zeros(nx,length(total_Time)+1); 
duh=zeros(nu,length(total_Time));
%error = zeros(nu,length(total_Time));
x0= X0*ones(nx,1);    % Assigning Initial condition
dxh(:,1) = x0-Tss;

%%
der = 0;              % derivative of error (time series)
if N_input == 2
    Kp = [22.44; 8];
    Ki = [2.155; 0.02];
    Kd = 0;% Kp*Tu/3;
else
    Kp = [22.44; 8 ;ones(N_nz,1)*0.005 ];
    Ki = [2.155; 0.02; ones(N_nz,1)*0.04];
    Kd = 0;% Kp*Tu/3;
end
aw = 10;

%%

m=9:2:9+(N_nz-1)*2;
midnoz = m(:,ceil(end/2));
F_d=B_d(:,nu+1:end);
B_d=B_d(:,1:nu);
%F_dd = ones(MatSize,1)*0.25905e-4;
aaa=0;
for k =1:1: T_sim
       
        
        % calculate error
        if N_input == 2 
            error2(k+aw) = -dxh(midnoz ,k);
        else
            xx=[dxh(midnoz,k); dxh(midnoz ,k); dxh(9,k)];
            for i=1:1:N_nz-1
                xx=[xx; dxh(9+i*2,k)];
            end
            error2(:,k+aw) =-xx;
        end
        
        % P term
            du_k = Kp.*error2(end);
        % I term
        if N_input == 2
            du_k = Ki.*sum(error2(1,k:k+aw),2) + du_k;
            error1(:,k)=sum(error2(1,k:k+aw),2);
        else
            du_k(1:2+N_nz) = Ki.*sum(error2(1:N_nz+2,k:k+aw),2) + du_k;
            error1(:,k)=sum(error2(1:N_nz+2,k:k+aw),2);
        end
          
        % limit PID controller output
            
            for i = 1:1:nu
                if du_k(i) > U_h_max-Us(i) && i<3
                    du_k(i) = U_h_max-Us(i);
                elseif du_k(i) < U_h_min-Us(i) && i<3
                    du_k(i) = U_h_min-Us(i);
                elseif du_k(i) > U_p_max-Us(i) && i>=3
                    du_k(i) = U_p_max-Us(i);
                elseif du_k(i) < U_p_min-Us(i)  && i>=3
                    du_k(i) = U_p_min-Us(i);
                end
            end
        % simulate process behaviour
        %du_k=uh(:,k)-Us;
        df_k=[Dist_Prof.Qin(:,k)-Qinss;Dist_Prof.QRD(:,k)-Qrdss;Dist_Prof.QN(:,k)-QNss;Dist_Prof.QN(:,k)-QOss;Dist_Prof.L3(:,k)-L3ss];

        dxh(:,k+1) = A_d*dxh(:,k) + B_d*du_k + F_d*df_k;%+F_dd;
        duh(:,k) = du_k;
        
    com=((((k-1))/T))*Ts;

        if floor(com*10) ~=aaa
            waitbar(floor(com*10)/10,f,'Simulating PI...');
            aaa=floor(com*10); 
        end   
        
end
xhlin=Tss+dxh;
uhlin=Us+duh(1:nu,:);

ttt=toc(t_complete);
disp(['Completion of simimulation takes ' num2str(ttt) ' seconds'])

%%

figure
hold on
ii =1:1:N_nz;

for i=ii
    plot(total_Time(1:T_sim),xhlin(7+i*2,1:T_sim))
end 

plot(total_Time(1:T_sim),T_ref*ones(length(total_Time(1:T_sim)),1),'black--'),
title('DFA Ink Temperature'),

xlabel('Time [s]'),ylabel('Temperature $^\circ C$')
figp
axis([0 T floor(min([xhlin(9,1:T_sim) xhlin(11,1:T_sim) xhlin(13,1:T_sim)]))-2 ceil(max([xhlin(9,1:T_sim) xhlin(11,1:T_sim) xhlin(13,1:T_sim) T_ref]))+2])

figure, plot(total_Time(1:T_sim),xhlin(1,1:T_sim),total_Time(1:T_sim),xhlin(5,1:T_sim))
hold on,plot(total_Time(1:T_sim),70*ones(length(total_Time(1:T_sim)),1),'black--'), 
title('Reservoir and Distributor Solid Left'),legend('$T^{R1}$','$T^{D1}$','$T_{Ref}$','location','best'),xlabel('Seconds [s]'),ylabel('Temperature $^\circ C$')
figp
axis([0 T floor(min([xhlin(1,1:T_sim) xhlin(4,1:T_sim) xhlin(5,1:T_sim) xhlin(7,1:T_sim)]))-2 ceil(max([xhlin(1,1:T_sim) xhlin(4,1:T_sim) xhlin(5,1:T_sim) xhlin(7,1:T_sim) T_ref]))+2])

if N_input == 1
figure, plot(total_Time(1:T_sim),uhlin(3,1:T_sim),total_Time(1:T_sim),uhlin(4,1:T_sim),total_Time(1:T_sim),uhlin(5,1:T_sim))
title('Heat Flux in Nozzles'),legend('HN1','HN2','HN3'),xlabel('Seconds [s]'),ylabel('Flux')
axis([0 T U_p_min U_p_max])
figp
end

figure, plot(total_Time(1:T_sim),uhlin(1,1:T_sim),total_Time(1:T_sim),uhlin(2,1:T_sim))
title('Heat Flux in solid blocks 1 and 4'),legend('$H^{1}$','$H^{4}$'),xlabel('Seconds [s]'),ylabel('Flux')
axis([0 T U_h_min U_h_max])
figp


%%
figureHandle=createfigure(total_Time(1:T_sim),xhlin(1,1:T_sim),xhlin(3,1:T_sim),xhlin(4,1:T_sim),...
     xhlin(5,1:T_sim),xhlin(6,1:T_sim),xhlin(7,1:T_sim),...
     xhlin(8,1:T_sim),xhlin(9,1:T_sim),xhlin(10,1:T_sim),...
     xhlin(11,1:T_sim),xhlin(12,1:T_sim),xhlin(13,1:T_sim),xhlin(14,1:T_sim),T);
