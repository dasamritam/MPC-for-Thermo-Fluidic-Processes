t_complete=tic;
%-------------------------------------------------------------------
% loads the state space matrices where the,
% control inputs := Q_N1, Q_N2, Q_N3, H1, H4
% disturbance := Q_in, Q_RD, Q_o1, Q_o2, Q_o3
%-------------------------------------------------------------------
A_th =@(L3,Qin,Qrd,Qn,Qo) Model_cl.A(L3,Qin,Qrd,Qn,Qo);
select = @(M,r,c) M(r,c);
%N_input = input('Enter 1 for 5 inputs and 2 for 2 input actuation ');
if N_input == 1
Bb =@(L3,Qin,Qrd,Qn,Qo) [select(Model_cl.B(L3,Qin,Qrd,Qn,Qo),1:MatSize,1)+select(Model_cl.B(L3,Qin,Qrd,Qn,Qo),1:MatSize,4) (select(Model_cl.B(L3,Qin,Qrd,Qn,Qo),1:MatSize,5)+select(Model_cl.B(L3,Qin,Qrd,Qn,Qo),1:MatSize,7))];    
    for i =1:1:N_nz
        eval(['Bb=@(L3,Qin,Qrd,Qn,Qo)[Bb(L3,Qin,Qrd,Qn,Qo) select(Model_cl.B(L3,Qin,Qrd,Qn,Qo),1:MatSize,' num2str(7+i*2) ')];']);
    end
    Bb=@(L3,Qin,Qrd,Qn,Qo)[Bb(L3,Qin,Qrd,Qn,Qo)  Model_cl.f(L3,Qin,Qrd,Qn,Qo)];
else
Bb =@(L3,Qin,Qrd,Qn,Qo) [select(Model_cl.B(L3,Qin,Qrd,Qn,Qo),1:MatSize,1)+select(Model_cl.B(L3,Qin,Qrd,Qn,Qo),1:MatSize,4) (select(Model_cl.B(L3,Qin,Qrd,Qn,Qo),1:MatSize,5)+select(Model_cl.B(L3,Qin,Qrd,Qn,Qo),1:MatSize,7))];    
end

%% Discretisation- 4th order Taylor series expansion
A_d =@(L3,Qin,Qrd,Qn,Qo) (1/24 * A_th(L3,Qin,Qrd,Qn,Qo)^4 *Ts^4 ) + (1/6 * A_th(L3,Qin,Qrd,Qn,Qo)^3 *Ts^3 )+ (1/2 * A_th(L3,Qin,Qrd,Qn,Qo)^2 *Ts^2 ) + ( A_th(L3,Qin,Qrd,Qn,Qo)*Ts ) + eye(size(A_th(L3,Qin,Qrd,Qn,Qo),1));
B_d_int =@(L3,Qin,Qrd,Qn,Qo) (1/24 * A_th(L3,Qin,Qrd,Qn,Qo)^3*Bb(L3,Qin,Qrd,Qn,Qo) *Ts^4 ) + (1/6 * A_th(L3,Qin,Qrd,Qn,Qo)^2*Bb(L3,Qin,Qrd,Qn,Qo) *Ts^3 )+ (1/2 * A_th(L3,Qin,Qrd,Qn,Qo)^1*Bb(L3,Qin,Qrd,Qn,Qo) *Ts^2 ) + ( Bb(L3,Qin,Qrd,Qn,Qo)*Ts );

if N_input ==1
    B_d =@(L3,Qin,Qrd,Qn,Qo) select(B_d_int(L3,Qin,Qrd,Qn,Qo),1:MatSize,1:2+N_nz);
    f_d =@(L3,Qin,Qrd,Qn,Qo) select(B_d_int(L3,Qin,Qrd,Qn,Qo),1:MatSize,2+N_nz+1);
else 
    B_d =@(L3,Qin,Qrd,Qn,Qo) select(B_d_int(L3,Qin,Qrd,Qn,Qo),1:MatSize,1:2);
    f_d =@(L3,Qin,Qrd,Qn,Qo) select(B_d_int(L3,Qin,Qrd,Qn,Qo),1:MatSize,3);
end
%% reference creation
nx = MatSize;
nu=2+N_nz;

total_Time = 0:Ts:T;
Rk = ones(N_nz+2,1)*T_ref;

%% Calling paramerter/distubance profile 

%[Dist_Prof,~] = createdisnozzles(N_nz,T,Ts,N_per,L0,Res,Ink);
% [Dist_Prof,~] = createdis(N_nz,T,Ts,N_per,0.01);
% Dist_Prof.L3 = ones(size(Dist_Prof.Qin))*Res.h*2;
%%
xh=zeros(nx,length(total_Time)+1); 
uh=zeros(nu,length(total_Time));
%error = zeros(nu,length(total_Time));
x0= X0*ones(nx,1);    % Assigning Initial condition
xh(:,1) = x0;

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
tic
m=9:2:9+(N_nz-1)*2;
midnoz = m(:,ceil(end/2));
%B_D_K = cell(1,1001);
aaa=0;
for k =1:1: T_sim
    % Asigning the parameters/disturbance for $N$ horizon to variables
        %Q_in  =  Dist_Prof.Qin(k);   
        %Q_RD  =  Dist_Prof.QRD(k);
        %QQ=[Q_in Q_RD Dist_Prof.QN(:,k)' Dist_Prof.QN(:,k)' Dist_Prof.L3(:,k)];

        %subvec = [Qin; Qrd ;Qn ;Qo ;L3];
tt=tic;
        A_D_K   = A_d(Dist_Prof.L3(k),Dist_Prof.Qin(k),Dist_Prof.QRD(k),Dist_Prof.QN(:,k),Dist_Prof.QN(:,k)); %double(subs(A_d,subvec,QQ'));
        B_D_K   = B_d(Dist_Prof.L3(k),Dist_Prof.Qin(k),Dist_Prof.QRD(k),Dist_Prof.QN(:,k),Dist_Prof.QN(:,k));;
        F_DD    =  f_d(Dist_Prof.L3(k),Dist_Prof.Qin(k),Dist_Prof.QRD(k),Dist_Prof.QN(:,k),Dist_Prof.QN(:,k));;

t(k)=toc(tt) ;       
        % calculate error
        if N_input == 2 
            error(k+aw) = Rk(1)-xh(midnoz ,k);
        else
            xx=[xh(midnoz,k); xh(midnoz ,k); xh(9,k)];
            for i=1:1:N_nz-1
                xx=[xx; xh(9+i*2,k)];
            end
            error(1:N_nz+2,k+aw) = Rk-xx;
        end
        % P term
            u_k = Kp.*error(end);
        % I term
        if N_input == 2
            u_k = Ki.*sum(error(1,k:k+aw),2) + u_k;
            error1(:,k)=sum(error(1,k:k+aw),2);
        else
            u_k = Ki.*sum(error(1:N_nz+2,k:k+aw),2) + u_k;
            error1(:,k)=sum(error(1:N_nz+2,k:k+aw),2);
        end
          
        % limit PID controller output
            for i = 1:1:size(B_D_K,2)
                if u_k(i) > U_h_max && i<3
                    u_k(i) = U_h_max;
                elseif u_k(i) < U_h_min && i<3
                    u_k(i) = U_h_min;
                elseif u_k(i) > U_p_max && i>=3
                    u_k(i) = U_p_max;
                elseif u_k(i) < U_p_min  && i>=3
                    u_k(i) = U_p_min;
                end
            end
        % simulate process behaviour
        
        xh(:,k+1) = A_D_K*xh(:,k) + B_D_K*u_k + F_DD;
        uh(:,k) = u_k;
     
        com=((((k-1))/T))*Ts;

        if floor(com*10) ~=aaa
            waitbar(floor(com*10)/10,f,'Simulating PI...');
            aaa=floor(com*10); 
        end   
        
end
ttt=toc(t_complete);
disp(['Completion of simimulation takes ' num2str(ttt) ' seconds'])
%%

figure(1)
hold on
ii =1:1:N_nz;

for i=ii
    plot(total_Time(1:T_sim),xh(7+i*2,1:T_sim))
end 

plot(total_Time(1:T_sim),T_ref*ones(length(total_Time(1:T_sim)),1),'black--'),
title('DFA Ink Temperature'),

xlabel('Time [s]'),ylabel('Temperature $^\circ C$')
figp
axis([0 T floor(min([xh(9,1:T_sim) xh(11,1:T_sim) xh(13,1:T_sim)]))-2 ceil(max([xh(9,1:T_sim) xh(11,1:T_sim) xh(13,1:T_sim) T_ref]))+2])

figure(2), plot(total_Time(1:T_sim),xh(1,1:T_sim),total_Time(1:T_sim),xh(5,1:T_sim))
hold on,plot(total_Time(1:T_sim),70*ones(length(total_Time(1:T_sim)),1),'black--'), 
title('Reservoir and Distributor Solid Left'),legend('$T^{R1}$','$T^{D1}$','$T_{Ref}$','location','best'),xlabel('Seconds [s]'),ylabel('Temperature $^\circ C$')
figp
axis([0 T floor(min([xh(1,1:T_sim) xh(4,1:T_sim) xh(5,1:T_sim) xh(7,1:T_sim)]))-2 ceil(max([xh(1,1:T_sim) xh(4,1:T_sim) xh(5,1:T_sim) xh(7,1:T_sim) T_ref]))+2])

if N_input == 1
figure(3), plot(total_Time(1:T_sim),uh(3,1:T_sim),total_Time(1:T_sim),uh(4,1:T_sim),total_Time(1:T_sim),uh(5,1:T_sim))
title('Heat Flux in Nozzles'),legend('HN1','HN2','HN3'),xlabel('Seconds [s]'),ylabel('Flux')
axis([0 T U_p_min U_p_max])
figp
end

figure(4), plot(total_Time(1:T_sim),uh(1,1:T_sim),total_Time(1:T_sim),uh(2,1:T_sim))
title('Heat Flux in solid blocks 1 and 4'),legend('$H^{1}$','$H^{4}$'),xlabel('Seconds [s]'),ylabel('Flux')
axis([0 T U_h_min U_h_max])
figp

%%
% %% Video
% delete('temperature_blocks_try_1')
% vobj = VideoWriter('temperature_blocks_try_1','Motion JPEG AVI');
% vobj.FrameRate = 1/Ts;  %just to be cinematic 
% vobj.Quality=100; %no compression
% open(vobj); 
% figure(5)
% for ind = 1:T_sim
%     
%     subplot(211)
%     bar([8:2:6+N_nz*2],xh([8:2:6+N_nz*2],ind),'b')
%     hold on
%     plot([6 8+N_nz*2],[T_ref T_ref],'r--')
%     title('Temperature error')
%     ylabel('T_{nozzle}-T_{ref}')
%     xlabel('Nozzle')
%     ylim([60 80])
%     hold off
%     
%     subplot(212)
%     bar([8:2:6+N_nz*2],xh([8:2:6+N_nz*2],ind),'b')
%     hold on
%     plot([6 8+N_nz*2],[T_ref T_ref],'r--')
%     title('Temperature error')
%     ylabel('T_{nozzle}-T_{ref}')
%     xlabel('Nozzle')
%     ylim([69.5 70.5])
%     
%     frame = getframe(5); %get image of whats displayed in the figure
%     writeVideo(vobj, frame);
%     clc 
%     complete_video = round(((((ind-1)*10)/T)*100)/100)
%     hold off
% end
% %close the object so its no longer tied up by matlab
%  close(vobj);
%  close(gcf) %close figure since we don't need it anymore

%%
figureHandle=createfigure(total_Time(1:T_sim),xh(1,1:T_sim),xh(3,1:T_sim),xh(4,1:T_sim),...
     xh(5,1:T_sim),xh(6,1:T_sim),xh(7,1:T_sim),...
     xh(8,1:T_sim),xh(9,1:T_sim),xh(10,1:T_sim),...
     xh(11,1:T_sim),xh(12,1:T_sim),xh(13,1:T_sim),xh(14,1:T_sim),T);
