t_complete=tic;
%-------------------------------------------------------------------
% loads the state space matrices where the,
% control inputs := Q_N1, Q_N2, Q_N3, H1, H4
% disturbance := Q_in, Q_RD, Q_o1, Q_o2, Q_o3
%-------------------------------------------------------------------
A_th = Lin.A ;

B_b=[Lin.B Lin.f];

%% Discretisation- 4th order Taylor series expansion
A_d = (1/24 * A_th^4 *Ts^4 ) + (1/6 * A_th^3 *Ts^3 )+ (1/2 * A_th^2 *Ts^2 ) + ( A_th*Ts ) + eye(size(A_th,1));
Bb = (1/24 * A_th^3*B_b *Ts^4 ) + (1/6 * A_th^2*B_b *Ts^3 )+ (1/2 * A_th^1*B_b *Ts^2 ) + ( B_b*Ts );

F_d=Bb(:,nu+1:end);
B_d=Bb(:,1:nu);
 

%% Tuning Parameters 
 
% Choosing State and input Weights
Q = diag ([10 10 10 10 10 10 10 10 repmat([10000 10],1,N_nz)] );

if N_input ==1
    R = diag ([.001 .001 ones(1,N_nz)]);
else 
    R = diag ([.001 .001]);
end
%% reference creation

nx = size(A_d,1);
nu=size(B_d,2);


total_Time = 0:Ts:T;
Ini_Ref = Tss([3 6 9:2:8+N_nz*2]);


Ref_Prof = repmat(Ini_Ref ,1,N_per+1);
Rk = reshape(Ref_Prof,(2+N_nz)*(N_per+1),1);
C_r= sparse(3:N_nz+2,9:2:9+(N_nz-1)*2,ones(1,N_nz),N_nz+2,nx); 
C_r(1,3)=1;
C_r(2,6)=1;
%% constraint matrix declaration

x_max = X_max*ones(nx,1);
x_min = X_min*ones(nx,1);

if N_input ==1
    u_max = [ones(2,1)*U_h_max;ones(N_nz,1)*U_p_max];
    u_min = [ones(2,1)*U_h_min;ones(N_nz,1)*U_p_min];
else
    u_max = ones(2,1)*U_h_max;
    u_min = ones(2,1)*U_h_min;
end
%% Calling paramerter/distubance profile 
%[Dist_Prof,~] = Disturbance_Profilepipe (Ts);
% [Dist_Prof,~] = createdis(N_nz,T,Ts,N_per,0.01);
% Dist_Prof.L3 = ones(size(Dist_Prof.Qin))*Res.h*2;
%[Dist_Prof,~] = createdisnozzles(N_nz,T,Ts,N_per,L0,Res,Ink);
%% Constraints matrices creation
% getEDMcalc -creates the constraint matices  Mcal, Dcal, Ecal
[ c,Mqi,Ei,bi,MN,bN,Mcal,Dcal,Ecal] = getEDMcalc(nx, nu, x_max, x_min, u_max, u_min, N_per);


%% loop starts 


% Pre allocating size of the variables
xh=zeros(nx,length(total_Time)+1); 
uh=zeros(nu,length(total_Time));
Xs=zeros(nx,length(total_Time));
Us=zeros(nu,length(total_Time));

x0= X0*ones(nx,1);    % Assigning Initial condition
xh(:,1) = x0;



%% Creation of Prediction Matrices
 Gamma = zeros(N_nz*N_per,nu*N_per);
 %C_reye=[]
% subvec = [Qin; Qrd ;Qn ;Qo ;L3];
phi=A_d;
for i=2:N_per
    phi=[phi;A_d^i];
end
gamma=[];
a=[];
for i=1:N_per
    for n=N_per:-1:1
        b=n-N_per+(i-1);
        if b>=0
            a=[a A_d^(n-N_per+(i-1))*B_d];
        else
            a=[a zeros(size(B_d,1),size(B_d,2))];
        end
    end
    %C_reye=blkdiag(C_reye,C_r);
    gamma=[gamma;a];
    a=[];
end

Phi=phi;
Gamma = gamma;
% Creating the inequality constraint matrices which is of the form L U_k <= c + W x_{0|k} + Zeta 
L = Mcal*Gamma+Ecal;
W = -Dcal-Mcal*Phi;

% Creating Stabilising Terminal set using DLQR
[K_mpc,P]=dlqr(A_d,B_d,Q,R);
% Creation of constraint weighing matrices 
[ Psi, Omega ] = QRPN2PsiOmega( Q,R,10*Q,N_per );


SHinv = zeros((nx+N_nz)*N_per+nx+N_nz,(nx+nu)*N_per+nx+nu);
SHinv(N_per*(nx)+1:N_per*nx+nx,N_per*(nx)+1:N_per*nx+nx) = A_d-eye(nx);
SHinv(N_per*(nx)+1:N_per*nx+nx,end-nu+1:end) = B_d;
SHinv((N_per+1)*nx+1:(N_per+1)*nx+(N_nz+2),1:nx) = C_r;
for i = 1:N_per
    SHinv((i-1)*nx+1:i*(nx),(i-1)*nx+1:i*(nx)) =  A_d;
    SHinv((i-1)*nx+1:i*(nx),(i)*nx+1:(i+1)*(nx)) = -eye(nx);
    SHinv((i-1)*nx+1:i*(nx),N_per*nx+nx+((i-1)*nu+1:(i)*(nu))) = B_d;
    SHinv((N_per+1)*nx+((i)*(N_nz+2)+1:(i+1)*(N_nz+2)),(i)*nx+1:(i+1)*(nx)) = C_r;
end

%% Finding G F Linv for the quad prog cost function
G = 2*(Psi+Gamma'*Omega*Gamma);
G1 = (G+G')/2;
[G2,~] = chol(G1,'lower');
Linv = G2\eye(size(G1,1));
% Built F Matrix
F = 2*Gamma'*Omega*Phi;
%%
warning('off','all')
opt1 = mpcqpsolverOptions;
options = optimoptions(@quadprog,'Display','off');
opt = optimoptions('lsqlin','Algorithm','interior-point','MaxIterations',300,'display','off');
aaa=0;    
for k =1:1: T_sim  %length(total_Time)
%% Solving the Linear equation using constraint least square 
Fn=[];
for i=1:N_per
    df_k=[Dist_Prof.Qin(:,k-1+i)-Qinss;Dist_Prof.QRD(:,k-1+i)-Qrdss;Dist_Prof.QN(:,k-1+i)-QNss;Dist_Prof.QN(:,k-1+i)-QOss;Dist_Prof.L3(:,k-1+i)-L3ss];
    Fn = [Fn;F_d*df_k];
end
Fn = [Fn;zeros(nx,1)];
%XsUs = lsqlin(SHinv,[Fn;Rk],[],[],[],[],lb,ub,[],opt);
XsUs=inv(SHinv)*[Fn;Rk];
% XsUs  = [ x_{0|k}^{ss}\\ X_k^{ss}\\ U_k^{ss}\\ u_{N|k}^{ss}]
% XsUs2 = [ X_k^{ss}\\ U_k^{ss}]
XsUs2=[zeros((nx+nu)*N_per,nx), eye((nx+nu)*N_per), zeros((nx+nu)*N_per,nu)]*XsUs;
 Xs(:,k) = XsUs2(1:nx);
 Us(:,k) = XsUs2((N_per)*nx+1:(N_per)*nx+nu);


%% Finding Y H for the quad prog cost function
% Build Y Matrix
Y = -2*[Gamma'*Omega  Psi]*XsUs2;
% H = 2*Gamma'*Omega;

%% Receeding horizon optimisation

    t_qp=tic;
    [u_qp1,fval(k),exitflag] = quadprog(G ,F*xh(:,k)+Y ,L,W*xh(:,k)+c,[],[],[],[],[],options);
    t_quad(k)=toc(t_qp);
%      tic
%     if k==1
%         iA0 = false(size(c+W*xh(:,k)));
%     end
%  
%         [u_qp,status(k),iA0,~] = mpcqpsolver(Linv,F*xh(:,k)+Y,-L,-c-W*xh(:,k),[],zeros(0,1),iA0,opt1);
%         
%     t_mpcqp(k)=toc;
    u_k = u_qp1(1:nu);  
    
%%% --------Discrete Statespace model -----------
     xh(:,k+1) = A_d*xh(:,k) + B_d*u_k + F_d*df_k;
     uh(:,k) = u_k;

     yh(:,k) = C_r*xh(:,k);
   

        com=((((k-1))/T))*Ts;

        if floor(com*10) ~=aaa
            waitbar(floor(com*10)/10,f,'Simulating MPC...');
            aaa=floor(com*10); 
        end
end
ttt=toc(t_complete);
disp(['Completion of simimulation takes ' num2str(ttt) ' seconds'])
% save(['var_MPCpipe_VD_' num2str(nm)],'xh','yh','uh','total_Time','N','Ts','R','Q')

%%

figure
hold on
ii =1:1:N_nz;

for i=ii
    plot(total_Time(1:T_sim),xh(7+i*2,1:T_sim))
end 

plot(total_Time(1:T_sim),T_ref*ones(length(total_Time(1:T_sim)),1),'black--'),
title('DFA Ink Temperature MPC'),

xlabel('Time [s]'),ylabel('Temperature $^\circ C$')
figp
axis([0 T floor(min([xh(8,1:T_sim) xh(10,1:T_sim) xh(12,1:T_sim)]))-2 ceil(max([xh(8,1:T_sim) xh(10,1:T_sim) xh(12,1:T_sim) T_ref]))+2])

figure, plot(total_Time(1:T_sim),xh(1,1:T_sim),total_Time(1:T_sim),xh(4,1:T_sim))
hold on,plot(total_Time(1:T_sim),70*ones(length(total_Time(1:T_sim)),1),'black--'), 
title('Reservoir and Distributor Solid Left'),legend('$T^{R1}$','$T^{D1}$','$T_{Ref}$','location','best'),xlabel('Seconds [s]'),ylabel('Temperature $^\circ C$')
figp
axis([0 T floor(min([xh(1,1:T_sim) xh(4,1:T_sim) xh(5,1:T_sim) xh(7,1:T_sim)]))-2 ceil(max([xh(1,1:T_sim) xh(4,1:T_sim) xh(5,1:T_sim) xh(7,1:T_sim) T_ref]))+2])

if N_input == 1
figure, plot(total_Time(1:T_sim),uh(3,1:T_sim),total_Time(1:T_sim),uh(4,1:T_sim),total_Time(1:T_sim),uh(5,1:T_sim))
title('Heat Flux in Nozzles'),legend('HN1','HN2','HN3'),xlabel('Seconds [s]'),ylabel('Flux')
axis([0 T U_p_min U_p_max])
figp
end

figure, plot(total_Time(1:T_sim),uh(1,1:T_sim),total_Time(1:T_sim),uh(2,1:T_sim))
title('Heat Flux in solid blocks 1 and 4'),legend('$H^{1}$','$H^{4}$'),xlabel('Seconds [s]'),ylabel('Flux')
axis([0 T U_h_min U_h_max])
figp

%%
figureHandle=createfigure(total_Time(1:T_sim),xh(1,1:T_sim),xh(3,1:T_sim),xh(4,1:T_sim),...
     xh(5,1:T_sim),xh(6,1:T_sim),xh(7,1:T_sim),...
     xh(8,1:T_sim),xh(9,1:T_sim),xh(10,1:T_sim),...
     xh(11,1:T_sim),xh(12,1:T_sim),xh(13,1:T_sim),xh(14,1:T_sim),T);


% %% Video
% delete('temperature_blocks_50_200')
% vobj = VideoWriter('temperature_blocks_50_200','Motion JPEG AVI');
% vobj.FrameRate = 10;  %just to be cinematic 
% vobj.Quality=100; %no compression
% open(vobj); 
% figure(5)
% for ind = 1:T_sim
%     subplot(211)
%     bar([8:2:6+N_nz*2],xh([8:2:6+N_nz*2],ind)-T_ref)
%     title('Temperature error')
%     ylabel('T_{nozzle}-T_{ref}')
%     xlabel('Nozzle')
%     ylim([-10 10])
%     subplot(212)
%     bar([8:2:6+N_nz*2],xh([8:2:6+N_nz*2],ind)-T_ref)
%     title('Temperature error')
%     ylabel('T_{nozzle}-T_{ref}')
%     xlabel('Nozzle')
%     ylim([-0.05 0.05])
%     frame = getframe(5); %get image of whats displayed in the figure
%     writeVideo(vobj, frame);
%     clc 
%     complete_video = round(((((ind-1)*10)/T)*100)/100)
%     
% end
% %close the object so its no longer tied up by matlab
%  close(vobj);
%  close(gcf) %close figure since we don't need it anymore