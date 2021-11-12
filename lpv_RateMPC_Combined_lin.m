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
Q = 1000*eye(N_nz);
if N_input ==1
    R = diag ([.0001 .0001 ones(1,N_nz)]);
else 
    R = diag ([.0001 .0001]);
end


%% reference creation

Fin_T = 1000;
nx = size(B_d,1);
nu=size(B_d,2);


total_Time = 0:Ts:T;
Ini_Ref = ones(N_nz,1)*T_ref;

Rk = repmat(Ini_Ref ,N_per,1);

C_r= sparse(1:N_nz,9:2:9+(N_nz-1)*2,ones(1,N_nz),N_nz,nx);
%% constraint matrix declaration

x_max = X_max*ones(nx,1);
x_min = X_min*ones(nx,1);

if N_input ==1
    u_max = [ones(2,1)*U_h_max;ones(N_nz,1)*U_p_max];
    u_min = [ones(2,1)*U_h_min;ones(N_nz,1)*U_p_min];
    del_u_max = del_U_max*[ones(2,1)*U_h_max;ones(N_nz,1)*U_p_max];
    del_u_min = -del_U_min*[ones(2,1)*U_h_max;ones(N_nz,1)*U_p_max];
else
    u_max = ones(2,1)*U_h_max;
    u_min = ones(2,1)*U_h_min;
    del_u_max = del_U_max*ones(2,1)*U_h_max;
    del_u_min = -del_U_min*ones(2,1)*U_h_max;
end

Ymax = Y_max*ones(N_nz,1);
Ymin = Y_min*ones(N_nz,1);

%% Calling paramerter/distubance profile 
%[Dist_Prof,Dkminus1] = createdis(N_nz,T,Ts,N_per,0.01);
%Dist_Prof.L3 = ones(size(Dist_Prof.Qin))*Res.h*2;
%[Dist_Prof,Dkminus1] = createdisnozzles(N_nz,T,Ts,N_per,L0,Res,Ink);
%% Constraints matrices creation


Mi          = [zeros(4*nu,N_nz); -eye(N_nz);eye(N_nz)];
M0          = [zeros(4*nu,N_nz); -eye(N_nz);eye(N_nz)];

E0i         = [-eye(nu);eye(nu);-eye(nu);eye(nu);zeros(2*N_nz,nu)];

Eoff        = [zeros(2*nu,nu);-eye(nu);eye(nu);zeros(2*N_nz,nu)];
Mcal = [];
Ecal = [];%zeros((4*nu+2*nr)*(N+1),N*nu);
for j = 1:N_per
    if j==1
    Mcal = blkdiag(Mi,Mcal);
    Ecal = blkdiag(Ecal,E0i);
    else 
    Mcal = blkdiag(Mi,Mcal);
    Ecal = blkdiag(Ecal,E0i);
    Eoffdia = repmat(Eoff,1,j-1);
    Ecal((4*nu+2*N_nz)*(j-1)+1:(4*nu+2*N_nz)*j,1:nu*(j-1)) = Eoffdia; 
    end
end
Mcal = [zeros(size(Mi,1),size(Mcal,2));Mcal];
Ecal = [Ecal;zeros(size(E0i,1),size(Ecal,2))];
Dcal = [M0 ; zeros((4*nu+2*N_nz)*(N_per),N_nz)];
Mcal= sparse(Mcal);
Ecal= sparse(Ecal);

%% loop starts 
tic
% Pre allocating size of the variables
xp=zeros(nx,length(total_Time)+1);
yh=zeros(N_nz,length(total_Time));
X=zeros(nx+N_nz,length(total_Time)+1);

du=zeros(nu,length(total_Time));
uh=zeros(nu,length(total_Time));

% Assigning Initial condition
x0= X0*ones(nx,1);    % Assigning Initial condition
x_m1 = x0;  % [66.1427842726780 54.3628428223182 66.1427842726780 66.2119752738876 56.3365030051215 66.2119752738876 66.1889908131009 65.9357310598404 66.1780823615753 65.9306764653475 66.1780823615753 65.9357310598404 66.1889908131009]';
xp(:,1) = x0;
zero_nxnr = zeros(nx,N_nz);
eye_nr = eye(N_nz);
zero_nrnx = zeros(N_nz,nx);
C_ext = [zero_nrnx eye_nr];
A_ext = [A_d       zero_nxnr;...
                     C_r*A_d   eye_nr];
         B_ext = [B_d ; C_r* B_d]; 
warning('off','all')
aaa=0;
%% Creation of Prediction Matrices
 Gamma = zeros(N_nz*N_per,nu*N_per);
 %C_reye=[]
% subvec = [Qin; Qrd ;Qn ;Qo ;L3];
phi=C_ext*A_ext;
for i=2:N_per
    phi=[phi;C_ext*A_ext^i];
end
gamma=[];
a=[];
for i=1:N_per
    for n=N_per:-1:1
        b=n-N_per+(i-1);
        if b>=0
            a=[a C_ext*A_ext^(n-N_per+(i-1))*B_ext];
        else
            a=[a zeros(size(C_ext*B_ext,1),size(C_ext*B_ext,2))];
        end
    end
    %C_reye=blkdiag(C_reye,C_r);
    gamma=[gamma;a];
    a=[];
end

Phi=phi;
Gamma = gamma;

%%  Rate constraint
opt = mpcqpsolverOptions;
opt.FeasibilityTol = 1.0e-3;
delxp(:,1)= xp(:,1)-x_m1;
    xp(:,1)=xp(:,1);
     u_kpre = zeros(nu,1);% u_max;%
     y0 = C_r*x0;
     yh(:,1)=y0;
     X(:,1) = [delxp(:,1);yh(:,1)];
     
     L = Mcal*Gamma+Ecal;
     W = - (Dcal*C_ext+Mcal*Phi);
     P=50000*Q;
     [ Psi, Omega ] = QRPN2PsiOmega( Q,R,P,N_per );
     
    G = 2*(Psi+Gamma'*Omega*Gamma);
    G1 = (G+G')/2;
    [G2,~] = chol(G1,'lower');
    Linv = G2\eye(size(G1,1));

    opt = mpcqpsolverOptions;
    options = optimoptions(@quadprog,'Display','off');

    % Built F Matrix
    F = 2*Gamma'*Omega;
%%
for k =1:1: T_sim%length(total_Time)
bi          = [ -del_u_min;
                     del_u_max;
                     -u_min + u_kpre;
                     u_max- u_kpre;
                     -Ymin;
                     Ymax];
bN          = [zeros(4*nu,1); -Ymin; Ymax];
c = repmat(bi,N_per,1);
c = [c;bN];
if k==1
    iA0 = false(size(c+W*X(:,1)));
end

    x_k = X(:,k);
        
% tic
%         [du_qp1,fval,exitflag] = quadprog(G,F*(Phi*x_k - Rk),L,W*x_k+c ,[],[],[],[],[],options);
% t_quad(k)=toc;
tic
        [du_qp,status(k),iA0,~] = mpcqpsolver(Linv,F*(Phi*x_k - Rk),-L,-c-W*x_k,[],zeros(0,1),iA0,opt);
        
 t_mpcqp(k)=toc;
    u_k = u_kpre +du_qp(1:nu); % u(k) = u(k-1)+ \Del u(k)
    uh(:,k) = u_k;
    du(:,k) = du_qp(1:nu);
    df_k=[Dist_Prof.Qin(:,k)-Qinss;Dist_Prof.QRD(:,k)-Qrdss;Dist_Prof.QN(:,k)-QNss;Dist_Prof.QN(:,k)-QOss;Dist_Prof.L3(:,k)-L3ss];

    %-------Discrete State Space model --------------
    xp(:,k+1) = A_d*xp(:,k) + B_d*u_k + F_d*df_k;

    yh(:,k+1) = C_r*xp(:,k+1);
    delxp(:,k+1) = xp(:,k+1)-xp(:,k);
    % ------Extended State Space
    X(:,k+1) = [delxp(:,k+1);yh(:,k+1)];
    u_kpre = uh(:,k);
    com=((((k-1))/T))*Ts;

        if floor(com*10) ~=aaa
            waitbar(floor(com*10)/10,f,'Simulating Rate MPC...');
            aaa=floor(com*10); 
        end
end
ttt=toc(t_complete);
disp(['Completion of simimulation takes ' num2str(ttt) ' seconds'])
% save(['var_RateMPCpipe_VD_N_' num2str(N) 'DeluTest2'],'xp','yh','uh','total_Time','N','Ts','R','Q')

%%

% figure(1), plot(total_Time(1:T_sim),yh(1,1:T_sim),total_Time(1:T_sim),yh(2,1:T_sim),total_Time(1:T_sim),yh(3,1:T_sim))
% title('DFA Ink Temperature'),legend('$T^{F2}$','$T^{F4}$','$T^{F6}$'),xlabel('Seconds [s]'),ylabel('Temperature $^\circ C$')
% hold on,plot(total_Time(1:T_sim),70*ones(length(total_Time(1:T_sim)),1),'black--'), 
% figp

% figure(1)
% subplot(121)
%     hold on
%     ii =1:1:N_nz;
%     for i=ii
%         plot(total_Time(1:T_sim),xp(6+i*2,1:T_sim))
%     end 
%     plot(total_Time(1:T_sim),T_ref*ones(length(total_Time(1:T_sim)),1),'black--'),
%     title('DFA Ink Temperature'),
%     xlabel('Time [s]'),ylabel('Temperature $^\circ C$')
%     figp
%     axis([0 T floor(min([xp(8,1:T_sim) xp(10,1:T_sim) xp(12,1:T_sim)]))-2 ceil(max([xp(8,1:T_sim) xp(10,1:T_sim) xp(12,1:T_sim) T_ref]))+2])
% subplot(122)
%     hold on
%     ii =1:1:N_nz;
%     for i=ii
%         plot(total_Time(1:T_sim),xp(6+i*2,1:T_sim))
%     end 
%     plot(total_Time(1:T_sim),T_ref*ones(length(total_Time(1:T_sim)),1),'black--'),
%     title('DFA Ink Temperature'),
%     xlabel('Time [s]'),ylabel('Temperature $^\circ C$')
%     figp
%     axis([0 T T_ref-0.02 T_ref+0.02])
figure(1)
hold on
ii =1:1:N_nz;

for i=ii
    plot(total_Time(1:T_sim),xp(7+i*2,1:T_sim))
end 

plot(total_Time(1:T_sim),T_ref*ones(length(total_Time(1:T_sim)),1),'black--'),
title('DFA Ink Temperature RMPC with 2 inputs'),

xlabel('Time [s]'),ylabel('Temperature $^\circ C$')
figp
axis([0 T floor(min([xp(9,1:T_sim) xp(11,1:T_sim) xp(13,1:T_sim)]))-2 ceil(max([xp(9,1:T_sim) xp(11,1:T_sim) xp(13,1:T_sim) T_ref]))+2])

figure(2), plot(total_Time(1:T_sim),uh(1,1:T_sim),total_Time(1:T_sim),uh(2,1:T_sim))
    title('Heat Flux in solid blocks 1 and 4'),legend('$u_{1}$','$u_{2}$'),xlabel('Seconds [s]'),ylabel('Heat Power $W$')
    figp
    
    figure(5), plot(total_Time(1:T_sim),du(1,1:T_sim),total_Time(1:T_sim),du(2,1:T_sim))
    title('Rate of Heat Flux in solid blocks 1 and 4'),legend('$\Delta u_1$','$\Delta u_2$'),xlabel('Seconds [s]'),ylabel('Flux')
    figp

if N_input ==1
    
    figure(3), plot(total_Time(1:T_sim),uh(3,1:T_sim),total_Time(1:T_sim),uh(4,1:T_sim),total_Time(1:T_sim),uh(5,1:T_sim))
    title('Heat Flux in Nozzles'),legend('$u_3$','$u_4$','$u_5$'),xlabel('Seconds [s]'),ylabel('Flux')
    axis([0 T 0 del_U_max])
    figp
    figure(4), plot(total_Time(1:T_sim),du(3,1:T_sim),total_Time(1:T_sim),du(4,1:T_sim),total_Time(1:T_sim),du(5,1:T_sim))
    title('Rate of Heat Flux in Nozzles'),legend('$\Delta u_3$','$\Delta u_4$','$\Delta u_5$'),xlabel('Seconds [s]'),ylabel('Heat Power')
    figp
end

 
%% 

 figureHandle=createfigure(total_Time(1:T_sim),xp(1,1:T_sim),xp(3,1:T_sim),xp(4,1:T_sim),...
     xp(5,1:T_sim),xp(6,1:T_sim),xp(7,1:T_sim),...
     xp(8,1:T_sim),xp(9,1:T_sim),xp(10,1:T_sim),...
     xp(11,1:T_sim),xp(12,1:T_sim),xp(13,1:T_sim),xp(14,1:T_sim),T);

 %%
 
 figure
 subplot(421)
 plot(1:T_sim,xp(1,1:T_sim)-xp(4,1:T_sim))
 title('Difference block 1 & 3')
 subplot(422)
 plot(1:T_sim,xp(5,1:T_sim)-xp(7,1:T_sim))
 title('Difference block 4 & 6')
 subplot(423)
 plot(1:T_sim,xp(8,1:T_sim)-xp(10,1:T_sim))
 title('Difference block 7 & 9')
 subplot(424)
 plot(1:T_sim,xp(8,1:T_sim)-xp(12,1:T_sim))
 title('Difference block 7 & 11')
 subplot(425)
 plot(1:T_sim,xp(8,1:T_sim)-xp(14,1:T_sim))
 title('Difference block 7 & 13')
 subplot(426)
 plot(1:T_sim,xp(10,1:T_sim)-xp(12,1:T_sim))
 title('Difference block 9 & 11')
 subplot(427)
 plot(1:T_sim,xp(10,1:T_sim)-xp(14,1:T_sim))
 title('Difference block 9 & 13')
 subplot(428)
 plot(1:T_sim,xp(12,1:T_sim)-xp(14,1:T_sim))
 title('Difference block 11 & 13')
 