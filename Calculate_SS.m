
tic
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
Bb =@(L3,Qin,Qrd,Qn,Qo) [select(Model_cl.B(L3,Qin,Qrd,Qn,Qo),1:MatSize,1)+select(Model_cl.B(L3,Qin,Qrd,Qn,Qo),1:MatSize,4) (select(Model_cl.B(L3,Qin,Qrd,Qn,Qo),1:MatSize,5)+select(Model_cl.B(L3,Qin,Qrd,Qn,Qo),1:MatSize,7)) Model_cl.f(L3,Qin,Qrd,Qn,Qo)];    
end


%% Tuning Parameters 
 
% Choosing State and input Weights
Q = diag ([10 10 10 10 10 10 10 10 repmat([10000 10],1,N_nz)] );

if N_input ==1
    R = diag ([.001 .001 ones(1,N_nz)]);
    nu=N_nz+2;
else 
    R = diag ([.001 .001]);
    nu=2;
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



total_Time = 0:Ts:T;
Ini_Ref = ones(N_nz,1)*T_ref;

% % % % % N=40;
Ref_Prof = repmat(Ini_Ref ,1,N_per+1);
Rk = reshape(Ref_Prof,N_nz*(N_per+1),1);%Rk=reshape(Ref_prof(:,k:k+N),3*(N+1),1);
C_r= sparse(1:N_nz,9:2:9+(N_nz-1)*2,ones(1,N_nz),N_nz,nx); 
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
[ c,Mqi,Ei,bi,MN,bN,Mcal,Dcal,Ecal] = getEDMcalc( nx, nu, x_max, x_min, u_max, u_min, N_per);
Mcal= sparse(Mcal);
Dcal= sparse(Dcal);
Ecal= sparse(Ecal);
c= sparse(c);


%% loop starts 
tic

% Pre allocating size of the variables

% A_DD = cell(1,length(total_Time));
% B_DD = cell(1,length(total_Time)); 
% F_DD = cell(1,length(total_Time));
% A_D_K = cell(1,N_per);
% B_D_K = cell(1,N_per);
% BigTerm = cell(N_per,1);

% x0= T_ref*ones(nx,1);    % Assigning Initial condition
% xh(:,1) = x0;

warning('off','all')
% k=1;
%for k =1:1: T_sim  %length(total_Time)
%% Phi Gamma A_D_K, B_D_K creation
% Asigning the parameters/disturbance for $N$ horizon to variables

QQ=[Qinss' Qrdss' QNss' QOss' L3ss'];


 %% Creation of Prediction Matrices
 Gamma = zeros(nx*N_per,nu*N_per);  %Pre allocation of size of Gamma matix

%      for i = 1:1:N_per
         subvec = [Qin; Qrd ;Qn ;Qo ;L3];
         d = QQ';
         
         A_D_K   = A_d(L3ss,Qinss,Qrdss,QNss,QOss); %double(subs(A_d,subvec,QQ'));
         B_D_K   = B_d(L3ss,Qinss,Qrdss,QNss,QOss);
         F_DD    =  f_d(L3ss,Qinss,Qrdss,QNss,QOss);

%      end

     Gamma = sparse(Gamma);

% formation of F(\theta1),...,F_{\thetaN}
Fn = [F_DD;ones(N_nz,1)*T_ref];


%% Creation of Target system  for reference tracking


SHinv = zeros((nx+N_nz),(nx+nu));
SHinv(1:nx,1:nx) = A_D_K-eye(nx);
SHinv(1:nx,end-nu+1:end) = B_D_K;
SHinv(nx+1:nx+N_nz,1:nx) = C_r;
% for i = 1:N_per
%     SHinv((i-1)*nx+1:i*(nx),(i-1)*nx+1:i*(nx)) =  A_D_K{i};
%     SHinv((i-1)*nx+1:i*(nx),(i)*nx+1:(i+1)*(nx)) = -eye(nx);
%     SHinv((i-1)*nx+1:i*(nx),N_per*nx+nx+((i-1)*nu+1:(i)*(nu))) = B_D_K{i};
%     SHinv((N_per+1)*nx+((i)*N_nz+1:(i+1)*(N_nz)),(i)*nx+1:(i+1)*(nx)) = C_r;
% end


%% Solving the Linear equation using constraint least square 

% constraints for the X^{ss}, U^{ss}
lb_x = x_min;
lb_u = u_min;
ub_u = u_max;
ub_x = x_max;

lb=[lb_x;lb_u];
ub=[ub_x;ub_u];
% Constrained Least square 
opt = optimoptions('lsqlin','Algorithm','interior-point','MaxIterations',300,'display','off');

XsUs = lsqlin(SHinv,[Fn],[],[],[],[],lb,ub,[],opt); 
% XsUs  = [ x_{0|k}^{ss}\\ X_k^{ss}\\ U_k^{ss}\\ u_{N|k}^{ss}]
% XsUs2 = [ X_k^{ss}\\ U_k^{ss}]

 Tss= XsUs(1:nx);
 Us = XsUs(nx+1:end);

ttt=toc;
disp(['Estimation of SS takes ' num2str(ttt) ' seconds'])


