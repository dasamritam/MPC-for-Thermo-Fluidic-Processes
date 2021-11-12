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
%% Tuning Parameters 
 
% Choosing State and input Weights
Q = diag ([10 10 10 10 10 10 10 10 repmat([10000 10],1,N_nz)] );

if N_input ==1
    R = diag ([.001 .001 ones(1,N_nz)]);
    nu = N_nz+2;
else 
    R = diag ([.001 .001]);
    nu = 2;
end

%% reference creation

nx = MatSize;%size(A_d,1);

total_Time = 0:Ts:T;
Ini_Ref = Tss([3 6 9:2:8+N_nz*2]);

% % % % % N=40;
Ref_Prof = repmat(Ini_Ref ,1,N_per+1);
Rk = reshape(Ref_Prof,(2+N_nz)*(N_per+1),1);%Rk=reshape(Ref_prof(:,k:k+N),3*(N+1),1);
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
[ c,Mqi,Ei,bi,MN,bN,Mcal,Dcal,Ecal] = getEDMcalc( nx, nu, x_max, x_min, u_max, u_min, N_per);
Mcal= sparse(Mcal);
Dcal= sparse(Dcal);
Ecal= sparse(Ecal);
c= sparse(c);


%% loop starts 


% Pre allocating size of the variables
xh=zeros(nx,length(total_Time)+1); 
uh=zeros(nu,length(total_Time));
Xs=zeros(nx,length(total_Time));
Us=zeros(nu,length(total_Time));
A_DD = cell(1,length(total_Time));
B_DD = cell(1,length(total_Time)); 
A_D_K = cell(1,N_per);
B_D_K = cell(1,N_per);
F_DD = cell(1,N_per);
BigTerm = cell(N_per,1);

x0= X0*ones(nx,1);    % Assigning Initial condition
xh(:,1) = x0;

warning('off','all')
opt1 = mpcqpsolverOptions;
options = optimoptions(@quadprog,'Display','off');
opt = optimoptions('lsqlin','Algorithm','interior-point','MaxIterations',300,'display','off');
aaa=0;

for k =1:1: T_sim  %length(total_Time)
%% Phi Gamma A_D_K, B_D_K creation
% Asigning the parameters/disturbance for $N$ horizon to variables
Q_in  =  Dist_Prof.Qin(k:k+N_per - 1);   
Q_RD  =  Dist_Prof.QRD(k:k+N_per - 1);

QQ=[Q_in' Q_RD' Dist_Prof.QN(:,k:k+N_per - 1)' Dist_Prof.QN(:,k:k+N_per - 1)' Dist_Prof.L3(:,k:k+N_per - 1)'];


 %% Creation of Prediction Matrices
 Gamma = zeros(nx*N_per,nu*N_per);  %Pre allocation of size of Gamma matix

 if k==1    % condition if the time instant k=1
     for i = 1:1:N_per
                 
         A_D_K{i}   = A_d(Dist_Prof.L3(k),Dist_Prof.Qin(k),Dist_Prof.QRD(k),Dist_Prof.QN(:,k),Dist_Prof.QN(:,k)); %double(subs(A_d,subvec,QQ'));
         B_D_K{i}   = B_d(Dist_Prof.L3(k),Dist_Prof.Qin(k),Dist_Prof.QRD(k),Dist_Prof.QN(:,k),Dist_Prof.QN(:,k));
         F_DD{i}    =  f_d(Dist_Prof.L3(k),Dist_Prof.Qin(k),Dist_Prof.QRD(k),Dist_Prof.QN(:,k),Dist_Prof.QN(:,k));

         % Formation of Phi, Gamma, Eta_del matrices for each time instance
         % 'k'
         if i==1
             BigTerm{i,1} = [A_D_K{i}  B_D_K{i}];
             BigTerm2{i,1} = [A_D_K{i} F_DD{i}];
         else
             BigTerm{i,1} = [A_D_K{i}*BigTerm{i-1,1}  B_D_K{i}];
             BigTerm2{i,1} = [A_D_K{i}*BigTerm2{i-1,1} F_DD{i}];
         end
         Phi ( (i-1)*nx+1:nx*i,1:nx ) = BigTerm{i,1}(:, 1:nx);
         Gamma(  (i-1)*nx+1:nx*i, 1:nu*i  ) =  BigTerm{i,1}(:, nx+1:end );  % A_DD^(i-j) * B_d;
         Eta_del((i-1)*nx+1:nx*i, 1:i) = BigTerm2{i,1}(:, (nx)+1:(nx)+ i  );
     end

 else
     % A_{i+N|k-1}=A_{i+N-1|k} for the purpose of constructing the Phi, Gamma, Eta_del
     % matrices
        
         A_D_K{1}=[];       % emptying the A_{i|k-1}
         B_D_K{1} = [];     % emptying the B_{i|k-1}
         F_DD{1}=[];        % emptying the F_{i|k-1}
         
         % shifting A_{i+1|k-1} = A_{i|k} ,..., A_{i+N|k-1}=A_{i+N-1|k} &
         % so on
         A_D_K = circshift(A_D_K,-1); 
         B_D_K = circshift(B_D_K,-1);
         F_DD = circshift(F_DD,-1);
         A_D_K{N_per} = A_d(Dist_Prof.L3(k),Dist_Prof.Qin(k),Dist_Prof.QRD(k),Dist_Prof.QN(:,k),Dist_Prof.QN(:,k)); %double(subs(A_d,subvec,QQ'));
         B_D_K{N_per} =  B_d(Dist_Prof.L3(k),Dist_Prof.Qin(k),Dist_Prof.QRD(k),Dist_Prof.QN(:,k),Dist_Prof.QN(:,k));
         F_DD{N_per} = f_d(Dist_Prof.L3(k),Dist_Prof.Qin(k),Dist_Prof.QRD(k),Dist_Prof.QN(:,k),Dist_Prof.QN(:,k));
     for i = 1:1:N_per
        if i==1
             BigTerm{i,1} = [A_D_K{i}  B_D_K{i}];
              BigTerm2{i,1} = [A_D_K{i} F_DD{i}];
         else
             BigTerm{i,1} = [A_D_K{i}*BigTerm{i-1,1}  B_D_K{i}];
              BigTerm2{i,1} = [A_D_K{i}*BigTerm2{i-1,1} F_DD{i}];
         end
         Phi ( (i-1)*nx+1:nx*i,1:nx ) = BigTerm{i,1}(:, 1:nx);
         Gamma(  (i-1)*nx+1:nx*i, 1:nu*i  ) =  BigTerm{i,1}(:, nx+1:end )  ;  % A_DD^(i-j) * B_d;
         Eta_del((i-1)*nx+1:nx*i, 1:i) = BigTerm2{i,1}(:, (nx)+1:(nx)+ i  );
     end
     
 end

% Formation of Affine_Eta_del  
Affine_Eta_del = Eta_del*ones(size(Eta_del,2),1);
% formation of F(\theta1),...,F_{\thetaN}
Fn = [-cell2mat(F_DD(:));zeros(nx,1)];

% Creating the inequality constraint matrices which is of the form L U_k <= c + W x_{0|k} + Zeta 
L = Mcal*Gamma+Ecal;
W = -Dcal-Mcal*Phi;

% Creating Stabilising Terminal set using DLQR
[K_mpc,P]=dlqr(A_D_K{N_per},B_D_K{N_per},Q,R);
% Creation of constraint weighing matrices 
[ Psi, Omega ] = QRPN2PsiOmega( Q,R,10*Q,N_per );


%% Creation of Target system  for reference tracking
% \begin{pmatrix}
% A_{0|k}    & -I_n    &       &  0   &  B_{0|k}    &        &  0     &  0\\
% & \ddots  &\ddots &      &         & \ddots & &    \\
% 0     &         &  A_{N-1|k}& -I_n & 0     &        & B_{N-1|k}& 0 \\
% 0     &         &  0&A_{N|k} -I_n & 0     &        &   & B_{N|k} \\
% C  &         &  0 &  0  &  0 &        &  0& 0\\
% & \ddots  &    &     &    & \ddots &   &   \\
% 0  &         &  C &  0  &  0 &        &  0& 0\\
% 0  &         &  0 &  C  &  0 &        &  0& 0\\
% \end{pmatrix} 
% \begin{pmatrix}
% x_{0|k}^{ss}\\ X_k^{ss}\\ U_k^{ss}\\ u_{N|k}^{ss}
% \end{pmatrix} = 
% \begin{pmatrix}
% -\Xi_k \\ 0_{n\times1} \\ \mathcal{R}_k
% \end{pmatrix}

SHinv = zeros((nx+N_nz)*N_per+nx+N_nz,(nx+nu)*N_per+nx+nu);
SHinv(N_per*(nx)+1:N_per*nx+nx,N_per*(nx)+1:N_per*nx+nx) = A_D_K{N_per}-eye(nx);
SHinv(N_per*(nx)+1:N_per*nx+nx,end-nu+1:end) = B_D_K{N_per};
SHinv((N_per+1)*nx+1:(N_per+1)*nx+(N_nz+2),1:nx) = C_r;
for i = 1:N_per
    SHinv((i-1)*nx+1:i*(nx),(i-1)*nx+1:i*(nx)) =  A_D_K{i};
    SHinv((i-1)*nx+1:i*(nx),(i)*nx+1:(i+1)*(nx)) = -eye(nx);
    SHinv((i-1)*nx+1:i*(nx),N_per*nx+nx+((i-1)*nu+1:(i)*(nu))) = B_D_K{i};
    SHinv((N_per+1)*nx+((i)*(N_nz+2)+1:(i+1)*(N_nz+2)),(i)*nx+1:(i+1)*(nx)) = C_r;
end


%% Solving the Linear equation using constraint least square 

% constraints for the X^{ss}, U^{ss}
lb_x = repmat(x_min,N_per+1,1);
lb_u = repmat(u_min,N_per+1,1);
ub_u = repmat(u_max,N_per+1,1);
ub_x = repmat(x_max,N_per+1,1);

lb=[lb_x;lb_u];
ub=[ub_x;ub_u];
% Constrained Least square 


%XsUs = lsqlin(SHinv,[Fn;Rk],[],[],[],[],lb,ub,[],opt); 
XsUs=inv(SHinv)*[Fn;Rk];
% XsUs  = [ x_{0|k}^{ss}\\ X_k^{ss}\\ U_k^{ss}\\ u_{N|k}^{ss}]
% XsUs2 = [ X_k^{ss}\\ U_k^{ss}]
XsUs2=[zeros((nx+nu)*N_per,nx), eye((nx+nu)*N_per), zeros((nx+nu)*N_per,nu)]*XsUs;
 Xs(:,k) = XsUs2(1:nx);
 Us(:,k) = XsUs2((N_per)*nx+1:(N_per)*nx+nu);


%% Finding G F Y H for the quad prog cost function
 
G = 2*(Psi+Gamma'*Omega*Gamma);
G1 = (G+G')/2;
[G2,~] = chol(G1,'lower');
Linv = G2\eye(size(G1,1));
% Built F Matrix
F = 2*Gamma'*Omega*Phi;
% Build Y Matrix
Y = -2*[Gamma'*Omega  Psi]*XsUs2;
H = 2*Gamma'*Omega*Affine_Eta_del;

%% Receeding horizon optimisation

    x_k = xh(:,k);
    options = optimoptions(@quadprog,'Display','off');
    t_qp=tic;
        [u_qp1,fval(k),exitflag] = quadprog(G ,F*xh(:,k)+Y+H ,L,W*xh(:,k)+c-Mcal*Affine_Eta_del,[],[],[],[],[],options);
    t_quad(k)=toc(t_qp);
%     tic
%     if k==1
%         iA0 = false(size(c+W*xh(:,k)));
%     end
%         [u_qp,status(k),iA0,~] = mpcqpsolver(Linv,F*xh(:,k)+Y+H,-L,-c-W*xh(:,k)+-Mcal*Affine_Eta_del,[],zeros(0,1),iA0,opt1);
%     t_mpcqp(k)=toc;
    u_k = u_qp1(1:nu);    
%%% --------Discrete Statespace model -----------
     xh(:,k+1) = A_D_K{1}*xh(:,k) + B_D_K{1}*u_k + F_DD{1};
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

figure(1)
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

figure(2), plot(total_Time(1:T_sim),xh(1,1:T_sim),total_Time(1:T_sim),xh(4,1:T_sim))
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