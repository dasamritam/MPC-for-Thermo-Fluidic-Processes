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
Q = 50000*eye(N_nz);
if N_input ==1
    R = diag ([.0001 .0001 0.01*ones(1,N_nz)]);
    nu=2+N_nz;
else 
    R = diag ([.0001 .0001]);
    nu=2;
end

%% reference creation

Fin_T = 1000;
nx = MatSize;



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
[Dist_Prof,Dkminus1] = createdisnozzles(N_nz,T,Ts,N_per,L0,Res,Ink);
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

A_D_K = cell(1,N_per);
B_D_K = cell(1,N_per);
A_ext = cell(1,N_per);
B_ext = cell(1,N_per);
fd = cell(1,N_per);
delfd = cell(N_per,1);
Eta_delDk = cell(N_per,1);
BigTerm = cell(N_per,1);
% Assigning Initial condition
x0= X0*ones(nx,1);    % Assigning Initial condition
x_m1 = x0;  % [66.1427842726780 54.3628428223182 66.1427842726780 66.2119752738876 56.3365030051215 66.2119752738876 66.1889908131009 65.9357310598404 66.1780823615753 65.9306764653475 66.1780823615753 65.9357310598404 66.1889908131009]';
xp(:,1) = x0;
% variables for extended state
zero_nxnr = zeros(nx,N_nz);
eye_nr = eye(N_nz);
zero_nrnx = zeros(N_nz,nx);
C_ext = [zero_nrnx eye_nr];

warning('off','all')
opt = mpcqpsolverOptions;
    options = optimoptions(@quadprog,'Display','off');
aaa=0;    
for k =1:1: T_sim%length(total_Time)
%% Phi Gamma A_D_K, B_D_K creation
% Asigning the parameters/disturbance for $N$ horizon to variables
Q_in  =  Dist_Prof.Qin(k:k+N_per - 1);   
Q_RD  =  Dist_Prof.QRD(k:k+N_per - 1);


%% Creation of Prediction Matrices
 Gamma = zeros(N_nz*N_per,nu*N_per);


 if k==1  || isempty(A_D_K{2})
     for i = 1:1:N_per
            

         A_D_K{i}   = A_d(Dist_Prof.L3(k),Dist_Prof.Qin(k),Dist_Prof.QRD(k),Dist_Prof.QN(:,k),Dist_Prof.QN(:,k)); %double(subs(A_d,subvec,QQ'));
         B_D_K{i}   = B_d(Dist_Prof.L3(k),Dist_Prof.Qin(k),Dist_Prof.QRD(k),Dist_Prof.QN(:,k),Dist_Prof.QN(:,k));
         fd{i}    =  f_d(Dist_Prof.L3(k),Dist_Prof.Qin(k),Dist_Prof.QRD(k),Dist_Prof.QN(:,k),Dist_Prof.QN(:,k));
         % Formation of Extended state A_ext, B_ext for A_ext x^ext + B_ext Del u
         A_ext{i} = [A_D_K{i}       zero_nxnr;...
                     C_r*A_D_K{i}   eye_nr];
         B_ext{i} = [B_D_K{i} ; C_r* B_D_K{i}];  
        
         if k==1 && i==1
            % array storing the k-1 th instant value of parameters for k =0  
             dminus1_ini = [Dkminus1.Qin ; Dkminus1.QRD; Dkminus1.QN;Dkminus1.QN;Dkminus1.L3]; 
             
            
             A_D0minus1 = A_d( Dkminus1.L3, Dkminus1.Qin, Dkminus1.QRD, Dkminus1.QN, Dkminus1.QN);
             B_D0minus1 = B_d( Dkminus1.L3, Dkminus1.Qin, Dkminus1.QRD, Dkminus1.QN, Dkminus1.QN);
             f_d0minus1 = f_d( Dkminus1.L3, Dkminus1.Qin, Dkminus1.QRD, Dkminus1.QN, Dkminus1.QN); 
               %Del A , Del B construction for term Del A x(k-1)+ Del B u(k-1)
             DelA{i} = [A_D_K{i}-A_D0minus1; C_r*(A_D_K{i}-A_D0minus1)];
             DelB{i} = [B_D_K{i}-B_D0minus1; C_r*(B_D_K{i}-B_D0minus1)];
              delfd{i,1} = [fd{i}-f_d0minus1; C_r*(fd{i}-f_d0minus1)];  
         else
             DelA{i} = [A_D_K{i}-A_D_K{i-1}; C_r*(A_D_K{i}-A_D_K{i-1})];
             DelB{i} = [B_D_K{i}-B_D_K{i-1}; C_r*(B_D_K{i}-B_D_K{i-1})];
             delfd{i} = [fd{i}-fd{i-1}; C_r*(fd{i}-fd{i-1})];  
         end
         
        if i==1
            BigTerm{i,1} = [A_ext{i}  B_ext{i} ];
            phi_gamForm{i,1} = C_ext*BigTerm{i,1} ;
            BigTerm_2{i,1} = [A_ext{i}  DelA{i} ];
            phi_LambdaForm{i,1} = C_ext*BigTerm_2{i,1};
            BigTerm_3{i,1} = [A_ext{i}  DelB{i} ];
            phi_UpsilonForm{i,1} = C_ext*BigTerm_3{i,1};
            
            BigTerm4{i,1} =  [A_ext{i}  delfd{i} ];
            phi_XiForm{i,1} = C_ext*BigTerm4{i,1};
        else
            BigTerm{i,1} = [A_ext{i}*BigTerm{i-1,1}  B_ext{i}];
            phi_gamForm{i,1} = C_ext*BigTerm{i,1};
            BigTerm_2{i,1} = [A_ext{i}*BigTerm_2{i-1,1}  DelA{i} ];
            phi_LambdaForm{i,1} = C_ext*BigTerm_2{i,1};
            BigTerm_3{i,1} = [A_ext{i}*BigTerm_3{i-1,1}  DelB{i} ];
            phi_UpsilonForm{i,1} = C_ext*BigTerm_3{i,1};
            
             BigTerm4{i,1} =  [A_ext{i}*BigTerm4{i-1,1}  delfd{i} ];
              phi_XiForm{i,1} = C_ext*BigTerm4{i,1};
        end
         Phi ( (i-1)*(N_nz)+1:N_nz*i, 1:(N_nz+nx) ) = phi_gamForm{i,1}(:, 1:(N_nz+nx));
         Gamma(  (i-1)*N_nz+1:N_nz*i, 1:nu*i  ) =  phi_gamForm{i,1}(:, (N_nz+nx)+1:(N_nz+nx)+ nu*i  );
         Lambda((i-1)*N_nz+1:N_nz*i, 1:(nx)*i) = phi_LambdaForm{i,1}(:, (N_nz+nx)+1:(N_nz+nx)+ nx*i  );
         Upsilon((i-1)*N_nz+1:N_nz*i, 1:nu*i) = phi_UpsilonForm{i,1}(:, (N_nz+nx)+1:(N_nz+nx)+ nu*i  );
         XI((i-1)*N_nz+1:N_nz*i, 1:i) = phi_XiForm{i,1}(:, (N_nz+nx)+1:(N_nz+nx)+ i  );

     end
     

 else
     
   
         % Utilising the concept of A_{i+1|k-1} = A_{i|k} ,...
     % A_{i+N|k-1}=A_{i+N-1|k} for the purpose of constructing the Phi, Gamma, Eta_del
     % for building Prediction matrices
%          A_D_K{k+N_per-1} = double(subs(A_d,subvec,d));
%          A_D_K{k+N_per-1} = double(subs(B_d,subvec,d));
         
         A_D_K{1}=[]; 
         B_D_K{1} = [];
         fd{1}=[];
         A_D_K = circshift(A_D_K,-1);
         B_D_K = circshift(B_D_K,-1);
         fd = circshift(fd,-1);
         fd{N_per} = f_d(Dist_Prof.L3(k),Dist_Prof.Qin(k),Dist_Prof.QRD(k),Dist_Prof.QN(:,k),Dist_Prof.QN(:,k)); %double(subs(A_d,subvec,QQ'));
         A_ext{1}=[]; 
         B_ext{1} = [];
         DelA{1} = [];
         DelB{1} = [];
         delfd{1} = [];
         A_ext = circshift(A_ext,-1);
         B_ext = circshift(B_ext,-1);
         DelA = circshift(DelA,-1);
         DelB = circshift(DelB,-1);
         delfd =circshift(delfd,-1);
         
         A_D_K{N_per} = A_d(Dist_Prof.L3(k),Dist_Prof.Qin(k),Dist_Prof.QRD(k),Dist_Prof.QN(:,k),Dist_Prof.QN(:,k)); %double(subs(A_d,subvec,QQ'));
         B_D_K{N_per} =  B_d(Dist_Prof.L3(k),Dist_Prof.Qin(k),Dist_Prof.QRD(k),Dist_Prof.QN(:,k),Dist_Prof.QN(:,k)); %double(subs(A_d,subvec,QQ'));

          A_ext{N_per} = [A_D_K{N_per}       zero_nxnr;...
                     C_r*A_D_K{N_per}   eye_nr];
         B_ext{N_per} = [B_D_K{N_per} ; C_r* B_D_K{N_per}]; 
         
         DelA{N_per} = [A_D_K{N_per}-A_D_K{N_per-1}; C_r*(A_D_K{N_per}-A_D_K{N_per-1})];
         DelB{N_per} = [B_D_K{N_per}-B_D_K{N_per-1}; C_r*(B_D_K{N_per}-B_D_K{N_per-1})];
         delfd{N_per} = [fd{N_per}-fd{N_per-1} ;C_r*(fd{N_per}-fd{N_per-1})];

     for i = 1:1:N_per
         if i==1
            BigTerm{i,1} = [A_ext{i}  B_ext{i} ];
            phi_gamForm{i,1} = C_ext*BigTerm{i,1} ;
            BigTerm_2{i,1} = [A_ext{i}  DelA{i} ];
            phi_LambdaForm{i,1} = C_ext*BigTerm_2{i,1};
            BigTerm_3{i,1} = [A_ext{i}  DelB{i} ];
            phi_UpsilonForm{i,1} = C_ext*BigTerm_3{i,1};
            BigTerm4{i,1} =  [A_ext{i}  delfd{i} ];
            phi_XiForm{i,1} = C_ext*BigTerm4{i,1};
         else
            BigTerm{i,1} = [A_ext{i}*BigTerm{i-1,1}  B_ext{i}];
            phi_gamForm{i,1} = C_ext*BigTerm{i,1};
            BigTerm_2{i,1} = [A_ext{i}*BigTerm_2{i-1,1}  DelA{i} ];
            phi_LambdaForm{i,1} = C_ext*BigTerm_2{i,1};
            BigTerm_3{i,1} = [A_ext{i}*BigTerm_3{i-1,1}  DelB{i} ];
            phi_UpsilonForm{i,1} = C_ext*BigTerm_3{i,1};
            BigTerm4{i,1} =  [A_ext{i}*BigTerm4{i-1,1}  delfd{i} ];
            phi_XiForm{i,1} = C_ext*BigTerm4{i,1};
         end
         Phi ( (i-1)*(N_nz)+1:N_nz*i, 1:(N_nz+nx) ) = phi_gamForm{i,1}(:, 1:(N_nz+nx));
         Gamma(  (i-1)*N_nz+1:N_nz*i, 1:nu*i  ) =  phi_gamForm{i,1}(:, (N_nz+nx)+1:(N_nz+nx)+ nu*i  );
         Lambda((i-1)*N_nz+1:N_nz*i, 1:(nx)*i) = phi_LambdaForm{i,1}(:, (N_nz+nx)+1:(N_nz+nx)+ nx*i  );
         Upsilon((i-1)*N_nz+1:N_nz*i, 1:nu*i) = phi_UpsilonForm{i,1}(:, (N_nz+nx)+1:(N_nz+nx)+ nu*i  );
         XI((i-1)*N_nz+1:N_nz*i, 1:i) = phi_XiForm{i,1}(:, (N_nz+nx)+1:(N_nz+nx)+ i  );
     end
     
 end

     
  Affine_Eta_del = XI*ones(size(XI,2),1);

 %%
 if k==1
    delxp(:,k)= xp(:,k)-x_m1;
    xp(:,k)=xp(:,k);
     u_kpre = zeros(nu,1);% u_max;%
     y0 = C_r*x0;
     yh(:,1)=y0;
     X(:,k) = [delxp(:,k);yh(:,k)];
 else

     u_kpre = uh(:,k-1);
 end
    %%  Rate constraint
bi          = [ -del_u_min;
                     del_u_max;
                     -u_min + u_kpre;
                     u_max- u_kpre;
                     -Ymin;
                     Ymax];
bN          = [zeros(4*nu,1); -Ymin; Ymax];
c = repmat(bi,N_per,1);
c = [c;bN];

L = Mcal*Gamma+Ecal;
W = - (Dcal*C_ext+Mcal*Phi);


 % [K_mpc,P]=dlqr(A_ext{N_per},B_ext{N_per},C_ext'*Q*C_ext,R);
P=50000*Q;%C_ext*P*C_ext';
% Creation of constraint weighing matrices 
[ Psi, Omega ] = QRPN2PsiOmega( Q,R,P,N_per );


G = 2*(Psi+Gamma'*Omega*Gamma);
G1 = (G+G')/2;
    [G2,~] = chol(G1,'lower');
    Linv = G2\eye(size(G1,1));
% Built F Matrix
F = 2*Gamma'*Omega;

%%

    x_k = X(:,k);
    if k==1
    Xpkminus1 = repmat(x_m1,N_per,1);
    Ukminus1 = repmat(u_kpre,N_per,1);
    iA0 = false(size(c+W*X(:,1)));
    else
        Xpkminus1 = repmat(xp(:,k-1),N_per,1);
    Ukminus1 = repmat(u_kpre,N_per,1);
    iA0=A0;
    end
    
 
tic    
      [du_qp1,fval,exitflag] = quadprog(G,F*(Phi*x_k + Lambda*Xpkminus1 + Upsilon*Ukminus1+Affine_Eta_del - Rk),L,W*x_k+c - Mcal*(Lambda*Xpkminus1 +Affine_Eta_del+ Upsilon*Ukminus1),[],[],[],[],[],options);
t_quad(k)=toc;
tic
        [du_qp,status(k),A0,~] = mpcqpsolver(Linv,F*(Phi*x_k + Lambda*Xpkminus1 + Upsilon*Ukminus1+Affine_Eta_del - Rk),-L,-(W*x_k+c - Mcal*(Lambda*Xpkminus1 +Affine_Eta_del+ Upsilon*Ukminus1)),[],zeros(0,1),iA0,opt);
        
 t_mpcqp(k)=toc;
    u_k = u_kpre +du_qp(1:nu); % u(k) = u(k-1)+ \Del u(k)
    uh(:,k) = u_k;
    du(:,k) = du_qp(1:nu);
    %-------Discrete State Space model --------------
    xp(:,k+1) = A_D_K{1}*xp(:,k) + B_D_K{1}*u_k + fd{1};

    yh(:,k+1) = C_r*xp(:,k+1);
    delxp(:,k+1) = xp(:,k+1)-xp(:,k);
    % ------Extended State Space
    X(:,k+1) = [delxp(:,k+1);yh(:,k+1)];
   
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
    axis([0 T 0 1])
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
 