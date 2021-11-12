% set parameters & symbols
nx      = 1;        % amount of states per ODE
nu      = 1;        % amount of inputs per ODE
ny      = 1;        % amount of sensor-outputs per ODE
nd      = 1;        % amount of disturbances per ODE
nz      = 1;        % amount of outputs per ODE




%% Adjacency matrix A for parallel interconnection
Ns=8:2:(N_nz+1)*2+6;
Nl=9:2:(N_nz)*2+7;
ns = size(Ns,2);
nl = size(Nl,2);

Adj=zeros(MatSize,MatSize);
Adj(1:7,1:7)=  [0 1 1 1 1 0 0;
                1 0 1 1 0 0 0;
                1 1 0 1 0 1 0;
                1 1 1 0 0 0 1;
                1 0 0 0 0 1 1;
                0 0 1 0 1 0 1;
                0 0 0 1 1 1 0];

for i=5:7
    if i == 5
        if mod(ns,2)==0
            Adj(i,Ns(1:(ns)/2)) = 1;
            %Adj(i,Ns((ns-1)/2+1)) = 1;
        else 
            Adj(i,Ns(1:(ns+1)/2)) = 1;
        end
    elseif i == 6
        Adj(i,Nl) = 1;
    else
        if mod(ns,2)==0
            Adj(i,Ns((ns+2)/2:ns)) = 1;
            %Adj(i,Ns((ns-1)/2+1)) = 1;
        else 
            Adj(i,Ns((ns+1)/2:ns)) = 1;
        end
    end
end
kkk =8:1:N_nz*2+8;
for i = kkk
    if mod(i,2)==1      % if i is odd
            Adj(i,6)=1;
            Adj(i,i-1)=1;
            Adj(i,i+1)=1;        
    else                % if i is even
       if i ==8
            Adj(i,5)=1;
            Adj(i,i+1)=1;
            Adj(i,i+2)=1;
       elseif i == N_nz*2+8
            Adj(i,7)=1;
            Adj(i,i-1)=1;
            Adj(i,i-2)=1;
       elseif mod(N_nz,2)==0 && i == kkk(ceil(end/2))
            Adj(i,i-2)=1;
            Adj(i,i-1)=1;
            Adj(i,5)=1;
            Adj(i,7)=1;
            Adj(i,i+1)=1;
            Adj(i,i+2)=1;
       elseif i<kkk(ceil(end/2))
            Adj(i,i-2)=1;
            Adj(i,i-1)=1;
            Adj(i,5)=1;
            Adj(i,i+1)=1;
            Adj(i,i+2)=1;
       elseif i>kkk(ceil(end/2))
           Adj(i,i-2)=1;
            Adj(i,i-1)=1;
            Adj(i,7)=1;
            Adj(i,i+1)=1;
            Adj(i,i+2)=1;
        end
    end
end

[row,col]=find(Adj);
sizeL=nnz(Adj);

%%  construct matrices
clc
clear Mi
ii=1;
q=1;
for i = 1:MatSize
   nv = sum(Adj(i,:)~=0);
   nw = sum(Adj(:,i)~=0);
   
   eval(['Mi.B' num2str(i) '_xd = zeros(nx,nd);']);
   eval(['Mi.B' num2str(i) '_xu = zeros(nx,nu);']);
   
   eval(['Mi.A' num2str(i) '_wx = ones(nv,1);']);
   eval(['Mi.A' num2str(i) '_wv = zeros(nw,nv);']);
   eval(['Mi.B' num2str(i) '_wd = zeros(nw,nd);']);
   eval(['Mi.B' num2str(i) '_wu = zeros(nw,nu);']);
   
   eval(['Mi.C' num2str(i) '_zx = zeros(nz,nx);']);
   eval(['Mi.C' num2str(i) '_zv = zeros(nz,nv);']);
   eval(['Mi.D' num2str(i) '_zd = zeros(nz,nd);']);
   eval(['Mi.D' num2str(i) '_zu = zeros(nz,nu);']);
   
   eval(['Mi.C' num2str(i) '_yx = 1;']);
   eval(['Mi.C' num2str(i) '_yv = zeros(ny,nv);']);
   eval(['Mi.D' num2str(i) '_yd = zeros(ny,nd);']);
   eval(['Mi.D' num2str(i) '_yu = zeros(ny,nu);']);
   
   eval(['Mi.F' num2str(i) '_xaff = zeros(nx,1);']);
   eval(['Mi.F' num2str(i) '_waff = zeros(nw,1);']);
   eval(['Mi.F' num2str(i) '_zaff = zeros(nz,1);']);
   eval(['Mi.F' num2str(i) '_yaff = zeros(ny,1);']);
   
   if i == 1  
        Mi.A1_xx =@(L3) - (Mod_Par.a_1 + Res.h_tca*Res.P*(Res.h-L3)+Mod_Par.h_tc1*Res.P*L3 + Mod_Par.ca_1)/(Mod_Par.m1*Mod_Par.c1); 
        Mi.A1_xv =@(L3)[Res.h_tca*Res.P*(Res.h-L3) Mod_Par.h_tc1*Res.P*L3  Mod_Par.ca_1  Mod_Par.a_1]/(Mod_Par.m1*Mod_Par.c1);  
        Mi.B1_xu = 1/(2*Mod_Par.m1*Mod_Par.c1);
   elseif i == 2  
        Mi.A2_xx =@(L3,Qin,Qrd) - (Res.h_tca *Res.P*(Res.h-L3) + Res.h_tca *Res.P*(Res.h-L3) + Air.rho*Air.c*Qin +Air.rho*Air.c*Qrd +2*Ink.lambda*Res.A_ch/L3)/(Mod_Par.b_2(L3));
        Mi.A2_xv =@(L3,Qin,Qrd) [Res.h_tca *Res.P*(Res.h-L3)  2*Ink.lambda*Res.A_ch/L3  Res.h_tca *Res.P*(Res.h-L3)]/Mod_Par.b_2(L3); 
        Mi.F2_xaff =@(L3,Qin,Qrd) Air.rho*Air.c*Qrd*T_airin/(Mod_Par.b_2(L3));
   elseif i == 3
        Mi.A3_xx =@(L3,Qin,Qrd) - (Mod_Par.h_tc1 *Res.P*L3 + Mod_Par.h_tc4 *Res.P*L3 + Ink.rho*Ink.c*Qin +Ink.rho*Ink.c*Qrd + 2*Air.lambda*Res.A_ch/(Res.h-L3))/(Mod_Par.b_3(L3));
        Mi.A3_xv =@(L3,Qin,Qrd) [Mod_Par.h_tc1 *Res.P*L3  2*Air.lambda*Res.A_ch/(Res.h-L3)  Mod_Par.h_tc4 *Res.P*L3    Ink.rho*Ink.c*Qrd]/Mod_Par.b_3(L3);
        Mi.B3_xd = zeros(nx,2);
        Mi.B3_wd = zeros(nw,2);
        Mi.D3_zd = zeros(nz,2);
        Mi.D3_yd = zeros(ny,2);
        Mi.F3_xaff =@(L3,Qin,Qrd) Ink.rho*Ink.c*Qin*T_in/(Mod_Par.b_3(L3));
   elseif i == 4
        Mi.A4_xx =@(L3) - (Mod_Par.a_4 + Res.h_tca*Res.P*(Res.h-L3) + Mod_Par.h_tc1*Res.P*L3 + Mod_Par.ca_4)/(Mod_Par.m4*Mod_Par.c4); 
        Mi.A4_xv =@(L3) [Mod_Par.ca_4  Res.h_tca*Res.P*(Res.h-L3) Mod_Par.h_tc1*Res.P*L3  Mod_Par.a_4]/(Mod_Par.m4*Mod_Par.c4); 
        Mi.B4_xu = 1/(2*Mod_Par.m4*Mod_Par.c4);
   elseif i == 5
        Mi.A5_xx = - (Mod_Par.a_1 + Mod_Par.h_tc5 *Mod_Par.A_cv5 + Mod_Par.ca_5)/(Mod_Par.m5*Mod_Par.c5);
        Mi.A5_xv = [Mod_Par.a_1 Mod_Par.h_tc5*Mod_Par.A_cv5  Mod_Par.ca_5]/(Mod_Par.m5*Mod_Par.c5);
           if mod(N_nz,2)==1      % if i is odd 
                for k = 1:(N_nz+1)/2
                    eval(['Mi.A5_xx = Mi.A5_xx - (Mod_Par.a_' num2str(row(k+2+ii)) '/(Mod_Par.m5*Mod_Par.c5));']);
                    eval(['Mi.A5_xv = [Mi.A5_xv  (Mod_Par.a_' num2str(row(k+2+ii)) '/(Mod_Par.m5*Mod_Par.c5))];']);
                end
           else
                for k = 1:(N_nz)/2
                    eval(['Mi.A5_xx = Mi.A5_xx - (Mod_Par.a_' num2str(row(k+2+ii)) '/(Mod_Par.m5*Mod_Par.c5));']);
                    eval(['Mi.A5_xv = [Mi.A5_xv  (Mod_Par.a_' num2str(row(k+2+ii)) '/(Mod_Par.m5*Mod_Par.c5))];']);
                end
                eval(['Mi.A5_xx = Mi.A5_xx - (Mod_Par.a_' num2str(row(k+3+ii)) '/(2*Mod_Par.m5*Mod_Par.c5));']);
                eval(['Mi.A5_xv = [Mi.A5_xv  (Mod_Par.a_' num2str(row(k+3+ii)) '/(2*Mod_Par.m5*Mod_Par.c5))];']);
           end
        Mi.B5_xu = 1/(2*Mod_Par.m5*Mod_Par.c5);
   elseif i == 6
        Mi.A6_xx =@(Qrd,Qn) - (Mod_Par.h_tc5 * Mod_Par.A_cv5 + Mod_Par.h_tc7 * Mod_Par.A_cv7 + Ink.rho*Ink.c*(Qrd+sum(Qn)))/(Mod_Par.b_6);
        Mi.A6_xv =@(Qrd,Qn) [Ink.rho*Ink.c*Qrd Mod_Par.h_tc5 * Mod_Par.A_cv5 Mod_Par.h_tc7 * Mod_Par.A_cv7 Ink.rho*Ink.c*Qn']/(Mod_Par.b_6);
   elseif i == 7
        Mi.A7_xx = - (Mod_Par.a_4 + Mod_Par.h_tc7 *Mod_Par.A_cv7 + Mod_Par.ca_7)/(Mod_Par.m7*Mod_Par.c7);
        Mi.A7_xv = [Mod_Par.a_4 Mod_Par.ca_7 Mod_Par.h_tc7*Mod_Par.A_cv7]/(Mod_Par.m7*Mod_Par.c7);
             if mod(N_nz,2)==1      % if i is odd 
                for k = 1:(N_nz+1)/2
                    eval(['Mi.A7_xx = Mi.A7_xx - (Mod_Par.a_' num2str(row(k+2+ii)) '/(Mod_Par.m7*Mod_Par.c7));']);
                    eval(['Mi.A7_xv = [Mi.A7_xv  (Mod_Par.a_' num2str(row(k+2+ii)) '/(Mod_Par.m7*Mod_Par.c7))];']);
                end
             else
                for k = 1:(N_nz)/2
                    eval(['Mi.A7_xx = Mi.A7_xx - (Mod_Par.a_' num2str(row(k+2+ii)) '/(Mod_Par.m7*Mod_Par.c7));']);
                    eval(['Mi.A7_xv = [Mi.A7_xv  (Mod_Par.a_' num2str(row(k+2+ii)) '/(Mod_Par.m7*Mod_Par.c7))];']);
                end
                eval(['Mi.A7_xx = Mi.A7_xx - (Mod_Par.a_' num2str(row(k+3+ii)) '/(2*Mod_Par.m7*Mod_Par.c7));']);
                eval(['Mi.A7_xv = [Mi.A7_xv  (Mod_Par.a_' num2str(row(k+3+ii)) '/(2*Mod_Par.m7*Mod_Par.c7))];']);
             end
        Mi.B7_xu = 1/(2*Mod_Par.m7*Mod_Par.c7); 
   elseif i == 8
        Mi.A8_xx = - (Mod_Par.a_8 + Mod_Par.h_tc8 *Mod_Par.A_cv8 + Mod_Par.ca_8)/(Mod_Par.m8*Mod_Par.c8);
        Mi.A8_xv = [Mod_Par.a_8 Mod_Par.h_tc8*Mod_Par.A_cv8 Mod_Par.ca_8]/(Mod_Par.m8*Mod_Par.c8);  
   elseif i == MatSize
            eval(['Mi.A' num2str(i) '_xx = - (Mod_Par.a_' num2str(i) ' + Mod_Par.h_tc' num2str(i) ' *Mod_Par.A_cv' num2str(i) ' + Mod_Par.ca_' num2str(i) ')/(Mod_Par.m' num2str(i) '*Mod_Par.c' num2str(i) ');']);
            eval(['Mi.A' num2str(i) '_xv = [Mod_Par.a_' num2str(i) ' Mod_Par.ca_' num2str(i) ' Mod_Par.h_tc' num2str(i) '*Mod_Par.A_cv' num2str(i) ']/(Mod_Par.m' num2str(i) '*Mod_Par.c' num2str(i) ');']);
   elseif i>8 && i<MatSize
        if mod(i,2)==1      % if i is odd
             eval(['Mi.A' num2str(i) '_xx=@(Qn,Qo) - (Mod_Par.h_tc' num2str(i-1) '*Mod_Par.A_cv' num2str(i-1) '+ Mod_Par.h_tc' num2str(i+1) '*Mod_Par.A_cv' num2str(i+1) ' + Ink.rho*Ink.c*(Qn(' num2str(q) ')+Qo(' num2str(q) ')))/(Mod_Par.b_' num2str(i) ');']);
             eval(['Mi.A' num2str(i) '_xv=@(Qn,Qo) [Ink.rho*Ink.c*Qn(' num2str(q) ') Mod_Par.h_tc' num2str(i-1) '*Mod_Par.A_cv' num2str(i-1) ' Mod_Par.h_tc' num2str(i+1) '*Mod_Par.A_cv' num2str(i+1) ' ]/Mod_Par.b_' num2str(i) ';']);
             eval(['Mi.B' num2str(i) '_xu=1/Mod_Par.b_' num2str(i) ';']);
             q=q+1;
        else                % if i is even 
            if mod(N_nz,2)==0      % if i is even
                m=8:2:8+(N_nz+1)*2-2;
                if i == m(ceil(end/2))
                    eval(['Mi.A' num2str(i) '_xv = [Mod_Par.a_' num2str(i) '/2 Mod_Par.a_' num2str(i) '/2 Mod_Par.ca_' num2str(i) ' Mod_Par.h_tc' num2str(i) '*Mod_Par.A_cv' num2str(i) ' Mod_Par.h_tc' num2str(i) '*Mod_Par.A_cv' num2str(i) ' Mod_Par.ca_' num2str(i) ']/(Mod_Par.m' num2str(i) '*Mod_Par.c' num2str(i) ');']);
                    eval(['Mi.A' num2str(i) '_xx=- (Mod_Par.a_' num2str(i) ' +2*Mod_Par.h_tc' num2str(i) '*Mod_Par.A_cv' num2str(i) ' + 2*Mod_Par.ca_' num2str(i) ')/(Mod_Par.m' num2str(i) '*Mod_Par.c' num2str(i) ');']);
                else
                    eval(['Mi.A' num2str(i) '_xx=- (Mod_Par.a_' num2str(i) ' +2*Mod_Par.h_tc' num2str(i) '*Mod_Par.A_cv' num2str(i) ' + 2*Mod_Par.ca_' num2str(i) ')/(Mod_Par.m' num2str(i) '*Mod_Par.c' num2str(i) ');']);
                    eval(['Mi.A' num2str(i) '_xv = [Mod_Par.a_' num2str(i) ' Mod_Par.ca_' num2str(i) ' Mod_Par.h_tc' num2str(i) '*Mod_Par.A_cv' num2str(i) ' Mod_Par.h_tc' num2str(i) '*Mod_Par.A_cv' num2str(i) ' Mod_Par.ca_' num2str(i) ']/(Mod_Par.m' num2str(i) '*Mod_Par.c' num2str(i) ');']);
                end 
            else
                eval(['Mi.A' num2str(i) '_xx=- (Mod_Par.a_' num2str(i) ' +2*Mod_Par.h_tc' num2str(i) '*Mod_Par.A_cv' num2str(i) ' + 2*Mod_Par.ca_' num2str(i) ')/(Mod_Par.m' num2str(i) '*Mod_Par.c' num2str(i) ');']);
                eval(['Mi.A' num2str(i) '_xv = [Mod_Par.a_' num2str(i) ' Mod_Par.ca_' num2str(i) ' Mod_Par.h_tc' num2str(i) '*Mod_Par.A_cv' num2str(i) ' Mod_Par.h_tc' num2str(i) '*Mod_Par.A_cv' num2str(i) ' Mod_Par.ca_' num2str(i) ']/(Mod_Par.m' num2str(i) '*Mod_Par.c' num2str(i) ');']);
            end
        end
   end
   ii=ii+nv;
end

clear i ii k kkk q