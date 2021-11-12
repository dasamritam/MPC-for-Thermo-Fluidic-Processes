function [Dist_Prof,Dkminus1]=createdisnozzles(N_nz,T,Ts,N,L0,Res,Ink)

Fmax=100;
TTs=1/Fmax;
Q_drop = 0.6e-12;

Dist_Prof.Qin = 0;
kkk=1;
Dist_Prof.L3(:,1)=L0;
Qin=0;
sys_l3 = c2d(ss(0,[-1/(Res.A_ch) 1/(Res.A_ch)],1,0),Ts);
AA=sys_l3.A;
BB=sys_l3.B;
for ii=0:TTs:T+N*TTs
    Pulse=randi([0, 1], [1,N_nz]);
    
    for iii= 0:1:TTs/Ts-1
        Dist_Prof.QN(:,kkk+iii) = Pulse'*Q_drop*Fmax;
        Dist_Prof.L3(:,kkk+1+iii) = AA*Dist_Prof.L3(:,kkk+iii)+BB*[ones(1,N_nz)*Dist_Prof.QN(:,kkk+iii);Dist_Prof.Qin(:,kkk+iii)];
        if Dist_Prof.L3(:,kkk+iii+1)<=Res.h_ch*0.2
           Dist_Prof.Qin(:,kkk+iii+1) = 15e-9;
        elseif  Dist_Prof.L3(:,kkk+iii+1)>=Res.h_ch*0.8 
            Dist_Prof.Qin(:,kkk+iii+1) = 0;
        else
            Dist_Prof.Qin(:,kkk+iii+1) = Dist_Prof.Qin(:,kkk+iii);
        end
    end
    kkk=kkk+TTs/Ts;
end

Dist_Prof.QRD = sum(Dist_Prof.QN);
Dkminus1.Qin = Dist_Prof.Qin(:,1);
Dkminus1.QRD = Dist_Prof.QRD(:,1);
Dkminus1.QN = Dist_Prof.QN(:,1);
Dkminus1.L3 = Dist_Prof.L3(:,1);
% figure
% plot(t,Dist_Prof.QRD)
% ylim([0 2e-10])
% 
% figure
% plot(t,Dist_Prof.L2(:,1:end-1))
