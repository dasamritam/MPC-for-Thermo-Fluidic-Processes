function [Dist_Prof,Dkminus1]=createdis(N_nz,T,Ts,N,f)
Off_in = 0.3e-10*Ts;
Amp_in = 0.3e-10*Ts;


t = (0:Ts:T+N);     % seconds
Dist_Prof.Qin = 0;

for i = 1:N_nz
    eval(['Dist_Prof.QN(' num2str(i) ',:)= Amp_in*satlin(sin(2*pi*f*t))+Off_in;']);
    eval(['Dist_Prof.Qin = Dist_Prof.Qin+Dist_Prof.QN(' num2str(i) ',:);']);
    eval(['Dkminus1.QN(' num2str(i) ')=Dist_Prof.QN(' num2str(i) ',1);'])
end
% Dist_Prof.Qin = Amp_in*satlin(sin(2*pi*f*t))+Off_in;
Dist_Prof.QRD = Dist_Prof.Qin;
Dkminus1.Qin = Dist_Prof.Qin(1,1);
Dkminus1.QRD = Dist_Prof.Qin(1,1);
end