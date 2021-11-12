clear all
close all
clc

syms L3 T2 T3 T4 T6 T1 Qin Qrd
syms h_c1 S1 L_res h_c4 airrho T_airin A_ch airc T_in inkrho inkc

f= (h_c1 *S1*(L_res-L3)*(T1-T2) + h_c4 *S1*(L_res-L3)*(T4-T2) + airrho*airc*Qin*(-T2) +airrho*airc*Qrd*(T_airin-T2) +h_c1*A_ch*(T3-T2))/(airrho*airc*(L_res-L3)*A_ch);
h = (h_c1 *S1*(L3)*(T1-T3) + h_c4 *S1*(L3)*(T4-T3) + inkrho*inkc*Qrd*(T6-T3) +inkrho*inkc*Qin*(T_in-T2) +h_c1*A_ch*(T2-T3))/(inkrho*inkc*(L3)*A_ch);
dif_f=simplify(diff(f,L3));
dif_h=simplify(diff(h,L3));