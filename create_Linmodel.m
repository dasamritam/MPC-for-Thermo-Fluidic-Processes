

A_non_lin =@(L3,Qin,Qrd,Qn,Qo) Model_cl.A(L3,Qin,Qrd,Qn,Qo);
f_non_lin =@(L3,Qin,Qrd,Qn,Qo) Model_cl.f(L3,Qin,Qrd,Qn,Qo);
select = @(M,r,c) M(r,c);
%N_input = input('Enter 1 for 5 inputs and 2 for 2 input actuation ');
if N_input == 1
B_non_lin =@(L3,Qin,Qrd,Qn,Qo) [select(Model_cl.B(L3,Qin,Qrd,Qn,Qo),1:MatSize,1)+select(Model_cl.B(L3,Qin,Qrd,Qn,Qo),1:MatSize,4) (select(Model_cl.B(L3,Qin,Qrd,Qn,Qo),1:MatSize,5)+select(Model_cl.B(L3,Qin,Qrd,Qn,Qo),1:MatSize,7))];    
    for i =1:1:N_nz
        eval(['B_non_lin=@(L3,Qin,Qrd,Qn,Qo)[B_non_lin(L3,Qin,Qrd,Qn,Qo) select(Model_cl.B(L3,Qin,Qrd,Qn,Qo),1:MatSize,' num2str(7+i*2) ')];']);
    end
    
else
B_non_lin =@(L3,Qin,Qrd,Qn,Qo) [select(Model_cl.B(L3,Qin,Qrd,Qn,Qo),1:MatSize,1)+select(Model_cl.B(L3,Qin,Qrd,Qn,Qo),1:MatSize,4) (select(Model_cl.B(L3,Qin,Qrd,Qn,Qo),1:MatSize,5)+select(Model_cl.B(L3,Qin,Qrd,Qn,Qo),1:MatSize,7))];    
end

TC = sym('TC',[MatSize 1]);
UC = sym('UC',[nu 1]);

fff=A_non_lin(L3,Qin,Qrd,Qn,Qo)*TC+B_non_lin(L3,Qin,Qrd,Qn,Qo)*UC+f_non_lin(L3,Qin,Qrd,Qn,Qo);

A_lin_eqn=simplify(jacobian(fff,TC));
B_lin_eqn=simplify(jacobian(fff,UC));
f_lin_eqn=simplify(jacobian(fff,[Qin; Qrd; Qn; Qo; L3]));

subvec = [TC; UC; Qin; Qrd ;Qn ;Qo ;L3];

subin  = [Tss; Us; Qinss; Qrdss; QNss; QOss; L3ss];

        Lin.A   = double(subs(A_lin_eqn,subvec,subin));
        Lin.B   = double(subs(B_lin_eqn,subvec,subin));
        Lin.f   = double(subs(f_lin_eqn,subvec,subin));
 