%% Creating L matrix which connects V = LW
%  The Elements of the matrix has a following structure
%       |v^{i,j}| = |w^{j,i}|  
%       |w^{i,j}| = |v^{j,i}| 
M = zeros(sizeL,sizeL);
for i=1:sizeL
    in1= intersect(find(row(i)==col), find(col(i)==row) );
    eval(['M(' num2str(i) ',' num2str(in1) ') = 1;']);
end

%% Block Diagonal structure of State Space matrices
for i=1:MatSize
    if i == 1 
        eval ( ['Diag_Axx =@(L3,Qin,Qrd,Qn,Qo) Mi.A' num2str(i) '_xx(L3);'])
        eval ( ['Diag_Axv =@(L3,Qin,Qrd,Qn,Qo) Mi.A' num2str(i) '_xv(L3);'])
        eval ( ['Diag_Bxd = Mi.B' num2str(i) '_xd;'])
        eval ( ['Diag_Bxu = Mi.B' num2str(i) '_xu;'])
        eval ( ['Diag_Awx = Mi.A' num2str(i) '_wx;'])
        eval ( ['Diag_Awv = Mi.A' num2str(i) '_wv;'])
        eval ( ['Diag_Bwd = Mi.B' num2str(i) '_wd;'])
        eval ( ['Diag_Bwu = Mi.B' num2str(i) '_wu;'])
        eval ( ['Diag_Czx = Mi.C' num2str(i) '_zx;'])
        eval ( ['Diag_Czv = Mi.C' num2str(i) '_zv;'])
        eval ( ['Diag_Dzd = Mi.D' num2str(i) '_zd;'])
        eval ( ['Diag_Dzu = Mi.D' num2str(i) '_zu;'])
        eval ( ['Diag_Cyx = Mi.C' num2str(i) '_yx;'])
        eval ( ['Diag_Cyv = Mi.C' num2str(i) '_yv;'])
        eval ( ['Diag_Dyd = Mi.D' num2str(i) '_yd;'])
        eval ( ['Diag_Dyu = Mi.D' num2str(i) '_yu;'])

        eval ( ['Diag_Fxaff =@(L3,Qin,Qrd) Mi.F' num2str(i) '_xaff;'])
        eval ( ['Diag_Fwaff = Mi.F' num2str(i) '_waff;'])
        eval ( ['Diag_Fzaff = Mi.F' num2str(i) '_zaff;'])
        eval ( ['Diag_Fyaff = Mi.F' num2str(i) '_yaff;'])
    elseif i == 4 
        eval ( ['Diag_Axx =@(L3,Qin,Qrd,Qn,Qo) blkdiag(Diag_Axx(L3,Qin,Qrd,Qn,Qo), Mi.A' num2str(i) '_xx(L3));'])
        eval ( ['Diag_Axv =@(L3,Qin,Qrd,Qn,Qo) blkdiag(Diag_Axv(L3,Qin,Qrd,Qn,Qo),Mi.A' num2str(i) '_xv(L3));'])
        eval ( ['Diag_Bxd = blkdiag(Diag_Bxd,Mi.B' num2str(i) '_xd);'])
        eval ( ['Diag_Bxu = blkdiag(Diag_Bxu,Mi.B' num2str(i) '_xu);'])
        eval ( ['Diag_Awx = blkdiag(Diag_Awx,Mi.A' num2str(i) '_wx);'])
        eval ( ['Diag_Awv = blkdiag(Diag_Awv,Mi.A' num2str(i) '_wv);'])
        eval ( ['Diag_Bwd = blkdiag(Diag_Bwd,Mi.B' num2str(i) '_wd);'])
        eval ( ['Diag_Bwu = blkdiag(Diag_Bwu,Mi.B' num2str(i) '_wu);'])
        eval ( ['Diag_Czx = blkdiag(Diag_Czx,Mi.C' num2str(i) '_zx);'])
        eval ( ['Diag_Czv = blkdiag(Diag_Czv,Mi.C' num2str(i) '_zv);'])
        eval ( ['Diag_Dzd = blkdiag(Diag_Dzd,Mi.D' num2str(i) '_zd);'])
        eval ( ['Diag_Dzu = blkdiag(Diag_Dzu,Mi.D' num2str(i) '_zu);'])
        eval ( ['Diag_Cyx = blkdiag(Diag_Cyx,Mi.C' num2str(i) '_yx);'])
        eval ( ['Diag_Cyv = blkdiag(Diag_Cyv,Mi.C' num2str(i) '_yv);'])
        eval ( ['Diag_Dyd = blkdiag(Diag_Dyd,Mi.D' num2str(i) '_yd);'])
        eval ( ['Diag_Dyu = blkdiag(Diag_Dyu,Mi.D' num2str(i) '_yu);'])
        % Stack as column vector 
        eval ( ['Diag_Fxaff =@(L3,Qin,Qrd) [Diag_Fxaff(L3,Qin,Qrd);Mi.F' num2str(i) '_xaff];'])
        eval ( ['Diag_Fwaff = [Diag_Fwaff;Mi.F' num2str(i) '_waff];'])
        eval ( ['Diag_Fzaff = [Diag_Fzaff;Mi.F' num2str(i) '_zaff];'])
        eval ( ['Diag_Fyaff = [Diag_Fyaff;Mi.F' num2str(i) '_yaff];'])
        
    elseif i == 2 || i == 3
        eval ( ['Diag_Axx =@(L3,Qin,Qrd,Qn,Qo) blkdiag(Diag_Axx(L3,Qin,Qrd,Qn,Qo), Mi.A' num2str(i) '_xx(L3,Qin,Qrd));'])
        eval ( ['Diag_Axv =@(L3,Qin,Qrd,Qn,Qo) blkdiag(Diag_Axv(L3,Qin,Qrd,Qn,Qo),Mi.A' num2str(i) '_xv(L3,Qin,Qrd));'])
        eval ( ['Diag_Bxd = blkdiag(Diag_Bxd,Mi.B' num2str(i) '_xd);'])
        eval ( ['Diag_Bxu = blkdiag(Diag_Bxu,Mi.B' num2str(i) '_xu);'])
        eval ( ['Diag_Awx = blkdiag(Diag_Awx,Mi.A' num2str(i) '_wx);'])
        eval ( ['Diag_Awv = blkdiag(Diag_Awv,Mi.A' num2str(i) '_wv);'])
        eval ( ['Diag_Bwd = blkdiag(Diag_Bwd,Mi.B' num2str(i) '_wd);'])
        eval ( ['Diag_Bwu = blkdiag(Diag_Bwu,Mi.B' num2str(i) '_wu);'])
        eval ( ['Diag_Czx = blkdiag(Diag_Czx,Mi.C' num2str(i) '_zx);'])
        eval ( ['Diag_Czv = blkdiag(Diag_Czv,Mi.C' num2str(i) '_zv);'])
        eval ( ['Diag_Dzd = blkdiag(Diag_Dzd,Mi.D' num2str(i) '_zd);'])
        eval ( ['Diag_Dzu = blkdiag(Diag_Dzu,Mi.D' num2str(i) '_zu);'])
        eval ( ['Diag_Cyx = blkdiag(Diag_Cyx,Mi.C' num2str(i) '_yx);'])
        eval ( ['Diag_Cyv = blkdiag(Diag_Cyv,Mi.C' num2str(i) '_yv);'])
        eval ( ['Diag_Dyd = blkdiag(Diag_Dyd,Mi.D' num2str(i) '_yd);'])
        eval ( ['Diag_Dyu = blkdiag(Diag_Dyu,Mi.D' num2str(i) '_yu);'])
        % Stack as column vector 
        eval ( ['Diag_Fxaff =@(L3,Qin,Qrd) [Diag_Fxaff(L3,Qin,Qrd);Mi.F' num2str(i) '_xaff(L3,Qin,Qrd)];'])
        eval ( ['Diag_Fwaff = [Diag_Fwaff;Mi.F' num2str(i) '_waff];'])
        eval ( ['Diag_Fzaff = [Diag_Fzaff;Mi.F' num2str(i) '_zaff];'])
        eval ( ['Diag_Fyaff = [Diag_Fyaff;Mi.F' num2str(i) '_yaff];'])
        
    elseif i == 5 || i == 7 
        eval ( ['Diag_Axx =@(L3,Qin,Qrd,Qn,Qo) blkdiag(Diag_Axx(L3,Qin,Qrd,Qn,Qo), Mi.A' num2str(i) '_xx);'])
        eval ( ['Diag_Axv =@(L3,Qin,Qrd,Qn,Qo) blkdiag(Diag_Axv(L3,Qin,Qrd,Qn,Qo),Mi.A' num2str(i) '_xv);'])
        eval ( ['Diag_Bxd = blkdiag(Diag_Bxd,Mi.B' num2str(i) '_xd);'])
        eval ( ['Diag_Bxu = blkdiag(Diag_Bxu,Mi.B' num2str(i) '_xu);'])
        eval ( ['Diag_Awx = blkdiag(Diag_Awx,Mi.A' num2str(i) '_wx);'])
        eval ( ['Diag_Awv = blkdiag(Diag_Awv,Mi.A' num2str(i) '_wv);'])
        eval ( ['Diag_Bwd = blkdiag(Diag_Bwd,Mi.B' num2str(i) '_wd);'])
        eval ( ['Diag_Bwu = blkdiag(Diag_Bwu,Mi.B' num2str(i) '_wu);'])
        eval ( ['Diag_Czx = blkdiag(Diag_Czx,Mi.C' num2str(i) '_zx);'])
        eval ( ['Diag_Czv = blkdiag(Diag_Czv,Mi.C' num2str(i) '_zv);'])
        eval ( ['Diag_Dzd = blkdiag(Diag_Dzd,Mi.D' num2str(i) '_zd);'])
        eval ( ['Diag_Dzu = blkdiag(Diag_Dzu,Mi.D' num2str(i) '_zu);'])
        eval ( ['Diag_Cyx = blkdiag(Diag_Cyx,Mi.C' num2str(i) '_yx);'])
        eval ( ['Diag_Cyv = blkdiag(Diag_Cyv,Mi.C' num2str(i) '_yv);'])
        eval ( ['Diag_Dyd = blkdiag(Diag_Dyd,Mi.D' num2str(i) '_yd);'])
        eval ( ['Diag_Dyu = blkdiag(Diag_Dyu,Mi.D' num2str(i) '_yu);'])
        % Stack as column vector 
        eval ( ['Diag_Fxaff =@(L3,Qin,Qrd) [Diag_Fxaff(L3,Qin,Qrd);Mi.F' num2str(i) '_xaff];'])
        eval ( ['Diag_Fwaff = [Diag_Fwaff;Mi.F' num2str(i) '_waff];'])
        eval ( ['Diag_Fzaff = [Diag_Fzaff;Mi.F' num2str(i) '_zaff];'])
        eval ( ['Diag_Fyaff = [Diag_Fyaff;Mi.F' num2str(i) '_yaff];'])
        
    elseif i == 6
        eval ( ['Diag_Axx =@(L3,Qin,Qrd,Qn,Qo) blkdiag(Diag_Axx(L3,Qin,Qrd,Qn,Qo), Mi.A' num2str(i) '_xx(Qrd,Qn));'])
        eval ( ['Diag_Axv =@(L3,Qin,Qrd,Qn,Qo) blkdiag(Diag_Axv(L3,Qin,Qrd,Qn,Qo),Mi.A' num2str(i) '_xv(Qrd,Qn));'])
        eval ( ['Diag_Bxd = blkdiag(Diag_Bxd,Mi.B' num2str(i) '_xd);'])
        eval ( ['Diag_Bxu = blkdiag(Diag_Bxu,Mi.B' num2str(i) '_xu);'])
        eval ( ['Diag_Awx = blkdiag(Diag_Awx,Mi.A' num2str(i) '_wx);'])
        eval ( ['Diag_Awv = blkdiag(Diag_Awv,Mi.A' num2str(i) '_wv);'])
        eval ( ['Diag_Bwd = blkdiag(Diag_Bwd,Mi.B' num2str(i) '_wd);'])
        eval ( ['Diag_Bwu = blkdiag(Diag_Bwu,Mi.B' num2str(i) '_wu);'])
        eval ( ['Diag_Czx = blkdiag(Diag_Czx,Mi.C' num2str(i) '_zx);'])
        eval ( ['Diag_Czv = blkdiag(Diag_Czv,Mi.C' num2str(i) '_zv);'])
        eval ( ['Diag_Dzd = blkdiag(Diag_Dzd,Mi.D' num2str(i) '_zd);'])
        eval ( ['Diag_Dzu = blkdiag(Diag_Dzu,Mi.D' num2str(i) '_zu);'])
        eval ( ['Diag_Cyx = blkdiag(Diag_Cyx,Mi.C' num2str(i) '_yx);'])
        eval ( ['Diag_Cyv = blkdiag(Diag_Cyv,Mi.C' num2str(i) '_yv);'])
        eval ( ['Diag_Dyd = blkdiag(Diag_Dyd,Mi.D' num2str(i) '_yd);'])
        eval ( ['Diag_Dyu = blkdiag(Diag_Dyu,Mi.D' num2str(i) '_yu);'])
        % Stack as column vector 
        eval ( ['Diag_Fxaff =@(L3,Qin,Qrd) [Diag_Fxaff(L3,Qin,Qrd);Mi.F' num2str(i) '_xaff];'])
        eval ( ['Diag_Fwaff = [Diag_Fwaff;Mi.F' num2str(i) '_waff];'])
        eval ( ['Diag_Fzaff = [Diag_Fzaff;Mi.F' num2str(i) '_zaff];'])
        eval ( ['Diag_Fyaff = [Diag_Fyaff;Mi.F' num2str(i) '_yaff];'])
        
    elseif i>=8 && mod(i,2)==0 % is even
        eval ( ['Diag_Axx =@(L3,Qin,Qrd,Qn,Qo) blkdiag(Diag_Axx(L3,Qin,Qrd,Qn,Qo), Mi.A' num2str(i) '_xx);'])
        eval ( ['Diag_Axv =@(L3,Qin,Qrd,Qn,Qo) blkdiag(Diag_Axv(L3,Qin,Qrd,Qn,Qo),Mi.A' num2str(i) '_xv);'])
        eval ( ['Diag_Bxd = blkdiag(Diag_Bxd,Mi.B' num2str(i) '_xd);'])
        eval ( ['Diag_Bxu = blkdiag(Diag_Bxu,Mi.B' num2str(i) '_xu);'])
        eval ( ['Diag_Awx = blkdiag(Diag_Awx,Mi.A' num2str(i) '_wx);'])
        eval ( ['Diag_Awv = blkdiag(Diag_Awv,Mi.A' num2str(i) '_wv);'])
        eval ( ['Diag_Bwd = blkdiag(Diag_Bwd,Mi.B' num2str(i) '_wd);'])
        eval ( ['Diag_Bwu = blkdiag(Diag_Bwu,Mi.B' num2str(i) '_wu);'])
        eval ( ['Diag_Czx = blkdiag(Diag_Czx,Mi.C' num2str(i) '_zx);'])
        eval ( ['Diag_Czv = blkdiag(Diag_Czv,Mi.C' num2str(i) '_zv);'])
        eval ( ['Diag_Dzd = blkdiag(Diag_Dzd,Mi.D' num2str(i) '_zd);'])
        eval ( ['Diag_Dzu = blkdiag(Diag_Dzu,Mi.D' num2str(i) '_zu);'])
        eval ( ['Diag_Cyx = blkdiag(Diag_Cyx,Mi.C' num2str(i) '_yx);'])
        eval ( ['Diag_Cyv = blkdiag(Diag_Cyv,Mi.C' num2str(i) '_yv);'])
        eval ( ['Diag_Dyd = blkdiag(Diag_Dyd,Mi.D' num2str(i) '_yd);'])
        eval ( ['Diag_Dyu = blkdiag(Diag_Dyu,Mi.D' num2str(i) '_yu);'])
        % Stack as column vector 
        eval ( ['Diag_Fxaff =@(L3,Qin,Qrd) [Diag_Fxaff(L3,Qin,Qrd);Mi.F' num2str(i) '_xaff];'])
        eval ( ['Diag_Fwaff = [Diag_Fwaff;Mi.F' num2str(i) '_waff];'])
        eval ( ['Diag_Fzaff = [Diag_Fzaff;Mi.F' num2str(i) '_zaff];'])
        eval ( ['Diag_Fyaff = [Diag_Fyaff;Mi.F' num2str(i) '_yaff];'])
        
    else
        eval ( ['Diag_Axx =@(L3,Qin,Qrd,Qn,Qo) blkdiag(Diag_Axx(L3,Qin,Qrd,Qn,Qo), Mi.A' num2str(i) '_xx(Qn,Qo));'])
        eval ( ['Diag_Axv =@(L3,Qin,Qrd,Qn,Qo) blkdiag(Diag_Axv(L3,Qin,Qrd,Qn,Qo),Mi.A' num2str(i) '_xv(Qn,Qo));'])
        eval ( ['Diag_Bxd = blkdiag(Diag_Bxd,Mi.B' num2str(i) '_xd);'])
        eval ( ['Diag_Bxu = blkdiag(Diag_Bxu,Mi.B' num2str(i) '_xu);'])
        eval ( ['Diag_Awx = blkdiag(Diag_Awx,Mi.A' num2str(i) '_wx);'])
        eval ( ['Diag_Awv = blkdiag(Diag_Awv,Mi.A' num2str(i) '_wv);'])
        eval ( ['Diag_Bwd = blkdiag(Diag_Bwd,Mi.B' num2str(i) '_wd);'])
        eval ( ['Diag_Bwu = blkdiag(Diag_Bwu,Mi.B' num2str(i) '_wu);'])
        eval ( ['Diag_Czx = blkdiag(Diag_Czx,Mi.C' num2str(i) '_zx);'])
        eval ( ['Diag_Czv = blkdiag(Diag_Czv,Mi.C' num2str(i) '_zv);'])
        eval ( ['Diag_Dzd = blkdiag(Diag_Dzd,Mi.D' num2str(i) '_zd);'])
        eval ( ['Diag_Dzu = blkdiag(Diag_Dzu,Mi.D' num2str(i) '_zu);'])
        eval ( ['Diag_Cyx = blkdiag(Diag_Cyx,Mi.C' num2str(i) '_yx);'])
        eval ( ['Diag_Cyv = blkdiag(Diag_Cyv,Mi.C' num2str(i) '_yv);'])
        eval ( ['Diag_Dyd = blkdiag(Diag_Dyd,Mi.D' num2str(i) '_yd);'])
        eval ( ['Diag_Dyu = blkdiag(Diag_Dyu,Mi.D' num2str(i) '_yu);'])
        % Stack as column vector 
        eval ( ['Diag_Fxaff =@(L3,Qin,Qrd) [Diag_Fxaff(L3,Qin,Qrd);Mi.F' num2str(i) '_xaff];'])
        eval ( ['Diag_Fwaff = [Diag_Fwaff;Mi.F' num2str(i) '_waff];'])
        eval ( ['Diag_Fzaff = [Diag_Fzaff;Mi.F' num2str(i) '_zaff];'])
        eval ( ['Diag_Fyaff = [Diag_Fyaff;Mi.F' num2str(i) '_yaff];'])
        
    end
    
end
%% Cascaded Diagonal interconnection system
%  Xdot = Diag_Axx X + Diag_Axv V + Diag_Bxu U_nr           Xdot = Diag_Axx X + Diag_Axv V + Diag_Bxu U
%   W   = Diag_Awx X + Diag_Awv V + Diag_Bwu U_nr   ==>      W   = Diag_Awx X + Diag_Awv V + Diag_Bwu U
%   Y   = Diag_Cyx X + Di   ag_Cyv V + Diag_Dyu U_nr         Y   = Diag_Cyx X + Diag_Cyv V + Diag_Dyu U
%

Model_cl.A =@(L3,Qin,Qrd,Qn,Qo) Diag_Axx(L3,Qin,Qrd,Qn,Qo) + Diag_Axv(L3,Qin,Qrd,Qn,Qo)*M*inv(eye(sizeL)- Diag_Awv*M)*Diag_Awx;
Model_cl.B =@(L3,Qin,Qrd,Qn,Qo) Diag_Bxu + Diag_Axv(L3,Qin,Qrd,Qn,Qo)*M*inv(eye(sizeL)- Diag_Awv*M)*Diag_Bwu;
Model_cl.Bz =@(L3,Qin,Qrd,Qn,Qo) Diag_Bxd + Diag_Axv*M*inv(eye(sizeL)- Diag_Awv*M)*Diag_Bwd;
Model_cl.C1 =  Diag_Cyx + Diag_Cyv*M*inv(eye(sizeL)- Diag_Awv*M)*Diag_Awx;
Model_cl.D1 =  Diag_Dyu + Diag_Cyv*M*inv(eye(sizeL)- Diag_Awv*M)*Diag_Bwu;
Model_cl.D1z = Diag_Dyd + Diag_Cyv*M*inv(eye(sizeL)- Diag_Awv*M)*Diag_Bwd;
Model_cl.C2 =  Diag_Czx + Diag_Czv*M*inv(eye(sizeL)- Diag_Awv*M)*Diag_Awx;
Model_cl.D2 =  Diag_Dzu + Diag_Czv*M*inv(eye(sizeL)- Diag_Awv*M)*Diag_Bwu;
Model_cl.D2z = Diag_Dzd + Diag_Czv*M*inv(eye(sizeL)- Diag_Awv*M)*Diag_Bwd;

Model_cl.f =@(L3,Qin,Qrd,Qn,Qo)  Diag_Fxaff(L3,Qin,Qrd) + Diag_Axv(L3,Qin,Qrd,Qn,Qo)*M*inv(eye(sizeL)- Diag_Awv*M)*Diag_Fwaff;
Model_cl.g1 = Diag_Fyaff + Diag_Cyv*M*inv(eye(sizeL)- Diag_Awv*M)*Diag_Fwaff;
Model_cl.g2 = Diag_Fzaff + Diag_Czv*M*inv(eye(sizeL)- Diag_Awv*M)*Diag_Fwaff;


%clear Diag_Axx Diag_Awv Diag_Awx Diag_Axv Diag_Bwd Diag_Bwu Diag_Bxd Diag_Bxu Diag_Cyv Diag_Cyx Diag_Czv Diag_Czx  Diag_Dyd Diag_Dyu Diag_Dzd Diag_Dzu Diag_Fwaff Diag_Fxaff Diag_Fyaff Diag_Fzaff i in1
% MIMO Transfer funtion of the form 
% xdot =  A x + B  u + Bz  d  + f
%  y   = C1 x + D1 u + D1z d  + g1
%  z   = C2 x + D2 u + D2z d  + g2



%save('var_closed_SS_sep12','A','B','Bz', 'C1','D1' ,'D1z', 'C2','D2' ,'D2z','f','g1','g2')

%%  uncomment the below commands to display the system in differential eqn
 
% pretty(A*x+B*u +  f)