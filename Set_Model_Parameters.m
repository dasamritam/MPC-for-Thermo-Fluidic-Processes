%% Set variables to model variables
    % mass
        Mod_Par.m1 = Res.m/2;       
        Mod_Par.m4 = Res.m/2;
        Mod_Par.m5 = Dis.m/2;       
        Mod_Par.m7 = Dis.m/2;
    % specific heat capacity
        Mod_Par.c1 = Res.c;         
        Mod_Par.c4 = Res.c;  
        Mod_Par.c5 = Dis.c;         
        Mod_Par.c7 = Dis.c;
        
    % conduction coefficients
        % area
            Mod_Par.A_cd1 = Glue1.Acd/2;       
            Mod_Par.A_cd4 = Glue1.Acd/2;
            Mod_Par.A_cd5 = Glue2.Acd/2;       
            Mod_Par.A_cd7 = Glue2.Acd/2;
            
            
        % thermal conductivity per layer of glue per subsystem
            Mod_Par.lambdar1 = Glue1.lambda;         
            Mod_Par.lambdar4 = Glue1.lambda;  
            Mod_Par.lambdar5 = Glue2.lambda;          
            Mod_Par.lambdar7 = Glue2.lambda; 
            
        % Glue thickness   
            Mod_Par.lr1 = Glue1.h;         
            Mod_Par.lr4 = Glue1.h; 
            Mod_Par.lr5 = Glue2.h;         
            Mod_Par.lr7 = Glue2.h;
            
    % Convection cofficients
        % Average Heat transfer
            Mod_Par.h_tc1 = Res.h_tc;      
            Mod_Par.h_tc4 = Res.h_tc;
            Mod_Par.h_tc5 = Dis.h_tc;
            Mod_Par.h_tc7 = Dis.h_tc;

        % Area of convection
            Mod_Par.A_cv1 = Res.A_cv;   
            Mod_Par.A_cv4 = Res.A_cv;
            Mod_Par.A_cv5 = Dis.A_cv;   
            Mod_Par.A_cv7 = Dis.A_cv;
            
        for i=1:2:(N_nz+1)*2-1
             eval(['Mod_Par.c' num2str(i+7) '= DFA.c;']);
             eval(['Mod_Par.lambdar' num2str(i+7) '= Glue2.lambda;']);
             eval(['Mod_Par.lr' num2str(i+7) '= Glue2.h;']);
             eval(['Mod_Par.h_tc' num2str(i+7) '= DFA.h_tc;']);
             eval(['Mod_Par.A_cv' num2str(i+7) '= DFA.A_cv;']);
            if i ==1
                eval(['Mod_Par.m' num2str(i+7) '= DFA.m_2;']);
                eval(['Mod_Par.A_cd' num2str(i+7) '= DFA.Acd2;']);
            elseif i== (N_nz+1)*2-1
                eval(['Mod_Par.m' num2str(i+7) '= DFA.m_2;']);
                eval(['Mod_Par.A_cd' num2str(i+7) '= DFA.Acd2;']);
            else
                eval(['Mod_Par.m' num2str(i+7) '= DFA.m_1;']);
                eval(['Mod_Par.A_cd' num2str(i+7) '= DFA.Acd1;']);
            end
        end

    % ink coefficients
            Mod_Par.x_3 =@(L3) L3;   
            Mod_Par.x_6 = Dis.h;
            
            for i=2:2:(N_nz)*2
                eval(['Mod_Par.x_' num2str(i+7) '= DFA.h;']);
            end
            
            Mod_Par.A_ch3 = Res.A_ch ;
            Mod_Par.A_ch6 = Dis.A_ch ;
            
            for i=2:2:(N_nz)*2
                eval(['Mod_Par.A_ch' num2str(i+7) '= DFA.A_ch;']);
            end
%% Set constants a_i b_i and ca_i
    % a_i
        Mod_Par.a_1 = Mod_Par.A_cd1*Mod_Par.lambdar1/Mod_Par.lr1;
        Mod_Par.a_4 = Mod_Par.A_cd4*Mod_Par.lambdar4/Mod_Par.lr4;
        Mod_Par.a_5 = Mod_Par.A_cd5*Mod_Par.lambdar5/Mod_Par.lr5;
        Mod_Par.a_7 = Mod_Par.A_cd7*Mod_Par.lambdar7/Mod_Par.lr7;
        
    % ca_i
     
        Mod_Par.ca_1 = (Res.h*Res.w-Res.h_ch*Res.w_ch)*Res.lambda*2/(Res.l/4);
        Mod_Par.ca_4 = (Res.h*Res.w-Res.h_ch*Res.w_ch)*Res.lambda*2/(Res.l/4);
        Mod_Par.ca_5 = (Dis.h*Dis.w-Dis.h_ch*Dis.w_ch)*Dis.lambda*2/(Dis.l/4);
        Mod_Par.ca_7 = (Dis.h*Dis.w-Dis.h_ch*Dis.w_ch)*Dis.lambda*2/(Dis.l/4);

        for i=1:2:(N_nz+1)*2-1
            eval(['Mod_Par.a_' num2str(i+7) '= Mod_Par.A_cd' num2str(i+7) '*Mod_Par.lambdar' num2str(i+7) '/Mod_Par.lr' num2str(i+7) ';']);
            
             if i== 1
              eval(['Mod_Par.ca_' num2str(i+7) '=(DFA.h*DFA.w-DFA.h_ch*DFA.d_ch)*DFA.lambda/((DFA.l_b+(DFA.l-Glue2.l_ch)/2)/2);']);
             elseif i== (N_nz+1)*2-1
              eval(['Mod_Par.ca_' num2str(i+7) '=(DFA.h*DFA.w-DFA.h_ch*DFA.d_ch)*DFA.lambda/((DFA.l_b+(DFA.l-Glue2.l_ch)/2)/2);']);
             else
              eval(['Mod_Par.ca_' num2str(i+7) '=(DFA.h*DFA.w-DFA.h_ch*DFA.d_ch)*DFA.lambda/(DFA.l_b/2);']); 
             end
        end
    
    % b_i
        Mod_Par.b_2 =@(L3) Air.rho*Air.c*(Res.h_ch-Mod_Par.x_3(L3))*Mod_Par.A_ch3;
        Mod_Par.b_3 =@(L3) Ink.rho*Ink.c*Mod_Par.x_3(L3)*Mod_Par.A_ch3;
        Mod_Par.b_6 = Ink.rho*Ink.c*Mod_Par.x_6*Mod_Par.A_ch6;
        
        for i=2:2:(N_nz)*2
            eval(['Mod_Par.b_' num2str(i+7) '= Ink.rho*Ink.c*Mod_Par.x_' num2str(i+7) '*Mod_Par.A_ch' num2str(i+7) ';']);
        end



        