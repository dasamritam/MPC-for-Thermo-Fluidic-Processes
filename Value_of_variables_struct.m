%% Parameters of the glue (Heat resistor)
    % shape of glue between Resevoir and Distributor Glue1
        Glue1.l = 70e-3;                                                    % [m] Lenght of the Glue1
        Glue1.w = 10e-3;                                                    % [m] Width of the Glue1
        Glue1.h = 20e-6;                                                    % [m] Height of the Glue1
     % Shape of the channel within Glue1   
        Glue1.l_ch = 50e-3;                                                 % [m] Lenght of the channel within Glue1
        Glue1.w_ch = 1.6e-3;                                                % [m] Width of the channel within Glue1
        Glue1.h_ch = Glue1.h;                                                  % [m] Height of the channel within Glue1
        
        Glue1.Acd        =   Glue1.l*Glue1.w-Glue1.l_ch*Glue1.w_ch;         % [m^2] area of heat conduction between Resevoir and Distributor
        
        Glue1.lambda    =   0.5;                                            % [W/mK]thermal conductivity of first layer of Glue1
    % shape of glue between Distributor and DFA Glue2
        Glue2.l = 45e-3;                                                    % [m] Lenght of the Glue2
        Glue2.w = 2.2e-3;                                                   % [m] Width of the Glue2
        Glue2.h = 5e-6;                                                     % [m] Height of the Glue2
    % Shape of the channel within Glue2     
        Glue2.l_ch = 44e-3;                                                 % [m] Lenght of the channel within Glue2
        Glue2.w_ch = 1.2e-3;                                                % [m] Width of the channel within Glue2
        Glue2.h_ch = Glue2.h;                                                  % [m] Height of the channel within Glue2
     
        Glue2.Acd        =   Glue2.l*Glue2.w-Glue2.l_ch*Glue2.w_ch;         % [m^2] area of heat conduction between Distributor and DFA
        
        Glue2.lambda	=   0.5;                                            % [W/mK]thermal conductivity of second layer of Glue2
%% parameter of the Ink 
        Ink.c               =   2000;                                       % [J/kgK]specific heat capacity of the ink in each block
        Ink.lambda          =   0.14;                                       % [W/mK]thermal conductivity of ink  0.14
        Ink.rho             =   1000;                                       % [kg/m^3]density of the flow
        
%% parameter of the Air 
        Air.c               =   1010;                                       % [J/kgK]specific heat capacity of the ink in each block
        Air.lambda          =   26e-3;                                       % [W/mK]thermal conductivity of ink  0.026
        Air.rho             =   1.225;                                       % [kg/m^3]density of the flow        
%% Parameters of the resevoir 
    % parameters
        Res.lambda      =   100;                                                % [W/mK]thermal conductivity of resevoir
        Res.c           =   710;                                                % [J/kgK]specific heat capacity of resevoir
        Res.rho         =   2260;                                               % [kg/m^3]density of resevoir (graphite)
        Res.Nu          =   7.541;                                              % Nusselt Number
    % block
        Res.l           =   70e-3;                                              % [m] Lenght of the resevoir
        Res.w           =   10e-3;                                              % [m] Width of the resevoir
        Res.h           =   9.9e-3;                                             % [m] Height of the resevoir 
    % channel    
        Res.l_ch        =   50e-3;                                              % [m] lenght of the channel within the resevoir
        Res.w_ch        =   1.6e-3;                                             % [m] width of the channel within the resevoir
        Res.h_ch        =   Res.h;                                              % [m] depth of the channel witing the resevoir
    % volume & mass  
        Res.V           =   Res.l*Res.w*Res.h-Res.l_ch*Res.w_ch*Res.h_ch;       % [m^3] volume of the resevoir 
        Res.m           =   Res.V*Res.rho;                                      % [kg] mass of the resevoir
    % cross section area of the channel in the resevoir
        Res.A_ch        =   Res.l_ch*Res.w_ch;                                  % [m^2]
    % area of heat convection between the resevoir and  the ink
        Res.A_cv        =   (Res.l_ch+Res.w_ch)*Res.h_ch;                       % [m^2]
        Res.P           =   Res.l_ch+Res.w_ch;
    % characteristic length of the resevoir    
        Res.D           =   4*Res.A_ch/(2*Res.l_ch+2*Res.w_ch);                 % [m]
    % average heat transfer coefficient in the Resevoir
        Res.h_tc        =   Res.Nu*Ink.lambda/Res.D;                             %[W/m^2K]
        Res.h_tca       =   Res.Nu*Air.lambda/Res.D;                             %[W/m^2K]
%% Parameters of the distributor
    % parameters
        Dis.lambda      =   100;                                                % [W/mK]   thermal conductivity of distributor
        Dis.c           =   710;                                                % [J/kgK]  specific heat capacity of distributor
        Dis.rho         =   2260;                                               % [kg/m^3] density of distributor (graphite)
        Dis.Nu          =   3.658;                                              % Nusselt Number
    % block
        Dis.l           =   70e-3;                                              % [m] Lenght of the distributor
        Dis.w           =   10e-3;                                              % [m] Width of the distributor
        Dis.h           =   10.5e-3;                                            % [m] Height of the distributor 
    % channel    
        Dis.l_ch        =   44e-3;                                              % [m] lenght of the channel within the distributor
        Dis.w_ch        =   1e-3;                                               % [m] width of the channel within the distributor
        Dis.h_ch        =   Dis.h;                                              % [m] depth of the channel witing the distributor
    % volume & mass  
        Dis.V        	=   Dis.l*Dis.w*Dis.h-Dis.l_ch*Dis.w_ch*Dis.h_ch;       % [m^3] volume of the distributor 
        Dis.m           =   Dis.V*Dis.rho;                                      % [kg] mass of the distributor
    % cross section area of the channel in the distributor
        Dis.A_ch        =   Dis.l_ch*Dis.w_ch;                                  % [m^2]
    % area of heat convection between the distributor and  the ink
        Dis.A_cv        =   (Dis.l_ch+Dis.w_ch)*Dis.h_ch;                       % [m^2]
    % characteristic length of the distributor   
        Dis.D           =   4*Dis.A_ch/(2*Dis.l_ch+2*Dis.w_ch);                 % [m]
    % average heat transfer coefficient in the Resevoir
        Dis.h_tc        =   Dis.Nu*Ink.lambda/Dis.D;                            %[W/m^2K]
%% Parameters of the DFA
    % parameters
        DFA.lambda      =   100;                                                % [W/mK]thermal conductivity of DFA
        DFA.c           =   788;                                                % [J/kgK]specific heat capacity of DFA
        DFA.rho         =   2260;                                               % [kg/m^3]density of DFA (silicon)
        DFA.Nu          =   3.658;                                              % Nusselt Number
    % block
        DFA.l           =   45e-3;                                              % [m] Lenght of the DFA
        DFA.w           =   1.6e-3;                                             % [m] Width of the DFA
        DFA.h           =   18.7e-3;                                            % [m] Height of the DFA
    % channel    
        DFA.d_ch        =   8e-5;                                               % [m] Diameter of the channel within the DFA
        DFA.h_ch        =   DFA.h;                                              % [m] depth of the channel witing the DFA
    % cross section area of the channel in the DFA
        DFA.A_ch        =   pi/4*DFA.d_ch^2;                                    % [m^2]
    % volume & mass  
        DFA.V           =   DFA.l*DFA.w*DFA.h-N_nz*(DFA.A_ch*DFA.h_ch);         % [m^3] volume of the DFA 
        DFA.m           =   DFA.V*DFA.rho;                                      % [kg] mass of the DFA
    % area of heat convection between the DFA and  the ink
        DFA.A_cv        =   (pi*DFA.d_ch/2)*DFA.h_ch;                           % [m^2]
    % characteristic length of the DFA    
        DFA.D           =   4*DFA.A_ch/(pi*DFA.d_ch);                           % [m]
    % average heat transfer coefficient in the Resevoir
        DFA.h_tc        =   DFA.Nu*Ink.lambda/DFA.D;                            %[W/m^2K]
                              
    % Distance & mass between the nozzles
        DFA.l_b     = (Glue2.l_ch)/(N_nz+1);                                      %[m]
        DFA.m_1     = (DFA.l_b*DFA.w*DFA.h-DFA.A_ch*DFA.h_ch)*DFA.rho;          % mass between 2 nozzles
        DFA.m_2     = (DFA.l_b*DFA.w*DFA.h+(DFA.l-Glue2.l_ch)/2*DFA.w*DFA.h-(DFA.A_ch/2)*DFA.h_ch)*DFA.rho;      % mass between 1 nozzle and the edge
            if DFA.l_b-DFA.d_ch  <=  20e-6                                      % When the distance between the nozzles becomes lesser then 20 mu the program stops
                disp('Area between nozzles is too small. Reduce the number of nozzles for this configuration')
                return  
            else
            end
     % heat area for each nozzle 
        DFA.Acd1    =   (Glue2.w-Glue2.w_ch)*DFA.l_b;
        DFA.Acd2    =   (Glue2.w-Glue2.w_ch)*DFA.l_b+Glue2.w*(Glue2.l-Glue2.l_ch)/2;
        
        