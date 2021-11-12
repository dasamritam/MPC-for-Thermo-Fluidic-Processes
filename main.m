clear all
close all
clc

T = 10;                                        % Simulation time in seconds
Ts = 5e-3;                                      % Ts is the sampling time
T_sim = T/Ts+1;


N_nz  = 15;                                     % Amount of nozzles
T_in  = 25;                                     % inlet temperature
T_airin = 20;                                   % inlet air Temperature
T_ref = 70;                                     % reference temperature

N_per = 10;                                     % prediction horizon MPC controllers

X_max = 85;                                     % maximal state constraint
X_min = 0;                                      % minimal state constraint

U_h_max = 20;                                   % Maximal input for heaters
U_h_min = 0;                                    % Mimimal input for heaters

U_p_max = 1e-2;                                 % Maximal input for Piezo elements
U_p_min = 0;                                    % Mimimal input for Piezo elements

del_U_max = 0.05;
del_U_min = 0.05;

Y_max = 85;
Y_min = 0;

X0 = 65;                                        % inital condition for all states

L0 = 0.003;                                    % initial Level of the Resevoir
N_input = 1;                                    % Enter 1 for Nz+2 inputs and 2 for 2 input actuation

MatSize= 8+2*N_nz;
v = sym('v',[1 MatSize]);
u = sym('u',[MatSize 1]);
x = sym('x',[MatSize 1]);
d = sym('d',[MatSize,1]);

syms  L3 Qin Qrd 
Qn = sym('Qn',[N_nz 1]);
Qo = sym('Qo',[N_nz 1]);

%% Create linear system based on the contructed model using spatial interconnections
    Value_of_variables_struct                                               % set variables 
    Set_Model_Parameters                                                    % set parameters
    Interconnection                                                         % makes interconnections between blockes
    Create_MIMO                                                             % creates the MIMO system
   
    [Dist_Prof,~] = createdisnozzles(N_nz,T,Ts,N_per,L0,Res,Ink);           % creates the flow variables and level of the resevoir
    Qinss=mean(Dist_Prof.Qin);
    Qrdss=mean(Dist_Prof.QRD);
    QNss=mean(Dist_Prof.QN,2);
    QOss=mean(Dist_Prof.QN,2);
    L3ss=mean(Dist_Prof.L3);
    
    Calculate_SS                                                            % Estimates the steady states & inputs using lsqlin function
    create_Linmodel                                                         % creates a linear model using the jacobian around the steady states & inputs
%% Simulation of different controllers. 
% uncomment the desired controller to run the simulation

f = waitbar(0,'Simulating ...');
% non-linear models
%           lpv_PI_Combined
%           lpv_RateMPC_Combined
           lpv_MPC_Combined

% linearized models
%           lpv_PI_combined_lin
%           lpv_RateMPC_Combined_lin
%           lpv_MPC_Combined_lin


% fast estimation for reference tracking mpc !!! only use when N_input = 1
%           lpv_MPC_Combined_FastSS_lin
%           lpv_MPC_Combined_FastSS

close(f)







