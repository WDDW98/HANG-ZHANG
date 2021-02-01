% -----------------------
%
% SM_C.m
% Scenario Modelling with Correlation 
%
% Functions include scenario creation, preview scenario, correlation,
% bridge response simulaton.
%
% Input:    V_L1+2_2008_cws.mat
%
% Hang Zhang
% 18/01/2021

clc, clear, close all

% --------------------------------------------
% **** User Inputs for Arrival Process ****
% --------------------------------------------

disp('Loading file...')
load('V_L1+2_2008_cws.mat') % Load Tennessee WIM data, when you load your own data, 
disp('File loaded')         % please change variable name as Veh


% Pick on hour traffic 
Hr = 16;                    % 16pm

% Define Probability Lane Change Car / Truck for 2/4/6 x Length
ParSet = 1;                 % 0 / 1 / 2 / 3 / 10 
                            % 0 for no lane change, 1 for gradually reduce Lane Changing Probabilities 
                            % 10 for immediately Lane Changing Probabilities 

% GAP TYPES
GapType = 1;                % 0 => Constant; 1 => Bailey Distrib 2 => Constant Car-Car

% Save / Display
SV = 0;                     % = 1 => save L1 and L2 Vehicles; = 0 do not save
Dsp = 0;                    % Dsp = 0 => don't display text at each lane change; 1 => display
Cnt = 0;                    % Count HGV / non HGV before and after   



% --------------------------------------------
% **** User Inputs for Creating Scenarios ****
% --------------------------------------------

% Do you need create scenarios

Create_S = 1;                         % 1 => yes, 0 => scenarios have been created

% Define the lenght of scenario
SN  = 150;                            % Any integers 

% Save scenarios  
SV1 = 1;                              % 1 => yes; 0 => no

% Preview scenarios 
PV  = 1;                              % 1 => yes; 0 => no    


% --------------------------------------
% **** User Inputs for correlation ****
% --------------------------------------

% Do you need apply correlation bewteen scenario 
Cor = 1;                               % 1 => yes; 0 => no    

% Bin width of scenario intensity 
Interval= 2;                           % 1 for 1 kN/m, 2 for 2 kN/m. 
                                       % It decides how many scenario in each bin.

up_control = 2;                        % 1 for add 1 kN/m to the top bin, 2 for add 2 kN/m to the top bin.
                                       % Some data has very less scenarios in top bins, 
                                       % this factor is for expanding to bin.

% Shape factor of gamma distribution 
u = 0.25;                              % 0.25 for add 0.25 to defult shape factor, it can be -0.5 to 0.5
                                       % A positive value of u will increase mean value and upper
                                       % tail value of distribtuion 

% -----------------------------------------------------------------------------
% **** Genaration of New Scenarios ****
% -----------------------------------------------------------------------------                                      
 
No_Noise = 5000;                       % how many scenario is going to generate

GapType = 1;                           % 0 = 1.5m constant gap, 
                                       % 1 = 18kph gap, 
                                       % 2 = car - car full stop gap
                                       
GVW_bin = 1;                           % GVW bin 1=1m,                                    
ub_GVW_bin = 32;                       % the last GVW bin ie. 32 means 32m to max
                                                                             
Max_No_HGV_change = 12;                % in one scenario, The max number of HGVs can be 
                                       % changed 
                                       
T_noise = 3;                           % how many times noise apply, 3 = 3*No_Noise
                                    
                               
% -----------------------------------------------------------------------------
% **** User Inputs for Simulating bridge response (total load for example) ****
% -----------------------------------------------------------------------------

% Simualtions days 
Workdays = 250*2;                       % 250 * simulation years

% Span and precision
Span = 1000;                            % FIXED Beam - 500 1000 2000m span
dx = 1;                                 % dx = 1m
    
% Total load influence line function
x = [0:dx:Span]';                       % x,y can be changed to any influence line function
y = ones(length(x),1);

% load intensity
LI = 1;                                 % 1 => load intensity; 0 => load  

% Return period 
Return = 1000;                          % 1000 years or 75 years                          

%------------------------------------------------------------------
[Veh_L1,Veh_L2] = HZ100_MS_Lite(ParSet,GapType,Dsp,Cnt,Hr,SV,Veh);
Veh=[];
disp('Arrival process applied')

if Create_S == 1 
    
[L1_SC, L2_SC, V_L2_car] = SC(SN,SV1);
disp('Scenario created')

elseif Create_S == 0
    
load('L1_SC.mat')
load('L2_SC.mat')
load('V_L2_car.mat')

end

if PV  == 1
Preview_SC (L1_SC,L2_SC);
end

%--------------------------------------------------------------------    
disp('Analyzing Correlation') 
[Inten_SM] = Weight_SC (L1_SC,L2_SC);
[cell_b]= Corr(Inten_SM,Interval,up_control);
[p_groupgamma] = Cor_Group_CERI(cell_b,u,Inten_SM);
disp('Correlation Analyzed')
%-------------------------------------------------------------------- 

disp('Create Noise')   
[Group_GVW] = Noise_GVW(Veh_L1,Veh_L2,ub_GVW_bin,GVW_bin);
[SC_L1_noise, SC_L2_noise] = noise_HZ2 (L1_SC,L2_SC,Group_GVW,No_Noise,GapType,Max_No_HGV_change);
disp('Noise Created')

%-------------------------------------------------------------------- 

[DaMax_SM_C] = SM_C_HZ(p_groupgamma,L1_SC,L2_SC,Inten_SM,V_L2_car,Veh_L1,Veh_L2,Span,dx,x,y,LI,Workdays,SC_L1_noise,SC_L2_noise,T_noise);
[BR] = Extrapolation(Return,DaMax_SM_C);
    