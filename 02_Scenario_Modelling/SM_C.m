% -----------------------
%
% SM_C.m
% Scenario Modelling with Correlation 
%
% Functions include scenario creation, preview scenario, correlation,
% bridge response simulaton.
%
% Input:    Veh_L1 Veh_L2
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

% -----------------------------------------------------------------
% **** User Inputs for Simulating bridge response (total load for example) ****
% -----------------------------------------------------------------

% Simualtions days 
Workdays = 250;

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
    
[L1_SC, L2_SC, V_L2_car] = SC(SN,SV1)
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
    
[Inten_SM] = Weight_SC (L1_SC,L2_SC);

[cell_b]= Corr(Inten_SM,Interval,up_control);

[p_groupgamma] = Cor_Group_CERI(cell_b,u,Inten_SM);

disp('Correlation applied')
%--------------------------------------------------------------------    

[DaMax_SM_C] = SM_C_HZ(p_groupgamma,L1_SC,L2_SC,Inten_SM,V_L2_car,Veh_L1,Veh_L2,Span,dx,x,y,LI,Workdays)

[BR] = Extrapolation(Return,DaMax_SM_C);
    