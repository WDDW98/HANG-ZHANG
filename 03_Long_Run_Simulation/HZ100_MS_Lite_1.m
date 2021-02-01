% -----------------------
%
% MQ100_MS_Lite_FRB.m
% Micro-Sim Lite - generate data with correct queue for full stop conidtion from 2 lane traffic free flowing traffic
%
%
% 2 lanes traffic (same direction)
% Create queue
% Inter-lane gap / length of current vehicle - Prob of lane change different cars & trucks
%
% MQuilligan
% 28.07.2017 R1
% 15.07.2017 R2
% 20.07.2017 AM update
% 10.01.2018 FRB Specific
% -----------------------




function [SC] = HZ100_MS_Lite_1(GapType,SC,CDF_Gap_Con18)



% -----------------------
% GAP TYPES

if GapType == 0
    ConstGap = 1.5;
    
    % Bailey Distribution Gap Distance - Congested, i.e. 18kph - load data
    
elseif GapType == 1
    load CDF_Gap_Con
    % Define CDF for later use
    CDF = CDF_Gap_Con18(:,2);
    
elseif GapType == 2
    load CDF_Gap
    % Define CDF for later use
    CDF = CDF_Gap(:,2); % Car - car
end

% -----------------------
% **** 2. Load Vehicle data
% -----------------------



% Add 5 new columns to Veh as more data added in this algorithm
SC(size(SC,1),59) = 0;




% -----------------------
% **** 3. Initialise variables
% -----------------------

NoVeh = size(SC,1);

% -----------------------
% Initialise Random Number Generator
rng(3,'twister');   % set seed for Random Number Generator to 1 - can produce predicable results


% -----------------------
% **** 4. Loop through each vehicle
% -----------------------



% -*-*-*-*-*-*-*-*-*-*-*-
% ### MAIN MICRO SIM LIGHT CODE ###
% -*-*-*-*-*-*-*-*-*-*-*-
for i = 1:NoVeh
    
    % -----------------------
    

            
            % Define Gap
            if GapType == 0
                Gap = ConstGap;
                
            elseif GapType == 1
                R = rand(1);
                [~, Ind] = min(abs(CDF-R));
                Gap = CDF_Gap_Con18(Ind,1);
                
            elseif GapType == 2
                R = rand(1);
                [~, Ind] = min(abs(CDF-R));
                Gap = CDF_Gap(Ind,1);  
            end
            % Record new gap
            SC(i,28) =  Gap;
            SC(i,56) =  Gap;

             
          
end

            
        end
        
        
 
    
    



