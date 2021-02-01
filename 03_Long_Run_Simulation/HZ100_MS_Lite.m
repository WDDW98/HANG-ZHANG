% -----------------------
%
% HZ100_MS_Lite.m
% Micro-Sim Lite - generate data with correct queue for full stop conidtion from 2 lane traffic free flowing traffic
%
%
% 2 lanes traffic (same direction)
% Create queue
% Inter-lane gap / length of current vehicle - Prob of lane change different cars & trucks
%
% MQuilligan
% Hang Zhang
% -----------------------





function [Veh_L1,Veh_L2] = HZ100_MS_Lite(ParSet,GapType,Dsp,Cnt,Hr,SV,Veh)



% -----------------------
% GAP TYPES

% Constant Gap Distance
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
Veh(size(Veh,1),59) = 0;
% 1a. Reduce Veh to that Specified Hour
Veh = Veh( Veh(:,6)==Hr, :);            % Select that hour only

% -----------------------
% Define Probability Lane Change Car / Truck for 2/4/6 x Length

% Trucks / HGVs - Cars / nonHGVs   [<1 1-4 >4]
if ParSet == 0                  % Parameter Set 0 - no lane changing allowed
    PLC_T = [0 0 0];
    PLC_C = [0 0 0];
elseif ParSet == 1              % Parameter Set 1 - gradually reduce Lane Changing Probabilities
    PLC_T = [0 0.2 0.4];
    PLC_C = [0 0.4 0.8];
elseif ParSet == 2
    PLC_T = [0 0.1 0.2];
    PLC_C = [0 0.45 0.90];
elseif ParSet == 3
    PLC_T = [0 0.05 0.1];
    PLC_C = [0 0.6 0.9];
elseif ParSet == 10
    PLC_T = [0 1 1];
    PLC_C = [0 1 1];
end


% -----------------------
% **** 3. Initialise variables
% -----------------------

NoVeh = size(Veh,1);

% -----------------------
% Initialise Random Number Generator
rng(3,'twister');   % set seed for Random Number Generator to 1 - can produce predicable results


% -----------------------
% **** 4. Loop through each vehicle
% -----------------------
tic


% -*-*-*-*-*-*-*-*-*-*-*-
% ### MAIN MICRO SIM LIGHT CODE ###
% -*-*-*-*-*-*-*-*-*-*-*-
for i = 1:NoVeh-1
    
    % -----------------------
    % Set or Reset Queue Distance L1 and L2 back to zero at end of day and Disp WkDay No
    if i == 1
        % Set Queue Distance L1 & L2 = 0
        QuDiL1 = 0;
        QuDiL2 = 0;
        
        % Display Wk Day No to Monitor Progress on Screen
        
        % disp(Veh(i+1,46))
        
    elseif Veh(i,46)-Veh(i-1,46) == 1
        
        % RESET Queue Distance L1 & L2 = 0
        QuDiL1 = 0;
        QuDiL2 = 0;
        
        % Display Wk Day No to Monitor Progress on Screen
        disp(Veh(i+1,46))
    end
    
    % -----------------------
    % Cl = 1 for CAR and 0 for TRUCK (incl bus)
    Cl = Veh(i,49);
    
    % Define Lane
    La = Veh(i,9);
    
    % -----------------------
    % Lane 1 first - need to find a way to not do this seperately
    if La == 1
        
        % -----------------------
        % Determine Difference in Queue Lengths & Record
        ILG = QuDiL1 - QuDiL2;                  % Determine Inter Lane Gap (Difference Queue Lengths)
        ILG_Ratio = ILG / Veh(i,30);            % Determine ILG / Veh Length Ration - use this to determine lane change
        
        Veh(i,55) = ILG;
        Veh(i,56) = ILG_Ratio;
        
        
        % -----------------------
        % Define Probability of Lane Change occuring
        
        % -----------------------
        % Define PLC if Ratios < 1
        if ILG_Ratio <= 1
            if Cl == 1      % Car
                PLC = PLC_C(1);
            elseif Cl == 0  % Truck
                PLC = PLC_T(1);
            end
            
            % -----------------------
            % Define PLC if Ratios > 1 and < 4
        elseif ILG_Ratio > 1 && ILG_Ratio <= 4
            if Cl == 1      % Car
                PLC = PLC_C(2);
            elseif Cl == 0  % Truck
                PLC = PLC_T(2);
            end
            
            % -----------------------
            % Define PLC if Ratios > 4
        elseif ILG_Ratio > 4
            if Cl == 1      % Car
                PLC = PLC_C(3);
            elseif Cl == 0  % Truck
                PLC = PLC_T(3);
            end
            
        end
        
        
        % -----------------------
        % Determine if Lane Change occurs
        R = rand(1);
        
        % Lane Change
        if R > PLC
            
            % NO Lane Change occurs - record
            LaCh = 0;
            Veh(i,57) = LaCh;
            Veh(i,58) = 1;          % New Lane = Old Lane
            
            % -----------------------
            % L1 New Queue Distance Length = OLD + GAP + ith Vehicle length
            
            % Define Gap
            if GapType == 0
                Gap = ConstGap;
                
            elseif GapType == 1
                R = rand(1);
                [A, Ind] = min(abs(CDF-R));
                Gap = CDF_Gap_Con18(Ind,1);
                
            elseif GapType == 2
                R = rand(1);
                [A Ind] = min(abs(CDF-R));
                Gap = CDF_Gap(Ind,1);  
            end
            % Record new gap
            Veh(i,28) =  Gap;
            
            % Determine new Queue Length and record
            QuDiL1 = QuDiL1 + Gap + Veh(i,30);
            Veh(i,59) =  QuDiL1;
            
            
        elseif R < PLC
            
            % LANE CHANGE OCCURING
            if Dsp == 1
                disp('Lane Change 1-2 Occuring')
            end
            
            % Record Lane Change
            LaCh = 1;
            Veh(i,57) = LaCh;
            Veh(i,58) = 2;
            
            
            % -----------------------
            % L2 New Queue Distance Length = OLD + GAP + ith Vehicle length
            
            % Define Gap
            if GapType == 0
                Gap = ConstGap;
                
            elseif GapType == 1
                R = rand(1);
                [~, Ind] = min(abs(CDF-R));
                Gap = CDF_Gap_Con18(Ind,1);
            elseif GapType == 2
                R = rand(1);
                [A Ind] = min(abs(CDF-R));
                Gap = CDF_Gap(Ind,1);
                
            end
            % Record new gap
            Veh(i,28) =  Gap;
            
            % Determine new Queue Length and record
            QuDiL2 = QuDiL2 + Gap + Veh(i,30);
            Veh(i,59) =  QuDiL2;
            
            
        end
        
        
        % -----------------------
        % LANE 2 repeated as above
        % -----------------------
    elseif La == 2
        
        % -----------------------
        % Determine Difference in Queue Lengths & Record
        ILG = QuDiL2 - QuDiL1;              % Determine Inter Lane Gap (Difference Queue Lengths)
        ILG_Ratio = ILG / Veh(i,30);      	% Determine ILG / Veh Length Ration - use this to determine lane change
        
        Veh(i,55) = ILG;
        Veh(i,56) = ILG_Ratio;
        
        % -----------------------
        % Define Probability of Lane Change occuring
        
        % -----------------------
        % Define PLC if Ratios < 1
        if ILG_Ratio <= 1
            if Cl == 1      % Car
                PLC = PLC_C(1);
            elseif Cl == 0  % Truck
                PLC = PLC_T(1);
            end
            
            % -----------------------
            % Define PLC if Ratios > 1 and < 4
        elseif ILG_Ratio > 1 && ILG_Ratio <= 4
            if Cl == 1      % Car
                PLC = PLC_C(2);
            elseif Cl == 0  % Truck
                PLC = PLC_T(2);
            end
            
            % -----------------------
            % Define PLC if Ratios > 4
        elseif ILG_Ratio > 4
            if Cl == 1      % Car
                PLC = PLC_C(3);
            elseif Cl == 0  % Truck
                PLC = PLC_T(3);
            end
            
        end
        
        
        % -----------------------
        % Determine if Lane Change occurs
        R = rand(1);
        
        % Lane Change
        if R > PLC
            
            % NO Lane Change occurs - record
            LaCh = 0;
            Veh(i,57) = LaCh;
            Veh(i,58) = 2;          % New Lane = Old Lane
            
            % -----------------------
            % L2 New Queue Distance Length = OLD + GAP + ith Vehicle length
            
            % Define Gap
            if GapType == 0
                Gap = ConstGap;
                
            elseif GapType == 1
                R = rand(1);
                [~, Ind] = min(abs(CDF-R));
                Gap = CDF_Gap_Con18(Ind,1);
            elseif GapType == 2
                R = rand(1);
                [A Ind] = min(abs(CDF-R));
                Gap = CDF_Gap(Ind,1);
                
            end
            % Record new gap
            Veh(i,28) =  Gap;
            
            
            % -----------------------
            % L2 New Queue Distance Length = OLD + GAP + ith Vehicle length
            QuDiL2 = QuDiL2 + Gap + Veh(i,30);
            Veh(i,59) =  QuDiL2;
            
            
        elseif R <= PLC
            
            % LANE CHANGE OCCURING
            if Dsp == 1
                disp('Lane Change 2-1 Occuring')
            end
            
            % Record Lane Change
            LaCh = 1;
            Veh(i,57) = LaCh;
            Veh(i,58) = 1;
            
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
                [A Ind] = min(abs(CDF-R));
                Gap = CDF_Gap(Ind,1);
                
            end
            % Record new gap
            Veh(i,28) =  Gap;
            
            % -----------------------
            % L1 Determine New Queue Distance Length = OLD + GAP + ith Vehicle length
            QuDiL1 = QuDiL1 + Gap + Veh(i,30);
            Veh(i,59) =  QuDiL1;
            
        end
        
        
    end
    
    
end

% -*-*-*-*-*-*-*-*-*-*-*-
% ### MAIN MICRO SIM LIGHT CODE ###
% -*-*-*-*-*-*-*-*-*-*-*-
toc


%error('Stop here for now - create histogram below');

if Cnt == 1
    
    % -----------------------
    % **** 4. Count HGV & Non HGV per lanes before and after MS Basic
    % -----------------------
    
    % Select Lane 1 Before + count non zeros = non HGV; remaining = HGV
    Veh_L1_B4 = Veh( Veh(:,9)==1, :);           % select data from column 9 old Lane
    HGVnon_1_B4 = nnz(Veh_L1_B4(:,49));
    HGV_1_B4 = size(Veh_L1_B4,1)-HGVnon_1_B4;
    
    % Select Lane 2 Before + count non zeros = non HGV; remaining = HGV
    Veh_L2_B4 = Veh( Veh(:,9)==2, :);           % select data from column 9 old Lane
    HGVnon_2_B4 = nnz(Veh_L2_B4(:,49));
    HGV_2_B4 = size(Veh_L2_B4,1)-HGVnon_2_B4;
    
    % Select Lane 1 AFTER + count non zeros = non HGV; remaining = HGV
    Veh_L1_A = Veh( Veh(:,58)==1, :);           % select data from column 9 old Lane
    HGVnon_1_A = nnz(Veh_L1_A(:,49));
    HGV_1_A = size(Veh_L1_A,1)-HGVnon_1_A;
    
    % Select Lane 2 AFTER + count non zeros = non HGV; remaining = HGV
    Veh_L2_A = Veh( Veh(:,58)==2, :);           % select data from column 9 old Lane
    HGVnon_2_A = nnz(Veh_L2_A(:,49));
    HGV_2_A = size(Veh_L2_A,1)-HGVnon_2_A;
    
    Summary = [ HGV_1_B4        HGV_1_A         HGV_2_B4        HGV_2_A;
        HGVnon_1_B4     HGVnon_1_A      HGVnon_2_B4     HGVnon_2_A];
    
    SummaryHGV = round(100*[  HGV_1_B4/(HGV_1_B4+HGV_2_B4)   HGV_2_B4/(HGV_1_B4+HGV_2_B4);
        HGV_1_A/(HGV_1_A+HGV_2_A)      HGV_2_A/(HGV_1_A+HGV_2_A) ])
    
    SummaryNonHGV = round(100*[   HGVnon_1_B4/(HGVnon_1_B4+HGVnon_2_B4)   HGVnon_2_B4/(HGVnon_1_B4+HGVnon_2_B4);
        HGVnon_1_A/(HGVnon_1_A+HGVnon_2_A)      HGVnon_2_A/(HGVnon_1_A+HGVnon_2_A) ])
    
    
    % -----------------------
    % **** Data on Lane Change ****
    % -----------------------
    
    disp('Generating Lane Change Data 2')
    
    % -----------------------
    % Select SOUTHBOUND Lanes 1 and 2
    Veh1 = Veh( Veh(:,9)==1 , :);         % Southbound L1
    Veh2 = Veh( Veh(:,9)==2 , :);         % Southbound L2
    
    
    % -----------------------
    % Select HGV and Non HGV Vehicles
    Veh1_HGV = Veh1( Veh1(:,49)==0 , :);            % Southbound L1 -
    Veh1_nonHGV = Veh1( Veh1(:,49)==1 , :);         % Southbound L1 -
    
    Veh2_HGV = Veh2( Veh2(:,49)==0 , :);            % Southbound L2 -
    Veh2_nonHGV = Veh2( Veh2(:,49)==1 , :);         % Southbound L2
    
    
    % -----------------------
    % Determine % of lane changes in each category
    
    L1_HGV_perCh = 100*sum((Veh1_HGV(:,57)))/length(Veh1_HGV);
    L2_HGV_perCh = 100*sum((Veh2_HGV(:,57)))/length(Veh2_HGV);
    
    L1_nonHGV_perCh = 100*sum((Veh1_nonHGV(:,57)))/length(Veh1_nonHGV);
    L2_nonHGV_perCh = 100*sum((Veh2_nonHGV(:,57)))/length(Veh2_nonHGV);
    
    % summary array
    Sum = [L1_HGV_perCh L2_HGV_perCh; L1_nonHGV_perCh L2_nonHGV_perCh]
    
    
    
end

Veh_L1 = Veh(Veh(:,58)==1,:);
Veh_L2 = Veh(Veh(:,58)==2,:);

if SV == 1 
save Veh_L1 Veh_L1
save Veh_L2 Veh_L2
end
    
end
