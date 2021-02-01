% -----------------------
%
% SM_C_HZ.m
%
% Sript run an total load influence line of 500m 1000m 2000m bridge
%
% 
%
% Hang Zhang
% 08/01/2021




function [DaMax_SM_C] = SM_C_HZ(p_groupgamma,L1_SC,L2_SC,Inten_SM,V_L2_car,Veh_L1,Veh_L2,Span,dx,x,y,LI,Workdays,SC_L1_noise,SC_L2_noise,T_noise)


        rng = [];
        p_group = p_groupgamma;
        
        for Noise = 1 : T_noise
        L1_SC = [L1_SC;SC_L1_noise];
        L2_SC = [L2_SC;SC_L2_noise];
        end
        
 
        No_Veh1 = hist(Veh_L1(:,46),unique(Veh_L1(:,46)))';
        No_Veh1(all(No_Veh1==0,2),:)=[];                    % clear 0 in vector
        No_Veh2 = hist(Veh_L2(:,46),unique(Veh_L2(:,46)))';
        No_Veh2(all(No_Veh2==0,2),:)=[];                    % clear 0 in vector
        
        V_L2_car(all(V_L2_car==0,2),:)=[];
        V_L2_car(all(V_L2_car(:,30)>20,2),:)=[];
        
    
    
    for k = 1 : length(L1_SC)
        
        Veh = L1_SC{k,1};
        Veh_S(:,(1:14)) = Veh(:,(14:27));
        Veh_S(:,(15:27)) = Veh(:,(33:45));
        Veh_S(:,28) = Veh(:,30);
        Veh_S(:,29) = Veh(:,56);
        Veh_S(:,30) = Veh(:,30);
        Veh_S(:,31) = Veh(:,13);
        L1_SC_S{k,1} = Veh_S;
        Veh_S =[];
    end
    
    
    
    
    for k = 1 : length(L2_SC)
        
        Veh = L2_SC{k,1};
        Veh_S(:,(1:14)) = Veh(:,(14:27));
        Veh_S(:,(15:27)) = Veh(:,(33:45));
        Veh_S(:,28) = Veh(:,30);
        Veh_S(:,29) = Veh(:,56);
        Veh_S(:,30) = Veh(:,30);
        Veh_S(:,31) = Veh(:,13);
        L2_SC_S{k,1} = Veh_S;
        Veh_S =[];
    end
    
    
    car_S(:,(1:14)) = V_L2_car(:,(14:27));
    car_S(:,(15:27)) = V_L2_car(:,(33:45));
    car_S(:,28) = V_L2_car(:,30);
    car_S(:,29) = V_L2_car(:,56);
    car_S(:,30) = V_L2_car(:,30);
    car_S(:,31) = V_L2_car(:,13);
    
    L1_SC = [];
    L2_SC = [];
    
    sB =  cell2mat(p_group(:,3));
    
 
 
    
    tic
    for   i = 1 : Workdays
         disp(i)
        
        %         Day_Veh = No_Veh(randperm(length(No_Veh),1));
        %         Day_Veh1 = No_Veh1(randperm(length(No_Veh1),1));
        % Day_Veh(i,1) = No_Veh(i,1);
        
        R_Day_Veh = randi(length(No_Veh1));
        
        Day_Veh1 = No_Veh1(R_Day_Veh);
        Day_Veh2 = No_Veh2(R_Day_Veh);
        
        
        Day_VehRH = [];
        Day_VehRH1 = [];
        
        %-----------------------------------------------------------------
        for h = 1: Day_Veh1+Day_Veh2
            
            
            if h ==  1
                
                r = randi(length(L1_SC_S));
                Senario = L1_SC_S{r,1};
                Senario1 = L2_SC_S{r,1};
                
                Int_Ind = (sum(Senario(:,31))+sum(Senario1(:,31)))/ (sum(Senario(:,29))+sum(Senario(:,30)));
                
                Day_VehRH = [Day_VehRH;Senario];
                Day_VehRH1 = [Day_VehRH1;Senario1];
                
                
                rng = [rng;Int_Ind];
                
            elseif  h ~= 1
                
                
                
                sA =  abs(sB-Int_Ind);
                
                index = find(sA==min(sA));
                
                
                
                RESULT1 = p_group{index(1,1),1};
                RESULT2 = p_group{index(1,1),2};
                
                V = 0:0.001:1; % define the numbers
                R = V(randi([1,numel(V)]));
                
                
                [~, Ind] = min(abs(RESULT2-R));
                Int =RESULT1(Ind);
                
                
                [~,index_1] = min(abs(Inten_SM-Int));
                % RESULT_1 = Inten_SM(index_1(1));
                MC_Int = index_1(1);
                
                r1 = MC_Int;
                r1 = r1(1,1);
                
                Senario = L1_SC_S{r1,1};
                Senario1 = L2_SC_S{r1,1};
                
                Int_Ind =  Int;
                
                rng = [rng;Int_Ind];
                
                if length(Day_VehRH)< Day_Veh1
                    
                    Day_VehRH = [Day_VehRH;Senario];
                    
                elseif length(Day_VehRH)>= Day_Veh1
                    
                    Day_VehRH = Day_VehRH (1:Day_Veh1,:);
                end
                
                %--------------------------------------------
                
                if length(Day_VehRH1)< Day_Veh2
                    
                    Day_VehRH1 = [Day_VehRH1;Senario1];
                    
                elseif length(Day_VehRH1)>= Day_Veh2
                    
                    Day_VehRH1 = Day_VehRH1 (1:Day_Veh2,:);
                end
                %--------------------------------------------
                if length(Day_VehRH)== Day_Veh1 && length(Day_VehRH1)== Day_Veh2
                    break
                end
                
                
            end
        end
        
        
        %-------------------------------------------------------------------------------------
        % insert car back
        %-------------------------------------------------------------------------------------
        for hh = 1 : length(Day_VehRH1)-1
            
            if Day_VehRH1(hh,1) == 0 && Day_VehRH1(hh+1,1) == 0 && (Day_VehRH1(hh,29)+ Day_VehRH1(hh+1,29))>20
                
                Car = car_S(randperm(length(car_S),1),:);
                Car(1,29)=0;
                ac =  Day_VehRH1(hh,29)+ Day_VehRH1(hh+1,29)-Car(1,28);
                
                Day_VehRH1(hh,:)=Car;
                Day_VehRH1(hh+1,29)=ac;
            end
        end
        
        %-------------------------------------------------------------------
        % L1
        %-------------------------------------------------------------------
        
        
        % Initialise Day_BrR array - ing a matrix - it help to speed up matlab if you can define the matrix at the start
        % as a zero matrix and then fill.
        % For one day bridge reponse
        Matrix = round(sum(Day_VehRH(:,28))+sum(Day_VehRH(:,29)))+Span+500;
        
        Day_BrR = zeros(Matrix,1) ;
        Day_BrR_S1 = zeros(Matrix,1) ;
        
        
        % -----------------------
        % Loop through one day vehicles
        for j=1:length(Day_VehRH(:,1))
            
            
            
            % vehicles in one day
            
            % Define Vehicle Properties
            
            F = [Day_VehRH(j,1:14)]; % vehicle axle load
            F(F==0)=[]; % clear 0 in vector
            A = [0 Day_VehRH(j,15:27)]; % vehicle axle space
            NoAx = length(F);
            
            
            
            % -----------------------
            % 3. Bridge Response
            % -----------------------
            % estimate length of Br = ?
            BrLe = round(1.2 * ( (sum(A)+Span) / dx) ); % Bridge Length - BrLe determine length of zero vector
            BrReBM = zeros(BrLe,1);                  % set up suitble scale of zero matrix
            BrReBM_S = zeros(BrLe,1 );
            
            for k=1:NoAx                              % NoAx times vector calculation
                
                
                % (sum(A(1:i)),1) is a column vector of Accumulating axle space
                % Axle space affect how much zero before vector b
                a = (round(sum(A(1:k))/dx));
                
                
                % axle load mutiply influence line function
                b = F(k)*y;
                
                
                BrReBM_S((a+1):length(b)+a,1) = BrReBM_S((a+1):length(b)+a,1)+b;
                
                
            end
            
            
            
            
            % One day bridge reponse
            
            VeLe = [0; Day_VehRH(:,28)];
            % Gap = Veh_L2_GAP(:,28);
            Gap = [0;Day_VehRH(:,29)];
            % (sum(h(1:j))+1.5) is a column vector of Accumulating vehicle space, 1.5m
            % is vehicle distance
            % Vehicle gap affect how much zero before vector b
            
            g =(round(sum(VeLe(1:j))+sum(Gap(1:j)))/dx);
            
            Day_BrR_S1((g+1):length(BrReBM_S)+g,1) = Day_BrR_S1((g+1):length(BrReBM_S)+g,1)+BrReBM_S;
            
            
            
        end
        
        %-------------------------------------------------------------------
        % L2
        %-------------------------------------------------------------------
        
        % Initialise Day_BrR array - ing a matrix - it help to speed up matlab if you can define the matrix at the start
        % as a zero matrix and then fill.
        % For one day bridge reponse
        Matrix = (round(sum(Day_VehRH1(:,28))+sum(Day_VehRH1(:,29)))+Span+500)/dx;
        
        Day_BrR = zeros(Matrix,1 ) ;
        Day_BrR_S2 = zeros(Matrix,1 ) ;
        
        
        % -----------------------
        % Loop through one day vehicles
        for jj=1:length(Day_VehRH1(:,1))
            
            
            
            % vehicles in one day
            
            % Define Vehicle Properties
            
            F = [Day_VehRH1(jj,1:14)]; % vehicle axle load
            F(F==0)=[]; % clear 0 in vector
            A = [0 Day_VehRH1(jj,15:27)]; % vehicle axle space
            NoAx = length(F);
            
            
            
            % -----------------------
            % 3. Bridge Response
            % -----------------------
            % estimate length of Br = ?
            BrLe = round(1.2 * ( (sum(A)+Span) / dx) ); % Bridge Length - BrLe determine length of zero vector
            BrReBM = zeros(BrLe,1 );                  % set up suitble scale of zero matrix
            BrReBM_S = zeros(BrLe,1);
            
            for kk=1:NoAx                              % NoAx times vector calculation
                
                
                % (sum(A(1:i)),1) is a column vector of Accumulating axle space
                % Axle space affect how much zero before vector b
                a = (round(sum(A(1:kk))/dx) );
                
                
                % axle load mutiply influence line function
                b = F(kk)*y;
                
                
                BrReBM_S((a+1):length(b)+a,1) = BrReBM_S((a+1):length(b)+a,1)+b;
            end
            
            
            
            VeLe = [0; Day_VehRH1(:,28)];
            % Gap = Veh_L2_GAP(:,28);
            Gap = [0;Day_VehRH1(:,29)];
            
            
            
            g = (round(sum(VeLe(1:jj))+sum(Gap(1:jj)))/dx);
            
            
            
            Day_BrR_S2((g+1):length(BrReBM_S)+g,1) = Day_BrR_S2((g+1):length(BrReBM_S)+g,1)+BrReBM_S;
            
            
            
            
        end
        
        %-----------------------------------------------------------
        % L1+L2
        %-----------------------------------------------------------
        
        if length(Day_BrR_S1) >= length(Day_BrR_S2)
            bb = zeros ((length(Day_BrR_S1))-(length(Day_BrR_S2)),1 );
            cc = [Day_BrR_S2 ;bb];
            aa = Day_BrR_S1 + cc;
            
        elseif length(Day_BrR_S1) < length(Day_BrR_S2)
            bb = zeros ((length(Day_BrR_S2))-(length(Day_BrR_S1)),1);
            cc = [Day_BrR_S1;bb];
            aa = Day_BrR_S2+ cc;
        end
        
        
        % -----------------------
        % 3. Max Daily Bridge Response
        % -----------------------
        
        cc=[];
        Day_BrR_S1=[];
        Day_BrR_S2=[];
        
        if LI == 1 
        Daily_M(i,1) = max(aa)./Span;   % load intensity 
        elseif LI == 0
        Daily_M(i,1) = max(aa);         % load 
        end
        
    end
    toc
    
    
   DaMax_SM_C = Daily_M;

    


end

