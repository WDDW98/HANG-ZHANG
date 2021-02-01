% -----------------------
%
%
%
% Noise.m
%
% script generates new scenarios 
%
%
%
%
% Hang Zhang
% 05/07/2020
% -----------------------
function[SC_L1_noise, SC_L2_noise] = noise_HZ2 (L1_SC,L2_SC,Group_GVW,No_Noise,GapType,Max_No_HGV_change)


load('CDF_Gap_Con.mat')


%-------------------------------
% Determine No SC
%------------------------------- 

% count number of HGV in scenario 

for i = 1 : length(L1_SC)
    
    SC1 = L1_SC{i};
    SC2 = L2_SC{i};
    
    if SC2(1,1)==0
        SC2(1,49)=1;
    elseif SC2(length(SC2(:,1)),1)==0
        SC2(length(SC2(:,1)),49)=1;
    end
    
    
     No_HGV_1 (i,1)= sum(SC1(:,49)==0);
     No_HGV_2 (i,1)= sum(SC2(:,49)==0);
     
end
No_HGV = No_HGV_1 + No_HGV_2;


L1_Veh= cell2mat(L1_SC);
L2_Veh= cell2mat(L2_SC);
Veh = [L1_Veh;L2_Veh];
Veh(Veh(:,2)==0,:)=[];
HGV = Veh(Veh(:,49)==0,:);
CAR = Veh(Veh(:,49)==1,:);
CAR(CAR(:,2)==0,:)=[];
CAR_L1 = CAR(CAR(:,58)==1,:);
CAR_L2 = CAR(CAR(:,58)==2,:);
  
  
% Determine bandwidth
[y,x] = ksdensity (No_HGV, 'npoint',100,'Bandwidth',1);
y1 = 1./sqrt(y);

yfit = y1(y1<=Max_No_HGV_change)';
loc = find(y1<=Max_No_HGV_change);
xfit = x(:,loc)';
f = fit(xfit,yfit,'poly2');

r1 = randperm(length(L1_SC),No_Noise);
% load ('r1.mat')
    
tic 

for i = 1 : No_Noise
    % disp(i)
    SAS = 0;
  
    SC1 = L1_SC{r1(i)};
    SC2 = L2_SC{r1(i)};
    SC2(SC2(:,2)==0,:)=[];
     
    if SC2(1,1)==0
        SC2(1,49)=1;
    elseif SC2(length(SC2(:,1)),1)==0
        SC2(length(SC2(:,1)),49)=1;
    end
    
    x_b = sum(SC1(:,49)==0) + sum(SC2(:,49)==0);
     record1(i,1)= x_b;
     [Add_HGV]= add_HGV(x_b,xfit,f,Max_No_HGV_change);
     Add_Nu_HGV = Add_HGV;
      record(i,1)= Add_Nu_HGV;
      % record3(i,1)= r;
      
    if Add_Nu_HGV >0
        
        HGV_add = HGV(randperm(length(HGV),Add_Nu_HGV),:);
        L1_HGV_add = HGV_add(HGV_add(:,58)==1,:);
        L2_HGV_add = HGV_add(HGV_add(:,58)==2,:);
        
        SC1_car = SC1(SC1(:,49)==1,:);
        SC2_car = SC2(SC2(:,49)==1,:);
        
        SC_car = [SC1_car;SC2_car];
        
        SC1_truck = SC1(SC1(:,49)==0,:);
        SC2_truck = SC2(SC2(:,49)==0,:);
        
        SC_truck = [SC1_truck;SC2_truck];  
        
        if sum(SC1_car(:,30)+SC1_car(:,28))>= sum(L1_HGV_add(:,30)+L1_HGV_add(:,28))
           team_car_L1 = cumsum(SC1_car(:,30)+SC1_car(:,28));
           [~, Ind] = min(abs(team_car_L1-sum(L1_HGV_add(:,30)+L1_HGV_add(:,28))));
           New_cars_L1 = SC1_car(Ind:length(team_car_L1),:);
           SAS = 1;
        end
        if sum(SC2_car(:,30)+SC2_car(:,28))>= sum(L2_HGV_add(:,30)+L2_HGV_add(:,28))
           team_car_L2 = cumsum(SC2_car(:,30)+SC2_car(:,28));
           [~, Ind] = min(abs(team_car_L2-sum(L2_HGV_add(:,30)+L2_HGV_add(:,28))));
           New_cars_L2 = SC2_car(Ind:length(team_car_L2),:);
           SAS = 1;
        end
        if sum(SC1_car(:,30)+SC1_car(:,28))< sum(L1_HGV_add(:,30)+L1_HGV_add(:,28))  
           New_cars_L1=[];
           SAS = 1; 
        end
        if sum(SC2_car(:,30)+SC2_car(:,28))< sum(L2_HGV_add(:,30)+L2_HGV_add(:,28))  
           New_cars_L2=[];
           SAS = 1; 
        end   
%         if sum(SC_car(:,30))>= sum(HGV_add(:,30))
%             
%             team_car = cumsum(SC_car(:,30));
%             [~, Ind] = min(abs(team_car-sum(HGV_add(:,30))));
%             New_cars = SC_car(1:Ind,:);
%             SAS = 1; 
%         elseif sum(SC_car(:,30))< sum(HGV_add(:,30))
%         % SC1_added = [SC1;HGV_add];
%         New_cars=[];
%         SAS = 1; 
%         end
  
    elseif Add_Nu_HGV ==0
        SC1_added = SC1;
        
    elseif Add_Nu_HGV <0
        
        No_HGV_1 = sum(SC1(:,49)==0);
        No_HGV_2 = sum(SC2(:,49)==0);
        
        No_HGV_12 = No_HGV_1+ No_HGV_2;
        HGV_exist_nu = No_HGV_12+Add_Nu_HGV;
        
        SC_T = [SC1;SC2];
        SC_T(SC_T(:,2)==0,:)=[];
        HGV_T = SC_T(SC_T(:,49)==0,:);
        
        if HGV_exist_nu<= 0
            
            Truck_minus_L1= SC1(SC1(:,49)==0,:); 
            SC1(SC1(:,49)==0,:)=[];
            Truck_minus_L2= SC2(SC2(:,49)==0,:); 
            SC2(SC2(:,49)==0,:)=[];
            
            
            CAR_L1_add = CAR_L1(randperm(length(CAR_L1),100),:); 
            CAR_L2_add = CAR_L2(randperm(length(CAR_L2),100),:); 
            
            if  ~isempty(Truck_minus_L1)
            add_car_L1 = cumsum(CAR_L1_add(:,30)+CAR_L1_add(:,28));
            [~, Ind] = min(abs(add_car_L1-sum(Truck_minus_L1(:,30)+Truck_minus_L1(:,28))));
            fill_cars_L1 = CAR_L1_add(1:Ind,:);
            elseif isempty(Truck_minus_L1)
            fill_cars_L1 = [];
            end
            
            if  ~isempty(Truck_minus_L2)   
            add_car_L2 = cumsum(CAR_L2_add(:,30)+CAR_L2_add(:,28));
            [~, Ind] = min(abs(add_car_L2-sum(Truck_minus_L2(:,28)+Truck_minus_L2(:,30))));
            fill_cars_L2 = CAR_L2_add(1:Ind,:);           
            elseif isempty(Truck_minus_L2)
            fill_cars_L2 = [];
            end    
            
            SAS = 2;   

           % fill_cars_L1 fill_cars_L2 SC1 SC2
           % fill_cars_L2=[] fill_cars_L1=[] New_cars_L1=[] New_cars_L1=[]
           
        elseif HGV_exist_nu> 0
            
            SC1(SC1(:,49)==0,:)=[];
            SC2(SC2(:,49)==0,:)=[];
            
            r=randperm( size(HGV_T,1) ); 
            HGV_tem=HGV_T(r, :);  
            
            HGV_add = HGV_tem(1:HGV_exist_nu,:);
            HGV_rest= HGV_tem(HGV_exist_nu+1: length(HGV_tem(:,1)),:);
            
            CAR_L1_add = CAR_L1(randperm(length(CAR_L1),100),:); 
            CAR_L2_add = CAR_L2(randperm(length(CAR_L2),100),:); 
            
            Truck_minus_L1 = HGV_rest(HGV_rest(:,58)==1,:); 
            Truck_minus_L2 = HGV_rest(HGV_rest(:,58)==2,:);
            
            if  ~isempty(Truck_minus_L1)
            add_car_L1 = cumsum(CAR_L1_add(:,30)+CAR_L1_add(:,28));
            [~, Ind] = min(abs(add_car_L1-sum(Truck_minus_L1(:,30)+Truck_minus_L1(:,28))));
            fill_cars_L1 = CAR_L1_add(1:Ind,:);
            elseif isempty(Truck_minus_L1)
            fill_cars_L1 = [];
            end
            
            if  ~isempty(Truck_minus_L2)   
            add_car_L2 = cumsum(CAR_L2_add(:,30)+CAR_L2_add(:,28));
            [~, Ind] = min(abs(add_car_L2-sum(Truck_minus_L2(:,30)+Truck_minus_L2(:,28))));
            fill_cars_L2 = CAR_L2_add(1:Ind,:);    
            elseif isempty(Truck_minus_L2)
            fill_cars_L2 = [];
            end            
          
            % SC1_added = [SC1;HGV_add];
            % fill_cars_L1 fill_cars_L2 HGV_add SC2  SC1
            SAS = 3;  
        end
        
    end
    
    % now we have added SC1 and SC2, next step is giving a GVW 
    if SAS == 1
    SC = [HGV_add;New_cars_L1;New_cars_L2;SC_truck]; 
    SC(SC(:,2)==0,:)=[];
    elseif  SAS == 0
    SC = [SC1_added;SC2];
    SC(SC(:,2)==0,:)=[];
    elseif  SAS == 2
   % SC = [ SC1; SC2];  
    SC = [fill_cars_L1; fill_cars_L2; SC1; SC2];
    SC(SC(:,2)==0,:)=[];
    
    elseif  SAS == 3
    SC = [fill_cars_L1; fill_cars_L2; HGV_add; SC2; SC1];
    SC(SC(:,2)==0,:)=[];
    
    end
    
         %--------------------------------------------
    SC_1_noise= SC(SC(:,58)==1,:);
    SC_2_noise= SC(SC(:,58)==2,:);
    
 
%     SM150_Veh_noise =  SC_L1_noise;
%     Inten_SM150_L1_noise = sum(SM150_Veh_noise(:,13))/(sum(SM150_Veh_noise(:,30))+sum(SM150_Veh_noise(:,56)));
%     
%     SM150_Veh_noise =  SC_L2_noise; 
%     Inten_SM150_L2_noise = sum(SM150_Veh_noise(:,13))/(sum(SM150_Veh_noise(:,30))+sum(SM150_Veh_noise(:,56)));
%     
% 
% Inten_SM150_noise_P1 = Inten_SM150_L1_noise + Inten_SM150_L2_noise
% %      if Inten_SM150_noise_P1<4
% %      break
% %      end
 %---------------------------------------------------  
    
    
    
    SC_L1 =  SC(SC(:,58)==1,:);
    SC_L2 =  SC(SC(:,58)==2,:);
    SC_L1(:,56)=  SC_L1(:,28);
    SC_L2(:,56)=  SC_L2(:,28);

    length_L1_SC =  sum(SC_L1(:,56)+ SC_L1(:,30));
    length_L2_SC =  sum(SC_L2(:,56)+ SC_L2(:,30));
    End_differ1 = length_L1_SC-length_L2_SC;
   %  (i,1)
%      if End_differ1>80
%          break
%      end
    
    
    for k = 1: length(SC(:,1))
        if round(SC(k,30)) < 32
        GVW_xy = Group_GVW(round(SC(k,30)),:);
        elseif round(SC(k,30)) >= 32
        GVW_xy = Group_GVW(32,:);  
        end
        GVW_x = GVW_xy{1,2};
        GVW_y = GVW_xy{1,3};
        r = rand(1);
        [~, Ind] = min(abs(GVW_y-r));
        New_GVW =GVW_x(Ind(randperm(length(Ind),1),1));
        axles_weight_odd = SC(k,14 :27);
        proportion = axles_weight_odd./sum(axles_weight_odd);
        axles_weight_new = New_GVW.*proportion;
        SC(k,13)= New_GVW;
        SC(k,14 :27)= axles_weight_new;
        
    end
    % now we have given new weight to Veh, next step is arrange vehicle properly
    r=randperm( size(SC,1) ); 
    SC=SC(r, :);                             
  
    
    [SC] =  HZ100_MS_Lite_1(GapType,SC,CDF_Gap_Con18);
    
    
        %--------------------------------------------
%     SC_L1_noise= SC(SC(:,58)==1,:);
%     SC_L2_noise= SC(SC(:,58)==2,:);
%     
%  
%     SM150_Veh_noise =  SC_L1_noise;
%     Inten_SM150_L1_noise = sum(SM150_Veh_noise(:,13))/(sum(SM150_Veh_noise(:,30))+sum(SM150_Veh_noise(:,56)));
%     
%     SM150_Veh_noise =  SC_L2_noise; 
%     Inten_SM150_L2_noise = sum(SM150_Veh_noise(:,13))/(sum(SM150_Veh_noise(:,30))+sum(SM150_Veh_noise(:,56)));
%     
% 
% Inten_SM150_noise_P1 = Inten_SM150_L1_noise + Inten_SM150_L2_noise
%      if Inten_SM150_noise_P1<4
%      break
%      end
 %---------------------------------------------------  
    
    % keep L1 SC and L2 SC same length 
    SC_L1 =  SC(SC(:,58)==1,:);
    SC_L2 =  SC(SC(:,58)==2,:);
    SC_L1(:,56)=  SC_L1(:,28);
    SC_L2(:,56)=  SC_L2(:,28);

     length_L1_SC =  sum(SC_L1(:,56)+ SC_L1(:,30));
     length_L2_SC =  sum(SC_L2(:,56)+ SC_L2(:,30));
     End_differ(i,1) = length_L1_SC-length_L2_SC;
     
    if length_L1_SC - length_L2_SC > 0 
        SC_L2(length(SC_L2(:,1)),56) = SC_L2(length(SC_L2(:,1)),56) +(length_L1_SC - length_L2_SC);
    elseif length_L1_SC - length_L2_SC < 0 
        SC_L1(length(SC_L1(:,1)),56) = SC_L1(length(SC_L1(:,1)),56) +(length_L2_SC - length_L1_SC);
    end

    SC_L1(:,28)=  SC_L1(:,56);
    SC_L2(:,28)=  SC_L2(:,56);
    % save new SC
    
    SC_L1_noise{i,1} = SC_L1;
    SC_L2_noise{i,1} = SC_L2;
    
    SC = [];
    SC_L1=[];
    SC_L2=[];
    Truck_minus_L1=[];
    Truck_minus_L2=[];
    HGV_add=[];
    New_cars_L1=[];
    New_cars_L2=[];
    SC_truck=[];
    
end
toc



end
% save r1 r1