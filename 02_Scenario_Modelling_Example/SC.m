
% -----------------------
%
% Traffic Senario .m
%
% sript extract traffic senario
% L1 & L2
%
% Hang Zhang
% 13/01/2018




function [L1_SC, L2_SC, V_L2_car] = SC(SN,SV1)

% ------------------------
% 1. Prepare Vairables
% ------------------------

disp('Loading data...')
load ('Veh_L1.mat')
load ('Veh_L2.mat')
disp('Data Loaded...')

Veh_Hr1 = Veh_L1;        
Veh_Hr2 = Veh_L2;
    

Veh_Hr1 (:,56)= Veh_Hr1 (:,28);
Veh_Hr2 (:,56)= Veh_Hr2 (:,28);

Veh_Hr1 (:,55)= cumsum((Veh_Hr1(:,30)+Veh_Hr1(:,56)));
Veh_Hr2 (:,55)= cumsum((Veh_Hr2(:,30)+Veh_Hr2(:,56)));



c=0;
d=0;


N0_car = 0;
N0_gap = 0;

No_SC_L1 = round((sum(Veh_Hr1(:,28))+ sum(Veh_Hr1 (:,30)))/SN);
No_SC_L2 = round((sum(Veh_Hr2(:,28))+ sum(Veh_Hr2 (:,30)))/SN);

No_SC = max(No_SC_L1,No_SC_L2)+1000;

car = zeros (No_SC,(length(Veh_Hr2(1,:)))) ;
gap = zeros (1,length(Veh_Hr2(1,:)));


    
    for i = 1: No_SC
        disp(i)
        b=0;
        k=0;
        
        
        for j = 1:length(Veh_Hr1(:,1))
            
            a = sum(  Veh_Hr1((1:j),30) ) + sum(  Veh_Hr1((1:j),56)   );
            
            if round(a)>= SN-10
                
                for k = 1: length(Veh_Hr2(:,1))
                    
                    b =  sum(  Veh_Hr2((1:k),30)   ) + sum(  Veh_Hr2((1:k),56) );
                    
                    if  a<=b && Veh_Hr2(k,16) ~= 0 && b-a >= Veh_Hr2(k,56)
                        
                        break
                    end
                    
                    if a<=b && b-a < Veh_Hr2(k,56)
                        
                        Veh_Hr2(k,56) = Veh_Hr2(k,56) -(b-a);
                        % Veh_Hr2(k,28)=Veh_Hr2(k,56);                       % replace last car to gap
                        L1_SC{i,1} = Veh_Hr1((1:j),:);
                        L2_SC{i,1} = Veh_Hr2((1:k),:);
                        
                        break
                        
                    end
                    
                    
                    
                    if  a<=b && Veh_Hr2(k,16) == 0 && b-a >= Veh_Hr2(k,56)
                        
                        g = (Veh_Hr2(k,30)) + (Veh_Hr2(k,56)); % last car in lane 2
                        f = a - (b-g);                         % distance in last scenario
                        car(i,:) = Veh_Hr2(k,:);               % record last car
                        

                        
                        Veh_Hr2(k,:) = gap;
                        Veh_Hr2(k,56)=f;                       % replace last car to gap
                        Veh_Hr2(k,28)=f;                       % replace last car to gap
                        
                        
                        L1_SC{i,1} = Veh_Hr1((1:j),:);
                        L2_SC{i,1} = Veh_Hr2((1:k),:);
                        
                        
                        
                        break
                        
                        
                        
                    end
                    
                    
                end
            end
            
            if round(a)>= SN-10
                
                
                if  a<b  && Veh_Hr2(k,16) == 0  && b-a >= Veh_Hr2(k,56)
                    d = b-a;
                    car1 = gap;
                    car1(1,56)=d;                       % replace first car to gap
                    
                    
                    
                    Veh_Hr1= Veh_Hr1( (j+1:length( Veh_Hr1(:,1))),:  );
                    Veh_Hr2= Veh_Hr2( (k+1:length( Veh_Hr2(:,1))),:  );
                    Veh_Hr2 = [car1;Veh_Hr2];
                    
                    break
                end
                    
                  %  ------------------------------------------------------    
                    
                   if a<=b && b-a < Veh_Hr2(k,56)
                        
            
                        
                    Veh_Hr1= Veh_Hr1( (j+1:length( Veh_Hr1(:,1))),:  );
                    Veh_Hr2= Veh_Hr2( (k+1:length( Veh_Hr2(:,1))),:  );
                    
                    
                    d = b-a;
                    car1 = gap;
                    car1(1,56)=d;                       % replace first car to gap
                    car1(1,28)=d;                       % replace first car to gap
                    Veh_Hr2 = [car1;Veh_Hr2];   
                    
                        break
                        
                   end
                %  ------------------------------------------------------     
                    
                    
                
            end
            
            
            
            
            
            
        end
        
        
        if k+1 > length(Veh_Hr2(:,1))
            break
        elseif j+1 > length(Veh_Hr1(:,1))
            break
        end
    end
    
    V_L2_car = car;
    
% ------------------------
% 2. Save Data
% ------------------------


   if SV1 == 1 
    save L1_SC L1_SC
    save L2_SC L2_SC
    save V_L2_car V_L2_car;
   end 

end


    

    
    

