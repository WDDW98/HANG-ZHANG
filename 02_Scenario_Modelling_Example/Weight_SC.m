% -----------------------
%
% Weight_SN_ME.m
%
% script works out intensity of scenario 
% 
% 
%
%
% Hang Zhang
% 18/01/2021
% -----------------------



function [Inten_SM] = Weight_SC (L1_SC,L2_SC)



SM150_Veh = [];

for i = 1 : length(L1_SC)
    
    SM150_Veh = L1_SC{i,1};
    
    Inten_SM150_L1(i,1) = sum(SM150_Veh(:,13))/(sum(SM150_Veh(:,30))+sum(SM150_Veh(:,56)));
    

end



for i = 1 : length(L2_SC)
    
    SM150_Veh = L2_SC{i,1};
    
    Inten_SM150_L2(i,1) = sum(SM150_Veh(:,13))/(sum(SM150_Veh(:,30))+sum(SM150_Veh(:,56)));
    

end


Inten_SM = Inten_SM150_L1 + Inten_SM150_L2;


end

