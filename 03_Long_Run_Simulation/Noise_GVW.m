% -----------------------
%
%
%
% Noise_GVW.m
%
% script...
%
%
%
%
% Hang Zhang
% 30/06/2020
% -----------------------

function [Group_GVW] = Noise_GVW(Veh_L1,Veh_L2,ub_GVW_bin,GVW_bin)


Veh_Normal = [Veh_L1;Veh_L2];
aa = sort(Veh_Normal(:,30));

lb = round(aa(1)-0.5);
ub = min(round(aa(length(aa)-1)+0.5),ub_GVW_bin);

a = [lb : GVW_bin : ub]';


for i = 1 : length(a)
    
    ii = a(i,1);
    
    if ii <= ub-1
    Veh = Veh_Normal( Veh_Normal(:,30)>=ii & Veh_Normal(:,30)<ii+1 ,:);
    Veh1{i,1} = Veh;
    
    elseif ii == ub
        
    Veh = Veh_Normal( Veh_Normal(:,30)>=ii ,:);
    Veh1{i,1} = Veh; 
    end
end   


for i = 1 : length(a)
   
   Veh2 = Veh1{i,1};
  [y_GVW,x_GVW] = ksdensity (Veh2(:,13), 'npoint',1000,'Function','cdf'); % 'Bandwidth',1
  
  Group_GVW(i,1) = {i};
  Group_GVW(i,2) = {x_GVW};
  Group_GVW(i,3) = {y_GVW};
end
end
