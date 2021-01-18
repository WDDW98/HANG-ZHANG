% -----------------------
%
%
%
% Correlation.m
%
% script classfies scenarios to bins 
%
%
% Hang Zhang
% 18/01/2021



function [cell_b]= Corr(Inten_SM,Interval,up_control)

a = [0:Interval:round(max(Inten_SM)-Interval+up_control)];

for i = 1 : length(a)
    
    
    Vector_b=[];
    
    
    for j = 1:length(Inten_SM)-1
        
        Inten_f = Inten_SM(j,1);
        Inten_b = Inten_SM(j+1,1);
        
        
        %---------------------------------------------------------------------
        if a(i) ~= round(max(Inten_SM)-Interval+up_control)
            if a(i)< Inten_f && Inten_f <a(i)+Interval
                
                Vector_b(j,1)= Inten_b;
                Vector_b(j,2)= j;
            end
        elseif a(i) == round(max(Inten_SM)-(Interval+up_control))
            
            if a(i)< Inten_f 
                
                Vector_b(j,1)= Inten_b;
                Vector_b(j,2)= j;
                
            end
            
        end
    end
    
    Vector_b(all(Vector_b==0,2),:)=[];
    
    if isempty(Vector_b)
        cell_b{i,1} = Vector_b;
        cell_b{i,2} = 0;
        continue
    end
    
    cell_b{i,1} = Vector_b;
    cell_b{i,2} = a(i);
    
    
    if a(i) == round(max(Inten_SM)-(Interval+up_control))
        
        break
    end
    
    
    
end

end



% save cell_b cell_b
