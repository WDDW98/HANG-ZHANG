% -----------------------
%
%
%
% Cor_Group.m
%
% script work out correlation between scenarios
%
%
%
%
% Hang Zhang
% 04/11/2019
% 18/18/2021
% -----------------------

function [p_groupgamma] = Cor_Group_CERI(cell_b,u,Inten_SM)

b_SC = (cell_b(:,1));
b_SC(cellfun(@isempty,b_SC))=[];
  
    figure
    for i = 1 : length(b_SC)
        
        a = b_SC{i,1};
  
        x2 = 0:0.001:max(Inten_SM+3);

        [phat, pci] = gamfit(a(:,1));
        y5 = gamcdf(x2,phat(1)+u,phat(2));
        plot( x2,y5,'--','LineWidth',1.5)
        hold on

        
        p_groupgamma{i,1} = x2;
        p_groupgamma{i,2} = y5;
        
    end
    
    grid on
    
    % legend('Sc. Bin 2-4 kN/m','Sc. Bin 4-6 kN/m','Sc. Bin 6-8 kN/m','Sc. Bin 8-10 kN/m','Sc. Bin 10-12 kN/m','Sc. Bin 12-14 kN/m','Sc. Bin 14-18 kN/m')
    
    xlabel ('Intensity (kN/m)')
    ylabel ('Cumulative Probability')
    xlim([0 18])
   %  title( { ['u=default + ', num2str(u(k)) ] } )
    
    b=cell2mat(cell_b(:,2));
    b(all(b==0,2),:)=[];
    b = num2cell(b);
    p_groupgamma(:,3) = b;
    
end
    
    
    
    
    
    
    
    
    
%     name_string = ['p_groupgamma' num2str(k) '=p_groupgamma'];
%     eval(name_string);
    
    



% 
% for i = 7 %:length(b_SC)
%     figure
%     a = b_SC{i,1};
%     
%     [y1, x1,bw] = ksdensity(a(:,1));
%     
%     
%     hold on
%     h = histogram(a(:,1),'Normalization','probability','BinWidth',bw*0.9)
%     hold on
%     x2 = 0:0.1:30;
%     
%  
%     
%     [phat1, pci1] = gamfit(a(:,1));
%     y7 = gampdf(x2,phat1(1)+0.25,phat1(2));
%     plot( x2,y7,'--','LineWidth',1.5)
%     
%     grid on
%     legend('Measured','Gamma')
%     xlabel ('Following Scenario Load Intensity (kN/m)')
%     ylabel ('Probability')
%     xlim([0 20])
%     % title( 'Front Intensity > 14kN/m' )
%     box on 
% end

