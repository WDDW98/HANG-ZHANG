 

function [Add_HGV]=add_HGV(x_b,xfit,f,Max_No_HGV_change)
 
 if x_b>=0 && x_b<min(xfit)
    y_b= Max_No_HGV_change;
 elseif x_b>=min(xfit) && x_b < max(xfit)
    y_b = f(x_b);
 elseif x_b >= max(xfit)
    y_b = Max_No_HGV_change;
 end
 
 pd = fitdist(0.2,'kernel','Kernel','triangle','Bandwidth',y_b/2);
 x = -y_b-2:.01:y_b+2;
 y = cdf(pd,x);
 r = rand(1);
 [~, Ind] = min(abs(y-r));
 Add_HGV =round(x(Ind));
 % plot(x,y)
end
 
%  

%  pd = fitdist(0,'kernel','Kernel','triangle','Bandwidth',1);
%  x = -y_b-5:.1:y_b+5;
%  y = pdf(pd,x);
%  r = rand(1);
%  [~, Ind] = min(abs(y-r));
%  Add_HGV =round(x(Ind));
%  figure 
%   plot(x,y)
%   xlabel('Number of Scenarios')
%     ylabel('Probability')
%     title('Example of Kernel Density Estaimator')
%    legend('Bandwidth=6')
%      grid on 