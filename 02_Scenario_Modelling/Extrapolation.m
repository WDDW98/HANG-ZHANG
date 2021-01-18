% -----------------------
%
% Extrapolation.m
%
% sript extrapolates load effects to 75 or 1000 years return period
%
% Hang Zhang
% 13/01/2021

% -----------------------
% 4. Probability Paper Plot
% -----------------------





function [BR] = Extrapolation(Return,DaMax_SM_C)

if Return == 75
    str = '75 Years Return Period ';
elseif Return == 1000
    str = '1000 Years Return Period ';
end





% Clear 0 in max daily bridge response
DaMax_SM_C(all(DaMax_SM_C==0,2),:)=[] ;
% 
% DaMax = Daily_M;

% Organise data
Re = sort(DaMax_SM_C);           % sort Daily Max 'Response' - low to high
N = length(Re);


% Define Plotting Points - Weibull Distribution
m=[1:1:N]';                  %

PlPt = m./(N+1);             % for a set of N ov
Fs = -log(-log(PlPt));



%---------------------------------------------------

%a standard extremal variate derived from the Gumbel distribution
n = Return*50*5;
Fr = 1-1/n;
G = -log(-log(Fr))

% Fit the first 30% data of measured L1
r = round(2*sqrt(length(Re)));
Re_ft = [Re(length(Re)-r:(length(Re)))];
Fs_ft = [Fs(length(Re)-r:(length(Re)))] ;
d = 1;
fit = polyfit(Re_ft,Fs_ft,d);

BR = (G-fit(2))/fit(1)
xfit = [Re_ft(1):0.01:BR];
yfit = fit(1)*xfit+fit(2);




    % Plot with y variate on y axis
figure
hold on
scatter(Re,Fs,   'r')


% plot trend line
plot (xfit, yfit,   'r')


% plot return period line
line([0,BR],[G, G])

% predict response in return year
line([BR,BR],[-2 G])



ylim([-2, 1.1*G])
xlim([round(min(DaMax_SM_C)-1),BR+1])

set(gca, 'XTick', [0:2:22])
text (2,G+0.5,str)



% ax = gca;
% ax.XAxis.Exponent = 1;                         % set to required exponent 10^3 10^6 10^9...

grid on

 %   title('Total Load Max Daily Bridge Response in L1 ');

ylabel('Standard Extremal Varuate');

% if InfType == 1
%     xlabel('Reaction (kN)');
% elseif InfType == 2
%     xlabel('Midspan Bending Moment (kNm)');
% elseif InfType == 3
xlabel('Load Intensity (kN/m)');
legend ('LSSM incl Correlation')

ax.XAxis.Exponent = 3;                          % set to required exponent
 
box on


 
