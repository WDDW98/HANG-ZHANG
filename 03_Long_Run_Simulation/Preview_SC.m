% -----------------------
%
% Preview_SC.m
%
% sript preview scenarios on 500m span
%
%
% Hang Zhang
% 18/01/2021





function [] = Preview_SC (L1_SC,L2_SC)


a = randi(length(L1_SC));
L1_SC = L1_SC((a-10:a),1);
dx = 0.01;


for k =1:10
    
    SC = L1_SC{k,1};
    SCC = [];
    
    for i = 1 : length(SC(:,1))
        
        
         Veh = zeros(round(SC(i,30)/dx),1);
         Gap = zeros(round(SC(i,56)/dx),1);
        
        for j = 1 : length (Veh)
            Veh(j,1) = SC(i,13);
        end
        
        SCC = [SCC;Veh;Gap];
        
    end
    if k>1 
       Daymax1 =  [zeros(length(T1000m{k-1,1}),1);SCC];
       T1000m{k,1}=Daymax1;
    elseif k==1
     T1000m{k,1}=SCC;
    end
    
end


L2_SC = L2_SC((a-10:a),1);

for k =1:10
    
    SC = L2_SC{k,1};
    SCC = [];
    
    for i = 1 : length(SC(:,1))
        
         Veh = zeros(round(SC(i,30)/dx),1);
         Gap = zeros(round(SC(i,56)/dx),1);
        
        
        for j = 1 : length (Veh)
            Veh(j,1) = SC(i,13);
        end
        
        SCC = [SCC;Veh;Gap];
        
        
    end
    if k>1 
       Daymax1 =  [zeros(length(T1000m2{k-1,1}),1);SCC];
       T1000m2{k,1}=Daymax1;
    elseif k==1
     T1000m2{k,1}=SCC;
    end
    
end




figure
subplot(2,1,1);
h1 = area(T1000m{1,1});
h1.FaceColor = [0 0.75 0.75];
hold on

h2 = area(T1000m{2,1});
h2.FaceColor = [0 0.95 0.50];
hold on

h3 =area(T1000m{3,1});
h3.FaceColor = [0 0.5 0.5];

hold on
h4 =area(T1000m{4,1});
h4.FaceColor = [0.9 0 0.2];

hold on
h5 =area(T1000m{5,1});
h5.FaceColor = [0.5 0.6 0.2];

ylabel('Weight (kN)');
xlim ([0,500/dx])
ylim ([0,500])

set(gca,'XTick', [0 10000 20000 30000 40000 50000])
set(gca,'xticklabels',{'0','100','200','300','400','500'});

set(gca, 'ytick', [0:250:500]);

title(['Vehicle Team in L1']);
grid on


subplot(2,1,2);
h1 = area(T1000m2{1,1});
h1.FaceColor = [0 0.75 0.75];
hold on

h2 = area(T1000m2{2,1});
h2.FaceColor = [0 0.95 0.50];
hold on

h3 =area(T1000m2{3,1});
h3.FaceColor = [0 0.5 0.5];

hold on
h4 =area(T1000m2{4,1});
h4.FaceColor = [0.9 0 0.2];

hold on
h5 =area(T1000m2{5,1});
h5.FaceColor = [0.5 0.6 0.2];


ylabel('Weight (kN)');
xlim ([0,500/dx])
ylim ([0,550])

set(gca,'XTick', [0 10000 20000 30000 40000 50000])
set(gca,'xticklabels',{'0','100','200','300','400','500'});

set(gca, 'ytick', [0:250:500]);
title(['Vehicle Team in L2']);

set(gca, 'ytick', [0:250:500]);
grid on

end


