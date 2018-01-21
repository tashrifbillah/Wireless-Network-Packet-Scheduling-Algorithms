%%%% EDF algorithm, it does EDF with arbitrary tie braking
%%%% Generates Deficit vs Time graph for EDF policy
%%%% QOS graph is plotted at 10^6 times slots

clc;
% clear all;
close all;

tic;

% File_Name= 'EDF_20_Links_10^6_All_3_Deadlines';

e= 3:6;
DL= zeros(1,length(e));

for r= 1:length(e)

T= 10^e(r);
lambda= [0.1*ones(1,2) 0.09 0.07*ones(1,5) 0.04*ones(1,5) 0.02*ones(1,7)];
% lambda= [ 0.26 0.19 0.15 0.14 0.09 0.07 0.04 0.02 0.02 .01 ];
% lambda= [0.4 0.35 0.24];
% lambda= [0.55 0.3 0.15];
% lambda= [0.9 0.9 0.9];

L= length(lambda);

p= ones(1,L)*0.80;
% p= [0.9 0.85 0.75];
% p= [0.9 0.9 0.8 0.8 0.7 0.7 0.7 0.7 0.65 0.65]; % Delivery ratio
% p= [0.9*ones(1,3) 0.8*ones(1,5) 0.75*ones(1,5) 0.65*ones(1,7) 0.6];

% tao= [3*ones(1,3) 1*ones(1,4) 2*ones(1,3)];
% tao= [ 4*ones(1,3) 1*ones(1,8) 3*ones(1,9)];
tao= 2*ones(1,L);
% tao= [5 3 4];

a= zeros(L,T);
s= a;
ED= a;
prior= a;
D1= a;
D= a;


for j=1:L
    a(j,: )= rand(1,T)>(1-lambda(j));
end

disp('Appropriateness of the data');
mean(a,2)
sum(mean(a,2))


% Stack of deadlines
z= zeros(sum(tao),T);

for t=1:T

    for j= 1:L        
        
        start= sum(tao(1:j-1))+1;
        last= sum(tao(1:j-1))+tao(j);
        jth_link_index= start:last;    
        
        temp= max(0,z(jth_link_index,max(1,t-1))-1);
        z(jth_link_index,t)= [0;temp(1:end-1)];
        
        if a(j,t)==1
            z(start,t)= tao(j);
        end
        
        if sum(z(jth_link_index,t))
            ED(j,t)= z(find(z(jth_link_index,t),1,'last')+start-1,t);
        end

    end  
    
    
    deadlines= ED(:,t);
    min_deadline= min(deadlines(deadlines>0));
    if (min_deadline)
        prior(find(ED( :,t)== min_deadline),t)= 1;
    end
    
    prior_ind= 0;    
    if sum(prior( :,t))
        prior_ind= find(prior( :,t));               
    end
    
    %%%%% Deficit calculation
    
    D(:,t)= D(:,max(1,t-1))+ a(:,t)*double(rand>1-p(j))'; 
    
    if sum(prior_ind)        
        chosen= datasample(prior_ind,1);
        s(chosen,t)= 1;
        D(chosen,t)= max(0,D(chosen,t)-1);
        
        start= sum(tao(1:chosen-1))+1;
        last= sum(tao(1:chosen-1))+tao(chosen);
        jth_link_index= start:last;
        z(find(z(jth_link_index,t),1,'last')+start-1,t)= 0;
    end   
        
    
    %%%% Drop calculation
    D1( :,t)= D1( :,max(1,t-1));
    
    for j= 1:L
            
            start= sum(tao(1:j-1))+1;
            last= sum(tao(1:j-1))+tao(j);
            jth_link_index= start:last;

            if z(last,t)==1 && s(j,t)==0
                D1(j,t)= D1(j,t)+1;
            end
        
    end
    
       
        
end


disp(['Total packet arrived ', num2str(sum(sum(a)))]);
disp(['Total packet scheduled ', num2str(sum(sum(s)))]);

disp('Delivery ratio- aggregate');
sum(sum(s))/sum(sum(a))

disp('Delivery ratio- linkwise');
EDF= sum(s,2)./sum(a,2)

disp('Difference of arrival and schedule');
sum(a,2)-sum(s,2)
disp('Total drop ');
sum(sum(a,2)-sum(s,2))

disp('Cumulative drop at the end');
D1(:,T)
sum(D1(:,T))

disp('Cumulative deficit at the end');
D(:,T)
DL(r)= sum(D(:,T));

end

mu= sum(s,2)/T;
[mu lambda']
mu-lambda'



%%

figure(1)
% plot(e, DL, 'k+', 'LineWidth',2); hold on;
plot(e, DL, 'LineWidth',2); hold on;
% ylim([1 500]);
xlabel('10^e');
ylabel('Deficit');
title('Deficit Comparison for EDF');
grid on;

figure(2);
plot(1:L, EDF, 'k+', 'LineWidth',2); hold on;
plot(1:L, p, 'LineWidth',2)
legend('EDF', 'Given');

plot(1:L,EDF, 'LineWidth',2); hold on;
xlabel('Links');
ylabel('Achieved delivery ratio');
title('QOS Comparison of EDF policy');
grid on;

toc;

h =  findobj('type','figure');
n = length(h);

mkdir(File_Name);
cd(File_Name);
for i=1:n
savefig(i,num2str(i));
end
cd ..



figure(8); hold on;
plot(1:20, backup_1,'k', 1:20, backup_3,'g', 'LineWidth',2);
legend('Given','LDF+EDF', 'EDF+LDF','Randomized','LDF','EDF');

figure(1); hold on;
plot(1:20, backup_1,'k', 1:20, backup_3,'g', 'LineWidth',2);
legend('Given','LDF+EDF', 'EDF+LDF','Randomized','LDF','EDF');



