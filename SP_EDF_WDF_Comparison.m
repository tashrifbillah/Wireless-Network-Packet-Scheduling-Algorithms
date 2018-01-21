%%%% Multiclass EDF algorithm, it does EDF with class weight based tie braking
%%%% Static priority, it schedules a link with packet in descending order of factor= tao*lambda 
%%%% Generates Weighted loss vs Time graph for both the policies
%%%% QOS graph is plotted at 10^6 times slots

clc;
clear all;
close all;

tic;

% File_Name= 'EDF_SP_10_Links_10^6_All_3_Deadlines';

% lambda= [0.1*ones(1,2) 0.09 0.07*ones(1,5) 0.04*ones(1,5) 0.02*ones(1,7)];
lambda= [ 0.26 0.19 0.15 0.14 0.09 0.07 0.04 0.02 0.02 .01 ];
% lambda= [0.4 0.35 0.24];
% lambda= [0.55 0.3 0.15];
% lambda= [0.9 0.9 0.9];

L= length(lambda);

% tao= [3*ones(1,3) 1*ones(1,4) 2*ones(1,3)];
% tao= [ 4*ones(1,3) 1*ones(1,8) 3*ones(1,9)];
tao= 2*ones(1,L);
% tao= [5 3 4];

temp_weight= rand(1,L);
w= sort(temp_weight)/sum(temp_weight); % Ascending lambda, ascending weight
% w= sort(temp_weight,'descend')/sum(temp_weight); % Ascending lambda, descending weight
w= round(w,4);
factor= 100*w.*lambda;
factor= round(factor,4);
[~, link_index]= sort(factor);
static_priority(link_index)= 1:10;

cost_sp= [ ];
cost_edf= [ ];
cost_w_deficit= [ ];

e= 3:6;
for r= 1:length(e)

T= 10^e(r);
a= zeros(L,T);
s= a;
ED= a;
prior= a;
D1= a;


for j=1:L
    a(j,: )= rand(1,T)>(1-lambda(j));
end

disp('Appropriateness of the data');
mean(a,2)
sum(mean(a,2))


%%%%%%%%% EDF Policy

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
    
   
    if sum(prior_ind)        
        temp_chosen= max(w(prior_ind));
        chosen= find(w==temp_chosen);
        s(chosen,t)= 1;
        
        
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



cost_edf= [cost_edf w*D1(:,T)/T];


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





%%%%%%%% Static priority policy
s= zeros(L,T);
ED= s;
prior= s;
D1= s;
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
            prior(j,t)= 1;
        end

    end   
    
    % Links with packets
    prior_ind= 0;    
    if sum(prior( :,t))
        prior_ind= find(prior( :,t));               
    end   
       
          
    if sum(prior_ind)        
        % chosen= datasample(prior_ind,1);
        temp_chosen= max(static_priority(prior_ind));
        chosen= find(static_priority==temp_chosen);
        s(chosen,t)= 1;
                
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

cost_sp= [cost_sp w*D1(:,T)/T];


disp(['Total packet arrived ', num2str(sum(sum(a)))]);
disp(['Total packet scheduled ', num2str(sum(sum(s)))]);

disp('Delivery ratio- aggregate');
sum(sum(s))/sum(sum(a))

disp('Delivery ratio- linkwise');
SP= sum(s,2)./sum(a,2)

disp('Difference of arrival and schedule');
sum(a,2)-sum(s,2)
disp('Total drop ');
sum(sum(a,2)-sum(s,2))

disp('Cumulative drop at the end');
D1(:,T)
sum(D1(:,T))





%%%%%%% Deficit updated according to weight of the links
s= zeros(L,T);
ED= s;
prior= s;
D1= s;
D= s;
z= zeros(sum(tao),T);
A= s;
A(:,1)= a(:,1);
R= s;


w= w*1/max(w);

for t= 1:T    
        
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
        
       
 %%%% Deficit queue update method
 
    D(:,t)= D(:,max(1,t-1))+ a(:,t)*double(rand>1-w(j))';         
    
      
    if sum(prior_ind)        
     
        links_pack= zeros(length(prior_ind),2);

        for k= 1:length(prior_ind)
            links_pack(k,: )= [D(prior_ind(k),t), prior_ind(k)]; % Table of [Deficit, Link Index]
        end
        
        maximum_deficit= max(links_pack( :,1));
        max_def_ind= find(links_pack( :,1)==maximum_deficit); % Find the links that have the largest deficit with packets in the input          
        chosen= links_pack(datasample(max_def_ind,1),2); % links_pack(datasample(def_ind,1),2) is the link chosen to schedule from the above table
        s(chosen,t)= 1;

        D(chosen,t)= max(0,D(chosen,t)-1);
        
        start= sum(tao(1:chosen-1))+1;
        last= sum(tao(1:chosen-1))+tao(chosen);
        jth_link_index= start:last;
        z(find(z(jth_link_index,t),1,'last')+start-1,t)= 0;
    
    end
    
    
    D1( :,t)= D1( :,max(1,t-1));
    
    for j= 1:L
            
            start= sum(tao(1:j-1))+1;
            last= sum(tao(1:j-1))+tao(j);
            jth_link_index= start:last;

            if z(last,t)==1 && s(j,t)==0
                D1(j,t)= D1(j,t)+1;
            end
        
    end
  

      
 
  if t>1
    for j= 1:L        
        A(j,t)= A(j,t-1)+a(j,t); % Cumulative arrival
        R(j,t)= R(j,t-1)+s(j,t); % Cumulative schedule
    end
  end
    
    
    
end

cost_w_deficit= [cost_w_deficit w*D1(:,T)/T];


disp(['Total packet arrived ', num2str(sum(sum(a)))]);
disp(['Total packet scheduled ', num2str(sum(sum(s)))]);

disp('Delivery ratio- aggregate');
sum(sum(s))/sum(sum(a))

disp('Delivery ratio- linkwise');
WDF= sum(s,2)./sum(a,2)

disp('Difference of arrival and schedule');
sum(a,2)-sum(s,2)
disp('Total drop ');
sum(sum(a,2)-sum(s,2))

disp('Cumulative drop at the end');
D1(:,T)
sum(D1(:,T))



end

% mu= sum(s,2)/T;
% [mu lambda']
% mu-lambda'



%%

figure(1)
plot(e, cost_sp,'-ro', e,cost_edf, '-go', e, cost_w_deficit, '-ko', 'LineWidth',2); hold on;
% ylim([1 500]);
xlabel('10^e');
ylabel('Cost due to packet drop');
legend('SP', 'EDF','WDF');
title('Weighted packet loss Comparison for EDF and SP and WDF');
grid on;

figure(2);
plot(1:L, SP, '-ro', 1:L, EDF, '-go', 1:L, WDF, '-ko', 'LineWidth',2); hold on;
legend('SP', 'EDF', 'WDF');
xlabel('Links');
ylabel('Achieved delivery ratio');
title('QOS Comparison of EDF and SP and WDF policy');
grid on;

toc;

% h =  findobj('type','figure');
% n = length(h);
% 
% mkdir(File_Name);
% cd(File_Name);
% for i=1:n
% savefig(i,num2str(i));
% end
% cd ..

