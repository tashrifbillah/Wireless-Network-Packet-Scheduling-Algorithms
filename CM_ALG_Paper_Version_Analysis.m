% Comparison between CM, and Reduced Cost
% Reduced cost is a policy when all the deadlines are one
% In the reduced case, we know that scheduling the highest
% weigth link with existing packet is the optimal policy
% Complete paper approach

clc;
clear all;
close all;

tic;

% L= 10;
% temp_lam= rand(1,L);
% temp_lam= temp_lam/sum(temp_lam);
% lambda= round(temp_lam,3);

% lambda= [ 0.23 0.20 0.16 0.14 0.09 0.07 0.04 0.02 0.02 .01 ];
% % lambda= linspace(0.01,1,20); lambda= lambda/sum(lambda); lambda= round(lambda,3);
% lambda= logspace(-2,0,10); lambda= lambda/sum(lambda); lambda= round(lambda,3);
% lambda= [(1-0.9)/9*ones(1,9) 0.9];

% A random combination of lambda that does not maintain an increasing or
% decreasing sequence
% lambda= 0.1*ones(1,10);
% lambda= [ 0.20 0.16 0.23 .01 0.14 0.09 0.07 0.04 0.02 0.02 ];
% lambda= 0.1*ones(1,10);
% lambda= [0.54 0.3 0.15];
% lambda= [0.4 0.3 0.29];
% lambda= [0.4 0.3 0.2 0.1];
% lambda= [0.3 0.25 0.20 0.15 0.1];
% lambda= 1/5*ones(1,5);
% lambda= [0.3 0.25 0.18 0.13 0.1 0.04];
% lambda= [0.5 0.4 0.1];
L= 10;
lambda= ones(1,L)*1/L;
% lambda= [ 0.7 0.3 ];
% lambda= [0.5 0.5];
% lambda= [0.8 0.8 0.8];

% tao= [3*ones(1,3) 1*ones(1,4) 2*ones(1,3)]; % 10 Links
% tao= [ 4*ones(1,3) 1*ones(1,8) 3*ones(1,9)]; % 20 Links
% tao= [5 3 4]; % 3 Links

L= length(lambda);
temp_weight= rand(1,L);
w= sort(temp_weight)/sum(temp_weight); % Ascending weight
w= sort(temp_weight,'descend')/sum(temp_weight); % Descending weight
w= round(w,4);

% L= length(lambda);
% w= 0.5.^(1:1:L); % 2.^(L:-1:1);
% w= (L:-1:1);
% w= w/sum(w);

T= 10^5;

tao_max= 2;

for repeat= 1:5
    
cost_CM= zeros(tao_max,1);
cost_red= cost_CM;
a= zeros(L,T);

for j=1:L
    a(j,: )= rand(1,T)>(1-lambda(j));
end

% disp('Appropriateness of the data');
% mean(a,2)
% sum(mean(a,2))    

for tao1= 1:tao_max

    tao1
    tao= tao1*ones(1,L);

%%%%%%%%% Reduced Policy %%%%%%%%%%%%
    if tao1==1

        s= zeros(L,T);
        D1= zeros(L,1);

        for t=1:T

            if sum(a(:,t))       
                ind= find(a(:,t),1);
                s(ind,t)= 1;    
            end

        end

        D1= sum(a,2)- sum(s,2);
        cost_red= w*D1/T./(1:tao_max)';

    end


    
%%%%%%%%%%%%% CM Policy %%%%%%%%%%%%%%%%

sm= zeros(L,T);
chain_res= [ ];
chain_link= [ ];
D2= zeros(L,1);

for t= 1:T   
    
    % Updating chain_link and chain_res on arrival

    if sum(a(:,t))
        
       for j= 1:L
           
            if a(j,t)      
               
               chain_link= sort([chain_link j+0.5]);
               place= find(chain_link==(j+0.5));                 
               chain_res= [chain_res(1:place-1) tao(j) chain_res(place:end)];
               chain_link(place)= j;

            end
       end
       
    end
    
  
    
    if sum(chain_res) 
         
       % Finding max weight packets      
       min_res= min(chain_res);
       max_res= max(chain_res);
              
       sched_res= [ ];
       sched_link= [ ];
       
       for i= min_res:max_res
            temp_ind= find(chain_res==i);
            sched_res= [sched_res chain_res(temp_ind)]; 
            sched_link= [sched_link chain_link(temp_ind)];
            
            temp= sortrows([sched_link;sched_res]');
            sched_link= temp(:,1)';
            sched_res= temp(:,2)';
            
            most= min(i,length(sched_link));            
            % Drop calculation
            for link= sched_link(most+1:end)
                D2(link)= D2(link)+1;
            end
           
            % Schedulable set
            sched_link= sched_link(1:most);
            sched_res= sched_res(1:most);
            
            
       end       
    
       
       % Finding packet to schedule
       [~,ind]= min(sched_res);
       sm(sched_link(ind),t)= 1;

       % Updating chain_link and chain_res
       sched_link(ind)= [];
       sched_res(ind)= [];

       chain_link= sched_link;
       chain_res= sched_res;            

       chain_res= chain_res-1;
    
    end
   
    
end   
   
cost_CM(tao1)= w*D2/T;
QOS(:,tao1)= sum(sm,2)./sum(a,2);


end % tao1 loop ends

disp('QOS as deadline increases from 1:tao_max');
disp(QOS);
per1= cost_CM./cost_red

repeat

disp('Weight transmission');
sum(QOS(:,2)'.*w.*lambda)/sum(w.*lambda)

end % Repeat 3 times and average the cost, loop ends


%%



disp('   CM     Reduced   ');
disp([cost_CM   cost_red]);


figure(1);
plot(1:tao_max, cost_CM, '-o', 1:tao_max, cost_red,'-s', 'linewidth',2); 
legend('CM Cost', 'Reduced Cost'); xlabel('Deadlines'); ylabel('Average Cost');
grid on;
title('Comparison of CM, and Reduced cost');
 
figure(2);
per1= cost_CM./cost_red;
plot(1:tao_max, per1,'-^','linewidth',2);
grid on;
xlabel('Deadlines'); ylabel('Ratio');
title('Ratio of CM to RED');


%%

toc;

% File_Name= 'EDF_SP_10_Links_10^6_All_3_Deadlines';
% 
% h =  findobj('type','figure');
% n = length(h);
% 
% mkdir(File_Name);
% cd(File_Name);
% for i=1:n
% savefig(i,num2str(i));
% end
% cd ..

