% Comparison between CM, and Reduced Cost
% Complete paper approach

clc;
clear all;
close all;

tic;

'The only correct version'

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
lambda= ones(1,L)*0.9;
% lambda= [ 0.7 0.3 ];
% lambda= [0.5 0.5];
% lambda= [0.8 0.8 0.8];


% lambda= [0.414 0.586]; % Two links, gamma_min
% lambda= [ 0.1736    0.4068    0.4196]; % Three links, gamma_min
% lambda= [0.0100    0.1837    0.3994    0.3868    0.0100    0.0100]; % Six links, gamma_min
% lambda= [0.1646    0.3966    0.4188    0.0100    0.0100]; % 5 links, gamma_min
% lambda= [0.1680    0.4014    0.4206    0.0100]; % 4 links, gamma_min
% lambda= [ 0.1573    0.3760    0.3967    0.0100    0.0100    0.0100    0.0100    0.0100    0.0100    0.0100]; % 10 links, gamma_min
% lambda=  [0.0100    0.1786    0.3832    0.3681    0.0100    0.0100 0.0100    0.0100    0.0100    0.0100]; % 10 links, gamma_min
% lambda= fliplr(lambda);
% L= length(lambda);

% tao= [3*ones(1,3) 1*ones(1,4) 2*ones(1,3)]; % 10 Links
% tao= [ 4*ones(1,3) 1*ones(1,8) 3*ones(1,9)]; % 20 Links
% tao= [5 3 4]; % 3 Links

L= length(lambda);
% temp_weight= rand(1,L);
% w= sort(temp_weight)/sum(temp_weight); % Ascending weight
% w= sort(temp_weight,'descend')/sum(temp_weight); % Descending weight
% w= round(w,4);

% L= length(lambda);
w= 0.8.^(1:1:L); % 2.^(L:-1:1);
% w= [ 10 7 4.9 3.43 2.40 ];
% w= [10 8 6.4 5.12 4.1 ];
% w= [ 10 0.1 0.05 0.03 0.02];
% w= (L:-1:1);
% w= w/sum(w);

T= 10^5;

tao_max= 2;
OUR= [ ];
PAPER= [ ];
CM= [ ];

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

for tao1= tao_max

    tao1
    tao= tao1*ones(1,L);

%%%%%%%%%%%% CM Policy %%%%%%%%%%%%%%%%

QOS= [ ];
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
   

QOS(:,tao1)= sum(sm,2)./sum(a,2); 
% disp('CM Weight transmission');
CM= [CM sum(QOS(:,2)'.*w.*lambda)/sum(w.*lambda)];







%%%%%%% Our Algorithm %%%%%%%%%

sm= zeros(L,T);
chain_res= [ ];
chain_link= [ ];
D2= zeros(L,1);
QOS= [ ];


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
    
       
       
              % RMG schedule
       % Finding packet to schedule
       [~,inde]= min(sched_res);
       indh= 1;

       if inde==indh
          ind= 1;       
       else
          packet_e= sched_link(inde);
          packet_h= sched_link(1);
    
          p= 0.99;
          
       if w(inde)/w(indh)<1 && w(inde)/w(indh)>p
          
          prob= rand>1-p; 
%         prob= rand>(1-w(packet_e)/w(packet_h));  
          


           if prob
               ind= inde;
           else
               ind= indh;
               if sched_res(inde)==1
                   sched_res(inde)= [ ];
                   sched_link(inde)= [ ];
               end
           end
           
       elseif w(inde)/w(indh)<p && w(inde)/w(indh)>0
           
       
           ind= indh;
           if sched_res(inde)==1
               sched_res(inde)= [ ];
               sched_link(inde)= [ ];
           end
           
       end    
           
       end
       
       sm(sched_link(ind),t)= 1;

       % Updating chain_link and chain_res
       sched_link(ind)= [];
       sched_res(ind)= [];       
           
       chain_link= sched_link;
       chain_res= sched_res;            

       
       chain_res= chain_res-1;
       
    end
    
end

QOS(:,tao1)= sum(sm,2)./sum(a,2);
% disp('Our Weight transmission');
OUR= [OUR sum(QOS(:,2)'.*w.*lambda)/sum(w.*lambda)];





%%%%%%%%%%%% Paper Algorithm %%%%%%%%%%%%%%%

sm= zeros(L,T);
chain_res= [ ];
chain_link= [ ];
D2= zeros(L,1);
QOS= [ ];

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
    
       
       
              % RMG schedule
       % Finding packet to schedule
       [~,inde]= min(sched_res);
       indh= 1;

       if inde==indh
          ind= 1;       
       else
          packet_e= sched_link(inde);
          packet_h= sched_link(1);
          
%           prob= rand>0.26; 
          prob= rand>(1-w(packet_e)/w(packet_h));  
              
           if prob
               ind= inde;
           else
               ind= indh;
               if sched_res(inde)==1
                   sched_res(inde)= [ ];
                   sched_link(inde)= [ ];
               end
           end

       end
       
       sm(sched_link(ind),t)= 1;

       % Updating chain_link and chain_res
       sched_link(ind)= [];
       sched_res(ind)= [];       
           
       chain_link= sched_link;
       chain_res= sched_res;            

       
       chain_res= chain_res-1;
       
    end
    
end

QOS(:,tao1)= sum(sm,2)./sum(a,2);
% disp('Paper Weight transmission');
PAPER= [PAPER sum(QOS(:,2)'.*w.*lambda)/sum(w.*lambda)];


end % tao1 loop ends


repeat

end % Repeat 3 times and average the cost, loop ends


N= length(CM);
plot(1:N,CM,'-o', 1:N,OUR,'--^', 1:N, PAPER, '-s', 'linewidth',2); grid on;
legend('CM', 'Probabilistic', 'RMG');
title('Ration of total packet weight transmission to total packet weight arrival');
xlabel('Iteration'); ylabel('Gain');



