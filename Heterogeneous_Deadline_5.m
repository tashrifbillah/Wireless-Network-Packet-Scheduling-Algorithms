%%%% Comparison of EDF_LDF....Only LDF with arbitrary tie braking
%%%% Improved drop calculation
%%%% Generates QOS graph at 10^e times slots
%%%% Generates Deficit vs Time graph
%%%% Also generates QOS graph at 10^6 times slots

clc;
clear all;
close all;

tic;

e= 3:5;
DL= zeros(1,length(e));
DW= zeros(1,length(e));
% File_Name= 'Heter_Deadline_20_Links_1_000_000_Time_Slots';

for r= 1:length(e)

T= 10^e(r);
% lambda= [0.1*ones(1,2) 0.09 0.07*ones(1,5) 0.04*ones(1,5) 0.02*ones(1,7)];
% lambda= [ 0.23 0.20 0.16 0.14 0.09 0.07 0.04 0.02 0.02 .01 ];
% lambda= ones(1,10)*0.7;
lambda= [ 0.20 0.16 .01 0.14 0.09 0.23  0.07 0.04 0.02 0.02 ]; % An arbitrary lambda
% lambda= [0.4 0.35 0.24];
% lambda= [0.54 0.3 0.15];
% lambda= [0.9 0.9 0.9];

L= length(lambda);
lambda= fliplr(lambda);

% p= [1.0000    0.9600    0.9344    0.9331    0.9144    0.9041    0.8802    0.8746    0.8716    0.8702]; % Required QOS for arbitrary lambda
% p= [1.0000    0.9980    0.9940    0.9902    0.9826    0.9698    0.9545    0.9329    0.9116    0.8893]; % Required QOS for the flipped lambda
% p= [1.0000    0.9540    0.9232    0.9035    0.8890    0.8810    0.8753    0.8723    0.8709    0.8694]; % Required QOS for the given lambda
p= ones(1,L)*0.80;
% p= [0.9 0.85 0.75];
% p= [0.9 0.9 0.8 0.8 0.7 0.7 0.7 0.7 0.65 0.65]; % Delivery ratio
% p= [0.9*ones(1,3) 0.8*ones(1,5) 0.75*ones(1,5) 0.65*ones(1,7)];

% tao= [3*ones(1,3) 1*ones(1,4) 2*ones(1,3)];
% tao= [4*ones(1,3) 2*ones(1,4) 3*ones(1,3)];
% tao= [5*ones(1,3) 2*ones(1,4) 3*ones(1,3)];
% tao= [ 4*ones(1,3) 1*ones(1,8) 3*ones(1,9)];
tao= 2*ones(1,L);
% tao= [4 2 3];

a= zeros(L,T);

for j=1:L
    a(j,: )= rand(1,T)>(1-lambda(j));
end

disp('Appropriateness of the data');
mean(a,2)
sum(mean(a,2))


% Traffic Pattern
% Traffic= zeros(L,T);
% for j= 1:L
%     for t= 1:T
%         Traffic(j,min(t:t+tao(j)-1,T))= a(j,t)+Traffic(j,min(t:t+tao(j)-1,T));
%     end
%     
% end 


%%%% Deficit queue based schedule
 
s= zeros(L,T);
ED= s;
D= s;
D1= s;
A= s;
R= s;
prior= s;
A(:,1)= a(:,1);

% Stack of deadlines
z= zeros(sum(tao),T);

% prior(:,1)= a(:,1);
% if sum(a(:,1))
%     s(datasample(find(a(:,1)),1),1)= 1;
% end
% 
% R(:,1)= s(:,1);


%%%%% EDF+LDF

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
 
    D(:,t)= D(:,max(1,t-1))+ a(:,t)*double(rand>1-p(j))';         
    
      
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
  

    %%%% Drop calculation   
 
  if t>1
    for j= 1:L        
        A(j,t)= A(j,t-1)+a(j,t); % Cumulative arrival
        R(j,t)= R(j,t-1)+s(j,t); % Cumulative schedule
        % D1(j,t)= A(j,t)-R(j,t); % Cumulative drop
    end
  end
    
    
    
end



disp(['Total packet arrived ', num2str(sum(sum(a)))]);
disp(['Total packet scheduled ', num2str(sum(sum(s)))]);

disp('Delivery ratio- aggregate');
sum(sum(s))/sum(sum(a))

disp('Delivery ratio- linkwise');
deficit_based= sum(s,2)./sum(a,2)


disp('Difference of arrival and schedule');
sum(a,2)-sum(s,2)
disp('Total drop ');
sum(sum(a,2)-sum(s,2))
disp('Cumulative deficit at the end');
DD= D(:,T)
disp(sum(D(:,T)))

disp('Cumulative drop at the end');
D1(:,T)
sum(D1(:,T))

DL(r)= sum(D(:,T));
figure(r)
plot(1:L, deficit_based, '+', 'LineWidth',2); hold on;
plot(1:L, p, 'LineWidth',2); hold on;

% legend('LDF', 'Given');

xlabel('Links');
ylabel('Achieved delivery ratio');
title(['QOS Comparison of LDF for time 10^' num2str(e(r))]);
grid on;


%%%% Take care of this end
% end

% figure(r+1)
% plot(e, DL, '+', 'LineWidth',2); hold on;
% plot(e, DL, 'LineWidth',2); hold on;
% 
% xlabel('10^e');
% ylabel('Deficit');
% title('Deficit Comparison for LDF');
% grid on;



% %%%% Saving the figures
% 
% h =  findobj('type','figure');
% n = length(h);
% mkdir(File_Name);
% cd(File_Name);
% for i=1:n
% savefig(i,num2str(i));
% end
% cd ..
% 
% toc;
% tic;


%%
%%%% Only LDF with arbitrary tie braking



% for r= 1:length(e)
% 
% T= 10^e(r);

%%%% Deficit queue based schedule
 
s= zeros(L,T);
ED= s;
D= s;
D1= s;
A= s;
R= s;
prior= s;
A(:,1)= a(:,1);

% Stack of deadlines
z= zeros(sum(tao),T);

% prior(:,1)= a(:,1);
% if sum(a(:,1))
%     s(datasample(find(a(:,1)),1),1)= 1;
% end
% 
% R(:,1)= s(:,1);


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
            prior(j,t)= 1;
        end

    end   
    
        
    prior_ind= 0;    
    if sum(prior( :,t))
        prior_ind= find(prior( :,t));               
    end 
        
       
 %%%% Deficit queue update method
 
    D(:,t)= D(:,max(1,t-1))+ a(:,t)*double(rand>1-p(j))';         
    
      
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
  

       
 
  if t>1
    for j= 1:L        
        A(j,t)= A(j,t-1)+a(j,t); % Cumulative arrival
        R(j,t)= R(j,t-1)+s(j,t); % Cumulative schedule
        % D1(j,t)= A(j,t)-R(j,t); % Cumulative drop
    end
  end
     

end


disp(['Total packet arrived ', num2str(sum(sum(a)))]);
disp(['Total packet scheduled ', num2str(sum(sum(s)))]);

disp('Delivery ratio- aggregate');
sum(sum(s))/sum(sum(a))

disp('Delivery ratio- linkwise');
deficit_based_without= sum(s,2)./sum(a,2)


disp('Difference of arrival and schedule');
sum(a,2)-sum(s,2)
disp('Total drop ');
sum(sum(a,2)-sum(s,2))
disp('Cumulative deficit at the end');
DD= D(:,T)
disp(sum(D(:,T)))

disp('Cumulative drop at the end');
D1(:,T)
sum(D1(:,T))

DW(r)= sum(D(:,T));
figure(r)
plot(1:L, deficit_based_without, '+', 'LineWidth',2); hold on;
% plot(1:L, p, 'LineWidth',2)

% legend('LDF','Given','LDF Without');

xlabel('Links');
ylabel('Achieved delivery ratio');
title(['QOS Comparison of LDF for time 10^' num2str(e(r))]);
grid on;

%%%% Take care of this end

end

% figure(r+1)
% plot(e, DW, 'k+', 'LineWidth',2); hold on;
% plot(e, DW, 'LineWidth',2); hold on;
% 
% xlabel('10^e');
% ylabel('Deficit');
% title('Deficit Comparison for LDF');
% grid on;


% end

%%

h =  findobj('type','figure');
n = length(h);
for i=1:n
    figure(i);
    legend('EDF-LDF','Given','LDF with Arbitrary');
end


figure(n+1);
plot(e,DL, e,DW, 'linewidth',2); grid on;
legend('EDF-LDF','LDF with arbitrary');
xlabel('10^e');
ylabel('Deficit');
title('Cumulative Deficit Comparison Among Policies');


figure(n+2);
plot(1:L,p,1:L,deficit_based,1:L,deficit_based_without,'linewidth',2);
legend('Given','EDF-LDF','LDF with Arbitrary');
grid on;
xlabel('Links');
ylabel('QOS');
title('Achieved QOS Comparison Among Policies');




%%% Saving the figures

h =  findobj('type','figure');
n = length(h);
mkdir(File_Name);
cd(File_Name);
for i=1:n
savefig(i,num2str(i));
end
cd ..

toc;



