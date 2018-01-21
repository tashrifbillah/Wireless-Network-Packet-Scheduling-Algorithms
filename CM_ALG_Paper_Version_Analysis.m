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

% cost_CM= cost_CM/3;
% cost_edf= cost_edf/3;
% cost_red= max(cost_edf)./(1:tao_max)';
% 
% disp('   CM     EDF    Reduced   ');
% disp([cost_CM   cost_edf   cost_red]);
% 
% 
% figure(1);
% plot(1:tao_max, cost_CM, '-^', 1:tao_max,cost_edf,'-o', 1:tao_max, cost_red,'-s', 'linewidth',2); 
% legend('CM Cost', 'EDF Cost','Reduced Cost'); xlabel('Deadlines'); ylabel('Average Cost');
% grid on;
% title('Comparison of CM, EDF, and Reduced cost');
% 
% 
% per1= cost_CM./cost_red;
% per2= cost_CM./cost_edf;
% per3= cost_edf./cost_red;
% figure(2);
% plot(1:tao_max, per1,'-^', 1:tao_max, per2,'-o', 1:tao_max, per3,'-s','linewidth',2);
% grid on;
% xlabel('Deadlines'); ylabel('Ratio');
% legend('CM to RED', 'CM to EDF', 'EDF to RED');
% 
% % CM algorithm correctness
% '     CM     EDF'
% [sum(sm,2)./sum(a,2) sum(se,2)./sum(a,2)]


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

lambda

'The way the proof is written'
'Cost analysis proof for all deadlines, induction over links'

xx= QOS(:,tao_max) - QOS(:,1);
zz= [ ];
for i= 1:L

     value= 0;
     for j=1:i-1

         value= value+(1-prod(1-lambda(1:j)))*w(j+1)*lambda(j+1);

     end

     zz(i,: )= [sum(xx(1:i)'.*w(1:i).*lambda(1:i)) (tao_max-1)/tao_max*value];
    

end
figure;
plot(1:length(zz),zz(:,1), '--o', 1:length(zz),zz(:,2), '-.o', 'linewidth',2);
grid on; legend('Sum(w.lambda.diff(QOS))','Probabilistic');
title('Proof evolution over links');
xlabel('Links');


'Inductive proof, only incremental cost'
mm= [];
for i= 1:L
    mm(i,: )= ([(QOS(i,tao_max)-QOS(i,1))*w(i)*lambda(i)    (tao_max-1)/tao_max*(1-prod(1-lambda(1:i-1)))*lambda(i)*w(i)]);

end
figure;
plot(1:length(mm),mm(:,1), '--o', 1:length(mm),mm(:,2), '-.o', 'linewidth',2);
grid on; legend('w.lambda.diff(QOS)','Probabilistic');
title('Incremental cost on both sides of the proof');
xlabel('Links');


[(0:L)' [1:tao_max; QOS]]

% [QOS(:,1)+0.8*0.2511 QOS(:,5)]
% [QOS(:,1)+0.5*0.2511 QOS(:,2)]


% (l1*l2)^tao*[2*(l1^tao+l2^tao)+3*(l1^(tao+1)+l2^tao)]



% soj= 0;
% i=1;
% while i<= T-tao_max
%     
%     if sum(a(:,i))       
%        backup= i+tao_max-1;
%        
%        j= i;
%        while j<=min(backup,T)
%            
%            if sum(a(:,j))
%                backup= j+tao_max-1;
%            end
%            
%            j= j+1;
%            
%        end
%        
%        soj= soj+(backup-i+1);
%        i= backup;
%         
%     end
%     
%     i= i+1;
%     
% end
% 
% 
% soj
% [soj/T 1-prod(1-lambda).^tao_max]
% 
% 
% ((QOS(:,1)+tao_max-1).*lambda'/tao_max)
% sum((QOS(:,1)+tao_max-1).*lambda'/tao_max)



% % Packet drop probability
% value= 3*lambda(1)*(1-lambda(1))^2*lambda(2)^3+ 9*(1-lambda(1))*lambda(1)^2*(1-lambda(2))*lambda(2)^2+ 3*(1-lambda(1))*lambda(1)^2*lambda(2)^3+...
%     3*lambda(1)^3*(1-lambda(2))^2*lambda(2)+ 3*lambda(1)^3*(1-lambda(2))*lambda(2)^2+lambda(1)^3*lambda(2)^3
% 
% 
% 
% i= 1:100;
% lam1i= lambda(1).^i;
% lam2i= lambda(2).^i;
% lamc1i= (1-lambda(1)).^i;
% lamc2i= (1-lambda(2)).^i;
% 
% comb= [ ];
% for j= 1:100
%     comb(j)= nchoosek(j+tao_max-2,tao_max-2);
% end
% 
% prob= comb.*(lamc1i.*lam2i + lamc2i.*lam1i);
% expect= prod(lambda)^tao_max*(1+sum(prob));
% 
% expect
% 1-expect*lambda(1)


% sum(w'.*lambda'.*[QOS(:,1)*(1-1/5)+1/5])/sum(w'.*lambda')

for i= 1:tao_max
    
    disp(['Tao ', num2str(i)]);
    sum(QOS(:,i)'.*lambda.*w)/sum(lambda.*w)

end

%%


%  'The way the proof is written'
%  'Cost analysis proof for all deadlines, induction over links'
%   
% for j=2:tao_max
%     'Deadline'
%     disp(j)
%     for i= 2:L
%         sum(w(1:i).*lambda(1:i)*(1-QOS(1:i,j)))-1/j*sum(w(1:i).*lambda(1:i)*(1-QOS(1:i,1)))
% 
%     end
% end
% 
% 
%  'Total cost analysis proof for deadline 5, induction over links'
% 
% for i= 1:L
% x(i,: )= [sum(w(1:i).*lambda(1:i)*(1-QOS(1:i,5))) 1/5*sum(w(1:i).*lambda(1:i)*(1-QOS(1:i,1)))];
% 
% end
% 
% 
% 
% 'Inductive proof, just the incremental changing factor'
% for i= 1:L
% disp([(QOS(i,5)-QOS(i,1))    4/5*(1-prod(1-lambda(1:i-1)))]);
% 
% end


