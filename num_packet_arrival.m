clc;
clear all;
close all;

T= 10^5;
% lambda= [0.4 0.3 0.29];
% lambda= [0.5 0.4 0.1];
% lambda= [0.4 0.3 0.15 0.09];
lambda= 0.1*ones(1,10);
% lambda= [ 0.23 0.20 0.16 0.14 0.09 0.07 0.04 0.02 0.02 .01 ];
% lambda= [ 0.5 ones(1,9)*0.5/9 ];
% lambda= [0.5 0.5];
% lambda= [ 0.20 0.16 0.23 .01 ];
% lambda= [0.1*ones(1,2) 0.09 0.07*ones(1,5) 0.04*ones(1,5) 0.02*ones(1,7)];


% L= 3; % length(lambda);
% lambda= rand(1,L);
% lambda= lambda/sum(lambda);
% lambda= round(lambda,3);


L= length(lambda);
a= [ ];
for j=1:L
    a(j,: )= rand(1,T)>(1-lambda(j));
end

tao_max= 3;

soj= 0;
i=1;
while i<= T-tao_max
% for i=1:T-tao_max
    
    if sum(a(:,i))
%        i
%        disp('sum of a(:,i)');
%        sum(a(:,i))
       backup= i+tao_max-1;
       
       j= i;
       while j<=min(backup, T)
           % for j= i+1:backup % i+tao_max-1
           if sum(a(:,j))
%                j
%                disp('sum of a(:,j)');
%                sum(a(:,j))
               backup= j+tao_max-1;
           end
           
           j= j+1;
           
       end
       
       soj= soj+(backup-i+1);
%        backup
       i= backup;
        
    end
    
    i= i+1;
    
end

soj
soj/T



num_arrivals= zeros(1,L);
z= [ ];
for j= 1:L

    z= nchoosek(1:L,j);
        
    for i= 1:nchoosek(L,j)
       arrive= z(i,: );
       noarrive= setdiff(1:L,arrive);
       num_arrivals(j)= num_arrivals(j)+ prod(lambda(arrive))*prod(1-lambda(noarrive)); 
        
    end


end

num_arrivals= [prod(1-lambda) num_arrivals];

atot= sum(a);
z= zeros(1,L);
for i= 0:L
    z(i+1)= length(find(atot==i));
end



disp('  #of Packets   Prob     Actual');
[(0:L)' num_arrivals' z'/T]


