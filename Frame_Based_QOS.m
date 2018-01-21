clc;
clear all;
close all;

T= 10^5;
% lambda= [0.4 0.3 0.29];
% lambda= [0.4 0.3 0.15 0.09]; % 0.04];
% lambda= 0.1*ones(1,10);
lambda= [ 0.23 0.20 0.16 0.14 0.09 0.07 0.04 0.02 0.02 .01 ];
% lambda= [0.1*ones(1,2) 0.09 0.07*ones(1,5) 0.04*ones(1,5) 0.02*ones(1,7)];

% lambda= fliplr(lambda);
L= length(lambda);
w= 0.5.^(L:-1:1);

a= [ ];
for j=1:L
    a(j,: )= rand(1,T)>(1-lambda(j));
end

tao_max= 3;

soj= 0;
i=1;
while i<= T-tao_max
    
    if sum(a(:,i))

       backup= i+tao_max-1;
       
       j= i;
       while j<=min(backup, T)
           
           if sum(a(:,j))

               backup= j+tao_max-1;
           end
           
           j= j+1;
           
       end
       
       soj= soj+(backup-i+1);
       i= backup;
        
    end
    
    i= i+1;
    
end

soj
soj/T


num_arrivals= ones(1,L);
for i= tao_max+1:L  % This are the links whose delivery ratio is being calculated

    
    
    temp= 0;
    for j= 0:tao_max-1  % This number of arrivals

        z= [ ];
        z= nchoosek(1:i-1,j); % In this indices of links

        for n= 1:nchoosek(i-1,j) % In this indices of links
           arrive= z(n,: );
           noarrive= setdiff(1:i-1,arrive);
           temp= temp+ prod(lambda(arrive))*prod(1-lambda(noarrive)); 

        end
        
    end
    
    num_arrivals(i)= temp;

end

[(1:L)' num_arrivals']

disp('Cost of the frame based policy');
sum((1- num_arrivals).*w)

% atot= sum(a);
% z= zeros(1,L);
% for i= 0:L
%     z(i+1)= length(find(atot==i));
% end
% 
% 
% 
% disp('  #of Packets   Prob     Actual');
% [(0:L)' num_arrivals' z'/T]


