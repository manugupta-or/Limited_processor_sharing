function [Squared_throughput] = Performance_index(n, p, d1,d2,q1,q2, theta, D, beta)


C0= 1;
C1 =1;

Total_States= (n+1)^2;

all_states = ones((n+1)^2, 2);

State_id =ones((n+1)^2, 1);
count =0;
for i = 0:n
        for j = 0:n
                count = count+1;
                all_states(count, 1) = i;
                all_states(count, 2) = j;
                State_id(count) = i+1 + j*(n+1);
        end
end
%%transition of type (1,0) to (0,1)

% disp(all_states)
% disp(State_id)
for i=1:Total_States
    ordered_states(State_id(i),1) = all_states(i,1);
    ordered_states(State_id(i),2) = all_states(i,2);
end
%disp(all_states +1)

TPM = zeros((n+1)^2);

count_route = 0;
X = 0
%Generating TPM
for i = 0:n
        for j = 0:n 
              q1_LPS_index = Whittle_indices_limited_PS(i,p,q1, theta, d1, D, beta);
              q2_LPS_index = Whittle_indices_limited_PS(j,p,q2, theta, d2, D, beta);
           if(q1_LPS_index >= q2_LPS_index && i< n)
              route =1; 
           elseif(q1_LPS_index <= q2_LPS_index && j< n)
                route=2;
           elseif(i==n && j <n)
               route = 2;
           elseif(i<n && j == n)   
               route =1;
           else
               
               route = 0;
               count_route = count_route +1;
               X = Unique_id(i,j,n);
               
           end
                if(i==n && j==n )
                    for count1 = 0:n
                        for count2 = 0:n
                        TPM(Unique_id(n,n,n), Unique_id(n-count1,n-count2,n)) = LPS_prob(n,count1,d1, q1)*LPS_prob(n,count2,d2, q2);
                        end
                    end    
               
                    
                elseif (route == 1)
                    %If arrival happens
                    for k=0:j
                        for count = 0:i+1
                        TPM(Unique_id(i,j,n), Unique_id(i+1-count,j-k,n)) = p*LPS_prob(i,count,d1, q1)*LPS_prob(j,k,d2, q2);
                        end
                    end
                    %If arrival doesn't happen
                    for k=0:j
                        for count = 0:i
                        TPM(Unique_id(i,j,n), Unique_id(i-count,j-k,n)) = TPM(Unique_id(i,j,n), Unique_id(i-count,j-k,n))+ (1-p)*LPS_prob(i,count,d1, q1)*LPS_prob(j,k,d2, q2);
                        end
                    end
                elseif(route == 2)
                    %If arrival happens
                    for k=0:i
                        for count = 0:j+1
                        TPM(Unique_id(i,j,n), Unique_id(i-k,j+1-count,n)) = p*LPS_prob(i,k,d1, q1)*LPS_prob(j,count,d2, q2);
                        end
                    end
                    %If arrival doesn't happen
                    for k=0:i
                        for count = 0:j
                        TPM(Unique_id(i,j,n), Unique_id(i-k,j-count,n)) = TPM(Unique_id(i,j,n), Unique_id(i-k,j-count,n))+ (1-p)*LPS_prob(i,k,d1, q1)*LPS_prob(j,count,d2, q2);
                        end
                    end

                    
                end
            
        end
end
disp('Total states')
disp(Total_States)


disp('ordered states')
disp(ordered_states)

disp('Transition probability matrix')

P = TPM
sum(P,2)


Pi = Stationary_distribution_new(P);

disp('staionary distribution')
disp(Pi)

E_number=0;
states = ones(Total_States,1);
throughput =0;
Squared_throughput = 0; 


for i=1:Total_States
   states(i) =  ordered_states(i,1) + ordered_states(i,2); 
   E_number  = E_number + states(i)*Pi(i);
   throughput = throughput + (mean_throughput(ordered_states(i,1), d1, q1) ... 
                    + mean_throughput(ordered_states(i,2), d2, q2))*Pi(i);
   Squared_throughput = Squared_throughput + (beta*states(i) + (1-beta)*(C_troughput(ordered_states(i,1), d1, q1, theta) ... 
                    + C_troughput(ordered_states(i,2), d2, q2,theta)))*Pi(i);
                
end

disp('Total expected number')
disp(E_number)

disp('number of times it goes to route')
disp(count_route)

disp('throughput')
disp(throughput)

%disp(Pi)

%disp(state)

Squared_throughput
%Stationary distribution



end

function c = C_troughput(n,d,q, theta)
% We have linear cost 
s = 0;
%     for i = 0:min(n,d)
%         s = s +i^2*LPS_prob(n,i,d, q);
%     end
%     c = s;

    for i = 0:min(n,d)
        s = s -(i*theta-i^2)*LPS_prob(n,i,d, q);
    end
    c = s;
    
end



function sum = mean_throughput(current_state, d, q)
sum = 0;
    for i = 0:d
        sum = sum + i^2*LPS_prob(current_state,i,d, q);
    end
    
end