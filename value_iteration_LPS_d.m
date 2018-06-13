% ----------------- Value iteration - Optimal ------------------------

function[r,VOpt, gOpt]=value_iteration_LPS_d(N1,N2,p,q1,q2,d1,d2, c1,c2)

disp('Value Iteration - Minimum Cost (= Best case)')

% mu11 = mu(1,1); % When state of env. is 1, departure rate for bandit 1
% mu12 = mu(1,2); % When state of env. is 1, departure rate for bandit 2
% mu21 = mu(2,1); % When state of env. is 2, departure rate for bandit 1
% mu22 = mu(2,2); % When state of env. is 2, departure rate for bandit 2
% 
% lambda11 = lambda(1,1); % When state of env. is 1, arrival rate for bandit 1
% lambda12 = lambda(1,2); % When state of env. is 1, arrival rate for bandit 2
% lambda21 = lambda(2,1); % When state of env. is 2, arrival rate for bandit 1
% lambda22 = lambda(2,2); % When state of env. is 2, arrival rate for bandit 2
% 
% q12=q(1); % Rate for environment going from state 1 to state 2
% q21=q(2); % Rate for environment going from state 2 to state 1

%For stability


if p>= q1+q2
    
       disp('Not stable') 
    return   
        
end

%gamma = max(lambda11,lambda21) + max(lambda12,lambda22) + max(mu11,mu21) + max(mu12,mu22) ; % Uniformization parameter

%Probabilities of going from one environment to the other in the time given
%by the exponential r.v. with rate gamma.

%For the stationary measure of the environments
%phi(a) is the stationary measure of the environment in state b
% M = q(1) + q(2) ; 
%     for z= 1:2
%         phi(z) = q(3-z) / M;
%     end


% coefficients for a linear cost function (b linear part)


%precgamma = 10^(-8)*gamma; % numerical precision

VOpt = zeros(N1+1,N2+1); % Initializing value function
r = zeros(N1+1,N2+1); % Initializing the array that keeps the policy
ConvergencePrecision = 10^-4; % The algorithm will stop after this precision is reached
MaxNumIterations = 80000;     %Maximum number of iterations if precision is not reached

iter = 0;

while (iter<MaxNumIterations)
    VOptold = VOpt;
    VOpt = zeros(N1+1,N2+1);
    MinDiff = realmax;
    MaxDiff = -realmax;
        for j = 0:N1 % j is the state of bandit 1
            for k = 0:N2 % k is the state of bandit 2
                %Adding the linear cost 
                VOpt(j+1,k+1) = cost(j,k,c1, c2);
                
                %Adding the cost of value function that do not depend on action 
                for b1 = 0:j
                    for b2 = 0:k
                        VOpt(j+1,k+1) = VOpt(j+1,k+1) + (1-p)*LPS_prob(j,b1,d1, q1) ...
                            * LPS_prob(k,b2,d2, q2)*VOptold(j+1-b1,k+1-b2);
                            
                    end
                end
                
               
                
                m1=0;
                m2=0;
               
                   
                % In what follows m1 refers to the transition when arrival
                % is dispatched to route to server 
                % m2 refers to the transition when arrival
                % is dispatched to server 2
                % m_block is the arrival transition and extra cost when the 
                % system decides to block an arriving user to the system
                for b1 = 0:j
                    for b2 = 0:k
                        m1 = m1 + LPS_prob(j,b1,d1, q1) ...
                            * LPS_prob(k,b2,d2, q2)*VOptold(min(j+1-b1,N1)+1,k+1-b2);
                        m2 = m2 + LPS_prob(j,b1,d1, q1) ...
                            * LPS_prob(k,b2,d2, q2)*VOptold(j+1-b1,min(k+1-b2,N2)+1);
                            
                    end
                end
                

                               
                if (m1>=m2)
                     r(j+1,k+1) = 2; % r determines the routing decision
                     VOpt(j+1,k+1) = VOpt(j+1,k+1) + p*m2;
                elseif (m2>m1) 
                     r(j+1,k+1) = 1;
                     VOpt(j+1,k+1) = VOpt(j+1,k+1)+p*m1;
                end
                    
             end
        end
        
    % save V preventively (for the case it gets stuck)
    if mod(iter,1000)==0
        save VMinCostSame VOpt;
    end

    if mod(iter,10)==0               %Every tenth iteration we update MinDiff and MaxDiff
        MinDiff = min(min(VOpt-VOptold))
        MaxDiff = max(max(VOpt-VOptold))
        Difference = MaxDiff-MinDiff
        if (MaxDiff-MinDiff)<ConvergencePrecision % algorithm has converged
            save rMinCostSame r;
            gOpt=(MaxDiff+MinDiff)/2;
            save gMinCostSame gOpt;
            
            fprintf('Value iteration converged in %3d iterations\n',iter+1);
            return
        end
    end
            
    if iter==(MaxNumIterations-1) % algorithm reached MaxNumIterations
        save rMinCostSame r;
        gOpt=(MaxDiff+MinDiff)/2;
        save gMinCostSame gOpt;
        %gOpt
        if MaxNumIterations >= 20000000
            display('Does not converge!')
            return
        end
            
        antw = 20000;
        if (antw >0 )
            MaxNumIterations=MaxNumIterations+antw;
        else 
            return
        end
    end
    
    iter=iter+1;
        
end
end
function c = cost(j,k,c1,c2)
% We have linear cost 
    c = (j)*c1 + (k)*c2;
    
end

function p = LPS_prob(current_state,departures,d, q)
% Computing the probability under LPS-d scheme
    if(departures <= min(current_state,d))
        if(d <= current_state)
            temp1 = (1-q/d)^(d-departures);
            temp2 = (q/d)^departures;
            p = nchoosek(d,departures)*temp1*temp2 ; 
        else
            temp1 = (1-q/current_state)^(current_state-departures);
            temp2 = (q/current_state)^departures;
             p = nchoosek(current_state,departures)*temp1*temp2;            
        end
    else
        p=0;
    end

end


