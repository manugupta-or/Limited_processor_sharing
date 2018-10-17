% ----------------- Value iteration - Optimal ------------------------

function[r,VOpt]=OptimalPolicySameEnvironments(N1,N2,lambda,mu,q,b)

disp('Value Iteration - Minimum Cost (= Best case)')
%speed scaling if i customers is i / i+1

mu11 = mu(1,1); % When state of env. is 1, departure rate for bandit 1
mu12 = mu(1,2); % When state of env. is 1, departure rate for bandit 2
mu21 = mu(2,1); % When state of env. is 2, departure rate for bandit 1
mu22 = mu(2,2); % When state of env. is 2, departure rate for bandit 2

lambda11 = lambda(1,1); % When state of env. is 1, arrival rate for bandit 1
lambda12 = lambda(1,2); % When state of env. is 1, arrival rate for bandit 2
lambda21 = lambda(2,1); % When state of env. is 2, arrival rate for bandit 1
lambda22 = lambda(2,2); % When state of env. is 2, arrival rate for bandit 2

q12=q(1); % Rate for environment going from state 1 to state 2
q21=q(2); % Rate for environment going from state 2 to state 1

%For stability
if rho(lambda,mu,q) >= 1
    
       disp('Not stable') 
    return   
        
end

gamma = max(lambda11,lambda21) + max(lambda12,lambda22) + max(mu11,mu21) + max(mu12,mu22) ; % Uniformization parameter

%Probabilities of going from one environment to the other in the time given
%by the exponential r.v. with rate gamma.

%For the stationary measure of the environments
%phi(a) is the stationary measure of the environment in state b
M = q(1) + q(2) ; 
    for z= 1:2
        phi(z) = q(3-z) / M;
    end


% coefficients for a linear cost function (b linear part)
b1 = b(1);
b2 = b(2);

precgamma = 10^(-8)*gamma; % numerical precision

VOpt = zeros(N1+1,N2+1); % Initializing value function
r = zeros(N1+1,N2+1); % Initializing the array that keeps the policy
ConvergencePrecision = 10^-6; % The algorithm will stop after this precision is reached
MaxNumIterations = 80000;     %Maximum number of iterations if precision is not reached

iter = 0;

while (iter<MaxNumIterations)
    VOptold = VOpt;
    VOpt = zeros(N1+1,N2+1);
    MinDiff = realmax;
    MaxDiff = -realmax;
        for j = 0:N1 % j is the state of bandit 1
            for k = 0:N2 % k is the state of bandit 2
                %I do a first part where I add costs that do not depend
                %de a, en los 4 posibles climas
                VOpt(j+1,k+1) = gamma*cost(j,k,b);
                m1=0;
                m2=0;
                for d = 1:2
                        
                        
                %agregar casos del borde !
                % f determinan a que envir pasamos
                VOpt(j+1,k+1) = VOpt(j+1,k+1) + phi(d) ...
                    * ( lambda(d,1)*VOptold(min(j+1,N1)+1,k+1)...
                    + lambda(d,2)*VOptold(j+1,min(k+1,N2)+1)...
                    + (max(lambda11,lambda21) + max(lambda12,lambda22) ...
                    - lambda(d,1) - lambda(d,2) ) *VOptold(j+1,k+1) ) ;
                            
                        
                    
                % In what follows m1 refers to the arrival transition when
                % the chosen action is to route to server 1
                % m2 refers to the arrival transition when the chosen action
                % is to route to server 2
                % m_block is the arrival transition and extra cost when the 
                % system decides to block an arriving user to the system

                            m1 = m1 + phi(d) * (mu(d,1)*(j/(j+1))*VOptold(max(j-1,0)+1,k+1)...
                                + (max(mu11,mu21) + max(mu12,mu22)-mu(d,1)*(j/(j+1)))*VOptold(j+1,k+1));
                            m2 = m2 + phi(d) * (mu(d,2)*(k/(k+1))*VOptold(j+1,max(k-1,0)+1)...
                                + (max(mu11,mu21) + max(mu12,mu22)-mu(d,2)*(k/(k+1)))*VOptold(j+1,k+1));
                        end
                               
                if (m1>=m2)
                     r(j+1,k+1) = 2;
                     VOpt(j+1,k+1) = VOpt(j+1,k+1) + m2;
                elseif (m2>m1) 
                     r(j+1,k+1) = 1;
                     VOpt(j+1,k+1) = VOpt(j+1,k+1)+m1;
                end
                    
                end
        end
        
       VOpt = VOpt/gamma;
    % save V preventively (for the case it gets stuck)
    if mod(iter,1000)==0
        save VMinCostSame VOpt;
    end

    if mod(iter,10)==0               %Every tenth iteration we update MinDiff and MaxDiff
        MinDiff = min(min(VOpt-VOptold));
        MaxDiff = max(max(VOpt-VOptold));
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
function c = cost(j,k,b)
% We have linear cost 
    c = (j)*b(1) + (k)*b(2);
    
end
function rh = rho(lambda,mu,q)
%For the stability.

%For the stationary measure of the environments
%phi(a,b) is the stationary measure of the environment of class a in state b
rh = 0;
M = q(1) + q(2) ; 
    for z= 1:2
        phi(z) = q(3-z) / M ;
    end

for a= 1:2
rh = rh + (phi(1) * lambda(1,a) + phi(2) * lambda (2,a)) / ...
    (phi(1) * mu(1,a) + phi(2) * mu(2,a)) ;
end

end

