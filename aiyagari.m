% =====================================================================
% ====    This program solves the Aiyagari self-insurance model   =====
% =====================================================================
clear;
clf;

tic

disp('Solving the Aiyagari self-insurance model by time iteration and the method of endogenous grid points')
disp('    ')
disp('Start by finding the equilibrium interest rate with a bisection method')


% ======================================================================
% ================        User Definition Area       ===================
% ======================================================================


Nstar = 2;                 % number of states for the efficiency shock.
alpha  = 0.4;             % elasticity of output to input of capital
beta   = 1.04^(-4);             % subjective discount factor
prob   = [.95 .05; .9 .1];  
                           % prob(i,j) = probability (A(t+1)=Aj | A(t) =
                           % Ai) - THIS NEEDS TO BE CONSISTENT WITH NSTATES
delta = .025;                % depreciation rate
Z = [1;0];                 % Effective labor shock values. (1=employed, 0=unemployed)
sigma = 2;                 % Coefficient of relative risk-aversion.
phi = 0;                   % Debt limit (for this program, the debt limit should be
                           % set equal to the natural debt limit: phi=-w*min(Z)/r). 

kmin =  phi+1e-3;          % minimum value of the capital grid
kmax =  10;                % maximum value of the capital grid 
N    = 2000;                % number of grid points
grid = 1;                  % choose how you want to space the grid:
                           %  1 = Equally spaced
                           %  2 = Logarithmic spacing (recommended)                           
                           %  3 = Chebyshev nodes 

% compute ergodic probability of the two states
Amat = zeros(Nstar+1,Nstar);
Amat(1:Nstar,1:Nstar) = eye(Nstar)-prob';
Amat(Nstar+1,:) = 1;
emat = eye(Nstar+1);
evec = emat(:,Nstar+1);

Pergodic = (inv(Amat'*Amat))*Amat'*evec;
                           
% ======================================================================
% =========   Nothing below needs to be changed by the user  ===========       
% ======================================================================


% Some useful numbers.
n = ones(1,Nstar); nn = ones(N,1); epsilon = 1e-5;

SSprob = prob^500;

% Following Aiyagari (1994), the initial interest rate is set as
% 1/beta-1-epsilon.
r = 1/beta-1-epsilon;               

% The grid is now computed. NOTE: With the method of endogenous grid
% points, the grid is given for 'kprime' and NOT for 'k'.
if grid == 1;                     % grid when equally spaced
kp = linspace(kmin,kmax,N)';

elseif grid == 2;                 % grid when logarithmically spaced
kptemp = linspace(0,log(kmax+1-kmin),N)';
kp = exp(kptemp)-1+kmin;

elseif grid == 3;                 % Chebyshec collocation nodes.
ZZ     = -cos((2*[1:N]'-1)*pi/(2*N)); 
kp = (ZZ+1)*(kmax-kmin)/2+kmin;
end

kpp = (1-delta)*kp*n;             % Initial guess for tomorrow's policy function (here zero net-investment).

Gamma = ones(N,1)*SSprob(1,:)./N; % Initial distribution (uniform)

% Start of the outer loop (for the interest rate).

s = 0;                               % Initial value of the counter (needed for the bisection method used
                                     % to find the equilibrium interest rate (see Aiyagari (1994))).

d1 = 1;                              % Initial value for the loop's convergence measure.
rcrit = 1e-6;                        % When d1<rcrit: STOP.
while d1 > rcrit

kn = ((r+delta)/alpha)^(1/(alpha-1));% By r implied capital labor ratio
w = (1-alpha)*kn^(alpha);            % Wage implied by the capital labor ratio.

s = s+1;                             % s is the # of iterations started.

mp = nn*w*Z'+(1+r)*kp*n;             % Cash on hand tomorrow.

% Start of inner loop (Household's problem).

d2 = 1;                         % Initial value for the loop's convergence measure.

m0 = 0;                         % Initial value for cash on hand.

while d2>1e-6                                
    EMUp = (mp-max(kpp,phi)).^(-sigma)*prob';                % Expected marginal utility tomorrow.
    m = (beta*(1+r)*EMUp).^(-1/sigma)+kp*n;                  % Cash on hand today.   
    d2 = max(max(abs(m0-m)./(1+abs(m))));                    % Update convergence measure.
    m0 = m;                                                  % Update initial value for cash on hand
    for i=1:length(Z)
         kpp(:,i) = interp1(m(:,i),kp,mp(:,i),'linear','extrap');    % Update the policy function by interpolation.
    end
end

% That's it: Inner loop completed -- household's problem is solved.

% Find the ergodic distribution. Here we wish to find the measure of
% agentes at each node of the grid kp, for each exogenous state (i.e. employed and unemployed)

k = (m-nn*w*Z')/(1+r);                                  % The grid of capital (N x Nstar)
d3 = 1;                                                 % Convergence measure for the distribution.

while d3>1e-6
    
for j=1:Nstar
G(:,j) = interp1(kp,Gamma(:,j),k(:,j));                 % G gives us the measure of agents at each k.
end

G(isnan(G)) = 0;                                        % Set to zero if k was outside of the domain (given by [kmin,kmax]).
Gamma1 = G*prob;                                        % Theoretically, Gamma1 is the updated Gamma, but due to the imperfect
Gamma1 = [SSprob(1,1)*Gamma1(:,1)/(sum(Gamma1(:,1))), ...
    SSprob(1,2)*Gamma1(:,2)/(sum(Gamma1(:,2)))];        % interpolation, it is needed to be normalized to sum to one.
d3 = max(max(abs(Gamma1-Gamma)./(1+Gamma)));            % Update the convergence measure.   
Gamma = Gamma1;                                         % And finally update Gamma as the normalized Gamma1.

end

meankn = sum(kp'*Gamma)/(sum(Gamma(:,1)));              % The average capital labor ratio.
rimplied = alpha*meankn^(alpha-1)-delta;                % and the implied interest rate 
wage = (1-alpha)*(meankn)^alpha;


% The Bisection Method update for the interest rate.

% The first iteration is treated differently. This is to get an initial
% bracket [r1,r2] where the equilibrium interest rate lies within.
% See Aiyagari (1994) for more on this procedure.

if s==1
    
    r1 = max(rimplied,0);
    r2 = r;
end

if rimplied>r
    r1 = r;
else
    r2 = r;
end
r3 = (r1+r2)/2;

if d1>rcrit
    disp(['Interest rate bracket' '  ' 'Conv. metric'])
    disp([r1 r2 d1])
end

d1 = abs(r3-r)./(1+abs(r));
r = r3;

% Bisection update completed.


end

% Program is done. Print out results 

disp('   ')
disp('Done!')
disp('   ')
disp('Market clearing interest rate is found at:')
disp(r)
disp('Using')
disp(s)
disp('iterations with the bisection method, the process took')
disp(toc)
disp('seconds.')
disp('   ')
disp('PARAMETER VALUES')
disp('   ')
disp('    alpha      beta      delta  ')
disp([alpha, beta, delta])
disp('   ')
disp('Values of Productivity')
disp(Z)
disp('Transition Probability Matrix for A')
disp(prob)
disp('   ')
disp('RESULTS')
disp('   ')
disp('    mean of K ')
disp(         sum(kp'*Gamma))
disp(' as compared to mean of deterministic steady-states: ')
disp(((1/beta-(1-delta))/(alpha))^(1/(alpha-1)))

plot(k,kp);
xlim([kmin kmax])
ylim([kmin kmax])
title('Incomplete Markets Model: POLICY FUNCTION for kprime');
hold on
plot(kp,kp,'k');
hold off

pause

plot(k,k*(1+r)+nn*w*Z'-kp*n);
xlim([kmin kmax])
title('Incomplete Markets Model: POLICY FUNCTION for C');

pause

plot(k,kp*n-(1-delta)*k);
xlim([kmin kmax])
title('Incomplete Markets Model: POLICY FUNCTION for I');

pause

Gamma1 = Gamma;
Gamma1(:,1) = Gamma(:,1)/Pergodic(1);
Gamma1(:,2) = Gamma(:,2)/Pergodic(2);

plot(kp,Gamma1);
xlim([kmin kmax])
title('Stationary DISTRIBUTION OF CAPITAL Conditional on state');

Gammatot = Gamma*Pergodic;

pause

plot(kp,Gammatot);
xlim([kmin kmax])
title('Stationary DISTRIBUTION OF CAPITAL');

pause

% And that's it.