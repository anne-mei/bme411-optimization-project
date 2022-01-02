% Implementation of Ant Colony Optimization

close all ;
clear all ;

tic

% previously randomly generated infection and testing states
% load a constant set while testing vaccination optimization
load('ru1.mat')
load('ru2.mat')

% ANT COLONY VARS 
MaxIteration= 5000; 
COLONY_SIZE = 500 ; % number of ants
evaporation = 0.5 ; % how much pheremone is evaporating in each iteration (percent)
Q           = 0.1 ; % information of total pheromone left by each ant

fitscore    = zeros(1,COLONY_SIZE) ;
ants        = {COLONY_SIZE}  ;
trails      = {COLONY_SIZE}  ;
tau         = (1/COLONY_SIZE)*ones(1,COLONY_SIZE)  ;
bestResults = {MaxIteration}  ;

% SIDERV0 VARS
M           = 365  ;
cr_ru3      = rand ;

% SETUP ANTS
for ant_i=1:COLONY_SIZE
    ru3_input = zeros(1,M+1);
    for i = 1:M+1
        if rand > cr_ru3
            ru3_input (i) = 1;
        end
    end
    h_rand  = randi([5,30]) ;
    h_input = h_rand / 100000 ;
    c_rand  = randi(12) ;
    c_input = c_rand / 1000 ;
    
    % create an initial trail for each ant in the colony
    trails{ant_i,:} = {h_input, c_input , ru3_input} ;
end

iteration = 1 ;
converge  = false ;

while (not(converge) & iteration < MaxIteration)
    % MOVE ANTS 
    if iteration == 1
        for ant_i=1:COLONY_SIZE
            h_input   = trails{ant_i}{1} ;
            c_input   = trails{ant_i}{2} ;
            ru3_input = trails{ant_i}{3} ;
            [ants{ant_i},fitscore(1,ant_i)] = siderv0_fn( ant_i, h_input, c_input , ru1, ru2, ru3_input );
            ants{ant_i}{3} = ant_i ;
        end
    else
        % determine current iterations' trails, dependent on tau (pheromone density)
        for ant_i=1:COLONY_SIZE
            curr_trail = ants{ant_i}{3};
            if rand > tau(1,ants{ant_i}{3}) 
                % follow other ant's path
                trail_rand  = randi(COLONY_SIZE);
                curr_trail  = trail_rand ;
                ants{ant_i}{3} = curr_trail ;    
            end
            % stay on own path
            h_input   = trails{curr_trail}{1} ;
            c_input   = trails{curr_trail}{2} ;
            ru3_input = trails{curr_trail}{3} ;
            [ants{ant_i},fitscore(1,ant_i)] = siderv0_fn( ant_i, h_input, c_input , ru1, ru2, ru3_input );
            ants{ant_i}{3} = curr_trail ;
        end
    end
    
    temp1 = vertcat(ants{:});
    curr_paths_temp = vertcat(temp1{:,3})' ;
    
    % normalize fitness scores
    max_fit = max(fitscore);
    min_fit = min(fitscore);
    norm_fit = 0.01 + (fitscore - min_fit) / ( max_fit - min_fit) ;
    fitscore = norm_fit ;
        
    % UPDATE TRAILS
    for ant_i=1:COLONY_SIZE
        pheremone_contribution = Q/fitscore(1,ant_i) ;
        tau(1,ants{ant_i}{3})  = tau(1,ants{ant_i}{3}) + pheremone_contribution ;
    end
    
    tau = (1 - evaporation)*tau ;
    
    % normalize tau scores
    max_tau = max(tau);
    min_tau = min(tau);
    norm_tau = 0.01 + (tau - min_tau) / ( max_tau - min_tau) ;
    tau = norm_tau ;
    
    % UPDATE BEST 
    [bestScore,bestScore_i] = min(fitscore) ;
    bestDeath = ants{bestScore_i}{1} ;
    bestCost  = ants{bestScore_i}{2} ;
    bestResults{iteration} = {bestScore, bestDeath, bestCost} ;  
    
    % break loop
    if all(curr_paths_temp == curr_paths_temp(1))
        converge = true ;
    end
    
    % ITERATE
    iteration = iteration + 1 ;
end

if (converge)
    % Print final output
    bestAnt = curr_paths_temp(1) ;
    fprintf('Final Output: \n') ;
    fprintf('Iteration %d\n',iteration) ;
    fprintf('h_bar       = %d\n',trails{bestAnt}{1}) ;
    fprintf('c3          = %d\n',trails{bestAnt}{2}) ;
    fprintf('ru3         = %d\n',sum(trails{bestAnt}{3})) ;
    fprintf('TotalDeaths = %d\n',bestDeath) ;
    fprintf('TotalCost   = %d\n',bestCost) ;
    fprintf('Fitness     = %d\n',0.02*bestDeath + bestCost) ;
end

toc