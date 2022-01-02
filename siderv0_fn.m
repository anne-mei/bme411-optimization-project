function [result, fitscore] = siderv0_fn( ant, h_bar_input, c_input, ru1_input, ru2_input,ru3_input );

    % SETUP FOR FORWARD-BACKWARD SWEEP METHOD
    test = -1       ; % Test variable; as long as variable is negative the while loops keeps repeating
    t0   = 0        ;  % Initial time
    tf   = 365      ;    % Final time
    delt = 0.00001  ;% Accepted tollerance
    M    = 365      ;    % Number of nodes
    t    = linspace(t0,tf,M+1); % Time variable where linspace creates M+1 equally spaced nodes between t0 and tf, including t0 and tf.
    h    = tf/M     ;% Spacing between nodes

    % MODEL PARAMETERS
    gamma_i = 1/14  ; % Recovery rate undetected: 14 days
    gamma_d = 1/14  ; % Recovery rate detected: 14 days
    gamma_a = 1/12.4;% Recovery rate threatend: 12.4 days
    beta    = 0.201 ;% Infection rate
    xi_i    = 0.053 ;% Rate infected undetected to acutely symtomatic
    xi_d    = 0.0053;% Rate infected detected to acutely symtomatic
    mu      = 0.0185;% Mortality rate
    mu_hat  = 5*mu  ;% Mortality rate when hospital capacity is exceeded
    nu      = 0.1   ;% Testing rate (no, slow, fast testing = 0, 0.05, 0.10)
    h_bar   = h_bar_input;% Hospital capacity rate (range from 0.00005 to 0.000292)
    psi     = 0.01  ;% Rate of vaccination

    % WEIGHT FACTORS
    c1 = 0 ;
    c2 = 0 ;
    c3 = 0 ;
    b1 = 1 ;
    b2 = 1 ;
    b3 = 1 ;
    
    % Weight on threatened population
    % Weigth on deceased population
    
    % INITIAL CONDITIONS MODEL
    x1 = zeros(1,M+1) ; % Susceptible
    x2 = zeros(1,M+1) ; % Infected - undetected
    x3 = zeros(1,M+1) ; % Infected - detected
    x4 = zeros(1,M+1) ; % Acutely symptomatic - Threatened
    x5 = zeros(1,M+1) ; % Recovered
    x6 = zeros(1,M+1) ; % Deceased
    x7 = zeros(1,M+1) ; % Vaccinated

    x1(1) = 1-0.00001 ;
    x2(1) = 0.00001 ;
    x3(1) = 0 ;
    x4(1) = 0 ;
    x5(1) = 0 ;
    x6(1) = 0 ;
    x7(1) = 0 ;
    
    % INITIAL GUESS FOR OPTIMAL CONTROL INPUT

    ru1 = ru1_input ;
    ru2 = ru2_input ;
    ru3 = ru3_input ;
    c1 = .02 ; c2 = .01 ; c3 = c_input ;
    u1 = c1*ones(1,M+1) ; % Control input for government intervention
    u2 = c2*ones(1,M+1) ; % Control input for testing
    u3 = c3*ones(1,M+1) ; % Control input for vaccinating
    u1 = u1.*ru1 ;
    u2 = u2.*ru2 ;
    u3 = u3.*ru3 ;
    sum_ru3 = sum(ru3);
    
    %cost of vaccination
    costvac = costvacfn(u1,u2,u3,h_bar);
    
    % INITIAL CONDITIONS ADOINT SYSTEM
    L1 = zeros(1,M+1) ;
    L2 = zeros(1,M+1) ;
    L3 = zeros(1,M+1) ;
    L4 = zeros(1,M+1) ;
    L5 = zeros(1,M+1) ;
    L6 = zeros(1,M+1) ;
    L7 = zeros(1,M+1) ;

    L1(M+1) = 0 ;
    L2(M+1) = 0 ;
    L3(M+1) = 0 ;
    L4(M+1) = 0 ;
    L5(M+1) = 0 ;
    L6(M+1) = c2;
    L7(M+1) = 0 ;

    loopcnt = 0 ; % Count number of loops
    while(test < 0)
        loopcnt = loopcnt + 1 ;
        oldu1 = u1 ;
        oldu2 = u2 ;
        oldu3 = u3 ;
        oldx1 = x1 ;
        oldx2 = x2 ;
        oldx3 = x3 ;
        oldx4 = x4 ;
        oldx5 = x5 ;
        oldx6 = x6 ;
        oldx7 = x7 ;

        oldL1 = L1 ;
        oldL2 = L2 ;
        oldL3 = L3 ;
        oldL4 = L4 ;
        oldL5 = L5 ;
        oldL6 = L6 ;
        oldL7 = L7 ;
        
        % SYSTEM DYNAMICCS
        for i=1:M
            % IMPACT HEALTHCARE CAPACITY ON MORTALITY RATE
            if x4(i) <= h_bar
                mu_bar = mu*x4(i) ;
            else
                mu_bar = mu*h_bar + mu_hat*(x4(i) - h_bar) ;
            end

            m11 = -beta*x1(i)*x2(i)*(1-u1(i)) - psi*u3(i)*x1(i);
            m12 = beta*x1(i)*x2(i)*(1-u1(i)) - gamma_i*x2(i) - xi_i*x2(i) - nu*x2(i)*(u2(i));
            m13 = nu*x2(i)*(u2(i))-gamma_d*x3(i)-xi_d*x3(i);
            m14 = xi_i*x2(i)+xi_d*x3(i)-gamma_a*x4(i)-mu_bar;
            m15 = gamma_i*x2(i) + gamma_d*x3(i) + gamma_a*x4(i);
            m16 = mu_bar;
            m17 = psi*u3(i)*x1(i);
            x1(i+1) = x1(i) + h*m11 ; %first order explicit Euler! seemed to work!
            x2(i+1) = x2(i) + h*m12 ;
            x3(i+1) = x3(i) + h*m13 ;
            x4(i+1) = x4(i) + h*m14 ;
            x5(i+1) = x5(i) + h*m15 ;
            x6(i+1) = x6(i) + h*m16 ;
            x7(i+1) = x7(i) + h*m17 ;
            test = 1 ;
        end
    end
%     % PLOT
%     figure(ant)
%     plot(x1)
%     hold on
%     plot(x2)
%     plot(x3)
%     plot(x4)
%     plot(x5)
%     plot(x6,'*')
%     plot(x7)
%     legend('S','Iu','ID','A', 'R','D','V')
%     hold off
%     title_str = append("Population ", string(ant), " Results") ;
%     save_str  = append("Population_", string(ant), ".png") ;
%     title( title_str ) ;
%     saveas( gcf, save_str);

    totdeaths= 100*sum(x6);
    totcost  = costvac + 100000*h_bar;
    
    % FITNESS SCORE / OBJECTIVE FUNCTION
    fitscore = c2*totdeaths + totcost ; %updated to multiply by c2
    
    result = {totdeaths, totcost} ;
%    print_output = [ sum_ru3 costvac totdeaths totcost ]
end