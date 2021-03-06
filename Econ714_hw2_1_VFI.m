% -------------------------------------------------------------------------
% Econ 714
% Homework 2
% Fall 2021 
% Javier Tasso
% Value Function Iteration
% -------------------------------------------------------------------------
clearvars
clc 
cd 'C:\Users\Javier\Dropbox\Aplicaciones\Overleaf\Econ714-PS2'
% -------------------------------------------------------------------------

% Set values of parameters 
alpha = 0.33;
beta = 0.97;
delta = 0.1;

%%
% -------------------------------------------------------------------------
% Question 1 - Find the steady state
% -------------------------------------------------------------------------

% Solve for the steady state level of capital 

    % First find an expression for consumption 
    % cons = max(0, x(1)^(alpha) * x(2)^(1-alpha) - delta * x(1))
    % we plug this expression in the denominator 

    % Now we reduce the thing to a system of two equations 
    
    options = optimset('Display','off');
    F = @(x) [(1-alpha) / max(0, x(1)^(alpha) * x(2)^(1-alpha) - delta * x(1)) ...
        * (x(1)/x(2))^(alpha) - x(2); ... 
        (1/(alpha*beta) - (1-delta)/alpha) * (x(1)/x(2))^(1-alpha) - 1];
    [x,~] = fsolve(F,[1;1],options);

    % Check that everything looks ok
    if x(1)<0 || x(2)<0 || x(2)>1
        
        disp('There may be corner solutions, the program stops here')
        stop

    end

    % Store the solution 
    cons_ss = max(0, x(1)^(alpha) * x(2)^(1-alpha) - delta * x(1)); 
    kapi_ss = x(1);
    labo_ss = x(2); 
    prod_ss = kapi_ss^(alpha) * labo_ss^(1-alpha); 
    valu_ss = (log(cons_ss) - (labo_ss)^2 / 2);

    % Check 
    if abs(prod_ss - (cons_ss + delta * kapi_ss)) > 10^(-7)

        disp('There is something wrong with the steady state values, the program stops here')
        stop 

    end

% Clean and keep the steady state 
clear F prod_ss options x 

% Export to csv 
values = round([cons_ss; kapi_ss; labo_ss], 4);

row_names = {'Consumption'; 'Capital'; 'Labor'};
  
vv = [row_names, num2cell(values)];
vvv = cell2table(vv, 'VariableNames', {'Variable', 'Value'});
  
writetable(vvv,'econ714_homework2_question1_output.csv')
clear vv vvv row_names values 

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%%
% -------------------------------------------------------------------------
% Question 2 - Value Function Iteration with a Fixed Grid
% -------------------------------------------------------------------------

% Define some things 
K_grid = linspace(0.7 * kapi_ss , 1.3 * kapi_ss, 251);
transition_matrix = [0.97, 0.03, 0; 0.01, 0.98, 0.01; 0, 0.03, 0.97]; 
V_prev = zeros(length(K_grid), 3); 
V = zeros(length(K_grid), 3);
V_aux = zeros(length(K_grid), 1); 
L = zeros(length(K_grid), 3); 
K = zeros(length(K_grid), 3); 
n_it = 0;
z = [-0.05, 0, 0.05]; 
sup_norm = 1; 

% Specify an initial guess 
    % Suppose labor equals its steady state level 
    % Do a grid search in this reduced problem 
    % V_prev will be our initial guess 
    % Do at most 1000 iterations, this is just for the initial guess 
    % It turns out it is not a very good guess. This can be improved 

    % While for the initial guess
    while sup_norm > 10^(-7) && n_it < 1000
        
        n_it = n_it + 1; 

        % Loop over technology states 
        for zz = 1:3
    
            z_value = z(zz); 
        
            % Loop over the current state of capital
            for ii = 1:length(K_grid)
                
                % Loop over capital in next period
                for jj = 1:length(K_grid)
        
                    cons = exp(z_value) * K_grid(ii)^(alpha) * labo_ss^(1-alpha) + ...
                        (1-delta) * K_grid(ii) - K_grid(jj);
        
                    if cons < 0
        
                       V_aux(jj,1) = -10^10;
        
                    else
        
                       V_aux(jj,1) = (1-beta) * (log(cons) - (labo_ss)^2/2) + ...
                           beta * V_prev(ii, :) * transpose(transition_matrix(zz, :)) ; 
        
                    end    
                       
                end
        
                [V(ii, zz) , ~] = max(V_aux);
                
            end
    
        end
    
        % Update the function and sup norm
        sup_norm = max(max(abs(V_prev - V)));
        V_prev = V; 

    end

% Now we start with value function iteration using this initial guess
V_initial_guess = V;     
sup_norm = 1;
sup_norm_ind = [1,1,1]; 
n_it = 0; 
clear cons ii jj zz z_value 
tol_level = 10^(-7);
tic
   
    % While of the algorithm 
    while sup_norm > tol_level && n_it < 10000


        % Loop over the current state of capital 
        for ii = 1:length(K_grid)

            % Loop over technology states 
            for zz = 1:3
                
                % Only do what follows if necessary! 
                if sup_norm_ind(zz) > tol_level
    
                    z_value = z(zz);
    
                    % Loop over capital next period
                    for jj = 1:length(K_grid)
    
                        % First compute l optimally given k and k_next
                        % The goal here is to get a grid for l
                        % For this use the intratemporal condition 
                        if intratemporal_condition(z_value, K_grid(ii), K_grid(jj), alpha, delta, 1) >= 0
                            
                            % Here I check that labor is not greater than 1, this
                            % is the intratemporal condition evaluated at l=1. If
                            % this number is positive, that means I am going to
                            % work all day
                            L(jj,zz) = 1; 
                
                        else
                
                            % Find the zero of the intratemporal condition 
                            L(jj, zz) = min(1,max(0, fzero(@(x) intratemporal_condition(z_value, K_grid(ii), ... 
                            K_grid(jj), alpha, delta, x), [10^(-7),1]))); 
                
                        end 
    
                    end
    
                    % Here I need to compute the value function choosing k' and
                    % doing the linear interpolation
                    K(ii,zz) = fminbnd(@(x) interpolated_value_function(x, K_grid, ...
                        z, V_prev, K_grid(ii), z(zz), alpha, L, delta, transition_matrix, ...
                        beta), min(K_grid), max(K_grid)); 
                    V(ii, zz) = - interpolated_value_function(K(ii,zz), K_grid, ...
                        z, V_prev, K_grid(ii), z(zz), alpha, L, delta, transition_matrix, ...
                        beta);             
    
                end
            
            end
        
        end

        % Update the function and sup norm
        sup_norm = max(max(abs(V-V_prev)));
        sup_norm_ind = max(abs(V-V_prev));
        n_it = n_it + 1;
        V_prev = V; 

    end

time_VFI_fixed_grid = toc; 

% Store results 
values = [transpose(K_grid), V, K, L, ones(length(K_grid),1) * n_it, ...
    ones(length(K_grid),1) * time_VFI_fixed_grid, ones(length(K_grid),1) * ...
    tol_level, V_initial_guess];

vvv = array2table(values, 'VariableNames', {'Kgrid', 'V1', 'V2', 'V3', 'K1', ...
    'K2', 'K3', 'L1', 'L2', 'L3', 'n_iter', 'time', 'tol_level', 'V1InitialGuess', ...
    'V2InitialGuess', 'V3InitialGuess'});
  
writetable(vvv,'econ714_homework2_question2_output.csv')
clear values vvv 
clearvars -except alpha beta cons_ss delta kapi_ss labo_ss valu_ss 

%%
% -------------------------------------------------------------------------
% Question 3 - Accelerator
% -------------------------------------------------------------------------

% See fixed tab to do this properly another day 
% Define some things 
K_grid = linspace(0.7 * kapi_ss , 1.3 * kapi_ss, 51);
transition_matrix = [0.97, 0.03, 0; 0.01, 0.98, 0.01; 0, 0.03, 0.97]; 
V_prev = zeros(length(K_grid), 3); 
V = zeros(length(K_grid), 3);
V_aux = zeros(length(K_grid), 1); 
cons =  zeros(length(K_grid), 3);
L = zeros(length(K_grid), 3); 
K = zeros(length(K_grid), 3); 
n_it = 0;
n_it_acc = 0; % To see when the accelerator starts 
z = [-0.05, 0, 0.05]; 
sup_norm = 1; 

% Specify an initial guess 
    % Suppose labor equals its steady state level 
    % Do a grid search in this reduced problem 
    % V_prev will be our initial guess 
    % Do at most 1000 iterations, this is just for the initial guess 
    % Same initial guess as before 

    % While for the initial guess
    while sup_norm > 10^(-7) && n_it < 1000
        
        n_it = n_it + 1; 

        % Loop over technology states 
        for zz = 1:3
    
            z_value = z(zz); 
        
            % Loop over the current state of capital
            for ii = 1:length(K_grid)
                
                % Loop over capital in next period
                for jj = 1:length(K_grid)
        
                    cons = exp(z_value) * K_grid(ii)^(alpha) * labo_ss^(1-alpha) + ...
                        (1-delta) * K_grid(ii) - K_grid(jj);
        
                    if cons < 0
        
                       V_aux(jj,1) = -10^10;
        
                    else
        
                       V_aux(jj,1) = (1-beta) * (log(cons) - (labo_ss)^2/2) + ...
                           beta * V_prev(ii, :) * transpose(transition_matrix(zz, :)) ; 
        
                    end    
                       
                end
        
                [V(ii, zz) , ~] = max(V_aux);
                
            end
    
        end
    
        % Update the function and sup norm
        sup_norm = max(max(abs(V_prev - V)));
        V_prev = V; 

    end

% Start the value function iteration using that initial guess
V_initial_guess = V;
sup_norm = 1;
sup_norm_ind = [1,1,1]; 
n_it = 0; 
clear cons ii jj zz z_value 
tol_level = 10^(-7);
tic
   
% Now start the value function iteration 

    % While of the algorithm 
    while sup_norm > tol_level && n_it < 10000

        if sup_norm > 10^2
            
            % Just in case everything explodes 
            disp('Value functions are not close to each other')
            stop

        end

        % If this is the case, I want to optimize 
        if n_it / 10 == floor(n_it/10) || n_it < 50

        % Loop over the current state of capital 
        for ii = 1:length(K_grid)

            % Loop over technology states 
            for zz = 1:3
                
                z_value = z(zz);

                % Only do what follows if necessary! 
                if sup_norm_ind(zz) > tol_level
                      
                    

                        % Loop over capital next period
                            % This is to get a grid of the optimal L 
                        for jj = 1:length(K_grid)
        
                            % First compute l optimally given k and k_next
                            % The goal here is to get a grid for l
                            % For this use the intratemporal condition 
                            if intratemporal_condition(z_value, K_grid(ii), K_grid(jj), alpha, delta, 1) >= 0
                                
                                % Here I check that labor is not greater than 1, this
                                % is the intratemporal condition evaluated at l=1. If
                                % this number is positive, that means I am going to
                                % work all day
                                L(jj,zz) = 1; 
                    
                            else
                    
                                % Find the zero of the intratemporal condition 
                                L(jj, zz) = min(1,max(0, fzero(@(x) ...
                                    intratemporal_condition(z_value, K_grid(ii), ... 
                                    K_grid(jj), alpha, delta, x), [10^(-7),1]))); 
                    
                            end 
        
                        end
                                       
                        % Here I need to compute the value function choosing k' and
                        % doing the linear interpolation
                        K(ii,zz) = fminbnd(@(x) interpolated_value_function(x, K_grid, ...
                            z, V_prev, K_grid(ii), z(zz), alpha, L, delta, transition_matrix, ...
                            beta), min(K_grid), max(K_grid)); 
                        V(ii, zz) = - interpolated_value_function(K(ii,zz), K_grid, ...
                            z, V_prev, K_grid(ii), z(zz), alpha, L, delta, transition_matrix, ...
                            beta); 

                        
                    
                   


                        

                end

                
     
            end
            

         
        end

         % Every other iteration
        else

            if n_it_acc == 0

                n_it_acc = n_it; 

            end
            
            % Assign to the policy function the number that's closer to
            % the grid 
            K_pol = K; 
            
            for zz = 1:length(z)
            for ii = 1:length(K_grid)
                
                [~, index] = min(abs(K(ii,zz)-K_grid)); 
                
                K_pol(ii,zz) = K_grid(index); 

            end
            end 

            % Create a matrix that is zero everywhere except when the
            % policy function matches the grid 
            Q = zeros(length(K_grid), length(K_grid));

            % row 1 
            for ii = 1:length(K_grid)
            for jj = 1:length(K_grid)

                if K_grid(ii) == K_pol(jj)

                    Q(ii,jj) = 1; 

                end

            end
            end

             
%%
            % need utility 
            output_1 = exp(z(1)) * transpose(K_grid).^(alpha) .* L(:,1).^(1-alpha);
            cons_1 = max(10^(-7), output_1 + (1-delta) * transpose(K_grid) - K(:,1));
            util_1 = (1-beta) * (log(cons_1) - L(:,1).^2/2);
            output_2 = exp(z(2)) * transpose(K_grid).^(alpha) .* L(:,2).^(1-alpha);
            cons_2 = max(10^(-7), output_2 + (1-delta) * transpose(K_grid) - K(:,2));
            util_2 = (1-beta) * (log(cons_2) - L(:,2).^2/2); 
            output_3 = exp(z(3)) * transpose(K_grid).^(alpha) .* L(:,3).^(1-alpha);
            cons_3 = max(10^(-7), output_3 + (1-delta) * transpose(K_grid) - K(:,3));
            util_3 = (1-beta) * (log(cons_3) - L(:,3).^2/2); 

            V_1_temp = zeros(length(K_grid) * length(z),1); 
            V_temporal = zeros(length(K_grid) * length(z),1); 

            % Generate a vector of all the value functions
            VV = [V(:,1); V(:,2); V(:,3)];
            eVV_temp = transition_matrix * transpose(V);

        
            
            eVV = [transpose(eVV_temp(1,:)); transpose(eVV_temp(2,:)); ...
                transpose(eVV_temp(3,:))];
            clear eVV_temp

            % Compute the function that I wnat to solve just to try 
            fun_temp = VV - [util_1; util_2; util_3] - beta * ...
                kron(ones(length(z),length(z)), Q) * eVV;

            clear fun_temp VV eVV_temp 

%%
            % I think this is the way to go, but it is too hard for me :( 
            %xxx = fsolve(@(x) x - [util_1; util_2; util_3] - beta * ...
               %kron(ones(length(z),length(z)), Q) * ); 
            

            stop 

            % Here I want to code the accelerator
            y = exp(z_value) * K_grid(ii)^(alpha) * L(ii,zz); 
            cons(ii,zz) = max(10^(-7), y + (1-delta) * K_grid(ii) - K(ii,zz)); 
            V_temp = V; 
            v_temp_1 = interpolated_value_function(K(ii,zz), K_grid, z, ...
                V_prev, K_grid(ii), -0.05, alpha, L, delta, transition_matrix, ...
                beta); 
            v_temp_2 = interpolated_value_function(K(ii,zz), K_grid, z, ...
                V_prev, K_grid(ii), 0, alpha, L, delta, transition_matrix, ...
                beta);
            v_temp_3 = interpolated_value_function(K(ii,zz), K_grid, z, ...
                V_prev, K_grid(ii), 0.05, alpha, L, delta, transition_matrix, ...
                beta);

            stop 
            v_bar_next = [v_temp_1, v_temp_2, v_temp_3] * ...
                transpose(transition_matrix(zz,:));
            V(ii,zz) = (1-beta) * (log(cons) - (L(ii,zz)^2)/2) + ...
                beta * V_temp(ii,:) * transpose(transition_matrix(zz,:));
            
    

                         

        end

        disp(V((length(K_grid)+1)/2, 2))

        % Update the function and sup norm
        sup_norm = max(max(abs(V-V_prev)));
        sup_norm_ind = max(abs(V-V_prev));
        n_it = n_it + 1;
        V_prev = V; 


         
     
    end

time_VFI_fixed_grid = toc; 

% Store results 
values = [transpose(K_grid), V, K, L, ones(length(K_grid),1) * n_it, ...
    ones(length(K_grid),1) * time_VFI_fixed_grid, ones(length(K_grid),1) * ...
    tol_level, V_initial_guess, ones(length(K_grid),1) * n_it_acc];

vvv = array2table(values, 'VariableNames', {'Kgrid', 'V1', 'V2', 'V3', 'K1', ...
    'K2', 'K3', 'L1', 'L2', 'L3', 'n_iter', 'time', 'tol_level', 'V1InitialGuess', ...
    'V2InitialGuess', 'V3InitialGuess', 'n_iter_acc'});
  
writetable(vvv,'econ714_homework2_question3_output.csv')
clear values vvv 
clearvars -except alpha beta cons_ss delta kapi_ss labo_ss valu_ss 

%%
% -------------------------------------------------------------------------
% Question 4 - Multigrid 
% -------------------------------------------------------------------------

% Grid 1 has 101 points
% Define some vectors and matrices to store things 
K_grid_1 = transpose(linspace(0.7 * kapi_ss , 1.3 * kapi_ss, 100 + 1));
transition_matrix = [0.97, 0.03, 0; 0.01, 0.98, 0.01; 0, 0.03, 0.97]; 
V_prev = zeros(length(K_grid_1), 3); 
V = zeros(length(K_grid_1), 3);
V_aux = zeros(length(K_grid_1), 1); 
L = zeros(length(K_grid_1), 3); 
K = zeros(length(K_grid_1), 3); 
n_it = 0;
z = [-0.05, 0, 0.05]; 
sup_norm = 1; 

% Specify an initial guess 
    % Suppose labor equals its steady state level 
    % Do a grid search in this reduced problem 
    % V_prev will be our initial guess 
    % Do at most 1000 iterations, this is just for the initial guess 
    % Same initial guess as before 

    % While for the initial guess
    while sup_norm > 10^(-7) && n_it < 1000
        
        n_it = n_it + 1; 

        % Loop over technology states 
        for zz = 1:3
    
            z_value = z(zz); 
        
            % Loop over the current state of capital
            for ii = 1:length(K_grid_1)
                
                % Loop over capital in next period
                for jj = 1:length(K_grid_1)
        
                    cons = exp(z_value) * K_grid_1(ii)^(alpha) * labo_ss^(1-alpha) + ...
                        (1-delta) * K_grid_1(ii) - K_grid_1(jj);
        
                    if cons < 0
        
                       V_aux(jj,1) = -10^10;
        
                    else
        
                       V_aux(jj,1) = (1-beta) * (log(cons) - (labo_ss)^2/2) + ...
                           beta * V_prev(ii, :) * transpose(transition_matrix(zz, :)) ; 
        
                    end    
                       
                end
        
                [V(ii, zz) , ~] = max(V_aux);
                
            end
    
        end
    
        % Update the function and sup norm
        sup_norm = max(max(abs(V_prev - V)));
        V_prev = V; 

    end

% Value function iteration in the first grid 
V_initial_guess = V;
sup_norm = 1;
sup_norm_ind = [1,1,1]; 
n_it = 0; 
clear cons ii jj zz z_value 
tol_level = 10^(-7);
tic; 

    while sup_norm > tol_level && n_it < 10000

        % Loop over the current state of technology 
        for zz = 1:3
            
            if sup_norm_ind(zz) > tol_level
            
                z_value = z(zz);
    
                % Loop over the current state of capital 
                for ii = 1:length(K_grid_1)
    
                    % Loop over the next state of capital 
    
                    for jj = 1:length(K_grid_1)
    
                        % First compute l optimally given k and k_next
                        % The goal here is to get a grid for l
                        % For this use the intratemporal condition 
                        if intratemporal_condition(z_value, K_grid_1(ii), K_grid_1(jj), ...
                                alpha, delta, 1) >= 0
                            
                            % Here I check that labor is not greater than 1, this
                            % is the intratemporal condition evaluated at l=1. If
                            % this number is positive, that means I am going to
                            % work all day
                            L(jj,zz) = 1; 
                
                        else
                
                            % Find the zero of the intratemporal condition 
                            L(jj, zz) = min(1,max(0, fzero(@(x) intratemporal_condition(z_value, K_grid_1(ii), ... 
                            K_grid_1(jj), alpha, delta, x), [10^(-7),1]))); 
                
                        end 
                             
                    end

                    % Here I need to compute the value function choosing k' and
                    % doing the linear interpolation
                    K(ii,zz) = fminbnd(@(x) interpolated_value_function(x, K_grid_1, ...
                        z, V_prev, K_grid_1(ii), z(zz), alpha, L, delta, transition_matrix, ...
                        beta), min(K_grid_1), max(K_grid_1)); 
                    V(ii, zz) = - interpolated_value_function(K(ii,zz), K_grid_1, ...
                        z, V_prev, K_grid_1(ii), z(zz), alpha, L, delta, transition_matrix, ...
                        beta);   
                      
                end
            
            end

        end
        
        % Update things
        n_it = n_it + 1;
        sup_norm = max(max(abs(V-V_prev)));
        sup_norm_ind = max(abs(V-V_prev));
        V_prev = V; 
        
    end

 
% Value function iteration on the second grid
K_grid_2 = transpose(linspace(0.7 * kapi_ss , 1.3 * kapi_ss, 1000 + 1));
V_prev_2 = zeros(length(K_grid_2), 3); 
V_2 = zeros(length(K_grid_2), 3);
L_2 = zeros(length(K_grid_2), 3); 
K_2 = zeros(length(K_grid_2), 3); 
n_it_2 = 0;
sup_norm = 1; 
sup_norm_ind = [1,1,1];

% Specify the initial guess interpolating what came from the previos grid
V_initial_guess_2 = zeros(length(K_grid_2),3);

% For each element on the new grid
for ii = 1:length(K_grid_2)
    
    K_temp = K_grid_1; 
    [~, i] = min(abs(K_grid_2(ii,1) - K_grid_1)); 
    K_temp(i) = 10^7;
    [~, j] = min(abs(K_grid_2(ii,1) - K_temp));
    clear K_temp
    i = min(i,j);
    j = i+1;

    % The new element is between i and j of the previous grid 

    % For each technology shock, complete V_initial_guess_2
    for zz = 1:3

        V_initial_guess_2(ii, zz) = V(i,zz) + ((V(j,zz) - V(i,zz)) ...
            / (K_grid_1(j) - K_grid_1(i))) * (K_grid_2(ii) - K_grid_1(i)); 

    end
     
end

% Done with the initial guess
V_prev = V_initial_guess_2; 

% Value function iteration in the second grid 
    while sup_norm > tol_level && n_it_2 < 10000

        % Loop over the current state of technology 
        for zz = 1:3
            
            if sup_norm_ind(zz) > tol_level
            
                z_value = z(zz);
    
                % Loop over the current state of capital 
                for ii = 1:length(K_grid_2)
    
                    % Loop over the next state of capital 
                    for jj = 1:length(K_grid_2)
    
                        % First compute l optimally given k and k_next
                        % The goal here is to get a grid for l
                        % For this use the intratemporal condition 
                        if intratemporal_condition(z_value, K_grid_2(ii), K_grid_2(jj), ...
                                alpha, delta, 1) >= 0
                            
                            % Here I check that labor is not greater than 1, this
                            % is the intratemporal condition evaluated at l=1. If
                            % this number is positive, that means I am going to
                            % work all day
                            L_2(jj,zz) = 1; 
                
                        else
                
                            % Find the zero of the intratemporal condition 
                            L_2(jj, zz) = min(1,max(0, fzero(@(x) intratemporal_condition(z_value, K_grid_2(ii), ... 
                            K_grid_2(jj), alpha, delta, x), [10^(-7),1]))); 
                
                        end 
      
                    end
                    
                    % Here I need to compute the value function choosing k' and
                    % doing the linear interpolation
                    K_2(ii,zz) = fminbnd(@(x) interpolated_value_function(x, K_grid_2, ...
                        z, V_prev, K_grid_2(ii), z(zz), alpha, L_2, delta, transition_matrix, ...
                        beta), min(K_grid_2), max(K_grid_2)); 
                    V_2(ii, zz) = - interpolated_value_function(K_2(ii,zz), K_grid_2, ...
                        z, V_prev, K_grid_2(ii), z(zz), alpha, L_2, delta, transition_matrix, ...
                        beta);   
    
                end
            
            end

        end
        
        % Update things
        n_it_2 = n_it_2 + 1;
        sup_norm = max(max(abs(V_2-V_prev)));
        sup_norm_ind = max(abs(V_2-V_prev));
        V_prev = V_2; 
        
    end

% Value function iteration in the third grid 
K_grid_3 = transpose(linspace(0.7 * kapi_ss , 1.3 * kapi_ss, 10000 + 1));
V_3 = zeros(length(K_grid_3), 3);
V_aux = zeros(length(K_grid_3), 1); 
L_3 = zeros(length(K_grid_3), 3); 
K_3 = zeros(length(K_grid_3), 3); 
n_it_3 = 0;
sup_norm = 1; 
sup_norm_ind = [1,1,1];

% Specify the initial guess interpolating what came from the previos grid
V_initial_guess_3 = zeros(length(K_grid_3),3);

% For each element on the new grid
for ii = 1:length(K_grid_3)
    
    K_temp = K_grid_2; 
    [~, i] = min(abs(K_grid_3(ii,1) - K_grid_2)); 
    K_temp(i) = 10^7;
    [~, j] = min(abs(K_grid_3(ii,1) - K_temp));
    clear K_temp
    i = min(i,j);
    j = i+1;

    % The new element is between i and j of the previous grid 

    % For each technology shock, complete V_initial_guess_3
    for zz = 1:3

        V_initial_guess_3(ii, zz) = V_2(i,zz) + ((V_2(j,zz) - V_2(i,zz)) ...
            / (K_grid_2(j) - K_grid_2(i))) * (K_grid_3(ii) - K_grid_2(i)); 

    end
     
end

% Done with the initial guess
V_prev = V_initial_guess_3; 

% Nos start the iteration 
    while sup_norm > tol_level && n_it_3 < 10000

        % Loop over the current state of technology 
        for zz = 1:3
            
            if sup_norm_ind(zz) > tol_level
            
                z_value = z(zz);
    
                % Loop over the current state of capital 
                for ii = 1:length(K_grid_3)
    
                    % Loop over the next state of capital 
    
                    for jj = 1:length(K_grid_3)
    
                        % First compute l optimally given k and k_next
                        % The goal here is to get a grid for l
                        % For this use the intratemporal condition 
                        if intratemporal_condition(z_value, K_grid_3(ii), K_grid_3(jj), ...
                                alpha, delta, 1) >= 0
                            
                            % Here I check that labor is not greater than 1, this
                            % is the intratemporal condition evaluated at l=1. If
                            % this number is positive, that means I am going to
                            % work all day
                            L_3(jj,zz) = 1; 
                
                        else
                
                            % Find the zero of the intratemporal condition 
                            L_3(jj, zz) = min(1,max(0, fzero(@(x) intratemporal_condition(z_value, K_grid_3(ii), ... 
                            K_grid_3(jj), alpha, delta, x), [10^(-7),1]))); 
                
                        end 
       
                    end
                    
                    % Here I need to compute the value function choosing k' and
                    % doing the linear interpolation
                    K_3(ii,zz) = fminbnd(@(x) interpolated_value_function(x, K_grid_3, ...
                        z, V_prev, K_grid_3(ii), z(zz), alpha, L_3, delta, transition_matrix, ...
                        beta), min(K_grid_3), max(K_grid_3)); 
                    V_3(ii, zz) = - interpolated_value_function(K_3(ii,zz), K_grid_3, ...
                        z, V_prev, K_grid_3(ii), z(zz), alpha, L_3, delta, transition_matrix, ...
                        beta);   
    
                end
            
            end

        end
        
        % Update things
        n_it_3 = n_it_3 + 1;
        sup_norm = max(max(abs(V_3-V_prev)));
        sup_norm_ind = max(abs(V_3-V_prev));
        V_prev = V_3; 
        % disp(n_it_3)
        
    end

% Store things 
time_multigrid = toc; % in seconds 

values = [K_grid_3, V_3, K_3, L_3, ones(length(K_grid_3),1) * n_it_3, ones(length(K_grid_3),1) ...
    * n_it_2, ones(length(K_grid_3),1) * n_it, ...
    ones(length(K_grid_3),1) * time_multigrid, ones(length(K_grid_3),1) * ...
    tol_level, V_initial_guess_3];

vvv = array2table(values, 'VariableNames', {'Kgrid', 'V1', 'V2', 'V3', 'K1', ...
    'K2', 'K3', 'L1', 'L2', 'L3', 'n_iter_3', 'n_iter_2', 'n_iter_1', 'time', 'tol_level', 'V1InitialGuess', ...
    'V2InitialGuess', 'V3InitialGuess'});
  
writetable(vvv,'econ714_homework2_question4_output.csv')
clear values vvv 
clearvars -except alpha beta cons_ss delta kapi_ss labo_ss valu_ss 

%% 
% -------------------------------------------------------------------------
% Plots
% -------------------------------------------------------------------------

% Import data
    % Set up the Import Options and import the data
    opts = delimitedTextImportOptions("NumVariables", 13);
    
    % Specify range and delimiter
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";
    
    % Specify column names and types
    opts.VariableNames = ["Kgrid", "V1", "V2", "V3", "K1", "K2", "K3", "L1", "L2", "L3", "n_iter", "time", "tol_level"];
    opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
    
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    
    % Import the data
    output_q2 = readtable("C:\Users\Javier\Dropbox\Aplicaciones\Overleaf\Econ714-PS2\econ714_homework2_question2_output.csv", opts);
    
    
    % Clear temporary variables
    clear opts

% Suppose we are in the ss, suddenly z=0.05 and after this z=0 forever
    
    % I can plot capital, labor, ouput and prices 
    T_max = 10+1;
    t = transpose(1:1:T_max); 
    k_trajectory = zeros(T_max, 1);
    l_trajectory = zeros(T_max, 1);

    % Extrat time and number of iter
    time_VFI_fixed_grid = output_q2 {1, 'time'};
    time_VFI_fixed_grid = round(time_VFI_fixed_grid / 60, 0);
    n_it_VFI_fixed_grid = output_q2 {1, 'n_iter'};
    tole_VFI_fixed_grid = output_q2 {1, 'tol_level'};

    % Extract the policy function when there is a positive shock
    K_grid = output_q2 {:, 'Kgrid'};
    policy_fun_pos = output_q2 {:, 'K3'}; 
    policy_fun_pos = [K_grid, policy_fun_pos];
    value_fun_pos = output_q2 {:, 'V3'};
    

    % First value of the trajectory 
    [~, i] = min(abs(kapi_ss - policy_fun_pos(:,1))); 
    k_trajectory(1,1) = policy_fun_pos(i,2);  
    clear m i 

    % Extract the policy function when z=0 
    policy_fun_neu = output_q2 {:, 'K2'};
    policy_fun_neu = [K_grid, policy_fun_neu];
    value_fun_neu = output_q2 {:, 'V2'};

    % Extract the policy function when z=-0.05 
    policy_fun_neg = output_q2 {:, 'K1'};
    policy_fun_neg = [K_grid, policy_fun_neg];
    value_fun_neg = output_q2 {:, 'V1'};
    
    % Fill the trajectory 
    for tt = 2:T_max 
        
        diff = abs(k_trajectory(tt-1,1) - policy_fun_neu(:,1));
        [~, i] = min(diff);
        diff(i) = 1000; 
        [m2, j] = min(diff);
        i = min(i,j);
        j = i+1; 
                 
        k_trajectory(tt,1) = policy_fun_neu(i, 2) + (policy_fun_neu(i,2) - policy_fun_neu(j,2)) ...
            / (policy_fun_neu(i,1) - policy_fun_neu(j,1)) * ...
            (k_trajectory(tt-1,1) - policy_fun_neu(i,1));

        if tt > 3 

            l_trajectory(tt-1,1) = min(1,max(0, fzero(@(x) intratemporal_condition(0, k_trajectory(tt-1,1), k_trajectory(tt,1), ...
                alpha, delta, x), [10^(-7),1])));

        else

            l_trajectory(tt-1,1) = min(1,max(0, fzero(@(x) intratemporal_condition(0.05, k_trajectory(tt-1,1), k_trajectory(tt,1), ...
                alpha, delta, x), [10^(-7),1])));

        end

   end

   clear j i tt

   % Production and consumption 
   c_trajectory = zeros(T_max, 1);

   y_trajectory = k_trajectory(1:T_max,:).^(alpha) .* l_trajectory(1:T_max,:).^(1-alpha);  
   y_trajectory(1,1) = exp(0.05) * k_trajectory(1,1)^(alpha) * l_trajectory(1,1)^(1-alpha); 
   c_trajectory(1:T_max-1,1) = max(0,y_trajectory(1:T_max-1,1) + (1-delta) * ...
       k_trajectory(1:T_max-1) - k_trajectory(2:T_max,1)); 

   k_trajectory = [kapi_ss; k_trajectory];
   l_trajectory = [labo_ss; l_trajectory];
   y_trajectory = [kapi_ss^(alpha) * labo_ss^(1-alpha); y_trajectory];
   c_trajectory = [cons_ss; c_trajectory]; 
    
   % plot(k_trajectory(1:T_max-1,:))
   % plot(l_trajectory(1:T_max-1,:))
   % plot(y_trajectory(1:T_max-1,:))
   % plot(c_trajectory(1:T_max-1,:))

   % Export figures 
   
   % Capital 
   figure(1)
   plot(t(1:T_max-1)-1, k_trajectory(1:T_max-1,:),'linewidth', 2)
   hold on
   yline(kapi_ss, 'linewidth', 1)
   xlabel('Time')
   ylabel('Capital')
   title('Capital over time')
   hold off
   saveas(gcf,'econ714_homework2_question2_plot_capital.png');
   close(figure(1))

   % Labor 
   figure(1)
   plot(t(1:T_max-1)-1, l_trajectory(1:T_max-1,:),'linewidth', 2)
   hold on
   yline(labo_ss, 'linewidth', 1)
   xlabel('Time')
   ylabel('Labor')
   title('Labor over time')
   hold off
   saveas(gcf,'econ714_homework2_question2_plot_labor.png');
   close(figure(1))

   % Output 
   figure(1)
   plot(t(1:T_max-1)-1, y_trajectory(1:T_max-1,:),'linewidth', 2)
   hold on
   yline(kapi_ss^(alpha) * labo_ss^(1-alpha), 'linewidth', 1)
   xlabel('Time')
   ylabel('Output')
   title('Output over time')
   hold off
   saveas(gcf,'econ714_homework2_question2_plot_output.png');
   close(figure(1))

   % Consumption 
   figure(1)
   plot(t(1:T_max-1)-1, c_trajectory(1:T_max-1,:),'linewidth', 2)
   hold on
   yline(cons_ss, 'linewidth', 1)
   xlabel('Time')
   ylabel('Consumption')
   title('Consumption over time')
   hold off
   saveas(gcf,'econ714_homework2_question2_plot_consumption.png');
   close(figure(1))

    
   % Policy function 
   figure(1)
   plot(K_grid, policy_fun_neg(:,2), 'linewidth', 2)
   hold on
   plot(K_grid, policy_fun_neu(:,2), 'linewidth', 2)
   plot(K_grid, policy_fun_pos(:,2), 'linewidth', 2)
   plot(K_grid, K_grid, 'linewidth', 1, 'color', 'black')
   xlabel('Capital')
   ylabel('Capital next period')
   title('Policy Functions')
   subtitle(['Iterations: ',num2str(n_it_VFI_fixed_grid), ...
      ' - Time (minutes): ', num2str(time_VFI_fixed_grid), ...
      ' - Tolerance level: ', num2str(tole_VFI_fixed_grid)])
   hold off
   legend({'z=-0.05','z=0','z=0.05'},'Location', 'southeast')
   saveas(gcf,'econ714_homework2_question2_plot_policy_function.png');
   close(figure(1))

   % Value function 
   figure(1)
   plot(K_grid, value_fun_neg(:,1), 'linewidth', 2)
   hold on
   plot(K_grid, value_fun_neu(:,1), 'linewidth', 2)
   plot(K_grid, value_fun_pos(:,1), 'linewidth', 2)
   xlabel('Capital')
   ylabel('Value')
   title('Value Functions')
   subtitle(['Iterations: ',num2str(n_it_VFI_fixed_grid), ...
      ' - Time (minutes): ', num2str(time_VFI_fixed_grid), ...
      ' - Tolerance level: ', num2str(tole_VFI_fixed_grid)])
   hold off
   legend({'z=-0.05','z=0','z=0.05'},'Location', 'southeast')
   saveas(gcf,'econ714_homework2_question2_plot_value_function.png');
   close(figure(1))

   old_K_grid = K_grid;

clearvars -except kapi_ss alpha beta cons_ss delta labo_ss valu_ss value_fun_neg ...
    value_fun_neu value_fun_pos old_K_grid

% Import data for question 4
    % Import data from text file
    % Script for importing data from the following text file:
    %
    %    filename: C:\Users\Javier\Dropbox\Aplicaciones\Overleaf\Econ714-PS2\econ714_homework2_question4_output.csv
    %
    % Auto-generated by MATLAB on 08-Dec-2021 17:50:56
    
    % Set up the Import Options and import the data
    opts = delimitedTextImportOptions("NumVariables", 16);
    
    % Specify range and delimiter
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";
    
    % Specify column names and types
    opts.VariableNames = ["Kgrid", "V1", "V2", "V3", "K1", "K2", "K3", "L1", "L2", "L3", "n_iter", "time", "tol_level", "V1InitialGuess", "V2InitialGuess", "V3InitialGuess"];
    opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
    
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    
    % Import the data
    output_q4 = readtable("C:\Users\Javier\Dropbox\Aplicaciones\Overleaf\Econ714-PS2\econ714_homework2_question4_output.csv", opts);
    
    
    % Clear temporary variables
    clear opts

    % Plot value and policy function 
    K_grid = output_q4 {:, 'Kgrid'};
    V_fun_1 = output_q4 {:, 'V1'};
    V_fun_2 = output_q4 {:, 'V2'};
    V_fun_3 = output_q4 {:, 'V3'};
    pol_fun_1 = output_q4 {:, 'K1'};
    pol_fun_2 = output_q4 {:, 'K2'};
    pol_fun_3 = output_q4 {:, 'K3'};
    time_min = output_q4 {1, 'time'};
    tole_le = output_q4 {1, 'tol_level'}; 

   % Policy function 
   figure(2)
   plot(K_grid, pol_fun_1, 'linewidth', 2)
   hold on
   plot(K_grid, pol_fun_2, 'linewidth', 2)
   plot(K_grid, pol_fun_3, 'linewidth', 2)
   plot(K_grid, K_grid, 'linewidth', 1, 'color', 'black')
   xlabel('Capital')
   ylabel('Capital next period')
   title('Policy Functions')
   subtitle(['Time (minutes): ', num2str(time_min), ...
      ' - Tolerance level: ', num2str(tole_le)])
   hold off
   legend({'z=-0.05','z=0','z=0.05'},'Location', 'southeast')
   saveas(gcf,'econ714_homework2_question4_plot_policy_function.png');
   close(figure(2))

   % Value function 
   figure(3)
   plot(K_grid, V_fun_1, 'linewidth', 2)
   hold on
   plot(K_grid, V_fun_2, 'linewidth', 2)
   plot(K_grid, V_fun_3, 'linewidth', 2)
   xlabel('Capital')
   ylabel('Value')
   title('Value Functions')
   subtitle(['Time (minutes): ', num2str(time_min), ...
      ' - Tolerance level: ', num2str(tole_le)])
   hold off
   legend({'z=-0.05','z=0','z=0.05'},'Location', 'southeast')
   saveas(gcf,'econ714_homework2_question4_plot_value_function.png');
   close(figure(3))

 
int_value_fun_q2 = zeros(length(K_grid),3); 
% K_grid_small = value_fun_neg(:,1); 
value_fun_q2 = [value_fun_neg, value_fun_neu, value_fun_pos];
clear value_fun_neg value_fun_neu value_fun_pos

% Rename some things
grid_k_complete = K_grid;
clear K_grid
K_grid = old_K_grid;
clear K_grid_small

 

for zz = 1:3

    for kk = 1:length(grid_k_complete)
    
        [~,i] = min(abs(grid_k_complete(kk) - K_grid));
        K_grid_temp = K_grid; 
        K_grid_temp(i) = 1000000;
        [~,j] = min(abs(grid_k_complete(kk) - K_grid_temp));
        
        i = min(i,j);
        int_value_fun_q2(kk,zz) = value_fun_q2(i,zz) + (value_fun_q2(i+1,zz) - value_fun_q2(i,zz)) / ...
            (K_grid(i+1) - K_grid(i)) * (grid_k_complete(kk) - K_grid(i));
              

    end

end

 

% Comparison with question 2
figure(4)
plot(grid_k_complete,abs(int_value_fun_q2(:,1) - V_fun_1), 'linewidth', 1)
hold on
plot(grid_k_complete,abs(int_value_fun_q2(:,2) - V_fun_2), 'linewidth', 1)
plot(grid_k_complete,abs(int_value_fun_q2(:,3) - V_fun_3), 'linewidth', 1)
title('Comparison with simple VFI')
% subtitle(['Time (seconds): ', num2str(time_cheby), ...
      % ' - Tolerance level: ', num2str(tol_level)])
xlabel('Capital')
ylabel('Difference (Absolute Value)')
legend({'z=-0.05','z=0','z=0.05'},'Location', 'southeast')
saveas(gcf,'econ714_homework2_question4_plot_comparison_q2.png');
close(figure(4))




 





  