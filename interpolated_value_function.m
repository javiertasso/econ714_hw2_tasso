function [vvv] = interpolated_value_function(k_next, K_grid, z_vector, ...
    V_prev, k_state, z_state, alpha, L, delta, transition_matrix, beta)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Identify the grid point that falls below k_next
lower_bound_index = max(sum(k_next>K_grid),1);

% The next value on the grid
upper_bound_index = min(lower_bound_index + 1,length(K_grid));

% Technology shock notation
if z_state == z_vector(1)

    zz = 1;

end

if z_state == z_vector(2)

    zz = 2;

end

if z_state == z_vector(3)

    zz = 3;

end

if lower_bound_index == upper_bound_index 
    
    cons = exp(z_state) * k_state^(alpha) * L(lower_bound_index, zz)^(1-alpha) + ...
    (1-delta) * k_state - K_grid(lower_bound_index); 

    if cons <= 0 

        vvv = -10^7;

    else

        vvv = (1-beta) * (log(cons) - L(lower_bound_index, zz)^2/2) + ...
        beta * V_prev(lower_bound_index, :) * transpose(transition_matrix(zz,:));

    end





else

    % Value function in the lower bound index 
    cons = exp(z_state) * k_state^(alpha) * L(lower_bound_index, zz)^(1-alpha) + ...
        (1-delta) * k_state - K_grid(lower_bound_index); 
    
    if cons <= 0 
    
        v_lower_bound = -10^7;
    
    else
    
        v_lower_bound = (1-beta) * (log(cons) - L(lower_bound_index, zz)^2/2) + ...
            beta * V_prev(lower_bound_index, :) * transpose(transition_matrix(zz,:));
    
    end
    
    % Value function in the upper bound index 
    cons = exp(z_state) * k_state^(alpha) * L(upper_bound_index, zz)^(1-alpha) + ...
        (1-delta) * k_state - K_grid(upper_bound_index); 
    
    if cons <= 0 
    
        v_upper_bound = -10^7;
    
    else
    
        v_upper_bound = (1-beta) * (log(cons) - L(upper_bound_index, zz)^2/2) + ...
            beta * V_prev(upper_bound_index, :) * transpose(transition_matrix(zz,:));
    
    end
    
    vvv = v_lower_bound + (v_upper_bound - v_lower_bound) / (K_grid(upper_bound_index) ...
        - K_grid(lower_bound_index)) * (k_next - K_grid(lower_bound_index)); 
    
    vvv = -vvv; 
    
    % vub = v_upper_bound;
    % vlb = v_lower_bound;

end






end