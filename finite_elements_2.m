function theta = finite_elements_2(k_current_state, z_current_state, ...
    n_elements_k, n_elements_z, beta, theta, tol_level, l_min, l_max, alpha, delta, z)
%UNTITLED3 Summary of this function goes here
%   Needs a value of the current state (k,z)

% Options and tolerance level
options=optimset('Display','Iter','TolFun',tol_level,'TolX',tol_level,...
    'MaxFunEvals', 10^5, 'MaxIter', 10^4);
theta = fsolve(@weighted_residuals, theta, options);

    function w_res = weighted_residuals(theta)

        w_res = zeros(n_elements_k * n_elements_z);
        
        % Suppose first that i = 1 and z = -0.05
        
            % Need: 
                % labor today 
                % consumption and output today
                % capital tomorrow 
                % possible values of z tomorrow (average them out)
                % labor tomorrow (on each case)
                % consumption and output tomorrow
                % depreciation rate and parameters 
            l = min(l_max,max(l_min, transpose(theta(:,1)) *  ...
                weight_fun_capital(k, k_elements)));
            y = exp(-0.05) * k^(alpha) * l^(1-alpha);
            c = (1-alpha) * y / (l^2); 
            kp = min(k_max, max(k_min, y + (1-delta) * k - c));  
            
            % This loop is over the three possible states tomorrow
            lp = zeros(length(z),1);
            yp = zeros(length(z),1);
            cp = zeros(length(z),1);

            for zz = 1:length(z)

                lp(zz,1) = min(l_max, max(l_min,transpose(theta(zz,1)) *  ...
                weight_fun_capital(kp, k_elements)));
                yp(zz,1) = exp(zz) * kp^(alpha) * lp(zz,1)^(1-alpha);
                cp(zz,1) = (1-alpha) * yp(zz,1) / (lp(zz,1))^2;

            end

            res(z=-0.05) = beta * c * transition_matrix(z_current_state, :)...
                * (1./cp .* (alpha * yp / kp + 1 - delta)) - 1; 

            
                 
       


    end


end