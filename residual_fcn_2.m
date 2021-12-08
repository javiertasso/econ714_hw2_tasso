function rho = residual_fcn_2(alpha,beta,delta,k_min,k_max,...
    rho,grid_k,T_k,Z,PI,n_polynomials,shock_num,M,l_min,l_max, tol_level)
% Solves for the coefficients associated to the Chebychev 
% polynomials. 
    % This function finds rho such that the residuals are zero 
    % That is why it starts with fsolve 

% Tolerance here initially was 10^(-15), I shouldn't change this
options=optimset('Display','Iter','TolFun',tol_level,'TolX',tol_level,...
    'MaxFunEvals', 10^5, 'MaxIter', 10^4);
rho = fsolve(@notime_iter,rho,options);

function res = notime_iter(rho)
    
    % Define some vectors to store things
    residual_section = zeros(n_polynomials*2,1);
    res = zeros(M*2,1);
        % const = (1-gamma)/theta; Only used when calculating the residual this
        % does not help here because we have a different utility function I
        % don't need that constant in this problem 
    
    % Note that res is initially a vector of zeros so rho will start at
    % zero I think that is the initial guess: 
    rho1 = rho(1:M,1);     % Coefficients for value fcn
    rho2 = rho(M+1:2*M,1); % Coefficients for labor
    
    % For each value of current the shock
    for z_index = 1:shock_num
        
        % Define some vector to store things 
        value = zeros(n_polynomials,1);
        g_l = zeros(n_polynomials,1);
        g_k = zeros(n_polynomials,1);
        g_c = zeros(n_polynomials,1);
        rho1_section = rho1(((z_index-1)*n_polynomials+1):z_index*n_polynomials);
        rho2_section = rho2(((z_index-1)*n_polynomials+1):z_index*n_polynomials);
        
        for k_index = 1:n_polynomials % Loop 1 over collocation point on k
            
            value(k_index) = dot(rho1_section,T_k(k_index,:)); % value fcn at each collocation points
            l = dot(rho2_section,T_k(k_index,:));   % labor at each collocation points
            k = grid_k(k_index);

             
            %disp(l)
            %disp(k)
     % In the original code here I have if l<0.1 and elseif l>0.9. Here I will
     % change it because in the steady state labor is greater than 0.9
            
            % Check that labor makes sense 
            if(l<l_min)
                l = l_min;
                disp('l break lower bound')
                 
                 
            elseif(l>l_max)
                l = l_max;
                disp('l break upper bound')
                 
            end

             
            g_l(k_index) = l;
            
            % Compute current output, current consumption and capital next period
            y = exp(Z(z_index))*k^(alpha)*l^(1-alpha); % Same production function 
            c = min((1-alpha) * y / (l^2), y); % This changed. Note that the value of consumption here is always positive
            kp = y+(1-delta)*k-c;
            
            % Check that capital next period makes sense 
            if( kp < k_min )
                kp = k_min + 0.01;
               disp('kp break lower bound')
                
            elseif((kp > y+(1-delta)*k-c) || (kp > k_max))
                kp = min(y+(1-delta)*k-c,k_max) - 0.01;
                disp('kp break upper bound')
            end

            % Capital next period may not be in the grid
            g_k(k_index) = kp;
    
            % Compute consumption again, now that we have capital next
            % period
            c = y+(1-delta)*k-kp;
            
            % This one can be negative, so check that 
            if(c < 0)
                disp('warning: c < 0')
                c = 10^(-10);
            end
            g_c(k_index) = c;
    
        end % Loop 1 over collocation point on k ends
    
        % scale k prime from [k_min,k_max] to [-1,1]
        g_k_scaled_down=(2*g_k-(k_min+k_max))/(k_max-k_min);
        
        % I didn't touch anything here
        % value of polynomials at each scaled k prime
        T_g_k=ones(n_polynomials,n_polynomials);
        T_g_k(:,2)=g_k_scaled_down;
        % Use recursive definition of Chebyshev polynomials here
        for i1=3:n_polynomials
            T_g_k(:,i1)=2*g_k_scaled_down.*T_g_k(:,i1-1)-T_g_k(:,i1-2);
        end
        
        % Calculate residual
        for k_index = 1:n_polynomials % Loop 2 over collocation point on k
            
            % Define vectors to store things 
            vp=zeros(shock_num,1);
            temp = zeros(shock_num,1);
            
            % Again for each value of the shock. We are going to average
            % this out, but we need it for the value function

            % What do I want to do in this loop? 
            for zp_index = 1:shock_num
    
                rho1_section = rho1(((zp_index-1)*n_polynomials+1):zp_index*n_polynomials);
                rho2_section = rho2(((zp_index-1)*n_polynomials+1):zp_index*n_polynomials);
                vp(zp_index) = dot(rho1_section,T_g_k(k_index,:)); % Here I'm using the coefficients to get an approximation for V in the next period and for every shock
                
                % Labor in the next period, to get consumption in the next
                % period which in turn is used to get capital in two
                % periods from now 
                lp = dot(rho2_section,T_g_k(k_index,:));     
                % In the original code here I have if l<0.1 and elseif l>0.9. Here I will
                % change it
                if( lp < l_min )
                    lp = l_min;
                    disp('lp break lower bound')
                elseif(lp > l_max)
                    lp = l_max;
                    disp('lp break upper bound')
                end
                

                % output and consumption in the next period, capital in two
                % periods from now
                yp = exp(Z(zp_index))*g_k(k_index)^(alpha)*lp^(1-alpha); % This stays the same
                cp = min((1-alpha)*yp/(lp^2), yp); % This changed 
                kpp = yp+(1-delta)*g_k(k_index)-cp; % This stays the same
    
                % Check that capital in two periods makes sense 
                if( kpp < k_min )
                    kpp = k_min + 0.01;
                    % disp('kpp break lower bound')
                elseif((kpp > yp+(1-delta)*g_k(k_index)-cp) || (kpp > k_max))
                    kpp = min(yp+(1-delta)*g_k(k_index)-cp,k_max) - 0.01;
                    % disp('kpp break upper bound')
                end
                
                % Recompute consumption tomorrow and check it makes sense
                cp = yp+(1-delta)*g_k(k_index)-kpp; % This stays the same 
                if(cp < 0)
                    disp('warning: cp < 0')
                    cp = 10^(-10);
                end             
                
                % What is Up? It is the period utility tomorrow
                % Up = (cp^nu*(1-lp)^(1-nu))^const;
                % Up = (log(cp) - (lp^2)/2); Not needed
    
                % What is Ucp? It is the marginal utility of consumption
                % tomorrow
                % Ucp = const*Up*nu/cp;
                Ucp = 1 / cp; 
    
                % What is Fkp? It is the marginal product of capital in the
                % next period 
                % This one stays the same 
                Fkp = alpha * exp(Z(zp_index)) * g_k(k_index)^(alpha-1) * lp^(1-alpha);
    
                % What is this temporal variable?
                % This could be the derivative of the value function wrt capital
                % Vp is the expected value, I know that becuase there is a dot
                % product when constructing it
                % I am not sure I understand what this temporal variable is
                % doing
                % temp(zp_index) = vp(zp_index)^(-const*(1-theta)) * Ucp * (Fkp+1-delta);
                % temp(zp_index) = vp(zp_index) * Ucp * (Fkp+1-delta);
                temp(zp_index) = Ucp * (Fkp + 1-delta);
            end
            
            % Right hand side of the euler equation
            % Why do we have a double dot product here? I don't get that 
            % euler_rhs = beta*dot(PI(z_index,:),vp.^(1-gamma)).^(1/theta-1)* ...
                         % dot(PI(z_index,:),temp);
            % euler_rhs = beta * dot(PI(z_index,:),vp) * dot(PI(z_index,:),temp);
            euler_rhs = beta * dot(PI(z_index,:), temp);
             
            % euler_rhs = beta * dot(PI(z_index,:), vp); 
    
            l = g_l(k_index);
            c = g_c(k_index);
    
            % Derivative with respect to consumption 
            % euler_lhs = const*(c^nu*(1-l)^(1-nu))^const*nu/c; 
            euler_lhs = 1 / c; % This was modified
            
            % Here I have to modify things because the objective function is
            % different. I will keep the (1-beta) to express everything in
            % units that make sense 
            % bellman_rhs = (1-beta)*(c^nu*(1-l)^(1-nu))^((1-gamma)/theta)...
                          % +beta*dot(PI(z_index,:),vp.^(1-gamma)).^(1/theta);
            % bellman_rhs = bellman_rhs^(theta/(1-gamma));
            bellman_rhs = (1-beta) * (log(c) - (l^2)/2) + beta * dot(PI(z_index,:), vp); 
            % bellman_rhs = (1-beta) * (log(c) - (l^2)/2) + beta * dot(PI(z_index,:), temp); 
            bellman_lhs = value(k_index);
            
            % Here compute the residuals 
            residual_section(k_index) = euler_rhs - euler_lhs;
            residual_section(n_polynomials+k_index) = bellman_rhs - bellman_lhs;
    
    
        end % Loop 2 over k ends
        
        res(((z_index-1)*n_polynomials*2+1):z_index*n_polynomials*2) = residual_section;
         
    end

end

end