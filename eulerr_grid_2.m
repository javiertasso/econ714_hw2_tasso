% This file computes the Euler Errors on the capital and exogenous stock
% grid

function [g_k,g_c,g_l,value_fcn,euler_error,max_error]= ...
         eulerr_grid_2(alpha,beta,delta,rho,Z, ...
         PI,k_min,k_max,grid_k_complete,shock_num,n_polynomials,grid_num,M)

% Not needed here
%const = (1-gamma)/theta;

grid_k_complete_scaled=(2*grid_k_complete-(k_min+k_max))/(k_max-k_min);

T_k_complete=ones(grid_num,n_polynomials);
T_k_complete(:,2)=grid_k_complete_scaled;
for i1=3:n_polynomials
    T_k_complete(:,i1)=2*grid_k_complete_scaled.*T_k_complete(:,i1-1)-T_k_complete(:,i1-2);
end     

rho1 = rho(1:M,1);     % Coefficients for value fcn
rho2 = rho(M+1:2*M,1); % coefficients for labor

euler_error=zeros(grid_num,shock_num);
value_fcn=zeros(grid_num,shock_num);
g_l=zeros(grid_num,shock_num);
g_c=zeros(grid_num,shock_num);
g_k=zeros(grid_num,shock_num);
for z_index = 1:shock_num
    for k_index = 1:grid_num % Loop 1 over collocation point on k

        rho1_section = rho1(((z_index-1)*n_polynomials+1):z_index*n_polynomials);
        rho2_section = rho2(((z_index-1)*n_polynomials+1):z_index*n_polynomials);
        value_fcn(k_index,z_index) = dot(rho1_section,T_k_complete(k_index,:)); % Value fcn at each collocation points
        g_l(k_index,z_index) = dot(rho2_section,T_k_complete(k_index,:));   % Labor at each collocation points

        y = exp(Z(z_index))*grid_k_complete(k_index)^alpha*g_l(k_index,z_index)^(1-alpha);
        g_c(k_index,z_index) = (1-alpha)*y/(g_l(k_index,z_index)^2);            
        g_k(k_index,z_index) = y+(1-delta)*grid_k_complete(k_index)-g_c(k_index,z_index);            

    end % Loop 1 over collocation point on k ends

    % Scale k prime from [k_min,k_max] to [-1,1]
    g_k_scaled_down=(2*g_k(:,z_index)-(k_min+k_max))/(k_max-k_min);
    % value of polynomials at each scaled k prime
    T_g_k=ones(grid_num,n_polynomials);
    T_g_k(:,2)=g_k_scaled_down;
    for i1=3:n_polynomials
        T_g_k(:,i1)=2*g_k_scaled_down.*T_g_k(:,i1-1)-T_g_k(:,i1-2);
    end     

    % Calculate residual        
    for k_index = 1:grid_num % Loop 2 over collocation point on k              
        vp = zeros(shock_num,1);
        temp = zeros(shock_num,1);
        for zp_index = 1:shock_num
            rho1_section = rho1(((zp_index-1)*n_polynomials+1):zp_index*n_polynomials);
            rho2_section = rho2(((zp_index-1)*n_polynomials+1):zp_index*n_polynomials);
            vp(zp_index) = dot(rho1_section,T_g_k(k_index,:));
            lp = dot(rho2_section,T_g_k(k_index,:));

            yp = exp(Z(zp_index))*g_k(k_index,z_index)^alpha*lp^(1-alpha);
            cp = (1-alpha)*yp/(lp^2);

            % Up = (cp^nu*(1-lp)^(1-nu))^const;
            Ucp = 1/cp;
            Fkp = alpha*exp(Z(zp_index))*g_k(k_index,z_index)^(alpha-1)*lp^(1-alpha);
            temp(zp_index) = Ucp*(Fkp+1-delta);
        end

        euler_rhs = beta * dot(PI(z_index,:),temp);

        % A = euler_rhs/(const*nu*(1-g_l(k_index,z_index))^((1-nu)*const));

        euler_error(k_index,z_index) = 1/g_c(k_index,z_index) - euler_rhs;

        if(log10(abs(euler_error(k_index,z_index))) < -17 )
            disp('euler_error')
            disp(euler_error(k_index,z_index))
            disp(log10(abs(euler_error(k_index,z_index))))
            euler_error(k_index,z_index) = 10^(-17);
        end

    end % Loop 2 over k ends

end % Loop over z ends

euler_error = log10(abs(euler_error));
max_error = max(euler_error,[],1);