function [psi] = weight_fun_capital(k, k_elements)
% weigth_fun_capital 
%   It requires a level of capital 
%   It requires the vector of elements 
%   It spits up the weigths 

psi = zeros(length(k_elements) - 1, 1); 

% First one
if k <= k_elements(2)

    psi(1) = (k - k_elements(1)) / (k_elements(2) - k_elements(1)); 

end

% Loop here
for ii = 2:(length(k_elements) - 1)

    if k <= k_elements(ii+1) && k_elements(ii) <= k
        
        psi(ii-1) = (k_elements(ii+1) - k) / (k_elements(ii+1) - k_elements(ii)); 
        psi(ii) = (k - k_elements(ii)) / (k_elements(ii+1) - k_elements(ii));
    
    end

end 

% Last one 
% if k >= k_elements(length(k_elements))

    % psi(9) = max(0, (k - k_elements(8) / (k_elements(9) - k_elements(8)))); 

% end


end