% --------------------------------------------------------------------
% Gustavo Grinsteins
% CU Boulder
% Mini-project
% One-dimensional model problem solver
% --------------------------------------------------------------------

% Return Values
% Nhat - shape functions
% Nhat_xi - shape functions derivative with respect to parent element
% Input Values
% k - polynomial degree
% q_n - number of quadrature points
% xi_q - quadrature point

function [Nhat,Nhat_xi] = Shape_Functions(k,q_n,xi_q)
    % Initializing shape function and its gradient
    Nhat  = zeros(k+1, q_n);
    Nhat_xi = zeros(k+1, q_n);
    if k == 1
        % Quadrature loop
        for q = 1 : q_n
             % Shape functions:
             Nhat(1, q) = 0.5 * (1 - xi_q(q));
             Nhat(2, q) = 0.5 * (1 + xi_q(q));
             % Derivative of shape functions:
             Nhat_xi(1, q) = - 0.5;
             Nhat_xi(2, q) =   0.5;
        end
    
    elseif k == 2
        % Quadrature loop
        for q = 1 : q_n
            % Shape functions:
            Nhat(1, q) = 0.5 * xi_q(q) * (xi_q(q)-1);
            Nhat(2, q) = -(xi_q(q)+1) * (xi_q(q)-1);
            Nhat(3, q) = 0.5 * xi_q(q) * (xi_q(q)+1);
            % Derivative of shape functions:
            Nhat_xi(1, q) = 0.5 * (2 * xi_q(q) - 1);
            Nhat_xi(2, q) = - 2 * xi_q(q);
            Nhat_xi(3, q) = 0.5 * (2 * xi_q(q) + 1);
        end

    elseif k == 3
        % Quadrature loop
        for q = 1 : q_n
            % Shape functions:
            Nhat(1, q) = - (9*xi_q(q)^3)/16 + (9*xi_q(q)^2)/16 + xi_q(q)/16 - 1/16;
            Nhat(2, q) = (27*xi_q(q)^3)/16 - (9*xi_q(q)^2)/16 - (27*xi_q(q))/16 + 9/16;
            Nhat(3, q) = (27*xi_q(q))/16 - (9*xi_q(q)^2)/16 - (27*xi_q(q)^3)/16 + 9/16;
            Nhat(4, q) = (9*xi_q(q)^3)/16 + (9*xi_q(q)^2)/16 - xi_q(q)/16 - 1/16;
            % Derivative of shape functions:
            Nhat_xi(1, q) = (9*xi_q(q))/8 - (27*xi_q(q)^2)/16 + 1/16;
            Nhat_xi(2, q) = (81*xi_q(q)^2)/16 - (9*xi_q(q))/8 - 27/16;
            Nhat_xi(3, q) = 27/16 - (81*xi_q(q)^2)/16 - (9*xi_q(q))/8;
            Nhat_xi(4, q) = (27*xi_q(q)^2)/16 + (9*xi_q(q))/8 - 1/16;
        end
        
    end

end 