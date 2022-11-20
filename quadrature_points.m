% --------------------------------------------------------------------
% Gustavo Grinsteins
% CU Boulder
% Mini-project
% One-dimensional model problem solver
% --------------------------------------------------------------------

%hard coded quadrature values up to k = 3. Assuming q_n = k+1

% Return Values
% xi_q - quadrature point
% w_q - quadrature weight
% Input Values
% k - polynomial degree

function [xi_q,w_q] = quadrature_points(k)
    if k == 1
        xi_q = [-1/sqrt(3),1/sqrt(3)];
        w_q = [1,1];
    elseif k == 2
        xi_q = [0,-sqrt(3/5),sqrt(3/5)];
        w_q = [8/9,5/9,5/9];
    else
        xi_q = [sqrt(3/7-2/7*sqrt(6/5)),-sqrt(3/7-2/7*sqrt(6/5)),sqrt(3/7+2/7*sqrt(6/5)),-sqrt(3/7+2/7*sqrt(6/5))];
        w_q = [(18+sqrt(30))/36,(18+sqrt(30))/36,(18-sqrt(30))/36,(18-sqrt(30))/36];
    end
end