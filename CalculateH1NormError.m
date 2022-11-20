% --------------------------------------------------------------------
% Gustavo Grinsteins
% CU Boulder
% Mini-project
% One-dimensional model problem solver
% --------------------------------------------------------------------

% Return Values
% h1_norm_error - Calculated h1_norm_error value
% Input Values
% u - analytical solution
% u_x - derivative of solution
% d - degress of freedom for 1D solution
% k - polynomial degree
% n_el - number of elements to be employed
% L - domain length

function h1_norm_error = CalculateH1NormError(u,u_x,d,k,n_el,L)

    %% Get Quadrature points
    [xi_q,w_q] = quadrature_points(k);
    %% Get Shape functions
    [N_hat,N_hat_xi] = Shape_Functions(k,length(xi_q),xi_q);
    h_e = L/n_el;
    x = 0:h_e:L;
    %Initialize value
    h1_error = 0;
    for e=1:n_el
        u_approx = 0;
        u_approx_xi = 0;
        for q = 1:length(xi_q)
            %Obtain approximations given k
            for a = 1:k+1
                A =  k*(e-1)+(a-1)+1;
                %Note: Output d already contains boundary conditions
                u_approx = u_approx + d(A)*N_hat(a,q);
                u_approx_xi = u_approx_xi + d(A)*N_hat_xi(a,q);
            end
            %u_approx = u_approx*(h_e/2)*w_q(q);
            u_approx_xi = u_approx_xi*(2/h_e);
            %Calculate parent element values
            xe_1 = x(e); 
            xe_xi_q = xe_1+h_e*(((xi_q(q))+1)/2);
            u_true = u(xe_xi_q);
            u_true_xi = u_x(xe_xi_q);
            %Calculate h1 norm error
            h1_error = (h1_error + ((u_true-u_approx)^2 + (L^2)*(u_true_xi-u_approx_xi)^2)*w_q(q)*(h_e/2));
        end
    end

    %Final Error Value
    h1_norm_error = sqrt(h1_error);

end
