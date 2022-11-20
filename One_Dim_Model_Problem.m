% --------------------------------------------------------------------
% Gustavo Grinsteins
% CU Boulder
% Mini-project
% One-dimensional model problem solver
% --------------------------------------------------------------------

% Return Values
% x - vector containing locations of the nodes
% d - finite element solutions at the nodes
% Input values
% k - polynomial degree for the approximation
% n_el - number of elements to be employed
% kappa - function handle for computing k at any point x
% f - function handle for computing f at any point x
% g_0 and g_L - boundary conditions
% L - domain length

function [x,d] = One_Dim_Model_Problem(k,n_el,kappa,f,g_0,g_L,L)
    %Initialize global system
    n = (k*n_el)-1; %DOF
    K = zeros(n,n);
    F = zeros(n,1);
    %element mesh size
    h_e = L/n_el;
    q_n = k+1; %quadrature points per element
    %% Get Quadrature points
    [xi_q,w_q] = quadrature_points(k);
    %% Get Shape functions
    [N_hat,N_hat_xi] = Shape_Functions(k,q_n,xi_q);
    %% FEM 1D Engine
    x = 0:h_e:L; %node locations
    for e=1:n_el
        %% Element Formation
        k_e = zeros(k+1,k+1);
        f_e = zeros(k+1,1);
        for q = 1:length(xi_q)
            for a = 1:k+1
                %Evaluate local quantities
                xe_1 = x(e); 
                xe_xi_q = xe_1+h_e*(((xi_q(q))+1)/2); 
                kappa_e_xe_xi = kappa(xe_xi_q );
                f_e_xe_xi = f(xe_xi_q);         
                for b = 1:k+1
                    k_e(a,b)=k_e(a,b)+kappa_e_xe_xi*N_hat_xi(a,q)*N_hat_xi(b,q)*w_q(q)*(2/h_e);
                end
                f_e(a) = f_e(a)+f_e_xe_xi*N_hat(a,q)*w_q(q)*(h_e/2);
            end
        end
        %% Element Assembly
        %Loop over rows of k^e
        for a = 1:k+1
            %Set Global Row IEN(a,e) = k*(e-1)+(a-1)
            A =  k*(e-1)+(a-1);
            if (A >= 1) && (A <= n)
                %Loop over columns of k^e
                for b = 1:k+1
                    %Set Global Columns
                    B =  k*(e-1)+(b-1);
                    if (B >= 1) && (B <= n)
                        K(A,B) = K(A,B) + k_e(a,b);
                    elseif (B == 0)
                        F(A) = F(A) - k_e(a,b)*g_0;                      
                    elseif B == n+1
                        F(A) = F(A) - k_e(a,b)*g_L; 
                    end
                end
                F(A) = F(A) + f_e(a);
            end
        end 
    end
    %Calculate and return values
    d_temp = K\F;
    d = [g_0, d_temp', g_L];
    x = 0:h_e/k:L;
end

