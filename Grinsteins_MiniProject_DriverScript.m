% --------------------------------------------------------------------
% Gustavo Grinsteins
% CU Boulder
% Mini-project
% One-dimensional model problem solver
% --------------------------------------------------------------------

%This driver script contains the logic that returns the responses for the
%mini-project problems

%% House Keeping
clc;
clear;
close all;
warning('off','all');

%% Problem 2 Code validation tests
%Graphical Validation of Problem 1
k = 1;
n_el = 3;
kappa = @(x) 1;
f = @(x) -x;
g_0 = 0;
g_L = 0;
L = 1;

[x,d] = One_Dim_Model_Problem(k,n_el,kappa,f,g_0,g_L,L);

fprintf('Rendering Figure 1 - P1 Approx k=1 n_el=3 \n')
figure(1)
plot(x,d), title('HW1 P1 Approx k=1 nel=3')
xlim([0 1])
hold on
fplot(@(x) (x.*(x.^2-1))/6,[0 1],'b')
legend('PDE Galerkin Approximation','PDE Analytical Solution')
hold off
grid on
%%%
k = 2;
n_el = 3;
kappa = @(x) 1;
f = @(x) -x;
g_0 = 0;
g_L = 0;
L = 1;

[x,d] = One_Dim_Model_Problem(k,n_el,kappa,f,g_0,g_L,L);

fprintf('Rendering Figure 2 - P1 Approx k=2 n_el=3 \n')
figure(2)
plot(x,d), title('HW1 P1 Approx k=2 nel=3')
xlim([0 1])
hold on
fplot(@(x) (x.*(x.^2-1))/6,[0 1],'b')
legend('PDE Galerkin Approximation','PDE Analytical Solution')
hold off
grid on
%%%
k = 3;
n_el = 5;
kappa = @(x) 1;
f = @(x) -x;
g_0 = 0;
g_L = 0;
L = 1;

[x,d] = One_Dim_Model_Problem(k,n_el,kappa,f,g_0,g_L,L);

fprintf('Rendering Figure 3 - P1 Approx k=3 n_el=3 \n')
figure(3)
plot(x,d), title('HW1 P1 Approx k=3 nel=3')
xlim([0 1])
hold on
fplot(@(x) (x.*(x.^2-1))/6,[0 1],'b')
legend('PDE Galerkin Approximation','PDE Analytical Solution')
hold off
grid on
%%%
%Calculate h1 norm error for increasing elements and compare
u = @(x) (x.*(x.^2-1))/6;
u_x = @(x) x^2/2 - 1/6;
n_el = 1:100;
errors_1 = zeros(1,length(n_el));
errors_2 = zeros(1,length(n_el));
errors_3 = zeros(1,length(n_el));
for i = 1:length(n_el)
    [~,d1] = One_Dim_Model_Problem(1,n_el(i),kappa,f,g_0,g_L,L);
    errors_1(i) = CalculateH1NormError(u,u_x,d1,1,n_el(i),L);
    [~,d2] = One_Dim_Model_Problem(2,n_el(i),kappa,f,g_0,g_L,L);
    errors_2(i) = CalculateH1NormError(u,u_x,d2,2,n_el(i),L);
    [~,d3] = One_Dim_Model_Problem(3,n_el(i),kappa,f,g_0,g_L,L);
    errors_3(i) = CalculateH1NormError(u,u_x,d3,3,n_el(i),L);
end
fprintf('Rendering Figure 4 - H1 Norm Comparison \n')
figure(4)
loglog(n_el,errors_1), title('HW1 Problem 1 Bubnov-Galerkin Approximation H1 Norm Error')
hold on
loglog(n_el,errors_2)
loglog(n_el,errors_3)
xlabel('Number of elements')
ylabel('H1 Norm Error')
legend('k = 1','k = 2','k = 3')
hold off
grid on
ConvRate1 = diff(log(errors_1))./diff(log(L./n_el));
ConvRate2 = diff(log(errors_2))./diff(log(L./n_el));
ConvRate3 = diff(log(errors_3))./diff(log(L./n_el));
fprintf('Convergence rate k = 1 -> %1.0f \n',ConvRate1(end))
fprintf('Convergence rate k = 2 -> %1.0f \n',ConvRate2(end))
fprintf('Convergence rate k = 3 -> %1.0f \n',ConvRate3(end))

%Manufactured solution
%Graphical Validation
kappa = @(x) 1;
f = @(x) (pi^2-1)*exp(-x)*sin(pi*x)+2*pi*exp(-x)*cos(pi*x);
g0 = 0;
gL = 0;
L = 1;
u = @(x) exp(-x)*sin(pi*x);
u_x = @(x) pi*exp(-x)*cos(pi*x) - exp(-x)*sin(pi*x);

[x,d] = One_Dim_Model_Problem(1,3,kappa,f,g_0,g_L,L);

fprintf('Rendering Figure 5 - Manufactured Solution Approx k=1 n_el=3 \n')
figure(5)
plot(x,d), title('Manufactured Solution Approx k=1 nel=3')
xlim([0 1])
hold on
fplot(@(x) exp(-x)*sin(pi*x),[0 1],'b')
legend('PDE Galerkin Approximation','PDE Analytical Solution')
hold off
grid on

[x,d] = One_Dim_Model_Problem(2,3,kappa,f,g_0,g_L,L);

fprintf('Rendering Figure 6 - Manufactured Solution Approx k=2 n_el=3 \n')
figure(6)
plot(x,d), title('Manufactured Solution Approx k=2 nel=3')
xlim([0 1])
hold on
fplot(@(x) exp(-x)*sin(pi*x),[0 1],'b')
legend('PDE Galerkin Approximation','PDE Analytical Solution')
hold off
grid on

[x,d] = One_Dim_Model_Problem(3,3,kappa,f,g_0,g_L,L);

fprintf('Rendering Figure 7 - Manufactured Solution Approx k=3 n_el=3 \n')
figure(7)
plot(x,d), title('Manufactured Solution Approx k=3 nel=3')
xlim([0 1])
hold on
fplot(@(x) exp(-x)*sin(pi*x),[0 1],'b')
legend('PDE Galerkin Approximation','PDE Analytical Solution')
hold off
grid on

%Calculting the h1 norm error
n_el = 1:100;
errors_1 = zeros(1,length(n_el));
errors_2 = zeros(1,length(n_el));
errors_3 = zeros(1,length(n_el));
for i = 1:length(n_el)
    [~,d1] = One_Dim_Model_Problem(1,n_el(i),kappa,f,g_0,g_L,L);
    errors_1(i) = CalculateH1NormError(u,u_x,d1,1,n_el(i),L);
    [~,d2] = One_Dim_Model_Problem(2,n_el(i),kappa,f,g_0,g_L,L);
    errors_2(i) = CalculateH1NormError(u,u_x,d2,2,n_el(i),L);
    [~,d3] = One_Dim_Model_Problem(3,n_el(i),kappa,f,g_0,g_L,L);
    errors_3(i) = CalculateH1NormError(u,u_x,d3,3,n_el(i),L);
end
fprintf('Rendering Figure 8 - H1 Norm Comparison for manufactured solution \n')
figure(8)
loglog(n_el,errors_1), title('u(x)=exp(-x)*sin(pi*x) Bubnov-Galerkin Approximation H1 Norm Error')
hold on
loglog(n_el,errors_2)
loglog(n_el,errors_3)
xlabel('Number of elements')
ylabel('H1 Norm Error')
legend('k = 1','k = 2','k = 3')
hold off
grid on
ConvRate1 = diff(log(errors_1))./diff(log(L./n_el));
ConvRate2 = diff(log(errors_2))./diff(log(L./n_el));
ConvRate3 = diff(log(errors_3))./diff(log(L./n_el));
fprintf('Convergence rate k = 1 -> %1.0f \n',ConvRate1(end))
fprintf('Convergence rate k = 2 -> %1.0f \n',ConvRate2(end))
fprintf('Convergence rate k = 3 -> %1.0f \n',ConvRate3(end))

%% Problem 3 - Cylindrical Rod Heating

kappa = @(x) 385; %Watts per meter per degree Celsius.
g_0 = 30; %Degrees Celsius
g_L = 30; %Degrees Celsius
L = 2; % Meters
h_max = 1500; % Watts per square meter
R = 0.025; % Meters
f = @(x) (2/(R))*h_max*exp(-100*((x/L)-0.5)^2);

%Calculate Max Temp using different mesh sizes and k degrees
n_el = 1:50;
maxT_1 = zeros(1,length(n_el));
maxT_2 = zeros(1,length(n_el));
maxT_3 = zeros(1,length(n_el));
for i = 1:length(n_el)
    [~,d1] = One_Dim_Model_Problem(1,n_el(i),kappa,f,g_0,g_L,L);
    maxT_1(i) = max(d1);
    [~,d2] = One_Dim_Model_Problem(2,n_el(i),kappa,f,g_0,g_L,L);
    maxT_2(i) = max(d2);
    [~,d3] = One_Dim_Model_Problem(3,n_el(i),kappa,f,g_0,g_L,L);
    maxT_3(i) = max(d3);
end
fprintf('Rendering Figure 9 - P3 Maximum temperature Convergence \n')
figure(9)
plot(n_el,maxT_1), title('Mini-Project P3: Maximum Temperature Convergence')
hold on
plot(n_el,maxT_2)
plot(n_el,maxT_3)
xlabel('Number of elements')
ylabel('Temperature (Celsius)')
legend('k = 1','k = 2','k = 3')
hold off
grid on
fprintf('P3 Max T Calculation %1.2f Degrees Celsius\n', maxT_3(end))

%% Problem 4 - Rod Tensile Test

%k = 1; 
%n_el = 3;
g_0 = 0; %meters
g_L = 0.00002; %Degrees Celsius
L = 0.1; % Meters
f = @(x) 0;
A_max = 0.0001; % square meters
A_min = 0.00002; %square meters
E = 200e9;%Pascals
kappa = @(x) A_max-(A_max-A_min)*exp(-50*((x/L)-0.5)^2); % Square Meters

%Calculate Max Temp using different mesh sizes and k degrees
n_el = 1:100;
maxSigma_1 = zeros(1,length(n_el));
maxSigma_2 = zeros(1,length(n_el));
maxSigma_3 = zeros(1,length(n_el));
for i = 1:length(n_el)
    [x1,d1] = One_Dim_Model_Problem(1,n_el(i),kappa,f,g_0,g_L,L);
    maxSigma_1(i) = max(E.*(diff(d1)./diff(x1)));
    [x2,d2] = One_Dim_Model_Problem(2,n_el(i),kappa,f,g_0,g_L,L);
    maxSigma_2(i) = max(E.*(diff(d2)./diff(x2)));
    [x3,d3] = One_Dim_Model_Problem(3,n_el(i),kappa,f,g_0,g_L,L);
    maxSigma_3(i) = max(E.*(diff(d3)./diff(x3)));
end
fprintf('Rendering Figure 10 - P4 Max Sigma Convergence \n')
figure(10)
plot(n_el,maxSigma_1), title('P4 Maximum Axial Stress Convergence ')
hold on
plot(n_el,maxSigma_2)
plot(n_el,maxSigma_3)
xlabel('Number of elements')
ylabel('Stress (Pascals)')
legend('k = 1','k = 2','k = 3')
hold off
grid on
fprintf('P4 Max Sigma Calculation %1.3f MPa\n', maxSigma_3(end)/(10^6))
