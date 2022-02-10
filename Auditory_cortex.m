clear all; clc; close all;

% Model of auditory primary cortex A1, following
% Loebel et al.(2007)
% https://doi.org/10.3389/neuro.01.1.1.015.2007

%% Parameters
columns = 5; %15;
cells_per_column = 4;%200;
N_E = 2;%100;
N_I = cells_per_column - N_E;

Tau_E = 0.001; % s
Tau_I = 0.001; % s

Tau_ref_E = 3*10^(-3); % s
Tau_ref_I = 3*10^(-3); % s

J_EE = [0.015,0.045,6,0.045,0.015]; % index goes from -2 to 2
J_IE = [0.0015,0.0035,0.5,0.0035,0.0015];
J_EI = -4; 
J_II = -0.5;

% synaptic resources
U = 0.5;
Tau_rec = 0.8; % s

%% Compact all parameters
params.N_E = N_E;
params.N_I = N_I;
params.columns = columns;
params.Tau_E = Tau_E;
params.Tau_I = Tau_I;
params.Taur_ref_E = Tau_ref_E;
params.Taur_ref_I = Tau_ref_I;
params.J_EE = J_EE;
params.J_EI = J_EI;
params.J_IE = J_IE;
params.J_II = J_II;
params.U = U;
params.Taur_rec = Tau_rec;

%% External input
% background activity
e_max = 10; % Hz
e_min = -10; % Hz

% stimulus
s = 0;

% Create background noise
e_E = zeros(N_E,columns); % Hz
e_I = zeros(N_I,columns); % Hz

rng('default');
for ii = 1:columns
    init_e_E = rand(N_E,1); %between (0,1)
    init_e_I = rand(N_I,1);
    % scale
    highest_eE = max(init_e_E);
    lowest_eE = min(init_e_E);
    highest_eI = max(init_e_I);
    lowest_eI = min(init_e_I);
    e_E_range = max(init_e_E) - min(init_e_E);
    e_I_range = max(init_e_I) - min(init_e_I);
    
    actual_eE = (1/e_E_range)*(init_e_E*(e_max-e_min)+(e_min*highest_eE-e_max*lowest_eE));
    actual_eI = (1/e_I_range)*(init_e_I*(e_max-e_min)+(e_min*highest_eI-e_max*lowest_eI));

    e_E(:,ii) = sort(actual_eE);
    e_I(:,ii) = sort(actual_eI);
end

%% Create neurons
E = rand(N_E,columns);
I = rand(N_I,columns);

%% Create synaptic resources
x = rand(size(E));
y = rand(size(I));

%% Get initial conditions
% set stimulus to 0 and let the system converge
X0 = [E(:); I(:); x(:); y(:)];
t_rest = 5*Tau_rec;
[t, X] = ode45(@(t,y) dXdt(t,y,params,0,e_E,e_I),[0,t_rest],X0);

[E,I,x,y] = reshape_X(X,[2,3,1], params);

E0 = E(:,:,end);
I0 = I(:,:,end);
x0 = x(:,:,end);
y0 = y(:,:,end);

%% Solve ODE
t_span = [0, 5*Tau_rec];
X0 = [E0(:); I0(:); x0(:); y0(:)];
[t, X] = ode45(@(t,y) dXdt(t,y,params,s,e_E,e_I),t_span,X0);

[E,I,x,y] = reshape_X(X,[1,3,2], params);


function dXdt = dXdt(t,X,params,s,e_E,e_I)
    [E,I,x,y] = reshape_X_single_time_point(X,params);
   
    N_E = params.N_E;
    N_I = params.N_I;
    columns = params.columns;
    Tau_E = params.Tau_E;
    Tau_I = params.Tau_I;
    Tau_ref_E = params.Taur_ref_E;
    Tau_ref_I = params.Taur_ref_I;
    J_EE = params.J_EE;
    J_EI = params.J_EI;
    J_IE = params.J_IE;
    J_II = params.J_II;
    U = params.U;
    Tau_rec = params.Taur_rec;
    
    in_col_exci = sum(U*x.*E, 1);
    between_col_exci = zeros(5,columns);
    for R = -2:2 
        aux = circshift(in_col_exci,-R,2);
        if R<0
            aux(:,1:-R)=0;
        elseif R>0
            aux(:,end-(R-1):end)=0;
        end
        between_col_exci(R+3,:) = aux;
    end
    between_col_exci = (J_EE/N_E)*between_col_exci;
    in_col_inhi = (J_EI/N_I)*sum(U*y.*I, 1);
    relu = max((between_col_exci + in_col_inhi + e_E + s),0);

    dEdt= (1/Tau_E)*(-E + (1-Tau_ref_E*E).*relu);

    %% dIdt
    in_col_exci = sum(E, 1);
    between_col_exci = zeros(5,columns);
    for R = -2:2 
        aux = circshift(in_col_exci,-R,2);
        if R<0
            aux(:,1:-R)=0;
        elseif R>0
            aux(:,end-(R-1):end)=0;
        end
        between_col_exci(R+3,:) = aux;
    end
    between_col_exci = (J_IE/N_E)*between_col_exci;
    in_col_inhi = (J_II/N_I)*sum(I, 1);
    relu = max((between_col_exci + in_col_inhi + e_I),0);

    dIdt= (1/Tau_I)*(-I + (1-Tau_ref_I*I).*relu);

    %% dxdt
    dxdt = (1-x)/Tau_rec - U*x.*E;
    %% dydt
    dydt = (1-y)/Tau_rec - U*y.*I;


    dXdt = [dEdt(:); dIdt(:); dxdt(:); dydt(:)];
end
function [E,I,x,y] = reshape_X_single_time_point(X,params)
% Gets a vector X of size (4*cells_per_column*columns, 1) and returns
% 4 arrays: E, I, x, y; each of them of size (N_i, columns), where i 
% is either E for excitatory cells or I for inhibitory cells.

    N_E = params.N_E;
    N_I = params.N_I;
    columns = params.columns; 
    % re-organize the vector X in its components
    cont = 1;
    E_col = X(1:N_E*columns);
    E = reshape(E_col,N_E,columns);
    cont = cont + size(E_col,1);
    
    I_col = X(cont:cont + N_I*columns-1);
    I = reshape(I_col,N_I,columns);
    cont = cont + size(I_col,1);
    
    x_col = X(cont:cont + N_E*columns-1);
    x = reshape(x_col,N_E,columns);
    cont = cont + size(x_col,1);

    y_col = X(cont: cont + N_I*columns-1);
    y = reshape(y_col,N_I,columns);
end
function [E,I,x,y] = reshape_X(X,order,params)
% Gets an array X of size (time_steps, 4*cells_per_column*columns) and returns
% 4 arrays: E, I, x, y; each of them of size (N_i, columns, time_steps), where i 
% is either E for excitatory cells or I for inhibitory cells.
% The order of the axes is given by the argument order:
% order = [2,3,1] -> (N_i, columns, time_steps). Useful when you need a
% 'photograph' of the variable at a single time_step.
%order = [1,3,2] -> (time_steps, columns, N_i. Useful for plotting the
%variable over time.

    N_E = params.N_E;
    N_I = params.N_I;
    columns = params.columns; 
    
    % Separate the vector X into E,I,x,y
    time_steps = size(X,1);
    cont = 1;
    
    E_col = X(:,cont:cont+N_E*columns-1);
    %E_col = X(1:time_steps*N_E*columns);
    cont = cont + size(E_col,2);
    I_col = X(:,cont:cont+N_I*columns-1);
    %I_col = X(1:time_steps*N_I*columns);
    cont = cont + size(I_col,2);
    x_col = X(:,cont:cont+N_E*columns-1);
    %x_col = X(1:time_steps*N_E*columns);
    cont = cont + size(x_col,2);
    y_col = X(:,cont:cont+N_I*columns-1);
    %y_col = X(1:time_steps*N_I*columns);

    % reshape each tensor in 3 dimensions, one for time, other for columns,
    % other for cells in the column
    E = reshape(E_col,time_steps,N_E,columns);
    I = reshape(I_col,time_steps,N_I,columns);
    x = reshape(x_col,time_steps,N_E,columns);
    y = reshape(y_col,time_steps,N_I,columns);

    % Re-organize the axes so that each tensor has the shape (cells_per_column,
    % columns, t)
    E = permute(E,order);
    I = permute(I,order);
    x = permute(x,order);
    y = permute(y,order);
end