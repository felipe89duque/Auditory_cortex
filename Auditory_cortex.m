clear variables; clc; close all;

% Model of auditory primary cortex A1, following
% Loebel et al.(2007)
% https://doi.org/10.3389/neuro.01.1.1.015.2007

%% Article Figures
% Use the function
% plot_article_figures(fig,params,s_params,e_E,e_I,initial_conditions,E0_non_zero),
% with fig equal the number of the figure in the article you want to plot.
% Use fig = -1 to plot all the figures
[params,s_params] = reset_parameters();
[initial_conditions,e_E,e_I,E0_non_zero] = get_initial_conditions_and_noise(params,10,0);
fig = 2;
plot_article_figures(fig,params,s_params,e_E,e_I,initial_conditions,E0_non_zero);

%% Extensions

%% Effect of random seed on the initial conditions
figure();
[params, ~] = reset_parameters();
seeds = 0:500;
naturally_active_exci_neurons = zeros(length(seeds),1);
for ii = 1:length(seeds)
    [initial_conditions, e_E,e_I,E0_non_zero] = get_initial_conditions_and_noise(params,seeds(ii),0);
    E0 = initial_conditions.E;
    I0 = initial_conditions.I;
    x0 = initial_conditions.x;
    y0 = initial_conditions.y;
    
    amount_of_active_neurons = sum(E0(:,1) > 0);
    naturally_active_exci_neurons(ii) = amount_of_active_neurons;
end
proportion_active_neurons = 100*naturally_active_exci_neurons / params.N_E;
plot(seeds,proportion_active_neurons,'k*');
hold on;
mean_active_neurons = mean(proportion_active_neurons);
std_active_neurons = std(proportion_active_neurons);
plot(seeds,mean_active_neurons*ones(length(seeds),1),'r','LineWidth',1.5);
plot(seeds,(mean_active_neurons - std_active_neurons)*ones(length(seeds),1),'r--','LineWidth',1.5);
plot(seeds,(mean_active_neurons + std_active_neurons)*ones(length(seeds),1),'r--','LineWidth',1.5);
xlabel('Random seed');
ylabel('% of excitatory neurons stimulated');
ylim([0,100]);
legend('',strcat(sprintf('mean = %.2f',mean_active_neurons),'%'),strcat('standard dev. = \pm',sprintf('%.2f',std_active_neurons),'%'));

%% Find if there is any seed for which a decreased background input lowrs the FTC (Figure 5D from article vs ours)
warning('off','all');
n_seeds = 500;
max_diffs = [0,0]; % max difference between increased and decrease noise, for either case (when increased > decreased and vise versa) 
article_seed = nan; % A seed that produces Figure 5D similar to the article's
not_article_seed = nan;
[params, s_params] = reset_parameters();
start_time = 0;
stimuli_durations = 0.05;

delta_e = 0.6;
e_max_inc = params.e_max + delta_e;
e_min_inc = params.e_min + delta_e;
e_max_dec = params.e_max - delta_e;
e_min_dec = params.e_min - delta_e;
e_max_vec = [e_max_inc,e_max_dec];
e_min_vec = [e_min_inc,e_min_dec];
for seed = 0:n_seeds
% Check if column 8 is more sensitive to stimulus when external noise is
% decreased
    inc_dec_sensitivity = zeros(1,2);
    for ii = 1:2
        params.e_min = e_min_vec(ii);
        params.e_max = e_max_vec(ii);
        [initial_conditions, e_E,e_I,E0_non_zero] = get_initial_conditions_and_noise(params,seed,0);
        E0 = initial_conditions.E;
        I0 = initial_conditions.I;
        x0 = initial_conditions.x;
        y0 = initial_conditions.y;
        for input = 0.2:0.2:9
            s = zeros(1,params.columns);
            s(8) = input;
            s = s_tuning_curve(s,s_params);
            [~,E,~,~,~] = solve_complex_stimuli(start_time,stimuli_durations,s,params,e_E,e_I,E0,I0,x0,y0,E0_non_zero);    
            E_col_avg = mean(E,3);
            E_col8 = E_col_avg(:,8);
            [~,E_indices]=findpeaks(E_col8,'MinPeakHeight', 15); %Completely arbitrary threshold O_O

            % Find if the stimulus presented caused a activation in column 8
            if length(E_indices)>= 1
               inc_dec_sensitivity(ii) = input;
               break
            end
        end
    end
    % If inc_dec_sensitivity(2) < inc_dec_sensitivity(1) then the result is
    % similar to the article's, otherwise it is the same as we got in our
    % Figure 5.
    if inc_dec_sensitivity(2) < inc_dec_sensitivity(1)
        if inc_dec_sensitivity(1) - inc_dec_sensitivity(2) > max_diffs(2)
            article_seed = seed;
            max_diffs(2) = inc_dec_sensitivity(1) - inc_dec_sensitivity(2);
        end
    else
        if inc_dec_sensitivity(2) - inc_dec_sensitivity(1) > max_diffs(1)
            not_article_seed = seed;
            max_diffs(1) = inc_dec_sensitivity(2) - inc_dec_sensitivity(1);
        end
    end
end

% Plot FTC for our results and for theirs
if ~isnan(article_seed)
    seeds = [not_article_seed,article_seed];
    plot_styles = ["k","k--"];
    plot_legends= ["Increased external noise", "Decreased external noise"];
    for seed_index = 1:length(seeds)
        subplot(1,2,seed_index);
        xlim([1,params.columns]);
        ylim([0,9]);
        xticks(2:2:params.columns);
        xlabel('Column index');
        ylabel('Input amplitude needed for PS in column 8 [Hz]');
        title(sprintf('seed = %d',seeds(seed_index)));
        hold on;
        for ii = 1:2
            params.e_min = e_min_vec(ii);
            params.e_max = e_max_vec(ii);
            [initial_conditions, e_E,e_I,E0_non_zero] = get_initial_conditions_and_noise(params,seeds(seed_index),0);
            E0 = initial_conditions.E;
            I0 = initial_conditions.I;
            x0 = initial_conditions.x;
            y0 = initial_conditions.y;
            FTC = frequency_tuning_curve(8,1,params,s_params,e_E,e_I,E0,I0,x0,y0,E0_non_zero);
            plot(FTC,plot_styles(ii),'DisplayName',plot_legends(ii));
            legend('-DynamicLegend');
            legend('show');
        end
    end
end
warning('on','all');

%% Effect of the seed in inter-column excitatory spread
figure();
[params, s_params] = reset_parameters();


start_time = 0.01; % s 
stimuli_durations = 0.05; % how long each stimulus lasts
s = zeros(length(stimuli_durations),params.columns);


input = 4;
seeds = 1:5;
for seed = seeds % arbitrary increasing stimulus
    [initial_conditions, e_E,e_I,E0_non_zero] = get_initial_conditions_and_noise(params,seed,0);
    E0 = initial_conditions.E;
    I0 = initial_conditions.I;
    x0 = initial_conditions.x;
    y0 = initial_conditions.y;
    s(1,:) = [0,0,0,0,0,0,0,input,0,0,0,0,0,0,0];
    s = s_tuning_curve(s,s_params);

    [t,E,I,x,y] = solve_complex_stimuli(start_time,stimuli_durations,s,params,e_E,e_I,E0,I0,x0,y0,E0_non_zero);
    % average activity per column
    E_col_avg = mean(E,3);
    
    t = 1000*t; % [ms]
    x_axis = 1:size(E_col_avg,2)+1; %columns (+1 so that the last column appears interpolated from last_column to last_column+1)
    y_axis = t;
    [X_axis,Y_axis] = meshgrid(x_axis,y_axis);
    subplot(length(seeds)+1,1,seed);
    grid off;
    surf(X_axis-0.5,Y_axis,cat(2,E_col_avg,zeros(size(t,1),1)),'edgecolor','none');
    xlim([0.5,15.5]);
    %ylim([0,start_time+stimuli_durations(1)]);
    xticks([]);
    title(sprintf('seed = %d',seed));
    colorbar;
    view(2);
    hold on;
end
subplot(length(seeds)+1,1,length(seeds)+1);
plot(1:params.columns,s);
xlim([0.5,15.5]);
xticks(1:params.columns);
yticks([2,4,6]);
xlabel('Column index');
ylabel('Input amp.');
%% Effect of reducing number of neurons per column
figure();
[params, s_params] = reset_parameters();
cells_per_column_vec = [100,50,25,10];
for ii = 1:length(cells_per_column_vec)
    
    params.cells_per_column = cells_per_column_vec(ii);
    params.N_E = floor(params.cells_per_column/2);
    params.N_I = params.cells_per_column - params.N_E;

    [initial_conditions, e_E,e_I,E0_non_zero] = get_initial_conditions_and_noise(params,'default',1);
    E0 = initial_conditions.E;
    I0 = initial_conditions.I;
    x0 = initial_conditions.x;
    y0 = initial_conditions.y;

    FTC = frequency_tuning_curve(8,1,params,s_params,e_E,e_I,E0,I0,x0,y0,E0_non_zero);
    plot(FTC,'DisplayName',sprintf('Cells per column = %d',params.cells_per_column));
    hold on;
end
legend()

%% Effect of changing N_E/N_I ratio
figure();
[params, s_params] = reset_parameters();
%N_E_vec = [0,0.2,0.5,0.7,1];
N_E_vec = [0.2,0.5,0.7];
for ii = 1:length(N_E_vec)
    params.N_E = floor(params.cells_per_column*N_E_vec(ii));
    params.N_I = params.cells_per_column - params.N_E;

    [initial_conditions, e_E,e_I,E0_non_zero] = get_initial_conditions_and_noise(params);
    E0 = initial_conditions.E;
    I0 = initial_conditions.I;
    x0 = initial_conditions.x;
    y0 = initial_conditions.y;

    FTC = frequency_tuning_curve(8,1,params,s_params,e_E,e_I,E0,I0,x0,y0,E0_non_zero);
    plot(FTC,'DisplayName',sprintf('N_E = %d',params.N_E));
    hold on;
end
legend()

%% Functions
function plot_article_figures(fig,params,s_params,e_E,e_I,X0,E0_non_zero)
E0 = X0.E;
I0 = X0.I;
x0 = X0.x;
y0 = X0.y;

if fig == -1
   for i = 2:8
       plot_article_figures(i,params,s_params,e_E,e_I,X0,E0_non_zero);
   end
end

%% Figure 2
if fig == 2
figure();
start_time = 0.1; % s 
stimuli_durations = [0.4]; % how long each stimulus lasts
s = zeros(max(size(stimuli_durations)),params.columns);
s(1,:) = [0,0,0,0,0,0,0,4,0,0,0,0,0,0,0];
s = s_tuning_curve(s,s_params);

[t,E,I,x,y] = solve_complex_stimuli(start_time,stimuli_durations,s,params,e_E,e_I,E0,I0,x0,y0,E0_non_zero);
% plot average activity per column
E_col_avg = mean(E,3);

t = 1000*(t-start_time); % scale to milliseconds, and offset so that the first stimulus is at t=0
max_amplitude = max(E_col_avg,[],'all');
for col =6:10
    subplot(3,5,col-5);
    plot(t,E_col_avg(:,col),'k','DisplayName',sprintf('col %d',col));
    ylim([0,max_amplitude]);
    
    subplot(3,5,col);
    plot(t,s(col)*heaviside(t),'k'); % The temporal envelope zeta used, is the heaviside function
    ylim([0,max(s)+1]);
    
end
subplot(3,5,[11,12,13,14,15]);
plot(6:10,s(6:10),'k');
xticks([6,7,8,9,10]);
xlabel('Column index');

%% Figure 3
elseif fig == 3
figure();
start_time = 0.025; % s 
stimuli_durations = [0.075]; % how long each stimulus lasts
s = zeros(max(size(stimuli_durations)),params.columns);

min_input = 2;
max_input = 8;
for ii = min_input:max_input % arbitrary increasing stimulus
    s(1,:) = [0,0,0,0,0,0,0,ii,0,0,0,0,0,0,0];
    s = s_tuning_curve(s,s_params);

    [t,E,I,x,y] = solve_complex_stimuli(start_time,stimuli_durations,s,params,e_E,e_I,E0,I0,x0,y0,E0_non_zero);
    % average activity per column
    E_col_avg = mean(E,3);
    
    x_axis = 1:size(E_col_avg,2)+1; %columns (+1 so that the last column appears interpolated from last_column to last_column+1)
    y_axis = t;
    [X_axis,Y_axis] = meshgrid(x_axis,y_axis);
    subplot(max_input+1,1,max_input-ii+1);
    grid off;
    surf(X_axis-0.5,Y_axis,cat(2,E_col_avg,zeros(size(t,1),1)),'edgecolor','none');
    xlim([0.5,15.5]);
    ylim([0,start_time+stimuli_durations(1)]);
    xticks([]);
    colorbar off;
    view(2);
    if (mod(ii,2) == 0) & (ii < 8)
        subplot(max_input+1,1,max_input+1);
        plot(1:params.columns,s);
        xlim([0.5,15.5]);
        xticks([1:params.columns]);
        yticks([2,4,6]);
    end
    hold on;
end

%% Figure 4
elseif fig == 4
figure();
start_time = 0; % s 
stimuli_durations = [0.05]; % how long each stimulus lasts

max_amp = 10;
step = 0.2;
plot_colors = ['k','r','g'];
color_cont = 1;
for tone_freq = 8:2:12
    t_peak_vec = zeros(1,1/step*(max_amp-1))*nan;
    cont = 1;
    for tone_amp = 1:step:max_amp
        s = zeros(max(size(stimuli_durations)),params.columns);
        s(1,tone_freq) = tone_amp;
        s = s_tuning_curve(s,s_params);
        [t,E,I,~,~] = solve_complex_stimuli(start_time,stimuli_durations,s,params,e_E,e_I,E0,I0,x0,y0,E0_non_zero);
        % average activity per column
        E_col_avg = mean(E,3);
        I_col_avg = mean(I,3);
        
        E_col8 = E_col_avg(:,8);
        I_col8 = I_col_avg(:,8);
        
        [~,E_indices]=findpeaks(E_col8,'MinPeakHeight', 15); %Completely arbitrary threshold O_O
        [~,I_indices]=findpeaks(I_col8,'MinPeakHeight', 15);         
        time_E_peak = t(E_indices);
        time_I_peak = t(I_indices);
        
        if ~isempty(time_E_peak)
            t_peak_vec(cont) = time_E_peak;
        end
        cont = cont+1;
    end
    hold on;
    subplot(3,1,3);
    plot(1:step:max_amp,1000*t_peak_vec,'Color',plot_colors(color_cont));
    color_cont = color_cont +1;
end
xlabel('Input Amplitude [Hz]');
ylabel('Latency to response [ms]');
legend('Stimulus at column 8','Stimulus at column 10','Stimulus at column 12');

ylim([0,40]);
xlim([0,10]);

subplot(3,1,1);
s = zeros(max(size(stimuli_durations)),params.columns);
s(1,8) = 2; % Weak input
s = s_tuning_curve(s,s_params);
[t,E,I,~,~] = solve_complex_stimuli(start_time,stimuli_durations,s,params,e_E,e_I,E0,I0,x0,y0,E0_non_zero);
% average activity per column
E_col_avg = mean(E,3);
I_col_avg = mean(I,3);
E_col8 = E_col_avg(:,8);
I_col8 = I_col_avg(:,8);
plot(1000*t,E_col8,'k','DisplayName','Excitatory population');
hold on;
plot(1000*t,I_col8,'k--','DisplayName','Inhibitory population');
legend();

subplot(3,1,2);
s = zeros(max(size(stimuli_durations)),params.columns);
s(1,8) = 8; % Strong input
s = s_tuning_curve(s,s_params);
[t,E,I,~,~] = solve_complex_stimuli(start_time,stimuli_durations,s,params,e_E,e_I,E0,I0,x0,y0,E0_non_zero);
% average activity per column
E_col_avg = mean(E,3);
I_col_avg = mean(I,3);
E_col8 = E_col_avg(:,8);
I_col8 = I_col_avg(:,8);
plot(1000*t,E_col8,'k','DisplayName','Excitatory population');
hold on;
plot(1000*t,I_col8,'k--','DisplayName','Inhibitory population');
legend();
xlabel('time [ms]');

%% Figure 5
elseif fig == 5
figure();
%Subplot A of the article:
[params,s_params] = reset_parameters();
subplot(2,2,1);
hold on;
FTC = frequency_tuning_curve(8,1,params,s_params,e_E,e_I,E0,I0,x0,y0,E0_non_zero);
plot(FTC,'k');
params.J_EE = [0,0,params.J_EE(3),0,0];
params.J_IE = [0,0,params.J_IE(3),0,0];
FTC = frequency_tuning_curve(8,1,params,s_params,e_E,e_I,E0,I0,x0,y0,E0_non_zero);
plot(FTC,'k--');
legend('with inter-column connections', 'without inter-column connections');
ylim([0,9]);
xlim([1,15]);
ylabel('Input amplitude [Hz]');


% Subplot B of the article:
[params,s_params] = reset_parameters();
subplot(2,2,3);
hold on;
delta_left_vec = [5,8,16,32];
line_style_vec = ['k','r','g','b'];
for ii = 1:size(delta_left_vec,2)
    s_params.delta_left = delta_left_vec(ii);
    FTC = frequency_tuning_curve(8,1,params,s_params,e_E,e_I,E0,I0,x0,y0,E0_non_zero);    
    plot(FTC,'Color',line_style_vec(ii),'DisplayName',strcat('\delta_{left}',sprintf('%d',delta_left_vec(ii)))); 
    legend('-DynamicLegend');
    legend('show');
end
ylim([0,9]);
xlim([1,15]);
xlabel('Column index');
ylabel('Input amplitude [Hz]');

% Subplot C of the article:
[params,s_params] = reset_parameters();
subplot(2,2,2);
hold on;
FTC = frequency_tuning_curve(8,1,params,s_params,e_E,e_I,E0,I0,x0,y0,E0_non_zero);
plot(FTC,'k');
params.J_II = params.J_II*2;
params.J_EI = params.J_EI*2;
FTC = frequency_tuning_curve(8,1,params,s_params,e_E,e_I,E0,I0,x0,y0,E0_non_zero);
plot(FTC,'k--');
legend('Normal settings', 'increased recurrent inh.');
ylim([0,9]);
xlim([1,15]);
% 
% Subplot D of the article:
[params,s_params] = reset_parameters();
delta_e = 1;
subplot(2,2,4);
hold on;
e_E_inc = e_E + delta_e/2;
e_I_inc = e_I + delta_e/2;
FTC = frequency_tuning_curve(8,1,params,s_params,e_E_inc,e_I_inc,E0,I0,x0,y0,E0_non_zero);
plot(FTC,'k');
e_E_dec = e_E - 2*delta_e;
e_I_dec = e_I - 2*delta_e;
FTC = frequency_tuning_curve(8,1,params,s_params,e_E_dec,e_I_dec,E0,I0,x0,y0,E0_non_zero);
plot(FTC,'k--');
legend('increased background input', 'decreased background input');
ylim([0,9]);
xlim([1,15]);
xlabel('Column index');

%% test difference in time between frequency_tuning_curve and the slow version
% warning('off','all');
% FTC1 = frequency_tuning_curve(8,params,s_params,e_E,e_I,E0,I0,x0,y0,E0_non_zero);
% FTC2 = frequency_tuning_curve2(8,1,params,s_params,e_E,e_I,E0,I0,x0,y0,E0_non_zero);
% subplot(1,2,1);
% plot(FTC1,'k');
% ylim([0,9]);
% xlim([1,15]);
% hold on;
% subplot(1,2,2);
% plot(FTC2,'k');
% ylim([0,9]);
% xlim([1,15]);
% warning('on','all');

%% Figure 6
elseif fig == 6
% Subplot A
fig6_part1 = figure();
start_time = 0.01; % s 
stimulus_time = 0.05;
end_time = start_time + 5*params.Tau_rec; % I want all simulations to end at the same time so that the plot looks nice
ISI_vec = [0.25,0.5,1,2,4];% inter stimulus vector
for ii = 1:length(ISI_vec)
    hold on;
    ISI = ISI_vec(ii)*params.Tau_rec;
    stimuli_durations = [stimulus_time,ISI,end_time-ISI]; % how long each stimulus lasts
    s = zeros(max(size(stimuli_durations)),params.columns);
    s(1,8) = 9; 
    s(2,8) = 0; 
    s(3,8) = 9; 
    s = s_tuning_curve(s,s_params);
    [t,E,I,x,y] = solve_complex_stimuli(start_time,stimuli_durations,s,params,e_E,e_I,E0,I0,x0,y0,E0_non_zero);
    E_col_avg = mean(E,3);
    E_col8 = E_col_avg(:,8);
    subplot(length(ISI_vec),1,ii);
    plot((t-start_time)/params.Tau_rec,E_col8,'k');
    xticks([0,0.5,1,2,4]);
    xlim([-params.Tau_rec/2,6*params.Tau_rec]);

end
superplot = axes(fig6_part1,'visible','off');
superplot.Title.Visible='on';
superplot.XLabel.Visible='on';
superplot.YLabel.Visible='on';
xlabel(superplot,'Inter-stimulus interval [\tau_{rec}]');
ylabel(superplot,'Network activity [Hz]');

% Subplot B
fig6_part2 = figure();
subplot(2,1,1);
ISI_vec = [0.25,0.5,1,2,4];% inter stimulus vector
amplitudes = 4:2:10;
P2P1_vec = zeros(max(amplitudes),length(ISI_vec))*nan;
for amp = amplitudes % Plot the ratio P2/P1 for several stimuli amplitudes
    for ii = 1:length(ISI_vec)
        ISI = ISI_vec(ii)*params.Tau_rec;
        stimuli_durations = [stimulus_time,ISI,stimulus_time]; % how long each stimulus lasts
        s = zeros(max(size(stimuli_durations)),params.columns);
        s(1,8) = amp; 
        s(2,8) = 0; 
        s(3,8) = amp; 
        s = s_tuning_curve(s,s_params);
        [t,E,I,x,y] = solve_complex_stimuli(start_time,stimuli_durations,s,params,e_E,e_I,E0,I0,x0,y0,E0_non_zero);
        E_col_avg = mean(E,3);
        E_col8 = E_col_avg(:,8);
        % Find both peaks and compute P2/P1
        [~,E_indices]=findpeaks(E_col8,'MinPeakHeight', 15); %Completely arbitrary threshold O_O
        if size(E_indices) < 2
            break
        end
        P2P1 = E_col8(E_indices(2))/E_col8(E_indices(1));
        P2P1_vec(amp,ii) = P2P1;        
    end
    hold on;
    plot(ISI_vec,P2P1_vec(amp,:),'DisplayName',strcat('Input Amplitude = ',sprintf('%d [Hz]',amp)));
end
legend();
xticks([0.5,1,2,4]);
xlabel('Inter-stimulus interval [\tau_{rec}]');
ylabel('Peak 2 / Peak 1 ratio');

% Subplot C
subplot(2,1,2);
cols = 4:8;
P2P1_vec = zeros(max(cols),length(ISI_vec))*nan;
for col = cols % Plot the ratio P2/P1 for response of column 8 to other columns being stimulated
    for ii = 1:length(ISI_vec)
        ISI = ISI_vec(ii)*params.Tau_rec;
        stimuli_durations = [stimulus_time,ISI,stimulus_time]; % how long each stimulus lasts
        s = zeros(max(size(stimuli_durations)),params.columns);
        s(1,col) = 20; 
        s(2,col) = 0; 
        s(3,col) = 20; 
        s = s_tuning_curve(s,s_params);
        [t,E,I,x,y] = solve_complex_stimuli(start_time,stimuli_durations,s,params,e_E,e_I,E0,I0,x0,y0,E0_non_zero);
        E_col_avg = mean(E,3);
        E_col8 = E_col_avg(:,8);
        % Find both peaks and compute P2/P1
        [~,E_indices]=findpeaks(E_col8,'MinPeakHeight', 15); %Completely arbitrary threshold O_O
        if size(E_indices) < 2
            P2P1_vec(col,:) = nan;       
            break
        end
        P2P1 = E_col8(E_indices(2))/E_col8(E_indices(1));
        P2P1_vec(col,ii) = P2P1;        
    end
    hold on;
    plot(ISI_vec,P2P1_vec(col,:),'DisplayName',strcat('Column stimulated = ',sprintf('%d [Hz]',col)));
end
legend();
xticks([0.5,1,2,4]);
xlabel('Inter-stimulus interval [\tau_{rec}]');
ylabel('Peak 2 / Peak 1 ratio (In column 8)');


%% Figure 7 (Issy Code)
elseif fig == 7
figure();
start_time = 0.025; % seconds to wait before first stimulus
time_at_end = 0.05; % seconds to wait after last stimulus

interstimulus_intervals = params.Tau_rec.*[1/8, 1/4, 1/2, 1, 2, 4]; % the sequence of interstimulus intervals to try

for i = 1:length(interstimulus_intervals)
    
    isi = interstimulus_intervals(i);
    stimuli_durations = [0.05, isi, 0.05, time_at_end]; % how long each stimulus lasts (i.e. there are 4 phases of stimulus). Each tone is provided for 0.05 seconds as per the caption in Figure 6A ("50 ms duration, amplitude of 10.."). The time between the tones is the ISI (interstimulus interval). At the end have a bit of time to see how the second tone affects activity.
    s = zeros(max(size(stimuli_durations)),params.columns); % empty stimulus matrix
    s(1,:) = [0,0,0,10,0,0,0,0,0,0,0,0,0,0,0]; % first phase is tone of amplitude 10 at the 4th column, 0 everywhere else (First tone, "Mask")
    s(2,:) = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]; % second phase is 0 everywhere (Interstimulus Interval)
    s(3,:) = [0,0,0,10,0,0,0,0,0,0,0,0,0,0,0]; % third phase is tone of amplitude 10 at the 4th column, 0 everywhere else (Second tone, "Probe")
    s(4,:) = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]; % fourth phase is 0 everywhere

    s = s_tuning_curve(s,s_params); % spread the pure tone across all the columns exponentially
 
    [t,E,I,x,y] = solve_complex_stimuli(start_time,stimuli_durations,s,params,e_E,e_I,E0,I0,x0,y0,E0_non_zero); % Run the ODE
    E_col_avg = mean(E,3); % average activity per column
    x_axis = 1:size(E_col_avg,2)+1; %columns (+1 so that the last column appears interpolated from last_column to last_column+1)
    y_axis = t;
    [X_axis,Y_axis] = meshgrid(x_axis,y_axis); % distribute the results evenly
    subplot(7,1,7-i); % fill in the (7-i)-th row of the subplot
    grid off;
    surf(X_axis-0.5,Y_axis,cat(2,E_col_avg,zeros(size(t,1),1)),'edgecolor','none'); % draw the surface (i.e. colour) over
    xlim([0.5,15.5]); 
    ylim([stimuli_durations(1) + stimuli_durations(2),start_time+sum(stimuli_durations)]); % only show the results in the last two periods (i.e. after the ISI)
    xticks([]);
    colorbar;
    caxis([0 70]) % max colour scaled to 70Hz as in paper
    colormap('jet')
    title(['ISI = ', num2str(isi), 's'])
    view(2); % show the color map in 2 dimensions
    %subplot(max_input+1,1,max_input+1);
    %plot(1:columns,s);
    %xlim([0.5,15.5]);
    %xticks([1:columns]);
    hold on;
end

subplot(7,1,7); % plot the initial tone response at the very bottom
grid off;
surf(X_axis-0.5,Y_axis,cat(2,E_col_avg,zeros(size(t,1),1)),'edgecolor','none'); 
xlim([0.5,15.5]);
ylim([0,start_time+stimuli_durations(1)+time_at_end]); % only show from the start to after the first tone to show how the A1 responds to the initial tone
xticks(1:params.columns);
xlabel('Column Index');
colorbar;
caxis([0 70])
colormap('jet')
view(2);


%% Figure 8
elseif fig == 8
figure();
[params,s_params] = reset_parameters();
% Response to the first input
FTC = frequency_tuning_curve(8,1,params,s_params,e_E,e_I,E0,I0,x0,y0,E0_non_zero);
plot(FTC,'k','LineWidth',2);
% Response to several inter stimulus intervals
ISI_vec = params.Tau_rec.*[0.5,1,2];% inter stimulus vector
hold on;
for ii = ISI_vec
    FTC = frequency_tuning_curve(8,2,params,s_params,e_E,e_I,E0,I0,x0,y0,E0_non_zero,ii);
    plot(FTC);
end
legend('Response to first input','ISI of \tau_{rec}/2','ISI of \tau_{rec}','ISI of 2*\tau_{rec}');
ylim([0,9]);
xlim([1,15]);
ylabel('Input amplitude [Hz]');
xlabel('Column index');

else 
    disp(strcat('No figure ',sprintf('%d',fig),' to plot'));
end
end
function [t_vec,E,I,x,y] = solve_complex_stimuli(start_time,dur_vec,s_matrix,params,e_E,e_I,E0,I0,x0,y0,E0_non_zero)
% Solve the auditory cortex system of ODEs for an arbitrary train of
% stimuli (s_matrix).
%
% s_matrix(onset_times, columns): Is a matrix with as many columns as
% auditory cortex columns, and as many rows as different stimulus onsets.
% dur_vec: vector of duration for each stimulus
   
    %% Solve ODE
    
    % if the starting time is different to 0, add a stimulus of 0 from t=0
    % to until the starting time
    dur_vec = dur_vec(:);
    if start_time ~= 0
        zero_stimulus = zeros(1,size(s_matrix,2));
        s_matrix = cat(1,zero_stimulus,s_matrix);
        dur_vec = cat(1,start_time,dur_vec);
    end
    
    time_sum = 0;
    E = [];
    I = [];
    x = [];
    y = [];
    t_vec = [];
    next_E0 = E0;
    next_I0 = I0;
    next_x0 = x0;
    next_y0 = y0;
    for jj = 1:size(dur_vec,1)
        t_span = [0, dur_vec(jj)];
        s = s_matrix(jj,:);
        X0 = [next_E0(:); next_I0(:); next_x0(:); next_y0(:)];
        [t, X] = ode45(@(t,X) dXdt(t,X,params,s,e_E,e_I,E0_non_zero),t_span,X0);

        t = t(:) + time_sum;
        time_sum = time_sum + dur_vec(jj);
        [E_jj,I_jj,x_jj,y_jj] = reshape_X(X,[1,3,2], params);% E(time,cols,N_E)
        
        E = cat(1,E,E_jj);
        I = cat(1,I,I_jj);
        x = cat(1,x,x_jj);
        y = cat(1,y,y_jj);
        t_vec = cat(1,t_vec,t);

        % the last E,I,x,y of each stimulus, are the initial conditions for
        % next one:
        [E0_jj,I0_jj,x0_jj,y0_jj] = reshape_X(X,[2,3,1], params);
        next_E0 = E0_jj(:,:,end);
        next_I0 = I0_jj(:,:,end);
        next_x0 = x0_jj(:,:,end);
        next_y0 = y0_jj(:,:,end);
           
    end
end
function dXdt = dXdt(t,X,params,s,e_E,e_I,E0_non_zero)
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
    Tau_rec = params.Tau_rec;
    
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
    relu = max((between_col_exci + in_col_inhi + e_E + s.*E0_non_zero),0);

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
function s = s_tuning_curve(s,s_params)
    columns = size(s,2);
    alpha = s_params.alpha;
    lambda_c = s_params.lambda_c;
    delta_left = s_params.delta_left;
    delta_right = s_params.delta_right;
    for tt = 1:size(s,1)
        A_vec = s(tt,:);      
        h = zeros(columns);
        for M = 1:size(s,2)
            for Q = 1:size(s,2)
                A = A_vec(M); 
                if Q > M
                    delta = delta_right;
                else
                    delta = delta_left;
                end
                lambda_s = ((A-alpha)/delta)*heaviside(A-alpha-eps) + lambda_c;
                h(M,Q) = A*exp(-abs(Q-M)/lambda_s);
            end
        end
        s(tt,:) = sum(h,1);
    end
end
function FTC_slow = frequency_tuning_curve_slow(col_to_tune, params,s_params,e_E,e_I,E0,I0,x0,y0,E0_non_zero)
    columns = params.columns;
    max_input = 9; % maximum input stimulus to simulate
    step_size = 0.2;
    start_time = 0; % s 
    stimuli_durations = 0.075; 

    FTC_slow = zeros(columns,1);
    for column = 1:columns
        for input = 0:step_size:max_input
            s = zeros(max(size(stimuli_durations)),columns);
            s(1,column) = input;
            s = s_tuning_curve(s,s_params);
            [t,E,I,x,y] = solve_complex_stimuli(start_time,stimuli_durations,s,params,e_E,e_I,E0,I0,x0,y0,E0_non_zero);    
            E_col_avg = mean(E,3);
            E_col_to_tune = E_col_avg(:,col_to_tune);
            [~,E_indices]=findpeaks(E_col_to_tune,'MinPeakHeight', 15); %Completely arbitrary threshold O_O

            % Find if the stimulus presented caused a activation in column 8
            if size(E_indices)>=1
               FTC_slow(column) = input;
               break
            else
               FTC_slow(column) = nan;
            end
        end
    end
end
function FTC = frequency_tuning_curve(col_to_tune,peak_number,params,s_params,e_E,e_I,E0,I0,x0,y0,E0_non_zero,optional_Tau_res)
% Computes the frequency tuning curve for a particular column, taking into account
% the known nature of a FCT (decreasing for columns < col_to_tune, and
% increasing for columns > col_to_tune)
    warning('off','all'); % To avoid getting a matlab warning when findpeaks function doesn't find any peaks
    columns = params.columns;
    max_input = 9; % maximum input stimulus to simulate
    step_size = 0.2;
    
    start_time = 0; % s 
    if nargin == 12 % This means that the optional Tau_res variable was given as input to the function
        rest_time = optional_Tau_res;
        stimuli_durations = ones(1,2*peak_number - 1)*0.075;
        stimuli_durations(2:2:end) = rest_time;
    else
        stimuli_durations = 0.075;
    end
    
    FTC = zeros(columns,1)*nan;   
    % First find the input threshold for the column to tune (should be the
    % lowest threshold of the entire system
    for input = step_size:step_size:max_input
        s = zeros(length(stimuli_durations),columns);
        s(1:2:end,col_to_tune) = input;
        s = s_tuning_curve(s,s_params);
        [~,E,~,~,~] = solve_complex_stimuli(start_time,stimuli_durations,s,params,e_E,e_I,E0,I0,x0,y0,E0_non_zero);    
        E_col_avg = mean(E,3);
        E_col_tuned = E_col_avg(:,col_to_tune);
        [~,E_indices]=findpeaks(E_col_tuned,'MinPeakHeight', 15); %Completely arbitrary threshold O_O

        % Find if the stimulus presented caused a activation in column to tune
        if length(E_indices)>=peak_number
           FTC(col_to_tune) = input;
           break
        end
    end
    
    % one loop from the column to tune down to column 1, the next one from
    % column to tune up to the end
    min_col_vec = [col_to_tune-1,col_to_tune+1];
    col_steps_vec = [-1,1];
    max_col_vec = [1,columns];
    
    for ii = 1:2
        x1 = col_to_tune;
        y1 = FTC(col_to_tune);
        [slope,intercept] = find_linear_params(x1,y1,x1-1,y1);
        for column = min_col_vec(ii):col_steps_vec(ii):max_col_vec(ii) % Sweep from col_to_tune to extreme columns   
            input = slope*column + intercept;
            % See if there is activation at this level, if there is, then start
            % decreasing amplitude to fin minimum, if not, start increasing
            % amplitude to find minimum activation amplitude.
            s = zeros(max(size(stimuli_durations)),columns);
            s(1:2:end,column) = input;
            s = s_tuning_curve(s,s_params);
            [t,E,I,x,y] = solve_complex_stimuli(start_time,stimuli_durations,s,params,e_E,e_I,E0,I0,x0,y0,E0_non_zero);    
            E_col_avg = mean(E,3);
            E_col_tuned = E_col_avg(:,col_to_tune);
            [~,E_indices]=findpeaks(E_col_tuned,'MinPeakHeight', 15); %Completely arbitrary threshold O_O
            % Find if the stimulus presented caused an activation in column to tune
            if length(E_indices)>= peak_number % Go downwards
               FTC(column) = input;
               input_stop = 0; % IMPORTANT: input_stop must be completely divisible by step i.e mod(input_stop,step) == 0
               step_sign = -1;
               logic_aux = true; 
            else % Go upwards
               FTC(column) = nan;
               input_stop = 12; % IMPORTANT: input_stop must be completely divisible by step i.e mod(input_stop,step) == 0
               step_sign = 1;
               logic_aux = false;
            end

            while xor(isnan(FTC(column)), logic_aux) && (input ~= input_stop)
                input = input + step_size*step_sign; 
                s = zeros(max(size(stimuli_durations)),columns);
                s(1:2:end,column) = input;
                s = s_tuning_curve(s,s_params);
                [t,E,I,x,y] = solve_complex_stimuli(start_time,stimuli_durations,s,params,e_E,e_I,E0,I0,x0,y0,E0_non_zero);    
                E_col_avg = mean(E,3);
                E_col_tuned = E_col_avg(:,col_to_tune);
                [~,E_indices]=findpeaks(E_col_tuned,'MinPeakHeight', 15); %Completely arbitrary threshold O_O

                % Find if the stimulus presented caused activation
                if length(E_indices) >= peak_number
                   FTC(column) = input;
                else
                    if logic_aux % If FTC is NaN and we are going downwards, then switch logic_aux to get out of the while loop
                        logic_aux = ~logic_aux;
                    end
                end            
            end
            x2 = column;
            y2 = FTC(column);
            x1 = column - col_steps_vec(ii);
            y1 = FTC(column - col_steps_vec(ii));
            [slope,intercept] = find_linear_params(x1,y1,x2,y2);
        end
    end
    warning('on','all'); % To avoid getting a matlab warning when findpeaks function doesn't find any peaks
end
function [slope, intercept] = find_linear_params(x1,y1,x2,y2)
    slope = (y2 - y1)/(x2 - x1);
    intercept = (y1*x2 - y2*x1)/(x2-x1);
end
function [initial_conditions,e_E,e_I,E0_non_zero] = get_initial_conditions_and_noise(params,seed,stimulate_all)
    %% 1.2. External input
    % Create background noise
    % seed: Random seed to generate external noise
    % stimulate_all: binary (0-1) variable that defines if all neurons
    % will receive stimulus, or onlty the naturally active ones
    
    rng(seed);
    init_e_E = rand(params.N_E,1); %between (0,1)
    init_e_I = rand(params.N_I,1);
    % scale
    highest_eE = max(init_e_E);
    lowest_eE = min(init_e_E);
    highest_eI = max(init_e_I);
    lowest_eI = min(init_e_I);
    e_E_range = max(init_e_E) - min(init_e_E);
    e_I_range = max(init_e_I) - min(init_e_I);

    actual_eE = (1/e_E_range)*(init_e_E*(params.e_max-params.e_min)+(params.e_min*highest_eE-params.e_max*lowest_eE));
    actual_eI = (1/e_I_range)*(init_e_I*(params.e_max-params.e_min)+(params.e_min*highest_eI-params.e_max*lowest_eI));

    e_E_col = sort(actual_eE);
    e_I_col = sort(actual_eI);
    
    e_E = repmat(e_E_col,1,params.columns);
    e_I = repmat(e_I_col,1,params.columns);

%     e_E = repmat(linspace(params.e_min,params.e_max,params.N_E)',1,params.columns);
%     e_I = repmat(linspace(params.e_min,params.e_max,params.N_I)',1,params.columns);

    %% 1.3. Create neurons
    E_ini = zeros(params.N_E,params.columns);
    I_ini = zeros(params.N_I,params.columns);
    x_ini = zeros(size(E_ini));
    y_ini = zeros(size(I_ini));

    %% 1.4. Get initial conditions
    % set stimulus to 0 and let the system converge
    X0 = [E_ini(:); I_ini(:); x_ini(:); y_ini(:)];
    t_rest = 5*params.Tau_rec;
    E0_non_zero = ones(size(E_ini));
    [t, X] = ode45(@(t,X) dXdt(t,X,params,0,e_E,e_I,E0_non_zero),[0,t_rest],X0,E0_non_zero);

    [E,I,x,y] = reshape_X(X,[2,3,1], params);

    E0 = E(:,:,end);
    I0 = I(:,:,end);
    x0 = x(:,:,end);
    y0 = y(:,:,end);

    initial_conditions.E = E0;
    initial_conditions.I = I0;
    initial_conditions.y = y0;
    initial_conditions.x = x0;
    % Decide if you want to stimulate only neurons that are naturally excited
    % or all neurons in the network.
    E0_non_zero = ones(size(initial_conditions.E)) - (E0 <= 0)*(1-stimulate_all);
end
function [params, s_params] = reset_parameters()
    columns = 15;
    cells_per_column = 200;
    N_E = 100;
    N_I = cells_per_column - N_E;

    Tau_E = 0.001; % [s]
    Tau_I = 0.001; % [s]

    Tau_ref_E = 3*10^(-3); % [s]
    Tau_ref_I = 3*10^(-3); % [s]

    J_EE = [0.015,0.045,6,0.045,0.015]; % index goes from -2 to 2
    J_IE = [0.0015,0.0035,0.5,0.0035,0.0015];
    J_EI = -4; 
    J_II = -0.5;

    % synaptic resources
    U = 0.5;
    Tau_rec = 0.8; % [s]

    % background activity
    e_max = 10; % [Hz]
    e_min = -10; % [Hz]

    % stimulus
    alpha = 2;
    lambda_c = 0.25;
    delta_left = 5;
    delta_right= 5;

    %% 1.1. Compress all parameters into a couple of variables
    params.cells_per_column = cells_per_column;
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
    params.Tau_rec = Tau_rec;
    params.e_min = e_min;
    params.e_max = e_max;

    s_params.alpha = alpha;
    s_params.lambda_c = lambda_c;
    s_params.delta_left = delta_left;
    s_params.delta_right= delta_right;
end