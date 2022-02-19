# Auditory_cortex
PBM project for modelling primary auditory cortex (A1) iso-frequency columns, following Loebel et al. (2007) model. (https://doi.org/10.3389/neuro.01.1.1.015.2007)

# How to use Auditory_cortex.m
The code is meant to be used for simulating an arbitrary number of columns, each containing an arbitrary (but equal among columns) number of excitatory and inhibitory neurons. The simulation outputs an activation level of each neuron (E or I, for excitatory and inhibitory neurons respectively) and its synaptic resources (x or y respectively) over time. Below, is a step by step guide to run your own simulations.

## 1. Input stimuli:
- Create a starting time variable (e.g. start_time = 0.1). This could be thought of as the delay (in seconds) before the first stimulus. 
- Create durations row vector, with as many elements as stimuli *sections*. One stimulus section, is the piece of time where there is no changes in the system's stimulus. For instance, dur = [0.2, 0.07] means at t = start_time the system is stimulated for 200ms, followed by a different stimulus of 70ms.  
- Create stimuli matrix. It should have as many rows as stimuli sections, and as many columns as the auditory cortex iso-frequency columns. the element *s(i,j)* should be the amplitude of the stimulus for column *j* in the cortex, at stimulus section *i*. e.g. using the duration vector *dur* from above, s = [0, 1, 0; 2, 1, 0] represents a system with 3 iso-frequency columns, stimulated at column 2 during the entire 270ms of simulation, with an amplitude of 1; and stimulated at column 1 for the last 70ms, with an amplitude of 2.
- Compute the tuning curves of each stimulus. As stated in equations 4-6 of the article, each column in the simulated cortex is sensitive to stimulus of the other columns as well. To apply these tuning curves on your stimuli matrix, call the function **s = s_tuning_curve(s, s_params)** where *s* is the stimuli matrix and *s_params* are the system's parameters (&alpha;, &delta;, &lambda;<sub>c</sub>, Explained in detail in the paper) 
## 2. Solve the ODEs:
Call the function **[t, E, I, x, y] = solve_complex_stimuli(start_time, dur, s, params, e_E, e_I)**. The output will consist of the following arrays:
- *t*: time vector, should go from 0 to (*start_time* + &Sigma;<sub>i</sub>*dur<sub>i</sub>*), and its size (number of elements) is automatically determined by the ODE solver (ode45).
- *E*, *I*: Activity levels for each excitatory and inhibitory neuron over time. Each variable has size (time x columns x cells per column). e.g. *E(10, 2, 50)* is the activity level of the 50th excitatory neuron from column 2, at the 10th time step of simulation.
- *x*, *y*: Synaptic resources for every neuron over time. They have the same dimensions as *E* and *I*.

### Example code for plotting a section of Figure 2 from the paper:
start_time = 0.1;  
dur = [0.4];  
s = [0,0,0,0,0,0,0,4,0,0,0,0,0,0,0];   
s = s_tuning_curve(s,s_params);  
[t,E,I,x,y] = solve_complex_stimuli(start_time, durations, s, params, e_E, e_I);  
E_col_avg = mean(E,3); % Average activity per column  
plot(t,E_col_avg(:,8));   
