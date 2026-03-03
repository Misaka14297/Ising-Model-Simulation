
L_values = [16, 32, 64]; % 系统尺寸
Tmin = 1.5;
Tmax = 3.5;
Tstep = 0.01; 
steps_per_T_base = 10000; % 低温下（T < 2.0K)增加至20000
averaging_runs = 10;
equilibration_steps = 5000; % 采样前的热平衡步骤
algorithms = {'metropolis', 'wolff'}; 

ising2d_mc_combined(L_values, Tmin, Tmax, Tstep, steps_per_T_base, averaging_runs, equilibration_steps, algorithms);