%% Code to run the Lotka–Volterra Model Class

clc; clear;

G = [log(2)/(25/60); log(2)/(40/60)];
K = [1e9; 1e9];
B = [1,-2];           % [b12 b21]
N0 = 1;

model = LV_model(G, K, B, N0);

%model.simulate_logistic_one_species_twofigs_repeat( ...
   % G(1), K(1), N0, 24, 0.2, 'solo_fighter_logistic', 30, 10,10 );

model = model.run_lotka_hours(24, 1, true);   % simulate to 100 h and plot

%model.plot_populations_equalratio();

%model.animate_populations_equalratio();


%clear classes; rehash;  % optional but helps reload the class

%obj.animate_equalratio_movingdot;  % uses defaults

% or with custom settings:
%model.animate_equalratio_twofigs_repeat('fightclub_two_pics_repeat', 30, 20, 10);
%model.animate_equalratio_movingdot('fightclub_equalratio_dot.mp4', 30, 100);
%model.animate_equalratio_movingdot_repeat('final_cut_mfc_4.mp4', 30, 100,100);
