%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Showcase the good convergence properties of the CTW model. To do this,
% we will generate simulated data from both a maximum entropy memoryless (iid)
% process and James' "even process", and estimate entropy rate with time series
% of varying length.
%
% Pedro Mediano, Jun 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Define simulation parameters
nb_trials = 10;                     % Number of repetitions for each simulation
t_vec = 100:300:3000;               % Length of each time series


%% Identically independently distributed (memoryless) process
lz_iid  = zeros([nb_trials, length(t_vec)]);
ctw_iid = zeros([nb_trials, length(t_vec)]);

for i=1:nb_trials
  for j=1:length(t_vec)
    T = t_vec(j);
    X = rand([1, T]) < 0.5;
    lz_iid(i,j)  = LZ76(X)*log2(T)/T;
    ctw_iid(i,j) = CTWEntropyRate(X);
  end
end


%% Even process (James et al., 2011)
lz_even  = zeros([nb_trials, length(t_vec)]);
ctw_even = zeros([nb_trials, length(t_vec)]);

for i=1:nb_trials
  for j=1:length(t_vec)
    T = t_vec(j);
    X = SimulateEvenProcess(T);
    lz_even(i,j)  = LZ76(X)*log2(T)/T;
    ctw_even(i,j) = CTWEntropyRate(X);
  end
end


%% Plot results
subplot(211); hold on;
errorbar(t_vec, mean(ctw_iid), std(ctw_iid));
errorbar(t_vec, mean(lz_iid), std(lz_iid));
legend({'CTW', 'LZ76'}, 'AutoUpdate', 'off');
plot([t_vec(1), t_vec(end)],  [1, 1], 'k--');
ylabel('Entropy rate');
title('Memoryless process');

subplot(212); hold on;
errorbar(t_vec, mean(ctw_even), std(ctw_even));
errorbar(t_vec, mean(lz_even), std(lz_even));
legend({'CTW', 'LZ76'}, 'AutoUpdate', 'off');
plot([t_vec(1), t_vec(end)],  [2/3, 2/3], 'k--');
ylabel('Entropy rate');
title('Even process');

