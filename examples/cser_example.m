%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Showcase the good convergence properties of state-space models. To do this,
% we will generate simulated data from random AR models of different orders,
% and estimate its entropy rate with data of varying length.
%
% Pedro Mediano, Feb 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run CSER once so MVGC is found and added to the path automatically
StateSpaceEntropyRate(randn([1,100]), 1);

%% Define simulation parameters
nb_trials = 5;                      % Number of repetitions for each simulation
t_vec = logspace(2, 5, 16);         % Length of each time series
q_vec = [0 3 5 9];                  % AR model order for simulated data


%% Run simulation
H = zeros([length(t_vec), nb_trials, length(q_vec)]);
for i=1:length(t_vec)
  disp(i);
  for j=1:nb_trials

    parfor (k=1:length(q_vec), 4)
      q = q_vec(k);
      T = floor(t_vec(i));

      % Simulate data
      if q > 0
        A = var_rand(1, q, 0.9, [], []);
        X = var_to_tsdata(A, 1, T);

      else
        X = randn([1, T]);
      end

      % Compute CSER
      H(i,j,k) = StateSpaceEntropyRate(X, 1);
    end
  end
end


%% Plot results
clf; hold on;
for k=1:length(q_vec)
  errorbar(t_vec, mean(H(:,:,k), 2), std(H(:,:,k), [], 2));
end
set(gca, 'Xscale', 'log');
true_h = 0.5*log(2*pi*exp(1));
plot([min(t_vec), max(t_vec)], [true_h true_h], 'k--');

