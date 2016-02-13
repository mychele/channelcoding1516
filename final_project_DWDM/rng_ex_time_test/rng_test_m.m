% test rng execution time
num_samples = 32816;
repetitions = 1000;

execution_times = zeros(repetitions, 1);

for i=1:100
    tic;
    vec = randi(2, num_samples, 1) - 1;
    execution_times(i) = toc;
end

disp(mean(execution_times)*1000);


execution_times = zeros(repetitions, 1);
sigma = 2;

for i=1:100
    tic;
    vec = sigma*randn(num_samples, 1);
    execution_times(i) = toc;
end

disp(mean(execution_times)*1000);