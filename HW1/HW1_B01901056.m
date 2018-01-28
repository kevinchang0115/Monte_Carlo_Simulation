% Author:      B01901056 ±i³Í´Q
% Date:        2016/3/4
% Purpose:     To examine the uniformity of the random number generator
% Description: 
%     This program contains both parts of the homework. After running the
%     program, we'll receive six graphs. There are five graphs for problem
%     (a), the first to the fifth run respectively, and one graph for 
%     problem (b), the means and standard deviations of the previous five
%     runs. It can be seen that the distribution of every interval for all
%     five runs are pretty uniform. It's also obvious that the mean of each
%     run is about 500 and the standard deviation is less than 36.
% Variables:
%     mean_array: the array of means of five distributions
%     std_array:  the array of standard deviations of five distributions
%     num_array:  the array of the numbers in the 20 intervals of [0,1]
%     interval:   the interval in which the generated number is
            
% Initialization
mean_array = zeros(1,5);
std_array = zeros(1,5);
% Run the random number generator for five times 
for i = 1:5   
    num_array = zeros(1,20);
    % generate 10000 numbers for each run
    for j = 1:10000
        % figure out which interval it is in and increase the corresponding
        % value in 'num_array' by 1
        interval = ceil(rand*20);    
        num_array(interval) = num_array(interval) + 1;
    end
    % plot the bar graph of each run
    subplot(2,3,i);
    bar([0.025:0.05:0.975], num_array);
    title(['Distribution ', num2str(i)]);
    xlabel('20 intervals in [0,1]');
    ylabel('number of occurrence');
    % record the mean and std of each run
    mean_array(i) = mean(num_array);
    std_array(i) = std(num_array);
end
% plot the graph of means and stds of 5 runs
subplot(2,3,6);
bar(mean_array);
errorbar([1:5], mean_array, std_array, '.');
title('Means and stds of 5 distributions');
xlabel('5 distributions')
ylabel('number of occurrence');