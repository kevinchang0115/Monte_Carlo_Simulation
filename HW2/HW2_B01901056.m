% Author:      B01901056 ±i³Í´Q
% Date:        2016/3/13
% Purpose:     To use the Monte Carlo model to compute the attenuation of a
%              collimated beam propagating in an absorbing medium
% Description: 
%     This program simulates the absorbance of photons for 5 times by 
%     sampling with uniform distribution. After running this program, we
%     will get 5 plots. Each plot contains a bar graph for the simulation 
%     of 10000 photons (blue) and a line graph for Beer-Lambert Law (red).
%     It's obvious that these two graphs match with each other.
% Variables:
%     mu:        absorption coefficient of the medium (cm^-1)
%     dz:        length of each depth interval (cm)
%     num:       the number of photons in simulation
%     run:       for how many times it runs the simulation
%     depth:     the depth of the medium (cm)
%     x_axis:    intervals in the medium
%     dep_array: the number of photons absorbed in each interval
%     z:         the depth which the photon reaches (cm)
%     interval:  the interval where 'z' belongs

% Initialization
mu = 10;
dz = 0.025;
num = 10000;
run = 5;
depth = 1;
x_axis = [dz/2:dz:depth-dz/2];
% Run the simulation for five times
for i = 1:run
    dep_array = zeros(1,depth/dz);    
    % simulate 10000 photons for each run
    for j = 1:num
        % sampling with uniform distribution
        z = -log(rand) / mu;
        % if the photon passes through the medium, it doesn't count in
        % 'dep_array'
        if (z>depth) continue; end
        interval = ceil(1/dz * z);
        dep_array(interval) = dep_array(interval) + 1;
    end
    % plot the graph of each run
    subplot(2,3,i);
    bar(x_axis, dep_array); hold on;
    p = plot(x_axis, num * dz * mu * exp(-mu*x_axis), 'r');
    set(p, 'linewidth', 2); hold off;
    title(['Simulation ', num2str(i)]);
    xlabel('depth (cm)'); ylabel('number of photons absorbed');
end
