% Author:      B01901056 ±i³Í´Q
% Date:        2016/3/20
% Purpose:     To use Monte Carlo simulation to compute the fraction of
%              light absorbed, reflected, transmiited
%              (isotropic, fixed weight)
% Description: 
%     This program simulates the absorbance and scattering of photons for 5
%     times by sampling travel distance(s) and scattering angle(phi &
%     theta) with uniform distribution. After running this program, we will
%     get 6 graphs. First five are bar graphs of the ratio of aborsorbance,
%     reflectance and transmittance for the simulation of 10000 photons and
%     the last one is of the average values. The solutions are all similar
%     to what the homework guide shows, R = 0.3616 and T = 0.3565.
% Variables:
%     ua:     absorption coefficient of the medium (cm^-1)
%     us:     scattering coefficient of the medium (cm^-1)
%     ut:     ua + ut (cm^-1)
%     num:    the number of photons in simulation
%     run:    for how many times it runs the simulation
%     depth:  the depth of the medium (cm)
%     a:      the number of absorbed photons
%     r:      the number of reflected photons
%     t:      the number of transmitted photons
%     A:      the ratio of absorbed photons to all photons
%     R:      the ratio of reflected photons to all photons
%     T:      the ratio of transmitted photons to all photons
%     s:      the distance a photon travels before absorbance or scattering
%             (cm)
%     x:      x position of the photon (cm) (not used)
%     y:      y position of the photon (cm) (not used)
%     z:      the depth which the photon reaches (cm)
%     cx:     cosin of the angle between the moving dirction and x-axis
%     cy:     cosin of the angle between the moving dirction and y-axis
%     cz:     cosin of the angle between the moving dirction and z-axis
%     cx_new: new cx
%     cy_new: new cy
%     cz_new: new cz
%     phi:    azimuthal angle of the scattering photon
%     theta:  polar angle of the scattering photon

% Initialization
ua = 10; us = 90;
num = 10000;
run = 5;
depth = 0.02;
ut = ua + us;
A = zeros(1,run);
R = zeros(1,run);
T = zeros(1,run);
% Run the simulation for five times
for i = 1:run
    a = 0; r = 0; t = 0;
    % simulate 10000 photons for each run
    for j = 1:num
        % x, y are not important in this question
        % x = 0; y = 0; 
        z = 0;
        % cx, cy are not important in this question
        % cx = 0; cy = 0;
        cz = 1;
        while (1)
            % update new position
            s = -log(rand) / ut;            
            % x = x + cx * s;
            % y = y + cy * s;
            z = z + cz * s;
            % reflection 
            if z > depth
                t = t + 1; break;
            % transmittance
            elseif z < 0
                r = r + 1; break;
            % absorbance
            elseif rand < ua / ut
                a = a + 1; break;
            % scattering
            else
                % update new direction
                phi = 2 * pi * rand;
                theta = acos(2 * rand - 1);
                % cx_new, cy_new are not important in this question
                % cx_new = (sin(theta)/sqrt(1-cz^2))*(cx*cz*cos(phi)-cy*sin(phi))+cx*cos(theta);
                % cy_new = (sin(theta)/sqrt(1-cz^2))*(cy*cz*cos(phi)+cx*sin(phi))+cy*cos(theta);
                cz_new = -sin(theta)*cos(phi)*sqrt(1-cz^2)+cz*cos(theta);
                % cx = cx_new; cy = cy_new;
                cz = cz_new;
            end
        end       
    end
    % calculate the ratio
    A(i) = a / num;
    R(i) = r / num;
    T(i) = t / num;
    % plot the graph
    subplot(2,3,i);
    y = [A(i);R(i);T(i)];
    bar(y);
    set(gca,'XTickLabel',{'Absorbance', 'Reflectance', 'Transmittance'})
    title(['Simulation ', num2str(i)]);
    xlabel('A.R.T'); ylabel('ratio');
    text(1:length(y), y, num2str(y), 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom');
end
% plot the graph
subplot(2,3,run+1);
y = [mean(A);mean(R);mean(T)];
bar(y)
set(gca,'XTickLabel',{'Absorbance', 'Reflectance', 'Transmittance'});
title('Average');
xlabel('A.R.T'); ylabel('ratio');
text(1:length(y), y, num2str(y), 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom');
