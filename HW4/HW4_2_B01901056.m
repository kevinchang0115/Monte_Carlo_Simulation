% Author:      B01901056 ±i³Í´Q
% Date:        2016/4/5
% Purpose:     To use Monte Carlo simulation to compute the fraction of
%              light absorbed, reflected, transmiited
%              (anisotropic, variable weight, "mismatched boundary")
% Description: 
%     This program is similar to HW3-2 except for using isotropic
%     scattering and mismatched boundary. After running this program, we
%     will get 6 graphs. First five are bar graphs of the ratio of
%     aborsorbance, reflectance and transmittance for the simulation of
%     10000 photons and the last one is of the average values. The
%     solutions are all similar to what the homework guide shows, R = 0.26.
% Variables:
%     ua:     absorption coefficient of the medium (cm^-1)
%     us:     scattering coefficient of the medium (cm^-1)
%     ut:     ua + us (cm^-1)
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
%     z:      the depth which the photon reaches (cm)
%     cz:     cosin of the angle between the moving dirction and z-axis
%     cz_new: new cz
%     phi:    azimuthal angle of the scattering photon
%     theta:  polar angle of the scattering photon
%     w:      the weight of a photon
%     min:    if weight less than min, play roulette
%     m:      if rand larger than 1/m, terminate it
%     n1:     the refractive index of tissue (n=1.5)
%     n2:     the refractive index of air (n=1)
%     Rs:     the function of the reflectance for s-polarized light
%     Rp:     the function of the reflectance for p-polarized light
%     Re:     the function of the reflectance for unpolarized light

% Initialization
clear;
ua = 10; us = 90;
ut = ua + us;
num = 10000;
run = 5;
depth = 20/ut;
min = 0.001;
m = 20;
A = zeros(1,run);
R = zeros(1,run);
T = zeros(1,run);
n1 = 1.5;
n2 = 1;
Rs = @(theta_i) ((n1*cos(theta_i)-n2*sqrt(1-(n1/n2*sin(theta_i)).^2))./(n1*cos(theta_i)+n2*sqrt(1-(n1/n2*sin(theta_i)).^2))).^2;
Rp = @(theta_i) ((n2*cos(theta_i)-n1*sqrt(1-(n1/n2*sin(theta_i)).^2))./(n2*cos(theta_i)+n1*sqrt(1-(n1/n2*sin(theta_i)).^2))).^2;
Re = @(theta_i) (abs(Rs(theta_i)+Rp(theta_i)))/2;
% Run the simulation for five times
for i = 1:run
    a = 0; r = 0; t = 0;
    % simulate 10000 photons for each run
    for j = 1:num
        z = 0;
        cz = 1;
        w = 1;
        while (1)
            % update new position
            s = -log(rand) / ut;
            z = z + cz * s;
            % transmittance
            if z > depth               
                if rand > Re(acos(abs(cz)))
                    t = t + w;
                    break;
                % internal reflection at the bottom of the medium
                else
                    z = 2*depth-z;
                    cz = -cz;
                end
            % reflection            
            elseif z < 0
                if rand > Re(acos(abs(cz)))                
                    r = r + w;
                    break;
                % internal reflection at the top of the medium
                else
                    z = -z;
                    cz = -cz;
                end
            end
            % absorbance
            % update w
            dw = w * (ua/ut);
            a = a + dw;
            w = w - dw;
            % if weight too low, roulette
            if w < min
                % terminate for some possibility
                if rand > 1/m
                    break;
                else
                    w = w * m;
                end
            end
            % scattering
            % update new direction
            phi = 2 * pi * rand;
            theta = acos(2 * rand - 1);               
            cz_new = -sin(theta)*cos(phi)*sqrt(1-cz^2)+cz*cos(theta);
            cz = cz_new;
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