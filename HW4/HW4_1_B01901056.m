% Author:      B01901056 ±i³Í´Q
% Date:        2016/4/5
% Purpose:     Calculate reflection coefficient of diffuse light and
%              compare the values with different index matching.
% Description: 
%     This program calculates the reflection coefficient against incident
%     angles from 0¢X to 90¢X from air into tissue. In addition, it
%     calculates the reflection coefficient of diffuse light, r_d, from
%     tissue to air or water respectively. We can figure out from the graph
%     that the index matching is essential to reflectance.
% Variables:
%     n_a:      the refractive index of air (n=1)
%     n_t:      the refractive index of tissue (n=1.4)
%     n_w:      the refractive index of water (n=1.33)
%     theta:    the incident angle (from 0 to 90 degrees) (radian)
%     theta_c1: the critical angle of light from tissue to air (radian)
%     theta_c2: the critical angle of light from tissue to water (radian)
%     x_axis:   x-axis of the graph
%     Rs:       the function of the reflectance for s-polarized light
%     Rp:       the function of the reflectance for p-polarized light
%     R:        the function of the reflectance for unpolarized light
%     r:        the reflection coefficient of light from air to tissue
%     F1:       the first function for integration of r_d
%     F2:       the second function for integration of r_d
%     r_d1:     the reflection coefficient of diffuse light from tissue to
%               air
%     r_d2:     the reflection coefficient of diffuse light from tissue to
%               water

% Initialization
clear;
n_a = 1;
n_t = 1.4;
n_w = 1.33;
theta = linspace(0,pi/2,90);
theta_c1 = asin(n_a/n_t);
theta_c2 = asin(n_w/n_t);
x_axis = linspace(0,90,90);

% (a) plot the reflection coefficient against incident angles (0„a~90„a)
% from air (n_a) into tissue (n_t).
% define Rs, Rp and R
Rs = @(n1,n2,theta_i) (abs((n1*cos(theta_i)-n2*sqrt(1-(n1/n2*sin(theta_i)).^2))./(n1*cos(theta_i)+n2*sqrt(1-(n1/n2*sin(theta_i)).^2)))).^2;
Rp = @(n1,n2,theta_i) (abs((n2*cos(theta_i)-n1*sqrt(1-(n1/n2*sin(theta_i)).^2))./(n2*cos(theta_i)+n1*sqrt(1-(n1/n2*sin(theta_i)).^2)))).^2;
R = @(n1,n2,theta_i) (Rs(n1,n2,theta_i)+Rp(n1,n2,theta_i))/2;
r = R(n_a,n_t,theta);
% plot
figure
plot(x_axis,r);
title('HW3 (A)');
xlabel('incident angle (¢X)'); ylabel('reflection coefficient');

% (b) calculate the reflection coefficient from tissue to an outside medium
% of air (n_a) and water (n_w)
% define functions for integration
F1 = @(n1,n2,theta_i) R(n1,n2,theta_i).*sin(2*theta_i);
F2 = @(n1,n2,theta_i) sin(2*theta_i);
% integration
r_d1 = integral(@(theta)F1(n_t,n_a,theta),0,theta_c1)+integral(@(theta)F2(n_t,n_a,theta),theta_c1,pi/2);
r_d2 = integral(@(theta)F1(n_t,n_w,theta),0,theta_c2)+integral(@(theta)F2(n_t,n_w,theta),theta_c2,pi/2);
% plot
figure
y = [r_d1;r_d2];
bar(y);
set(gca,'XTickLabel',{'Air', 'Water'});
title('HW3 (B)');
xlabel('outside medium'); ylabel('reflection coefficient');
text(1:length(y), y, num2str(y), 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom');