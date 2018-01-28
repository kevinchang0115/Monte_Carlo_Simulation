% Author:      B01901056 ±i³Í´Q
% Date:        2016/5/2
% Purpose:     To use Monte Carlo simulation to determine photon absorption
%              distribution for laser beams with a uniform distribution.

% Description: 
%     This program uses variable weight photons and anisotropic scattering
%     (Henyey-Greenstein phase function) to simulate photon absorption
%     distribution for a collimated beam that has a uniform distribution
%     with a radius of 0.5 mm and irradiance of 1 W/cm^2. After running
%     this program, we will get a subplot of six graphs. First five are
%     five runs of the simulation of 10000 photons, and the last one is the
%     graph of average values.

% Variables:
%     ua:     absorption coefficient of the medium (cm^-1)
%     us:     scattering coefficient of the medium (cm^-1)
%     ut:     ua + us (cm^-1)
%     num:    the number of photons in simulation
%     run:    for how many times it runs the simulation
%     depth:  the depth of the medium (cm)
%     width:  the width of the grid (cm) 
%     s:      the distance a photon travels before absorbance or scattering
%             (cm)
%     x:      x position of the photon (cm)
%     y:      y position of the photon (cm)
%     z:      the depth which the photon reaches (cm)
%     cx:     cosin of the angle between the moving dirction and x-axis
%     cy:     cosin of the angle between the moving dirction and y-axis
%     cz:     cosin of the angle between the moving dirction and z-axis
%     cx_new: new cx
%     cy_new: new cy
%     cz_new: new cz
%     phi:    azimuthal angle of the scattering photon
%     theta:  polar angle of the scattering photon
%     w:      the weight of a photon
%     g:      anisotropy factor
%     min:    if weight less than min, play roulette
%     m:      if rand larger than 1/m, terminate it
%     n1:     the refractive index of tissue (n=1.5)
%     n2:     the refractive index of air (n=1)
%     Rs:     the function of the reflectance for s-polarized light
%     Rp:     the function of the reflectance for p-polarized light
%     Re:     the function of the reflectance for unpolarized light
%     dr:     the grid unit in radius (cm)
%     dz:     the grid unit in depth (cm)
%     dv:     the unit volume of the grid (1/cm^3)
%     rnum:   the grid num in radius
%     dnum:   the grid num in depth
%     Agrid:  scattered photon number of the grid
%     Fgrid:  fluence rate of the grid (1/cm^2)
%     scatter:whether it's scattered photons or not
%     ir:     the index of the grid in radius
%     iz:     the index of the grid in depth
%     b_radius: the radius of the beam with uniform distribution (cm)
%     irradiance: the irradiance of the beam (W/cm^2)
%     total_power: the total power of the beam (W)
%     photon_power: the power of each photon (W)

% Initialization
clear;
ua = 6; us = 414;
ut = ua + us;
num = 10000;
run = 5;
depth = 0.15;
width = 0.3;
dr = 0.01;
dz = 0.01;
rnum = width/dr;
dnum = depth/dz;
min = 0.001;
m = 20;
g = 0.91;
Agrid = zeros(run,dnum,rnum);
Fgrid = zeros(run,dnum,rnum);
n1 = 1.37;
n2 = 1;
b_radius = 0.05;
irradiance = 1;
total_power = irradiance*b_radius^2*pi;
photon_power = total_power/num;
Rs = @(theta_i) ((n1*cos(theta_i)-n2*sqrt(1-(n1/n2*sin(theta_i)).^2))./(n1*cos(theta_i)+n2*sqrt(1-(n1/n2*sin(theta_i)).^2))).^2;
Rp = @(theta_i) ((n2*cos(theta_i)-n1*sqrt(1-(n1/n2*sin(theta_i)).^2))./(n2*cos(theta_i)+n1*sqrt(1-(n1/n2*sin(theta_i)).^2))).^2;
Re = @(theta_i) (abs(Rs(theta_i)+Rp(theta_i)))/2;
% Run the simulation for five times
for i = 1:run
    % simulate 10000 photons for each run
    for j = 1:num
        % collimated flat-top beam
        r = b_radius*sqrt(rand);        
        angle = 2*pi*rand;
        x = r*cos(angle);
        y = r*sin(angle);
        z = 0;
        cx = 0; cy = 0; cz = 1;
        w = 1;
        scatter = false;
        while (1)
            % update new position
            s = -log(rand) / ut;
            x = x + cx * s;
            y = y + cy * s;
            z = z + cz * s;
            % transmittance
            if z > depth               
                if rand > Re(acos(abs(cz)))
                    break;
                % internal reflection at the bottom of the medium
                else
                    z = 2*depth-z;
                    cz = -cz;
                end
            % reflection            
            elseif z < 0
                if rand > Re(acos(abs(cz)))                
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
            if scatter
                radius = sqrt(x^2+y^2);
                ir = floor(radius/dr);
                iz = floor(z/dz);
                if ir < rnum
                    Agrid(i,iz+1,ir+1) = Agrid(i,iz+1,ir+1) + dw;
                end
            else
                scatter = true;
            end
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
            theta = acos((1+g^2-((1-g^2)/(1-g+2*g*rand))^2)/(2*g));
            if abs(cz) > 0.99999
                cx = sin(theta)*cos(phi);
                cy = sin(theta)*sin(phi);
                cz = cos(theta)*cz/abs(cz);
            else
                cx_new = (sin(theta)/sqrt(1-cz^2))*(cx*cz*cos(phi)-cy*sin(phi))+cx*cos(theta);
                cy_new = (sin(theta)/sqrt(1-cz^2))*(cy*cz*cos(phi)+cx*sin(phi))+cy*cos(theta);           
                cz_new = -sin(theta)*cos(phi)*sqrt(1-cz^2)+cz*cos(theta);
                cx = cx_new;
                cy = cy_new;
                cz = cz_new;
            end
        end       
    end
end
% plot the heat map
for i = 1:run
    for j = 1:dnum
        for k = 1:rnum
            dv = (2*k-1)*pi*dz*dr^2;
            Fgrid(i,j,k) = photon_power*Agrid(i,j,k)/dv/num/ua;
        end
    end
end
figure
for k = 1:run
    subplot(2,3,k);
    imagesc(linspace(0,3,30), linspace(0,1.5,15), squeeze(Fgrid(i,:,:)));
    ylabel(colorbar,'fluence rate (W*cm^-2)');
    title(['Heat map ', num2str(k)]);
    xlabel('r (mm)'); ylabel('z (mm)');
end
subplot(2,3,6);
imagesc(linspace(0,3,30), linspace(0,1.5,15), squeeze(mean(Fgrid,1)));
ylabel(colorbar,'fluence rate (W*cm^-2)');
title('Average heat map');
xlabel('r (mm)'); ylabel('z (mm)');