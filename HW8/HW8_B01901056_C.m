% Author:      B01901056 ±i³Í´Q
% Date:        2016/5/16
% Purpose:     To use Monte Carlo simulation to determine photon absorption
%              distribution and the impulse response for the fluence rate
%              of scattered photons in a double-layer tissue model.

% Description: 
%     This program uses variable weight photons and anisotropic scattering
%     (Henyey-Greenstein phase function) to simulate photon absorption
%     distribution and the impulse response for the fluence rate. After
%     running this program, we will get three plots. The first is the
%     average absorbance, reflectance and transmittance. The second and the
%     third are the subplots of six graphs. Each for the absorption
%     distribution and the impulse response for the fluence rate. First
%     five graphs of each are five runs of the simulation of 10000 photons,
%     and the last one is the graph of average values.

% Variables:
%     ua1:    absorption coefficient of the medium1 (cm^-1)
%     us1:    scattering coefficient of the medium1 (cm^-1)
%     ut1:    ua1 + us1 (cm^-1)
%     ua2:    absorption coefficient of the medium2 (cm^-1)
%     us2:    scattering coefficient of the medium2 (cm^-1)
%     ut2:    ua2 + us2 (cm^-1)
%     num:    the number of photons in simulation
%     run:    for how many times it runs the simulation
%     t1:     the depth of the medium1 (cm)
%     t2:     the depth of the medium2 (cm)      
%     depth:  t1 + t2 (cm)
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
%     g1:     anisotropy factor of medium1
%     g2:     anisotropy factor of medium2
%     min:    if weight less than min, play roulette
%     m:      if rand larger than 1/m, terminate it
%     n0:     the refractive index of the air (n=1)
%     n1:     the refractive index of medium1 (n=1.4)
%     n2:     the refractive index of medium2 (n=1.4)
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
%     tissue1:whether the photon is in medium1 (yes:TRUE, no:FALSE)
%     tissue2:whether the photon is in medium2 (yes:TRUE, no:FALSE)

% Initialization
clear;
ua1 = 37; us1 = 480;
ua2 = 2.2; us2 = 220;
ut1 = ua1 + us1;
ut2 = ua2 + us2;
num = 10000;
run = 5;
t1 = 0.005;
t2 = 0.2;
depth = t1+t2;
width = 0.3;
dr = 0.01;
min = 0.001;
m = 20;
g1 = 0.79;
g2 = 0.79;
n0 = 1;
n1 = 1.4;
n2 = 1.4;
nd = 1.4;
Rs = @(n1,n2,theta_i) ((n1*cos(theta_i)-n2*sqrt(1-(n1/n2*sin(theta_i)).^2))./(n1*cos(theta_i)+n2*sqrt(1-(n1/n2*sin(theta_i)).^2))).^2;
Rp = @(n1,n2,theta_i) ((n2*cos(theta_i)-n1*sqrt(1-(n1/n2*sin(theta_i)).^2))./(n2*cos(theta_i)+n1*sqrt(1-(n1/n2*sin(theta_i)).^2))).^2;
Re = @(n1,n2,theta_i) (abs(Rs(n1,n2,theta_i)+Rp(n1,n2,theta_i)))/2;
A = zeros(1,run);
R = zeros(1,run);
T = zeros(1,run);
P = zeros(4,4,run);
L = 10*dr;
h = [0,0.1,0.2,0.4];
D = 5*dr;
FOV = [15,30,60,180];
% Run the simulation for five times
for i = 1:run
    % simulate 10000 photons for each run
    for j = 1:num
        % collimated flat-top beam
        x = 0; y = 0; z = 0;
        cx = 0; cy = 0; cz = 1;
        w = 1;
        scatter = false;
        tissue1 = true;
        tissue2 = false;
        s = -log(rand) / ut1;
        while (1)
            % update new position
            x = x + cx * s;
            y = y + cy * s;
            z = z + cz * s;
            % if in tissue, not hit boundary, absorbance and scattering
            if ((tissue1 && z<=t1 && z>=0) || (tissue2 && z<=depth && z>=t1))
                % "Absorbance"
                % update w
                if tissue1
                	dw = w * (ua1/ut1);
                else
                    dw = w * (ua2/ut2);                
                end             
                A(1,i) = A(1,i) + dw;
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
                % "Scattering"
                % update new direction and path
                phi = 2 * pi * rand;
                if tissue1
                    theta = acos((1+g1^2-((1-g1^2)/(1-g1+2*g1*rand))^2)/(2*g1));
                    s = -log(rand) / ut1;
                else
                    theta = acos((1+g2^2-((1-g2^2)/(1-g2+2*g2*rand))^2)/(2*g2));    
                    s = -log(rand) / ut2;
                end
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
            % hit boundary
            else
                % in tissue1
                if tissue1
                    % hit boundary of tissue1 and tissue2
                    if z > t1
                        % calculate the path passing the boundary
                        s = (z-t1) / cz;
                        % move to boundary
                        x = x - cx * s;
                        y = y - cy * s;
                        z = t1;
                        % update new path
                        s = s * ut1 / ut2;
                        % update in which tissue
                        tissue1 = false;
                        tissue2 = true;
                    % hit boundary of tissue1 and tissue0
                    else
                        % passing from tissue1 to tissue0 ("Reflection")
                        if rand > Re(n1,n0,acos(abs(cz)))
                            % update in which tissue              
                            R(1,i) = R(1,i) + w;
                            % calculate the transmitted angle
                            cz = sqrt(1-(n1/n0*sqrt(1-cz^2))^2);
                            % calculate the x,y position when reaching the probe
                            for m = 1:4
                                x_new = x + h(m)/cz*cx;
                                y_new = y + h(m)/cz*cy;
                                r = zeros(1,5);
                                for k = 1:4
                                    r = sqrt((x_new-L*dr)^2+y_new^2);
                                    if (r < D && acos(cz) < (FOV(k)/180*pi))
                                        if (rand > Re(n0,nd,acos(cz)) || h(m) == 0)
                                            P(m,k,i) = P(m,k,i) + w;
                                        end
                                    end   
                                end
                            end
                            break;
                        % internal reflection in tissue1
                        else
                            % calculate the path passing the boundary
                            s = z / cz;                            
                            % move to boundary
                            x = x - cx * s;
                            y = y - cy * s;
                            z = 0;
                            % update new direction
                            cz = -cz;
                        end                    
                    end
                % in tissue2
                else
                   % hit boundary of tissue1 and tissue2
                   if z < t1
                        % calculate the path passing the boundary
                        s = (z-t1) / cz;
                        % move to boundary
                        x = x - cx * s;
                        y = y - cy * s;
                        z = t1;
                        % update new path
                        s = s * ut2 / ut1;
                        % update in which tissue
                        tissue1 = true;
                        tissue2 = false;
                    % hit boundary of tissue2 and tissue0
                    else
                        % passing from tissue2 to tissue0 ("Transmittance")
                        if rand > Re(n2,n0,acos(abs(cz)))
                            % update in which tissue          
                            T(1,i)  = T(1,i) + w;
                            break;
                        % internal reflection in tissue2
                        else
                            % calculate the path passing the boundary
                            s = (z-depth) / cz;                            
                            % move to boundary
                            x = x - cx * s;
                            y = y - cy * s;
                            z = depth;
                            % update new direction
                            cz = -cz;
                        end                    
                    end                   
               end
            end
        end   
    end
end
% plot the graph
P = P/num*100;
y = [mean(P(1,1,:)) mean(P(1,2,:)) mean(P(1,3,:)) mean(P(1,4,:));
     mean(P(2,1,:)) mean(P(2,2,:)) mean(P(2,3,:)) mean(P(2,4,:));
     mean(P(3,1,:)) mean(P(3,2,:)) mean(P(3,3,:)) mean(P(3,4,:));
     mean(P(4,1,:)) mean(P(4,2,:)) mean(P(4,3,:)) mean(P(4,4,:))];
bar3(y);
set(gca,'XTickLabel',FOV, 'YTickLabel',h);
title('Average');
xlabel('FOV (¢X)'); ylabel('h (cm)'); zlabel('Percentage (%)');