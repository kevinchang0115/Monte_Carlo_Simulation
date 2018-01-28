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
dr = 0.0025;
dz = 0.0025;
rnum = width/dr;
dnum = depth/dz;
min = 0.001;
m = 20;
g1 = 0.79;
g2 = 0.79;
Agrid = zeros(run,dnum,rnum);
Fgrid = zeros(run,dnum,rnum);
n0 = 1;
n1 = 1.4;
n2 = 1.4;
Rs = @(n1,n2,theta_i) ((n1*cos(theta_i)-n2*sqrt(1-(n1/n2*sin(theta_i)).^2))./(n1*cos(theta_i)+n2*sqrt(1-(n1/n2*sin(theta_i)).^2))).^2;
Rp = @(n1,n2,theta_i) ((n2*cos(theta_i)-n1*sqrt(1-(n1/n2*sin(theta_i)).^2))./(n2*cos(theta_i)+n1*sqrt(1-(n1/n2*sin(theta_i)).^2))).^2;
Re = @(n1,n2,theta_i) (abs(Rs(n1,n2,theta_i)+Rp(n1,n2,theta_i)))/2;
A = zeros(1,run);
R = zeros(1,run);
T = zeros(1,run);
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
% plot the A.R.T graph
A = A/num; R = R/num; T = T/num;
y = [mean(A);mean(R);mean(T)];
bar(y);
set(gca,'XTickLabel',{'Absorbance', 'Reflectance', 'Transmittance'});
title('Average');
xlabel('A.R.T'); ylabel('ratio');
text(1:length(y), y, num2str(y), 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom');

% plot the heat map of absorption distribution of scattered photons
figure
for k = 1:run
    subplot(2,3,k);
    imagesc(linspace(0,width*10,rnum), linspace(0,depth*10,dnum), squeeze(Agrid(i,:,:)));
    ylabel(colorbar,'scattered photons (1/cm^3)');
    title(['Simulation ', num2str(k)]);
    xlabel('r (mm)'); ylabel('z (mm)');
end
subplot(2,3,6);
imagesc(linspace(0,width*10,rnum), linspace(0,depth*10,dnum), squeeze(mean(Agrid,1)));
ylabel(colorbar,'scattered photons (1/cm^3)');
title('Average');
xlabel('r (mm)'); ylabel('z (mm)');

% plot the heat map of impulse response for the fluence rate of scattered photons
for i = 1:run
    for j = 1:dnum
        for k = 1:rnum
            dv = (2*k-1)*pi*dz*dr^2;
            if dnum <= t1/dz
                Fgrid(i,j,k) = Agrid(i,j,k)/dv/num/ua1;
            else
                Fgrid(i,j,k) = Agrid(i,j,k)/dv/num/ua2;                
            end
        end
    end
end
figure
for k = 1:run
    subplot(2,3,k);
    imagesc(linspace(0,width*10,rnum), linspace(0,depth*10,dnum), squeeze(Fgrid(i,:,:)));
    ylabel(colorbar,'fluence rate (1/cm^2)');
    title(['Simulation', num2str(k)]);
    xlabel('r (mm)'); ylabel('z (mm)');
end
subplot(2,3,6);
imagesc(linspace(0,width*10,rnum), linspace(0,depth*10,dnum), squeeze(mean(Fgrid,1)));
ylabel(colorbar,'fluence rate (1/cm^2)');
title('Average');
xlabel('r (mm)'); ylabel('z (mm)');