%%% Neuharth and Mittelstaedt, "Temporal variations in plume flux: %%%
%%% Characterizing pulsations from tilted plume conduits in a      %%%
%%% rheologically complex mantle"                                  %%%
%%% Script to calculate and plot the critical viscosity ratio      %%%
%%% for a series of plume radii and wavenumbers (Fig. 9b)          %%%

clear, clf
%figure(1),clf;

% Plume parameters
minr = 10e3; % minimum radius (m)
maxr = 100e3; % maximum radius (m)
theta = 67; % Plume tilt (degrees).
% Plume conduit length above transition zone in our model (m).
um_length = 660e3-35e3; % Depth to 660-km transition - lithosphere thickness


% Other parameters
l = 0; % Index for critical vr
p = 0; % Index for critical vr
max_vr = 1000; % Maximum viscosity ratio where this is 10x actual vr.
drho = 30;  % delta rho (rho_mantle - rho_plume, kg/m^3) 
g = 9.81;   % gravity (m/s^2)
mantle_visc = 1.85e20; % Mantle viscosity (Pas) at timestep 0, 300 km depth.


% Loop through radius
for j = minr:5e2:maxr
    l = l+1; % Index for critical vr
    p=0;
 % Loop through wavenumbers
 for t = 1e-6:5.5e-7:1e-4
     p= p+1; % Index for critical vr

  % Loop through viscosity ratios.
  for i = 0:1:max_vr

   x(i+1) = i/10; % Viscosity ratio for determining critical ratio. 
   h = j*2;  % plume diameter (m)
   k = t;   % Wavelength (1/m). ~1e-5 to 7.85e-5 for our models.

   vr = i/10; % Viscosity ratio
   plume_visc = mantle_visc./vr; % From ratio determine plume viscosity (Pa s).


   % Determining growth rate from Houseman and Molnar, 1997, Eq. 10.
   % Where drho is density change, h is plume diameter, and k is the 
   % wavelength.
   q1 = (drho*g*h)/(2.*mantle_visc); % (1/s)
   q2_num = (vr^2*(sinh(k*h)*cosh(k*h) - k*h) + vr*(sinh(k*h)^2-k^2*h^2));
   q2_den = (vr^2*(k*h*cosh(k*h)^2+k^3*h^3) + 2*vr*k*h*sinh(k*h)*cosh(k*h) ...
          + k*h*sinh(k*h)^2 - k^3*h^3);

   q = q1*(q2_num/q2_den); % Growth rate in 1/s.

   % Convert rate to time, and convert to Myr.
   instability_time(i+1) = (q^-1)/(60*60*24*365)/1e6; % (myr)

   % Stokes velocity with viscosity ratio dependence, from table 7
   % in Meriaux et al., 2011
   ks = 0.13*(vr)^0.62;  
   vs_meriaux = ((ks*drho*g*(h/2)^2)/mantle_visc); % (m/s)

   % Velocity of rising tilted conduit. Eq. from Skilbeck and Whitehead, 1978. 
   % Calculate aspect ratio of plume length to radius. 
   ar = (um_length/cosd(theta))/(h/2);
   vs_tiltedconduit = log(ar).*((drho.*g.*(h./2).^2)./(8.*mantle_visc)) ...
                   .*(3 - cosd(2.*(90-theta)) ); % (m/s)

   % substitute pipe flow velocity for "rise rate" - below is the average
   % velocity through the pipe (Q/(pi*R^2) and then accounting for the angle
   % of the plume.
   % tilt = 67 degrees
   vs_pipe1 = cosd(theta).*(drho.*g.*(h/2).^2)./(8.*plume_visc); % (m/s)
   % tilt = 33.5 degrees
   vs_pipe2 = cosd(theta/2).*(drho.*g.*(h/2).^2)./(8.*plume_visc); % (m/s) 
   % tilt = 0 degrees
   vs_pipe3 = cosd(0).*(drho.*g.*(h/2).^2)./(8.*plume_visc); % (m/s)

   % Calculate rise time for three separate velocities based on
   % plume tilt.
   rise_time1(i+1) = um_length/vs_pipe1/(60*60*24*365)/1e6; % (myr)
   rise_time2(i+1) = um_length/vs_pipe2/(60*60*24*365)/1e6; % (myr)
   rise_time3(i+1) = um_length/vs_pipe3/(60*60*24*365)/1e6; % (myr)
  end % End VR loop iteration


  % Now find when the instability time first exceeds the rise time,
  % indicating where instability formation is no longer possible 
  % (critical ratio). For given wavenumber and plume radius.
  % Calculate for each plume tilt.
  %
  % Note: This assumes rise time is initially greater than instability
  % time. This is the case with vs_meriaux or vs_pipe, but for the 
  % vs_tiltedconduit the opposite is true, in which case the critical ratio
  % indicates where instability formation becomes possible.
  for i=1:length(instability_time)
      cross(p,l) = max_vr/10;
      if instability_time(i) > rise_time1(i)
          cross(p,l) = x(i);
          break
      end        
  end

  for i=1:length(instability_time)
      cross2(p,l) = max_vr/10;
      if instability_time(i) > rise_time2(i)
          cross2(p,l) = x(i);
          break
      end        
  end

  for i=1:length(instability_time)
      cross3(p,l) = max_vr/10;
      if instability_time(i) > rise_time3(i)
          cross3(p,l) = x(i);
          break
      end        
  end

 end % End wavenumber loop iteration
end % End plume radius loop iteration

x = (minr:500:maxr)/1000; % radius (km)
y = 1e-6:5.5e-7:1e-4; % wavenumber (1/m)
levels = [1, 5, 10, 20, 40, 60, 100];
x2 = [80 40 40 80];

% box for modeled plumes. Min from r1_tz41. max from r2_4cm
y2 = [1.1785e-5 1.1785e-5 6.0283e-5 6.0283e-5];

c1 = [229/255,   229/255,   229/255];
c2 = [225/255,   123/255,   123/255];
c3 = [107/255,   185/255,   190/255];
c4 = [142/255,   110/255, 169/255];
h2 = fill(x2,y2,'k','LineWidth',0.5);
hold on
set(h2, 'edgealpha',0.5);
set(h2,'facealpha',0.1,'facecolor','k');
[C, h] = contour(x,y,cross,levels,'k','ShowText','on','LineWidth',1.75,'LabelSpacing',175); 
clabel(C,h,'Color','k','FontSize',12,'FontWeight','bold')
[C, h] = contour(x,y,cross2,levels,':','Color',c4,'LineWidth',0.75,'LabelSpacing',150); 
clabel(C,h,'Color',c4,'FontSize',10)
[C, h] = contour(x,y,cross3,levels,'--','Color',c3,'LineWidth',0.75,'LabelSpacing',150); 
clabel(C,h,'Color',c3,'FontSize',10)
ylim([1e-6 7e-5])
yticks([1e-5 2e-5 3e-5 4e-5 5e-5 6e-5 7e-5]) 
xticks([20 40 60 80 100])
title('Critical viscosity ratio')
xlabel('Plume radius (km)', 'FontWeight','bold');
ylabel('Wavenumber (m^-^1)','FontWeight','bold');

hold on

% White box for text %
x3 = [77 43 43 77];
y3 = [3.35e-5 3.35e-5 4.35e-5 4.35e-5];
h2 = fill(x3,y3,'k','LineWidth',0.5);
hold on
set(h2, 'edgealpha',0);
set(h2,'facealpha',1,'facecolor',c1);
text(45,3.85e-5,'Modelled plumes','FontWeight','bold','Color','k','FontSize',14);

% Scatter point from other figure
xp = 70; % (km)
yp = 3.1416e-5; % (1/m)
scatter(xp,yp,75,'k','filled')
set(gca,'FontSize',12)

saveas(gcf,'7b.pdf')