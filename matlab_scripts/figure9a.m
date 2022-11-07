%%% Neuharth and Mittelstaedt, "Temporal variations in plume flux: %%%
%%% Characterizing pulsations from tilted plume conduits in a      %%%
%%% rheologically complex mantle"                                  %%%
%%% Script to calculate and plot the rise and instability time for %%%
%%% a plume, and determine the critical viscosity ratio for a      %%%
%%% single set of parameters (Fig. 9a)                             %%%

clear, clf

% PARAMETERS

% Plume parameters
h = 140e3;  % plume diameter (m)
theta = 67; % plume tilt (degrees).
% Plume conduit length above transition zone in our model (m).
um_length = 660e3-35e3; % Depth to 660-km transition - lithosphere thickness

% Other parameters
v = 0.08; % plate speed (m/yr)
per = 2.5e6;   % average periodicity of instabilities in our models (yr)
mantle_visc = 1.85e20; % Mantle viscosity (Pa s) at timestep 0, 300 km depth.
drho = 30;  % delta rho (rho_mantle - rho_plume, kg/m^3) 
g = 9.81;   % gravity (m/s^2)


% i represents the 10*(viscosity ratio).
for i = 1:1:2000

x(i+1) = i/10; % Viscosity ratio for determining critical ratio and plotting.  
k = (2*pi)/(per*v);   % Wavelength (1/m). ~1e-5 to 7.85e-5 for our models.

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
instability_time(i+1) = (q^-1)/(60*60*24*365)/1e6; % (Myr)


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
% of the plume
vs_pipe = cosd(theta).*(drho.*g.*(h/2).^2)./(8.*plume_visc); % (m/s)

% Here we use the pipe flow, as it is the most restricive and allows us
% to have both tilt and viscosity ratio dependence in the rise time.
vs = vs_pipe; 
rise_time(i+1) = um_length/vs/(60*60*24*365)/1e6; % (Myr)



end

% Now find when the instability time first exceeds the rise time,
% indicating where instability formation is no longer possible 
% (critical ratio).
% Note: This assumes rise time is initially greater than instability
% time. This is the case with vs_meriaux or vs_pipe, but for the 
% vs_tiltedconduit the opposite is true, in which case the critical ratio
% indicates where instability formation becomes possible.
for i=1:length(instability_time)
    cross = 200;
    if instability_time(i) > rise_time(i)
        cross = x(i);
        yy = instability_time(i);
        break
    end
        
end

% removes any zeros.
instability_time(instability_time==0)=[];
x(x==0)=[];
rise_time(rise_time==0)=[];

% Plot
c1 = [225/255,   123/255,   123/255];
figure(1),clf;
cx = [cross cross];
cy = [0 5];
plot(x,rise_time,'k','LineWidth',2)
hold on
plot(x,instability_time,'k-.','LineWidth',2)
xlabel('Viscosity ratio','FontWeight','bold');
ylabel('Time (Myr)','FontWeight','bold');
xlim([0 75]);
ylim([0 5]);
set(gca,'FontSize',12)

tzpx = [0, 0, cross, cross];
tzpy = [0, 5, 5, 0];
umpx = [cross, cross, 75, 75];
patch(umpx, tzpy, c1, 'FaceAlpha', 0.25)
text(44,3.2, ['Plume radius: ', num2str(h/2/1000), ' km'], 'FontSize', 11);
text(44,2.9, ['Plume tilt: ', num2str(theta), char(176)], 'FontSize', 11);
text(44,2.6, ['Wavenumber: ', num2str(k)], 'FontSize', 11);
text(44,2.3, ['Critical viscosity ratio: ', num2str(cross)], 'FontSize', 11);
text(44,3.9, 'No instability', 'FontSize', 14, 'FontWeight', 'bold');
text(44,3.6, 'formation', 'FontSize', 14, 'FontWeight', 'bold');
text(20,3.9, 'Instability', 'FontSize', 14, 'FontWeight', 'bold');
text(20,3.6, 'formation', 'FontSize', 14, 'FontWeight', 'bold');
text(20,3.3, 'possible', 'FontSize', 14, 'FontWeight', 'bold');
legend('Rise time', 'Instability time')

% Scatter point for critical VR
scatter(cross,yy,75,'k','filled')
saveas(gcf,'7a.pdf')


