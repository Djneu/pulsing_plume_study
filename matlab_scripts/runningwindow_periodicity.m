% Script to calculate periodicities at given intervals 
% within a specific window.
% To run this, load the volume flux .mat
clear

name = 'vf_r2_tz31_v8_t550_p1.mat'
load(['../numerics/volume_flux/', name])

% Time UMP begins at
umpt = 202;

% Window size
interval = 25;
% Alpha
signif = 0.01;

%Find out when the max volume flux occurs
[M,I] = max(final410)

% Maximum signal we consider.
maxsig = 10; 
xy = maxsig; 
%Set max year as the year max volume flux occurs
maxyear = final410(I(2),1)
periodicity = [0, 0, 0];

% These two runs had a maximum volume flux during a TZP pulse and not plume
% arrival.
if contains('vf_r2_tz31_p10.mat',name)
    maxyear = 140;
end

if contains('vf_r2_tz11.mat',name)
    maxyear = 118.3;
end

% Make sure window doesn't take data before start of UMP
llow = ceil((umpt - maxyear + interval/2)/10)*10

% Make sure we don't take a point where the window has no data.
hhigh = floor((400-maxyear-interval/2)/10)*10

for yy = llow:10:hhigh

    signif = 0.01;
    %Find out when the max volume flux occurs
    [M,I] = max(final410);
    xy = maxsig; 


    %Make new matrix that only keeps data from 5 Myr before max volume flux,
    %and 200 Myr after max volume flux. (Generally max is when plume first
    %passes the chosen depth.)
    for i = 1:length(final410)
      if final410(i,1) >= (maxyear+(yy-interval/2)) && final410(i,1) <= (maxyear+(yy+interval/2))
         time200(i,1) = final410(i,1);
         time200(i,2) = final410(i,2);
      else
         time200(i,:) = 0;
      end
    end

    %Get rid of all the zero elements.
    time200 = time200(any(time200,2),:);

    % If we have nothing for some reason, go to next interval.
    if(isempty(time200) == 1)
      break
    end

    %Set first element as 0 so we can track over the total data we kept
    time200(:,1) = time200(:,1)-time200(1,1);

 % Calculate periodogram
 [P,f,alpha,M] = lomb(time200(:,2),time200(:,1));
 % Convert frequency to period.
 period = 1./f;
 
 % Remove any periods above our maximum considered signal.
 for ii=1:length(period)
   if period(ii,1) >= maxsig
       P(ii,1) = NaN;
       period(ii,1) = NaN;
       alpha(ii,1) = NaN;
   end
 end

  period(any(isnan(period), 2), :) = [];
  P(any(isnan(P), 2), :) = [];
  alpha(any(isnan(alpha), 2), :) = [];
 
  % Keep only significant periods.
  for iii=1:length(alpha)
      if alpha(iii,1) > signif
          sig_period(iii,1) = NaN;
          sig_P(iii,1) = NaN;
      else
          sig_period(iii,1) = period(iii);
          sig_P(iii,1) = P(iii);
      end
  end

  sig_period(any(isnan(sig_period), 2), :) = [];
  sig_P(any(isnan(sig_P), 2), :) = []; 
 
  Pyy = sig_P/max(sig_P);
  mp = sig_period.*Pyy;
  meann = sum(mp)/sum(Pyy);
  
  
TF = islocalmax(P); 
PyyP = P(TF);
periodp(:,1) = period(TF);

a = [signif];
z = -log(1-(1-a).^(1/M));

for jj=1:length(periodp)
    if PyyP(jj) < z
        PyyP(jj) = NaN;
        periodp(jj,1) = NaN;
    end
end

periodp(any(isnan(periodp), 2), :) = [];
PyyP(any(isnan(PyyP), 2), :) = []; 


if length(periodp) == 0
   peak_time(1,:) = 0;
else
   peak_time(1,:) = periodp;
end

yy
% Add all significant signals to final matrix.
for i=1:length(peak_time)
    temp_per(i,1) = peak_time(i); 
    temp_per(i,2) = yy;
    temp_per(i,3) = yy+maxyear;
    periodicity = [periodicity; temp_per(i,:)];
end

% To plot a specific example
if yy == 50

% subplot(2,1,1),
% plot(time200(:,1),time200(:,2), 'k', 'LineWidth', 2)
% title('Time vs. Volume flux 4 cm/yr (depth 100 km)');
% xlabel('Time (Myr)');
% ylabel('Volume flux (m^3/yr)');
% xlim([0 400]);
% 

%fig = figure
%fig.PaperUnits = 'inches';
%fig.PaperPosition = [0 0 6 4.43]; 
% subplot(2,1,2)
% plot(period,P, 'k', 'LineWidth', 1)
% hold on
% scatter(periodp,PyyP,'r*')
% title('Lomb-normalized Periodogram')
% ylabel('Power')
% xlabel('Period (Myr)')
% xlim([0 maxsig]);
% xticks(0:10:maxsig)
% grid on;
%     styles = {':','-.','--'};
%     La = length(a);
%     hold on;
%     for i=1:La
%         line([0,maxsig],[z(i),z(i)],'Color','k','LineStyle',styles{ceil(i*3/La)});
%         text(maxsig-6,z(i)+0.03*max(P),strcat('\alpha = ',num2str(a(i))),'fontsize',10);
%     end
end

clearvars -except yy final410 periodicity maxyear interval maxsig
end

for i=1:length(periodicity)
    if periodicity(i,1) == 0
        periodicity(i,2) = 0;
        periodicity(i,3) = 0;
    end
end

periodicity = periodicity(any(periodicity,2),:);




