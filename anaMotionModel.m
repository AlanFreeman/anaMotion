%%% This code simulates changes in light-dark latency differences with stimulus speed measured in physiological experiments. 
%%% It uses Gaussian functions to simulate different stimulus speeds, a Naka-Rushton function to simulate neural blur and
%%% a spike threshold changes to simulate the fast decay o light-dark latency differences with velocity. 

clear all; close all;

%%% Neuronal blur. 
ONnre=1.4;                          %%% ON Naka-Rushton exponent
ONnrlum50=0.36;                     %%% ON Naka-Rushton half-response luminance; 
ONRmax=80;                          %%% ON Naka-Rushton Rmax 
OFFnre=2.3;                         %%% OFF Naka-Rushton exponent
OFFnrlum50=0.5;                     %%% OFF Naka-Rushton half-response luminance; 
OFFRmax=100;                        %%% OFF Naka-Rushton Rmax 

spk_th_factor = [1 2 2.5 3 3.5];    %%% spike threshold value for each velocity. Same value for ON and OFF
half_max = 0.5;                     %%% Half-max based on the response to lights (used to calculate the L-D latency difference)

%%% Time of stimulus motion
timemov=200;                        %%% time of motion in units of 10 ms (2,000 ms for 5 deg/s, 700 for 15 deg/s, 350 for 30 deg/s, 180 for 60 deg/s)
timewindow=timemov*2;
velocities=[5 10 16.15 30 60];      %%% in deg/s
n_vel=length(velocities);

%%% Show psths after applying neuronal blur to ON and OFF visual responses 
figure(1)
set(1,'Position',[2048 696 560 420])
clf
for vel=1:n_vel
    
    %%% generate psths for different stimulus velocities using Gaussian functions
    stimdur(vel)=timemov/velocities(vel);                           
    stim = fspecial('gaussian',timewindow, stimdur(vel));
    stim=stim(timewindow/2,:);
    stim=stim/max(stim);
    
    %%% apply neuronal blur to the responses
    stiml = ONRmax*((stim.^ONnre)./(ONnrlum50.^ONnre+stim.^ONnre));
    stimd = OFFRmax*((stim.^OFFnre)./(OFFnrlum50.^OFFnre+stim.^OFFnre));
    
    %%% plot the responses to light (red lines) and dark stimuli (blue lines)
    subplot(n_vel,1,vel);
    plot(stiml,'r'); 
    hold on;
    plot(stimd,'b'); 
    hold off;
    
    if vel == 1
        title('Responses after applying neuronal blur')
    end
    
    if vel == 5
        xlabel('Time (10s of ms)')
    end

end

%%% Show psths after applying neuronal blur and spike threshold to ON and OFF visual responses 
figure(2)
set(2,'Position',[2626 693 560 420])
clf
for vel=1:n_vel
    
    %%% generate psths for different stimulus velocities using gaussian functions
    stimdur(vel)=timemov/velocities(vel);                           
    stim = fspecial('gaussian',timewindow, stimdur(vel));
    stim=stim(timewindow/2,:);
    stim=stim/max(stim);
    
    %%% apply neuronal blur to the responses
    stiml2 = (ONRmax*((stim.^ONnre)./(ONnrlum50.^ONnre+stim.^ONnre))).^spk_th_factor(vel); %%%put the gaussians through a nonlinearity function
    stimd2 = (OFFRmax*((stim.^OFFnre)./(OFFnrlum50.^OFFnre+stim.^OFFnre))).^spk_th_factor(vel);

    %%% calculate latency at half maximum response to light stimuli
    [~,Id(vel)] = min(abs(stimd2 - half_max*max(stiml2)));
    [~,Il(vel)] = min(abs(stiml2 - half_max*max(stiml2)));
    
    %%% calculate the latency difference between lights and darks. Then, divide by 100 to get difference in seconds (each bin in time is 10 ms then divided by 1000 ms to get seconds)
    tempint(vel) = ((Il(vel)-Id(vel))/100)*-1;                  
    
    %%% plot the responses to light (red lines) and dark stimuli (blue lines)
    subplot(n_vel,1,vel);
    plot(stiml2./max(stimd2),'r');
    hold on;
    plot(stimd2./max(stimd2),'b'); 
    xlim([100 300])
    set(gca,'TickDir','out')
    box off
    %ylabel('Response')
    
    if vel == 1
        title('Responses after applying neuronal blur and spike threshold')
    end
        
    if vel == 5
        xlabel('Time (10s of ms)')
    end
    
end

%%% Show the rise time of visual responses and latency at half-maximum response
figure(3)
set(3,'Position',[3206 697 560 420])
clf
for vel=1:n_vel
    
    %%% generate psths for different stimulus velocities using Gaussian functions
    stimdur(vel)=timemov/velocities(vel);
    stim = fspecial('gaussian',timewindow, stimdur(vel));
    stim=stim(timewindow/2,:);
    stim=stim/max(stim);
    
    %%% apply neuronal blur to the responses    
    stiml2 = (ONRmax*(stim.^ONnre)./(ONnrlum50.^ONnre+stim.^ONnre)).^spk_th_factor(vel); %%%put the gaussians through a nonlinearity function
    stimd2 = (OFFRmax*(stim.^OFFnre)./(OFFnrlum50.^OFFnre+stim.^OFFnre)).^spk_th_factor(vel);

    %%% plot the responses to light (red lines) and dark stimuli (blue lines)
    subplot(n_vel,1,vel);
    plot(stiml2(100:200)./max(stimd2), 'r');
    hold on;
    plot(stimd2(100:200)./max(stimd2), 'b')
    hold on;
    plot(0:2:100,half_max*(max(stiml2)./max(stimd2)):half_max*(max(stiml2)./max(stimd2)),'.k','MarkerSize',2);
    box off
    xlim([0 100])
    set(gca,'TickDir','out');
    set(gca,'XTick',[]);
    if vel == 1
        title('Rise Time')
    end
    if vel == 5
        xlabel('Time (10s of ms)')
        set(gca,'XTick',[0:10:100]);
    end
end

%%% Show Luminance Response Function for ON and OFF 
figure(4)
set(4,'Position',[2047 180 560 420])
clf
lum = [0:0.001:1];

ON_LRF = ONRmax*(lum.^ONnre)./(ONnrlum50.^ONnre+lum.^ONnre);
OFF_LRF = OFFRmax*(lum.^OFFnre)./(OFFnrlum50.^OFFnre+lum.^OFFnre);

plot(lum,ON_LRF, 'r');
hold on;
plot(lum,OFF_LRF, 'b')
box off
xlim([0 1])
set(gca,'TickDir','out');
xlabel('Luminance (%)')
ylabel('Response')

%%% Plot Model and data from latency differences
figure(5)
set(5,'Position',[2625 180 560 420])
clf
avg_data = [ 0.07254;-0.00042782;-0.0030334;-0.00057827;-0.0025167];
semilogx(velocities,tempint, '-','Color','k') 
hold on 
semilogx(velocities,avg_data, 'o','Color','b')
xlim([4 70]); 
ylim([-0.02 0.1]);
box off
set(gca,'TickDir','out');
set(gca,'YTick',[-0.05 0 0.05 0.1]);
set(gca,'XTick',[5 10 20 30 60]);
% hold on;plot(5:2:80,0:0,'.k');
xlabel('Velocity (deg/s)')
ylabel('Light-dark latency difference (s)')
