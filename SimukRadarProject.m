
clc
close all
clear all

load('DataCube_ELE691_Simuk.mat');
x=Radar_Data;
y=Radar.LFM;
%%%%%%%%%%%%%%%%%%%%%%%%--Matched Filtering--%%%%%%%%%%%%%%%%%%%%%%%%%%%
j=conj((fliplr(y)));
for i=1:256
     PulseCompressedData(i,:)=conv(x(i,:),j);  
end                                                             
PulseCompressedData=PulseCompressedData(:,[999:1:9998]); 


%%%%%%%%%%%%%%%%%%%%%%%%%--Hamming Window--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=rectwin(256);
for n=1:9000

    WindowAppliedData(:,n)=PulseCompressedData(:,n).*P; %%Window applied
    
end


%%%%%%%%%%%%%%%%%%%%%%%%--Doppler Processing--%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:9000

    DopplerProcessedData(:,n)=fftshift(fft(WindowAppliedData(:,n))); 
    
end

Z=DopplerProcessedData;
%%%%%%%%%%%%%%%%%%%%%%%%--Calculations--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PRF=1/Radar.PRI;

minVelocity=-PRF/2*Radar.lambda;
maxVelocity=PRF/2*Radar.lambda;
velocity=minVelocity:(maxVelocity-minVelocity)/255:maxVelocity;

c=3*10^8;

minRange=c*Radar.Pulse_Width/2;
maxRange=c*Radar.PRI/2;
range=(minRange:(maxRange-minRange)/8999:maxRange);

%%%%%%%%%%%%%%%%%%%%%%%--Garphs Plot--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(range,20*log10(abs(PulseCompressedData(1,:))))
title('Matched filter output for a single pulse')
xlabel('Range(m)')
ylabel('Signal power db')



figure
plot(velocity,20*log10(abs(Z(:,668))));
title('Doppler response of max target')
xlabel('Velocity(m/sec)')
ylabel('Magnitude')   
             

figure
plot(range,20*log10(abs(Z(153,:))))
title('Range response of max target')
xlabel('Range(m)')
ylabel('Magnitude')



figure
mesh(range,velocity,20*log10(abs(Z)))
title('Surface plot of Range vs Doppler')
xlabel('Range(m)')
ylabel('Velocity (m/sec)')
zlabel('Signal Power')





