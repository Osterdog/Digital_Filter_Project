clc;
clear all;
close all;

load 'AssignmentData15.mat'; %loads the signal data into matlab workspace
t = 0:1/Fs:49452/Fs; % creates a time vector from the provided sampling frequency

figure(1)
subplot(2,1,1);
plot(t,xc);
title('Signal wave form from Lab3Data.mat');
ylabel('Clean Signal');
xlabel('time');
subplot(2,1,2);
plot(t,xn);
ylabel('Noisy Signal');
xlabel('time');
%% FFT of clean signal %%
Lc = length(xc);
NFFTc = 2^nextpow2(Lc);
Xfftc = fft(xc,NFFTc)/Lc;
mag_Xfftc = abs(Xfftc(1:NFFTc/2+1));
fc = Fs/2*linspace(0,1,NFFTc/2+1);
% FFT of noisy signal %
Ln = length(xn);
NFFTn = 2^nextpow2(Ln);
Xfftn = fft(xn,NFFTn)/Ln;
mag_Xfftn = abs(Xfftn(1:NFFTn/2+1));
fn = Fs/2*linspace(0,1,NFFTn/2+1);

% Single-sided amp plot
figure(2)
subplot(2,1,1);
plot(fc,mag_Xfftc);
title('Single-sided amp the given spectrum');
ylabel('Amplitude |xc|');
xlabel('Frequency (Hz)');
subplot(2,1,2);
plot(fn,mag_Xfftn);
axis([0 6000 0 8e-3])
title('Single-sided amp the given spectrum');
ylabel('Amplitude |xn|');
xlabel('Frequency (Hz)');
%% Filter Noisy Signal %%
cut1 = 3999/(Fs/2); %lower end frequency of stop band
cut2 = 4000.28/(Fs/2); %higher end frequency of stop band

[b,a] = ellip(3,2,23.1,[cut1 cut2],'stop');
freqz(b,a) %frequency response of filter

y = filtfilt(b,a,xn);
%FFT of filter:
Lstop = length(y);
NFFTstop = 2^nextpow2(Lstop);
Xfftstop = fft(y,NFFTstop)/Lstop;
mag_Xfftstop = abs(Xfftstop(1:NFFTstop/2+1));
ystop = Fs/2*linspace(0,1,NFFTstop/2+1);

%Single-sided amp plot
figure(3)
subplot(3,1,1);
plot(fc,mag_Xfftc);
title('Single-sided amp the given spectrum');
ylabel('|xc|');
xlabel('Frequency (Hz)');
subplot(3,1,2);
plot(fn,mag_Xfftn);
axis([0 6000 0 8e-3])
title('Single-sided amp the given spectrum');
ylabel('|xn|');
xlabel('Frequency (Hz)');
subplot(3,1,3);
plot(ystop,mag_Xfftstop);
axis([0 6000 0 8e-3])
title('Single-sided amp the given spectrum');
ylabel('|filtered xn|');
xlabel('Frequency (Hz)');
%% Singled sided amplitude spectrum comparing filtered and clean signals
figure(4)
plot(ystop,mag_Xfftstop,'r');
axis([0 6000 0 8e-3])
title('Single-sided amp the given spectrum');
hold on
plot(fc,mag_Xfftc,'b');
title('Single-sided amp the given spectrum');
ylabel('|xc| and |filtered|');
xlabel('Frequency (Hz)');
%% Sound System Specs %%
lowend = 38/(Fs/2); %Hz
cross1 = 260/(Fs/2); %Hz
cross2 = 2600/(Fs/2); %Hz
highend = 20000/(Fs/2); %Hz
%% Subwoofer %%
[d,c] = cheby1(3,2,[lowend cross1],'bandpass');
freqz(d,c); %frequency response of cheby1 bandpass for subwoofer

sub = filtfilt(d,c,y);
%calculating single-sided fft of subwoofer
Lsub = length(sub); 
NFFTsub = 2^nextpow2(Lsub);
Xfftsub = fft(y,NFFTsub)/Lsub;
mag_Xfftsub = abs(Xfftsub(1:NFFTsub/2+1));
ysub = Fs/2*linspace(0,1,NFFTsub/2+1);

figure(5) %single sided amp spectrum of subwoofer to filtered
subplot(2,1,1);
plot(ysub,mag_Xfftsub);
ylabel('|Subwoofer|');
xlabel('Frequency (Hz)');
subplot(2,1,2);
plot(ystop,mag_Xfftstop);
title('Single-sided amp the given spectrum');
ylabel('|filter|');
xlabel('Frequency (Hz)');

figure(6) %waveform comparison subwoofer and filtered
plot(t,y,'b');
hold on
plot(t,sub,'r');
xlabel('Time');
ylabel('Signal')
%% Woofer %%
[f,e] = cheby1(3,2,[cross1 cross2],'bandpass');
woof = filtfilt(f,e,y);
freqz(f,e) %frequency response of cheby1 bandpass

figure(7)
plot(t,y,'b'); %waveform comparison woofer and filtered
xlabel('Time (s)');
ylabel('Signal');
hold on
plot(t,woof,'r');
%% Tweeter %%
[h,g] = cheby1(6,2,cross2,'high');
tweet = filtfilt(h,g,y);
freqz(h,g) %frequency response of cheby1 filter

figure(8) %waveform comparison tweet and filtered
plot(t,y,'b');
xlabel('Time (s)');
ylabel('Signal');
hold on
plot(t,tweet,'r');
%% Comparing the three Speakers %%
figure(9)
plot(t,y,'m');
title('Comparison of Waveforms');
xlabel('Time (s)');
ylabel('Signal');
hold on
plot(t,sub,'k');
hold on
plot(t,woof,'b');
hold on
plot(t,tweet,'r');
%Find the sum%
system = sub+woof+tweet;

figure(10)
plot(t,y,'b')
title('Filtered Waveform vs Simulation Waveform Comparison');
xlabel('Time (s)');
ylabel('Signal');
hold on
plot(t,system,'r');
%% Compare with clean %%
figure(11)
plot(t,xc,'b');
title('Clean Waveform vs System Waveform comparison');
xlabel('Time (s)');
ylabel('Signal');
hold on
plot(t,system,'r');
%% SS Amp spectrum b/w clean and system %%
Lsystem = length(system);
NFFTsystem = 2^nextpow2(Lsystem);
Xfftsystem = fft(y,NFFTsystem)/Lsystem;
mag_Xfftsystem = abs(Xfftsystem(1:NFFTsystem/2+1));
ysystem = Fs/2*linspace(0,1,NFFTsystem/2+1);

figure(12)
plot(fc,mag_Xfftc,'r');
title('Single-Sided Amplitude given the Spectrum');
ylabel('|xc| and |system|');
xlabel('Frequency (Hz)');
hold on
plot(ysystem,mag_Xfftsystem,'b');