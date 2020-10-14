clear
clc

SF = 10 ;
BW = 125e3 ;
fc = 915e6 ;
Power = 14 ;

message = "Hello World!" ;

%% Sampling
Fs = 10e6 ;
Fc = 921.5e6 ;
%% Transmit Signal
signalIQ = LoRa_Tx(message,BW,SF,Power,Fs,Fc - fc) ;
Sxx = 10*log10(rms(signalIQ).^2) ;
disp(['Transmit Power   = ' num2str(Sxx) ' dBm'])
%% Plots
figure(1)
spectrogram(signalIQ,500,0,500,Fs,'yaxis','centered')
figure(2)
obw(signalIQ,Fs) ;
%% Received Signal
message_out = LoRa_Rx(signalIQ,BW,SF,2,Fs,Fc - fc) ;
%% Message Out
disp(['Message Received = ' char(message_out)])

