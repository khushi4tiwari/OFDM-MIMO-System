%%%%%%% Adaptive Modulation Simulation %%%%%%%%%%
clear all;
close all;
clc;
%Given Average SNR
SNR=10:1:30;

N=length(SNR);
N_bits=zeros(1,length(SNR));
eff=zeros(1,length(SNR));
%No. of Symbols
n=100;

%% VEHICULAR %%

%Channel A (Vehicular) Parameters for a 6-tap model
r1_r = sqrt(1).*randn(n*N,1)/sqrt(2);
r1_i = sqrt(1).*randn(n*N,1)/sqrt(2);
r1 = complex(r1_r,r1_i);

r2_r = sqrt(0.79).*randn(n*N,1)/sqrt(2);
r2_i = sqrt(0.79).*randn(n*N,1)/sqrt(2);
r2 = complex(r2_r,r2_i);

r3_r = sqrt(0.126).*randn(n*N,1)/sqrt(2);
r3_i = sqrt(0.126).*randn(n*N,1)/sqrt(2);
r3 = complex(r3_r,r3_i);

r4_r = sqrt(0.1).*randn(n*N,1)/sqrt(2);
r4_i = sqrt(0.1).*randn(n*N,1)/sqrt(2);
r4 = complex(r4_r,r4_i);

r5_r = sqrt(0.032).*randn(n*N,1)/sqrt(2);
r5_i = sqrt(0.032).*randn(n*N,1)/sqrt(2);
r5 = complex(r5_r,r5_i);

r6_r = sqrt(0.01).*randn(n*N,1)/sqrt(2);
r6_i = sqrt(0.01).*randn(n*N,1)/sqrt(2);
r6 = complex(r6_r,r6_i);

%Delay
t1=1;
t2=30;
t3=70;
t4=110;
t5=170;
t6=250;



%%% Magnitude Response of the frequency selective channel
gamma=zeros(64,1);

for i=1:length(SNR)
    for t=1:n
        
 h=zeros(64,1);
        
 h(t1)=r1(i*t);
 h(t2)=r2(i*t);
 h(t3)=r3(i*t);
 h(t4)=r4(i*t);
 h(t5)=r5(i*t);
 h(t6)=r6(i*t);
       
%Frequency Response of the channel
 h_fft=abs(fft(h,64));  
 count=0;
        
 No=1/(10^(SNR(i)/10));

%%%%% Calculating Spectral Efficiency       
        for k=9:1:56
 %gamma=(|H(f)|^2)*SNR
 gamma(k)=10*log10((h_fft(k)^2)*(1/No));
            
           if ((gamma(k)>=7) && (gamma(k)<10))             %BPSK
                count=count+1;
            elseif ((gamma(k)>=10) && (gamma(k)<15.25))     %QPSK
                count=count+2;
            elseif ((gamma(k)>=15.25) && (gamma(k)<18))     %8-PSK
                count=count+3;
            elseif ((gamma(k)>=18) && (gamma(k)<23.5))      %16-QAM
                count=count+4;
            elseif (gamma(k)>=23.5)                        %64-QAM     
                count=count+6;    
            end
        end
        N_bits(i)=count;
        eff(i)=eff(i)+N_bits(i);
    end
    eff(i)=eff(i)/n;
end
eff_vehicular=eff/48;

plot(SNR,eff_vehicular,'-r');
grid on; hold on;

%% PEDESTRIAN %%        
         
%Channel A Parameters for a 6-tap model
r1_r = sqrt(1).*randn(n*N,1)/sqrt(2);
r1_i = sqrt(1).*randn(n*N,1)/sqrt(2);
r1 = complex(r1_r,r1_i);

r1_r = sqrt(0.107).*randn(n*N,1)/sqrt(2);
r2_i = sqrt(0.107).*randn(n*N,1)/sqrt(2);
r2 = complex(r1_r,r2_i);

r3_r = sqrt(0.012).*randn(n*N,1)/sqrt(2);
r3_i = sqrt(0.012).*randn(n*N,1)/sqrt(2);
r3 = complex(r3_r,r3_i);

r4_r = sqrt(0.00525).*randn(n*N,1)/sqrt(2);
r4_i = sqrt(0.00525).*randn(n*N,1)/sqrt(2);
r4 = complex(r4_r,r4_i);

t1=1;
t2=110;
t3=190;
t4=410;


%%% Magnitude Response of the frequency selective channel
gamma=zeros(64,1);

for i=1:1:length(SNR)
    for t=1:1:n
        
        h=zeros(64,1);
        
        h(t1)=r1(i*t);
        h(t2)=r2(i*t);
        h(t3)=r3(i*t);
        h(t4)=r4(i*t);
       
 %Frequency Response of the channel       
 h_fft=abs(fft(h,64));   
        
 count=0;       
 No=1/(10^(SNR(i)/10));
        
 for k=9:1:56
 gamma(k)=10*log10((h_fft(k)^2)*(1/No));
            if ((gamma(k)>=7) && (gamma(k)<10))             %BPSK
                count=count+1;
            elseif ((gamma(k)>=10) && (gamma(k)<15.25))     %QPSK
                count=count+2;
            elseif ((gamma(k)>=15.25) && (gamma(k)<18))     
                count=count+3;
            elseif ((gamma(k)>=18) && (gamma(k)<23.5))      %16-QAM
                count=count+4;
            elseif (gamma(k)>=23.5)                        %64-QAM     
                count=count+6;    
            end
  end
        N_bits(i)=count;
        eff(i)=eff(i)+N_bits(i);
   end
    eff(i)=eff(i)/n;
end
eff_pedestrian=eff/48;


%%% Plot of Spectral Efficiency for Vehicular and Pedestrian Channel Profile
plot(SNR,eff_pedestrian,'-b');
title('Adaptive modulation - Empirical Spectral eff for the channel models');
xlabel('SNR  (dB)');
ylabel('Avg. Spectral eff');
legend('Vehicular','Pedestrian');
grid on; hold off;
