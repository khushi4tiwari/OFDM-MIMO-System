clear all;
close all;
clc;

SNR=10:1:30;
N=length(SNR);
N_bits=zeros(1,length(SNR));
eff=zeros(1,length(SNR));
n=100;

%% VEHICULAR %%

% Generating 4 channels for 2*2 MIMO system ---Channel Parameters for a 6-tap model

r11_r = sqrt(1).*randn(n*N,1)/sqrt(2);
r11_i = sqrt(1).*randn(n*N,1)/sqrt(2);
r11 = complex(r11_r,r11_i);

r12_r = sqrt(0.79).*randn(n*N,1)/sqrt(2);%Channel 1
r12_i = sqrt(0.79).*randn(n*N,1)/sqrt(2);
r12 = complex(r12_r,r12_i);
r13_r = sqrt(0.126).*randn(n*N,1)/sqrt(2);
r13_i = sqrt(0.126).*randn(n*N,1)/sqrt(2);
r13 = complex(r13_r,r13_i);
r14_r = sqrt(0.1).*randn(n*N,1)/sqrt(2);
r14_i = sqrt(0.1).*randn(n*N,1)/sqrt(2);
r14 = complex(r14_r,r14_i);
r15_r = sqrt(0.032).*randn(n*N,1)/sqrt(2);
r15_i = sqrt(0.032).*randn(n*N,1)/sqrt(2);
r15 = complex(r15_r,r15_i);
r16_r = sqrt(0.01).*randn(n*N,1)/sqrt(2);
r16_i = sqrt(0.01).*randn(n*N,1)/sqrt(2);
r16 = complex(r16_r,r16_i);
tau11=1;         
tau12=30;
tau13=70;
tau14=110;
tau15=170;
tau16=250;


r21_r = sqrt(1).*randn(n*N,1)/sqrt(2);%Channel 2
r21_i = sqrt(1).*randn(n*N,1)/sqrt(2);
r21 = complex(r21_r,r21_i);
r22_r = sqrt(0.79).*randn(n*N,1)/sqrt(2);
r22_i = sqrt(0.79).*randn(n*N,1)/sqrt(2);
r22 = complex(r22_r,r22_i);
r23_r = sqrt(0.126).*randn(n*N,1)/sqrt(2);
r23_i = sqrt(0.126).*randn(n*N,1)/sqrt(2);
r23 = complex(r23_r,r23_i);
r24_r = sqrt(0.1).*randn(n*N,1)/sqrt(2);
r24_i = sqrt(0.1).*randn(n*N,1)/sqrt(2);
r24 = complex(r24_r,r24_i);
r25_r = sqrt(0.032).*randn(n*N,1)/sqrt(2);
r25_i = sqrt(0.032).*randn(n*N,1)/sqrt(2);
r25 = complex(r25_r,r25_i);
r26_r = sqrt(0.032).*randn(n*N,1)/sqrt(2);
r26_i = sqrt(0.032).*randn(n*N,1)/sqrt(2);
r26 = complex(r26_r,r26_i);
tau21=1;
tau22=30;
tau23=70;
tau24=110;
tau25=170;
tau26=250;

r31_r = sqrt(1).*randn(n*N,1)/sqrt(2); %Channel 3
r31_i = sqrt(1).*randn(n*N,1)/sqrt(2);
r31 = complex(r31_r,r31_i);
r32_r = sqrt(0.79).*randn(n*N,1)/sqrt(2);
r32_i = sqrt(0.79).*randn(n*N,1)/sqrt(2);
r32 = complex(r32_r,r32_i);
r33_r = sqrt(0.126).*randn(n*N,1)/sqrt(2);
r33_i = sqrt(0.126).*randn(n*N,1)/sqrt(2);
r33 = complex(r33_r,r33_i);
r34_r = sqrt(0.1).*randn(n*N,1)/sqrt(2);
r34_i = sqrt(0.1).*randn(n*N,1)/sqrt(2);
r34 = complex(r34_r,r34_i);
r35_r = sqrt(0.032).*randn(n*N,1)/sqrt(2);
r35_i = sqrt(0.032).*randn(n*N,1)/sqrt(2);
r35 = complex(r35_r,r35_i);
r36_r = sqrt(0.032).*randn(n*N,1)/sqrt(2);
r36_i = sqrt(0.032).*randn(n*N,1)/sqrt(2);
r36 = complex(r36_r,r36_i);
tau31=1;
tau32=30;
tau33=70;
tau34=110;
tau35=170;
tau36=250;


r41_r = sqrt(1).*randn(n*N,1)/sqrt(2); % Channel 4
r41_i = sqrt(1).*randn(n*N,1)/sqrt(2);
r41 = complex(r41_r,r41_i);
r42_r = sqrt(0.79).*randn(n*N,1)/sqrt(2);
r42_i = sqrt(0.79).*randn(n*N,1)/sqrt(2);
r42 = complex(r42_r,r42_i);
r43_r = sqrt(0.126).*randn(n*N,1)/sqrt(2);
r43_i = sqrt(0.126).*randn(n*N,1)/sqrt(2);
r43 = complex(r43_r,r43_i);
r44_r = sqrt(0.1).*randn(n*N,1)/sqrt(2);
r44_i = sqrt(0.1).*randn(n*N,1)/sqrt(2);
r44 = complex(r44_r,r44_i);
r45_r = sqrt(0.032).*randn(n*N,1)/sqrt(2);
r45_i = sqrt(0.032).*randn(n*N,1)/sqrt(2);
r45 = complex(r45_r,r45_i);
r46_r = sqrt(0.032).*randn(n*N,1)/sqrt(2);
r46_i = sqrt(0.032).*randn(n*N,1)/sqrt(2);
r46 = complex(r46_r,r46_i);

tau41=1;
tau42=30;
tau43=70;
tau44=110;
tau45=170;
tau46=250;

%%% Magnitude Response of the frequency selective channel
mag_res=zeros(64,1);

for i=1:1:length(SNR)
    for t=1:1:n
        
        h11=zeros(64,1);
        h12=zeros(64,1);
        h21=zeros(64,1);
        h22=zeros(64,1);
        
        
        
        h11(tau11)=r11(i*t);
        h11(tau12)=r12(i*t);
        h11(tau13)=r13(i*t);
        h11(tau14)=r14(i*t);
        h11(tau15)=r15(i*t);
        h11(tau16)=r16(i*t);
        
        h12(tau21)=r21(i*t);
        h12(tau22)=r22(i*t);
        h12(tau23)=r23(i*t);
        h12(tau24)=r24(i*t);
        h12(tau25)=r25(i*t);
        h12(tau26)=r26(i*t);
        
        h21(tau31)=r31(i*t);
        h21(tau32)=r32(i*t);
        h21(tau33)=r33(i*t);
        h21(tau34)=r34(i*t);
        h21(tau35)=r35(i*t);
        h21(tau36)=r36(i*t);
        
        h22(tau41)=r41(i*t);
        h22(tau42)=r42(i*t);
        h22(tau43)=r43(i*t);
        h22(tau44)=r44(i*t);
        h22(tau45)=r45(i*t);
        h22(tau46)=r46(i*t);
        
        h11_fft=abs(fft(h11,64));   %Frequency Response of the channel
        h12_fft=abs(fft(h12,64));
        h21_fft=abs(fft(h21,64));
        h22_fft=abs(fft(h22,64));
        
        
        
        %%%Calculating svd
        H=zeros(1,128);
        j=0;
        for k=1:64
        H((1+(2*j)):(2+(2*j)))=svd([h11_fft(k) h12_fft(k);h21_fft(k) h22_fft(k)]);
        j=j+1;
        end
        m=1;
        Hmax=zeros(1,64);
        for k=1:2:128
        Hmax(m)= max(H(k:k+1));
        m=m+1;
        end
        
        
        count=0;
        
        No=1/(10^(SNR(i)/10));
        
        for r=9:1:56
            mag_res(r)=10*log10((Hmax(r)^2)*(1/No));
            if ((mag_res(r)>=7) && (mag_res(r)<10))             %BPSK
                count=count+1;
            elseif ((mag_res(r)>=10) && (mag_res(r)<15.25))     %QPSK
                count=count+2;
            elseif ((mag_res(r)>=15.25) && (mag_res(r)<18))     %8-PSK
                count=count+3;
            elseif ((mag_res(r)>=18) && (mag_res(r)<23.5))      %16-QAM
                count=count+4;
            elseif (mag_res(r)>=23.5)                        %64-QAM     
                count=count+6;    
            end
        end
        N_bits(i)=count;
        eff(i)=eff(i)+N_bits(i);
    end
    eff(i)=eff(i)/n;
end
b_vehicular=eff/48;

plot(SNR,b_vehicular,'-b*');
grid on; hold on;


%% Pedestrian %%
         
%Channel Parameters for a 6-tap model

r11_r = sqrt(1).*randn(n*N,1)/sqrt(2);
r11_i = sqrt(1).*randn(n*N,1)/sqrt(2);
r11 = complex(r11_r,r11_i);
r21_r = sqrt(1).*randn(n*N,1)/sqrt(2);
r21_i = sqrt(1).*randn(n*N,1)/sqrt(2);
r21 = complex(r21_r,r21_i);
r31_r = sqrt(1).*randn(n*N,1)/sqrt(2);
r31_i = sqrt(1).*randn(n*N,1)/sqrt(2);
r31 = complex(r31_r,r31_i);
r41_r = sqrt(1).*randn(n*N,1)/sqrt(2);
r41_i = sqrt(1).*randn(n*N,1)/sqrt(2);
r41 = complex(r41_r,r41_i);

r12_r = sqrt(0.107).*randn(n*N,1)/sqrt(2);
r12_i = sqrt(0.107).*randn(n*N,1)/sqrt(2);
r12 = complex(r12_r,r12_i);
r22_r = sqrt(0.107).*randn(n*N,1)/sqrt(2);
r22_i = sqrt(0.107).*randn(n*N,1)/sqrt(2);
r22 = complex(r22_r,r22_i);
r32_r = sqrt(0.107).*randn(n*N,1)/sqrt(2);
r32_i = sqrt(0.107).*randn(n*N,1)/sqrt(2);
r32 = complex(r32_r,r32_i);
r42_r = sqrt(0.107).*randn(n*N,1)/sqrt(2);
r42_i = sqrt(0.107).*randn(n*N,1)/sqrt(2);
r42 = complex(r42_r,r42_i);

r13_r = sqrt(0.012).*randn(n*N,1)/sqrt(2);
r13_i = sqrt(0.012).*randn(n*N,1)/sqrt(2);
r13 = complex(r13_r,r13_i);
r23_r = sqrt(0.012).*randn(n*N,1)/sqrt(2);
r23_i = sqrt(0.012).*randn(n*N,1)/sqrt(2);
r23 = complex(r23_r,r23_i);
r33_r = sqrt(0.012).*randn(n*N,1)/sqrt(2);
r33_i = sqrt(0.012).*randn(n*N,1)/sqrt(2);
r33 = complex(r33_r,r33_i);
r43_r = sqrt(0.012).*randn(n*N,1)/sqrt(2);
r43_i = sqrt(0.012).*randn(n*N,1)/sqrt(2);
r43 = complex(r43_r,r43_i);

r14_r = sqrt(0.00525).*randn(n*N,1)/sqrt(2);
r14_i = sqrt(0.00525).*randn(n*N,1)/sqrt(2);
r14 = complex(r14_r,r14_i);
r24_r = sqrt(0.00525).*randn(n*N,1)/sqrt(2);
r24_i = sqrt(0.00525).*randn(n*N,1)/sqrt(2);
r24 = complex(r24_r,r24_i);
r34_r = sqrt(0.00525).*randn(n*N,1)/sqrt(2);
r34_i = sqrt(0.00525).*randn(n*N,1)/sqrt(2);
r34 = complex(r34_r,r34_i);
r44_r = sqrt(0.00525).*randn(n*N,1)/sqrt(2);
r44_i = sqrt(0.00525).*randn(n*N,1)/sqrt(2);
r44 = complex(r44_r,r44_i);


tau11=1;         
tau12=110;
tau13=190;
tau14=410;


tau21=1;
tau22=110;
tau23=190;
tau24=410;


tau31=1;
tau32=110;
tau33=190;
tau34=410;


tau41=1;
tau42=110;
tau43=190;
tau44=410;



%%% Magnitude Response of the frequency selective channel
mag_res=zeros(64,1);

for i=1:1:length(SNR)
    for t=1:1:n
        
       h11=zeros(64,1);
        h12=zeros(64,1);
        h21=zeros(64,1);
        h22=zeros(64,1);
        
        
        
        h11(tau11)=r11(i*t);
        h11(tau12)=r12(i*t);
        h11(tau13)=r13(i*t);
        h11(tau14)=r14(i*t);
       
        
        h12(tau21)=r21(i*t);
        h12(tau22)=r22(i*t);
        h12(tau23)=r23(i*t);
        h12(tau24)=r24(i*t);
        
        
        h21(tau31)=r31(i*t);
        h21(tau32)=r32(i*t);
        h21(tau33)=r33(i*t);
        h21(tau34)=r34(i*t);
        
        
        h22(tau41)=r41(i*t);
        h22(tau42)=r42(i*t);
        h22(tau43)=r43(i*t);
        h22(tau44)=r44(i*t);
       
        
        h11_fft=abs(fft(h11,64));   %Frequency Response of the channel
        h12_fft=abs(fft(h12,64));
        h21_fft=abs(fft(h21,64));
        h22_fft=abs(fft(h22,64));
        
        
        
        %%%Calculating svd
        H=zeros(1,128);
        j=0;
        for k=1:64
        H((1+(2*j)):(2+(2*j)))=svd([h11_fft(k) h12_fft(k);h21_fft(k) h22_fft(k)]);
        j=j+1;
        end
        m=1;
        Hmax=zeros(1,64);
        for k=1:2:128
        Hmax(m)= max(H(k:k+1));
        m=m+1;
        end
        
        count=0;
        
        No=1/(10^(SNR(i)/10));
        
        for r=9:1:56
            mag_res(r)=10*log10((Hmax(r)^2)*(1/No));
            if ((mag_res(r)>=7) && (mag_res(r)<10))             %BPSK
                count=count+1;
            elseif ((mag_res(r)>=10) && (mag_res(r)<15.25))     %QPSK
                count=count+2;
            elseif ((mag_res(r)>=15.25) && (mag_res(r)<18))     
                count=count+3;
            elseif ((mag_res(r)>=18) && (mag_res(r)<23.5))     %16-QAM
                count=count+4;
            elseif (mag_res(r)>=23.5)                          %64-QAM     
                count=count+6;    
            end
        end
        N_bits(i)=count;
        eff(i)=eff(i)+N_bits(i);
    end
    eff(i)=eff(i)/n;
end
b_pedestrian=eff/48;



plot(SNR,b_pedestrian,'-r^');
axis([10 30 2 6.4]);
title('Spectral Efficiency for 2*2 MIMO OFDM system for different channel models');
xlabel('SNR  (dB)');
ylabel('Spectral Eficiency');
legend('Vehicular','Pedestrian');
grid on; hold off;
