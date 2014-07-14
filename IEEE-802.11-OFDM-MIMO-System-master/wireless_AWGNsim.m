 %%% AWGN Simulation%%%
clc;
clear all;
close all;

% OFDM parameters for IEEE 802.11a
Nfft=64;
Nused=52;
Nref=4;
Ndata=48;
Nright=5;
Nleft =6;
Ndc=1;
Ncp=16; 
Nnull=12;
Ppos=13;

%S for the AWGN Channel simulation
SNR=3:2:11;
%% --------------------BPSK--------------------%%

   %%%---Steps involved for Transmitter:
%1)Generate random input sequence and map the signals
%2)Generate 64 samples with Nused=52,Nleft=6,Nright=5,Ndc=1
%3)Converting this 64 point samples from Serial to Parallel and taking IFFT
%4)Generate 16 cyclic prefix and append to 64 point making total 80 samples
%5)Add AWGN to the Sample
% -----This is the OFDM transmitted signal
l=1;
error=zeros(1,length(SNR));

for i=3:2:11
    Tx=[];
    init=[];
    rx_sym=[];
    Count_error=0;
    
    for blocks=1:100
        ip=randint(1,Ndata,2);
        init = [init ip];
        
        signal=pskmod(ip,2);
        pilot=13:Ppos:Nused;
        sym=setxor(1:Nused,pilot);
        tx(sym)=signal;
        tx(pilot)=1;
        sym_left=zeros(1,6);%Nleft=6
        sym_right=zeros(1,5);%Nright=5
        n=length(tx)/2;
        b=1:n;
        c=(n+1):length(tx);
        sym1=zeros(1,n);
        sym2=zeros(1,n);
        sym1=tx(b);
        sym2=tx(c);
        tx_seq = zeros(1,64);
        tx_seq=[sym_left sym1 Ndc sym2 sym_right];
        
        S_to_P=reshape(tx_seq,64,1);%Serial to Parallel
        tx_sym_IFFT=sqrt(64)*ifft(S_to_P,64);
        j=Nfft-Ncp;
        
        CP=tx_sym_IFFT(j+1:64,1); % generate cyclic prefix %
        tx_IFFT_CP=[CP
            tx_sym_IFFT]; %appending CP %
        
        P_to_S=reshape(tx_IFFT_CP,1,80);%Parallel to serial
        Tx=[Tx P_to_S];
        
                
        % adding noise %
        
        N0 = 1/(10.^(i/10));
        n = sqrt(N0)*(randn(1,80)+(1i*randn(1,80)))/sqrt(2);
        P_to_S = P_to_S + n;
       
               
        %%%---BPSK Receiver
%----Steps involved in Receiver structure
  % 1)Removing 16 cyclic prefix Nleft,Nright and Ndc
  % 2)Taking FFT of 64 point sample
  % 3)Demodulate the signal and compare with the input transmitted
  %  signal and calculate the BER
          
         ifft_sy = P_to_S(Ncp+1:end); 
         sertopar_sym_rx = reshape(ifft_sy,64,1);
         rx_fft = (sqrt(64))*fft(sertopar_sym_rx,64);
         PtoS_rx = reshape(rx_fft,1,Nfft);
         
         rx_Nleft =  PtoS_rx(Nleft+1:end);
         diff = Nfft-Nright-Nleft;
         rx_Nright = rx_Nleft(1:diff);  %removing Nleft and Nright
         
         y = rx_Nright(1:(Nused/2));    %removing Ndc
         z = rx_Nright((Nused/2)+2:end);
         rx = [y z];
       %%%Signal demodulation   
         rx_pilot = rx(sym);
         demod=[];
         op = pskdemod(rx_pilot,2);
         demod = [demod op];
         Count_error = Count_error + biterr(op,ip);
    end
    error(l) = Count_error/(Ndata*100);
    l = l+1;
    
end

figure(1);
semilogy(SNR,error,'-r^');
hold on;
grid on;
title('BER plots for modulation techniques with AWGN');
xlabel('SNR (dB)');
ylabel('Bit Error Rate');
                            
  

%% ---------------------------QPSK------------------------------%%

l=1;
error=zeros(1,length(SNR));
M_QPSK=4;


for i=3:2:11
    Tx=[];
    init=[];
    rx_sym=[];
    Count_error=0;
    
    for blocks=1:100
        ip=randint(1,Ndata,M_QPSK);
        init = [init ip];
        
        signal=pskmod(ip,M_QPSK);
        pilot=13:Ppos:Nused;
        sym=setxor(1:Nused,pilot);
        tx(sym)=signal;
        tx(pilot)=1;
        sym_left=zeros(1,Nleft);
        sym_right=zeros(1,Nright);
        n=length(tx)/2;
        b=1:n;
        c=(n+1):length(tx);
        sym1=zeros(1,n);
        sym2=zeros(1,n);
        sym1=tx(b);
        sym2=tx(c);
        tx_seq = zeros(1,64);
        tx_seq=[sym_left sym1 Ndc sym2 sym_right];
        
        S_to_P=reshape(tx_seq,64,1);
        tx_sym_IFFT=sqrt(Nfft)*ifft(S_to_P,Nfft);
        j=Nfft-Ncp;
        
        CP=tx_sym_IFFT(j+1:Nfft,1); % generate cyclic prefix %
        tx_IFFT_CP=[CP
            tx_sym_IFFT]; %appending CP %
        
        P_to_S=reshape(tx_IFFT_CP,1,80);
        Tx=[Tx P_to_S];
        
                
        % adding noise %
        
        N0 = 1/(10.^(i/10));
        n = sqrt(N0)*(randn(1,length(tx_IFFT_CP))+(1i*randn(1,length(tx_IFFT_CP))))/sqrt(2);
        P_to_S = P_to_S + n;
       
                
         %Receiver
                
         ifft_sy = P_to_S(Ncp+1:end); 
         sertopar_sym_rx = reshape(ifft_sy,Nfft,1);
         rx_fft = (sqrt(Nfft))*fft(sertopar_sym_rx,Nfft);
         PtoS_rx = reshape(rx_fft,1,Nfft);
         
         rx_Nleft =  PtoS_rx(Nleft+1:end);
         diff = Nfft-Nright-Nleft;
         rx_Nright = rx_Nleft(1:diff);  %removing Nleft and Nright
         
         y = rx_Nright(1:(Nused/2));    %Removing Ndc
         z = rx_Nright((Nused/2)+2:end);
         rx = [y z];
         
         rx_pilot = rx(sym);
         demod=[];
         op = pskdemod(rx_pilot,M_QPSK);
         demod = [demod op];
         Count_error = Count_error + biterr(op,ip);
    end
    error(l) = Count_error/(Ndata*200);
    l = l+1;
    
end

figure(1);
semilogy(SNR,error,'-b^');
hold on;
grid on;
title('BER plots for modulation techniques with AWGN');
xlabel('SNR (dB)');
ylabel('Bit Error Rate');

%% ------------------------ 16-QAM-------------------------------%%
        
                
l=1;
error=zeros(1,length(SNR));
M_QAM=16;


for i=3:2:11
    Tx=[];
    init=[];
    rx_sym=[];
    Count_error=0;
    
    for blocks=1:100
        ip=randint(1,Ndata,M_QAM);
        init = [init ip];
        
        signal=qammod(ip,M_QAM)/sqrt(42);  
        pilot=13:Ppos:Nused;
        sym=setxor(1:Nused,pilot);
        tx(sym)=signal;
        tx(pilot)=1;
        sym_left=zeros(1,Nleft);
        sym_right=zeros(1,Nright);
        n=length(tx)/2;
        b=1:n;
        c=(n+1):length(tx);
        sym1=zeros(1,n);
        sym2=zeros(1,n);
        sym1=tx(b);
        sym2=tx(c);
        tx_seq = zeros(1,64);
        tx_seq=[sym_left sym1 Ndc sym2 sym_right];
        
        S_to_P=reshape(tx_seq,64,1);
        tx_sym_IFFT=sqrt(Nfft)*ifft(S_to_P,Nfft);
        j=Nfft-Ncp;
        
        CP=tx_sym_IFFT(j+1:Nfft,1); % generate cyclic prefix %
        tx_IFFT_CP=[CP
            tx_sym_IFFT]; %appending CP %
        
        P_to_S=reshape(tx_IFFT_CP,1,80);
        Tx=[Tx P_to_S];
        
                
        % adding noise %
        
        N0 = 1/(10.^(i/10));
        n = sqrt(N0)*(randn(1,80)+(1i*randn(1,80)))/sqrt(2);
        P_to_S = P_to_S + n;
       
               
         %Receiver
               
         ifft_sy = P_to_S(Ncp+1:end); 
         sertopar_sym_rx = reshape(ifft_sy,Nfft,1);
         rx_fft = (sqrt(Nfft))*fft(sertopar_sym_rx,Nfft);
         PtoS_rx = reshape(rx_fft,1,Nfft);
         
         rx_Nleft =  PtoS_rx(Nleft+1:end);
         diff = Nfft-Nright-Nleft;
         rx_Nright = rx_Nleft(1:diff);  %removing Nleft and Nright
         
         y = rx_Nright(1:(Nused/2));
         z = rx_Nright((Nused/2)+2:end);
         rx = [y z];
         
         rx_pilot = rx(sym);
         rx_pilot = rx_pilot*sqrt(42);
         demod=[];
         op = qamdemod(rx_pilot,M_QAM);
         
         demod = [demod op];
         Count_error = Count_error + biterr(op,ip);
    end
    error(l) = Count_error/(Ndata*400);
    l = l+1;
    
end

figure(1);
semilogy(SNR,error,'-m^');
hold on;
grid on;
title('BER plots for modulation techniques with AWGN');
legend('BPSK','QPSK','16-QAM');
xlabel('SNR  (dB)');
ylabel('Bit Error Rate');
