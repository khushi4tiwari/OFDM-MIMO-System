 

clear all;
clc;
% Vehicular%

% OFDM parameters
Nfft=64;
Nused=52;
Nref=4;
Ndata=48;
Nright=5;
Nleft =6;
Ndc=1;
Ncp=16; %CP length 
Nnull=12;
Ppos=13;


SNR=3:2:11;
%% --------------BPSK-----%% 
                            

l=1;
error=zeros(1,length(SNR));
M_BPSK=2;

% channel parameters for 6-tap model

t1=1;
t2=30;
t3=70;
t4=110;
t5=170;
t6=250;

for i=3:2:11
    Tx=[];
    init=[];
    rx_sym=[];
    Count_error=0;
    
    for blocks=1:80
        ip=randint(1,Ndata,M_BPSK);
        
        modsignal=pskmod(ip,M_BPSK);
        pilot=13:Ppos:Nused;
        sym=setxor(1:Nused,pilot);
        tx_sym(sym)=modsignal;
        tx_sym(pilot)=1;
        Nleft_sym=zeros(1,Nleft);
        Nright_sym=zeros(1,Nright);
        n=length(tx_sym)/2;
        b=1:n;
        c=(n+1):length(tx_sym);
        sym1=zeros(1,n);
        sym2=zeros(1,n);
        sym1=tx_sym(b);
        sym2=tx_sym(c);
        tx_seq=[Nleft_sym sym1 Ndc sym2 Nright_sym];
        
        sertopar_sym=reshape(tx_seq,length(tx_seq),1);
        tx_sym_IFFT=sqrt(Nfft)*ifft(sertopar_sym,Nfft);
        j=Nfft-Ncp;
        
        CP=tx_sym_IFFT(j+1:Nfft,1); % generate cyclic prefix %
        tx_IFFT_CP=[CP
            tx_sym_IFFT]; %appending CP %
        
        partoser_sym=reshape(tx_IFFT_CP,1,80);
        Tx=[Tx partoser_sym];
        
        %Generating Channel model using Channel A
        p1real=sqrt(1).*randn(1,1)/sqrt(2);
        p1imag=sqrt(1).*randn(1,1)/sqrt(2);
        p1=complex(p1real,p1imag);
        
        p2real=sqrt(0.8).*randn(1,1)/sqrt(2);
        p2imag=sqrt(0.8).*randn(1,1)/sqrt(2);
        p2=complex(p2real,p2imag);
        
        p3real=sqrt(0.125).*randn(1,1)/sqrt(2);
        p3imag=sqrt(0.125).*randn(1,1)/sqrt(2);
        p3=complex(p3real,p3imag);
        
        p4real=sqrt(0.1).*randn(1,1)/sqrt(2);
        p4imag=sqrt(0.1).*randn(1,1)/sqrt(2);
        p4=complex(p4real,p4imag);
        
        p5real=sqrt(0.031).*randn(1,1)/sqrt(2);
        p5imag=sqrt(0.031).*randn(1,1)/sqrt(2);
        p5=complex(p5real,p5imag);
        
        p6real=sqrt(0.01).*randn(1,1)/sqrt(2);
        p6imag=sqrt(0.01).*randn(1,1)/sqrt(2);
        p6=complex(p6real,p6imag);
       
        z=zeros(1,7);
        z(t1)=p1;
        z(t2)=p2;
        z(t3)=p3;
        z(t4)=p4;
        z(t5)=p5;
        z(t6)=p6;
             
        z_fft=(1/sqrt(Nfft))*fft(z,Nfft);
        tx_sym_freqchannel=conv(z,partoser_sym);
        d=tx_sym_freqchannel;
        
        % adding noise %
        tx_sym_freqchannel = tx_sym_freqchannel(1:Nfft+Ncp);
        N0 = 1/(10.^(i/10));
        n = sqrt(N0)*(randn(1,length(tx_IFFT_CP))+(1i*randn(1,length(tx_IFFT_CP))))/sqrt(2);
        tx_sym_freqchannelAWGN = tx_sym_freqchannel + n;
        
       
        %Receiver 
         tx_sym_freqchannelAWGN =  tx_sym_freqchannelAWGN(Ncp + 1:end);
         tx_symfreqchanneleq = ifft(fft( tx_sym_freqchannelAWGN)./z_fft);
         
         sertopar_sym_rx = reshape(tx_symfreqchanneleq,Nfft,1);
         rx_fft = (1/sqrt(Nfft))*fft(sertopar_sym_rx,Nfft);
         partoser_sym_rx = reshape(rx_fft,1,Nfft);
         
         rx_Nleft =  partoser_sym_rx(Nleft+1:end);
         diff = Nfft-Nright-Nleft;
         rx_Nright = rx_Nleft(1:diff);  %removing Nleft and Nright
         
         y = rx_Nright(1:(Nused/2));
         z = rx_Nright((Nused/2)+2:end);
         rx = [y z];
         
         rx_pilot = rx(sym);
         demod=[];
         op = pskdemod(rx_pilot,M_BPSK);
         
         demod = [demod op];
         Count_error = Count_error + biterr(op,ip);
    end
    error(l) = Count_error/(Ndata*1000);
    l = l+1;
    
end


figure(1);
semilogy(SNR,error,'-r^');
hold on;
grid on;

  
%% QPSK %%%


                           
l=1;
error=zeros(1,length(SNR));
M_QPSK=4;

% channel parameters for 6-tap model

t1=1;
t2=30;
t3=70;
t4=110;
t5=170;
t6=250;

for i=3:2:11
    Tx=[];
    init=[];
    rx_sym=[];
    Count_error=0;
    
    for blocks=1:1000
        ip=randint(1,Ndata,M_QPSK);
        
        modsignal=pskmod(ip,M_QPSK);
        pilot=13:Ppos:Nused;
        sym=setxor(1:Nused,pilot);
        tx_sym(sym)=modsignal;
        tx_sym(pilot)=1;
        Nleft_sym=zeros(1,Nleft);
        Nright_sym=zeros(1,Nright);
        n=length(tx_sym)/2;
        b=1:n;
        c=(n+1):length(tx_sym);
        sym1=zeros(1,n);
        sym2=zeros(1,n);
        sym1=tx_sym(b);
        sym2=tx_sym(c);
        tx_seq=[Nleft_sym sym1 Ndc sym2 Nright_sym];
        
        sertopar_sym=reshape(tx_seq,length(tx_seq),1);
        tx_sym_IFFT=sqrt(Nfft)*ifft(sertopar_sym,Nfft);
        j=Nfft-Ncp;
        
        CP=tx_sym_IFFT(j+1:Nfft,1); % generate cyclic prefix %
        tx_IFFT_CP=[CP
            tx_sym_IFFT]; %appending CP %
        
        partoser_sym=reshape(tx_IFFT_CP,1,length(tx_IFFT_CP));
        Tx=[Tx partoser_sym];
        
        p1real=sqrt(1).*randn(1,1)/sqrt(2);
        p1imag=sqrt(1).*randn(1,1)/sqrt(2);
        p1=complex(p1real,p1imag);
        
        p2real=sqrt(0.8).*randn(1,1)/sqrt(2);
        p2imag=sqrt(0.8).*randn(1,1)/sqrt(2);
        p2=complex(p2real,p2imag);
        
        p3real=sqrt(0.125).*randn(1,1)/sqrt(2);
        p3imag=sqrt(0.125).*randn(1,1)/sqrt(2);
        p3=complex(p3real,p3imag);
        
        p4real=sqrt(0.1).*randn(1,1)/sqrt(2);
        p4imag=sqrt(0.1).*randn(1,1)/sqrt(2);
        p4=complex(p4real,p4imag);
        
        p5real=sqrt(0.031).*randn(1,1)/sqrt(2);
        p5imag=sqrt(0.031).*randn(1,1)/sqrt(2);
        p5=complex(p5real,p5imag);
        
        p6real=sqrt(0.01).*randn(1,1)/sqrt(2);
        p6imag=sqrt(0.01).*randn(1,1)/sqrt(2);
        p6=complex(p6real,p6imag);
       
        z=zeros(1,7);
        z(t1)=p1;
        z(t2)=p2;
        z(t3)=p3;
        z(t4)=p4;
        z(t5)=p5;
        z(t6)=p6;
             
        z_fft=(1/sqrt(Nfft))*fft(z,Nfft);
        tx_sym_freqchannel=conv(z,partoser_sym);
        d=tx_sym_freqchannel;
        
        % adding noise %
        tx_sym_freqchannel = tx_sym_freqchannel(1:Nfft+Ncp);
        N0 = 1/(10.^(i/10));
        n = sqrt(N0)*(randn(1,length(tx_IFFT_CP))+(1i*randn(1,length(tx_IFFT_CP))))/sqrt(2);
        tx_sym_freqchannelAWGN = tx_sym_freqchannel + n;
        
        %Receiver
        
         tx_sym_freqchannelAWGN =  tx_sym_freqchannelAWGN(Ncp + 1:end);
         tx_symfreqchanneleq = ifft(fft( tx_sym_freqchannelAWGN)./z_fft);
         
         sertopar_sym_rx = reshape(tx_symfreqchanneleq,Nfft,1);
         rx_fft = (1/sqrt(Nfft))*fft(sertopar_sym_rx,Nfft);
         partoser_sym_rx = reshape(rx_fft,1,Nfft);
         
         rx_Nleft =  partoser_sym_rx(Nleft+1:end);
         diff = Nfft-Nright-Nleft;
         rx_Nright = rx_Nleft(1:diff);  %removing Nleft and Nright
         
         y = rx_Nright(1:(Nused/2));
         z = rx_Nright((Nused/2)+2:end);
         rx = [y z];
         
         rx_pilot = rx(sym);
         demod=[];
         op = pskdemod(rx_pilot,M_QPSK);
        
         demod = [demod op];
         Count_error = Count_error + biterr(op,ip);
    end
    error(l) = Count_error/(Ndata*2000);
    l = l+1;
    
end


figure(1);
semilogy(SNR,error,'-b^');

grid on;
hold on;
legend('BPSK','QPSK');


                %%%%%%%%%%%%%%16-QAM%%%%%%%%%%%%%%%%%%
         
       
                
l=1;
error=zeros(1,length(SNR));
M_QAM=16;

% channel parameters for 6-tap model

t1=1;
t2=30;
t3=70;
t4=110;
t5=170;
t6=250;

for i=3:2:11
    Tx=[];
    init=[];
    rx_sym=[];
    Count_error=0;
    
    for blocks=1:1000
        ip=randint(1,Ndata,M_QAM);
        
        modsignal=qammod(ip,M_QAM)/sqrt(10); %normalise by sqrt(10)
        pilot=13:Ppos:Nused;
        sym=setxor(1:Nused,pilot);
        tx_sym(sym)=modsignal;
        tx_sym(pilot)=1;
        Nleft_sym=zeros(1,Nleft);
        Nright_sym=zeros(1,Nright);
        n=length(tx_sym)/2;
        b=1:n;
        c=(n+1):length(tx_sym);
        sym1=zeros(1,n);
        sym2=zeros(1,n);
        sym1=tx_sym(b);
        sym2=tx_sym(c);
        tx_seq=[Nleft_sym sym1 Ndc sym2 Nright_sym];
        
        sertopar_sym=reshape(tx_seq,length(tx_seq),1);
        tx_sym_IFFT=sqrt(Nfft)*ifft(sertopar_sym,Nfft);
        j=Nfft-Ncp;
        
        CP=tx_sym_IFFT(j+1:Nfft,1); % generate cyclic prefix %
        tx_IFFT_CP=[CP
            tx_sym_IFFT]; %appending CP %
        
        partoser_sym=reshape(tx_IFFT_CP,1,length(tx_IFFT_CP));
        Tx=[Tx partoser_sym];
        
        p1real=sqrt(1).*randn(1,1)/sqrt(2);
        p1imag=sqrt(1).*randn(1,1)/sqrt(2);
        p1=complex(p1real,p1imag);
        
        p2real=sqrt(0.8).*randn(1,1)/sqrt(2);
        p2imag=sqrt(0.8).*randn(1,1)/sqrt(2);
        p2=complex(p2real,p2imag);
        
        p3real=sqrt(0.125).*randn(1,1)/sqrt(2);
        p3imag=sqrt(0.125).*randn(1,1)/sqrt(2);
        p3=complex(p3real,p3imag);
        
        p4real=sqrt(0.1).*randn(1,1)/sqrt(2);
        p4imag=sqrt(0.1).*randn(1,1)/sqrt(2);
        p4=complex(p4real,p4imag);
        
        p5real=sqrt(0.031).*randn(1,1)/sqrt(2);
        p5imag=sqrt(0.031).*randn(1,1)/sqrt(2);
        p5=complex(p5real,p5imag);
        
        p6real=sqrt(0.01).*randn(1,1)/sqrt(2);
        p6imag=sqrt(0.01).*randn(1,1)/sqrt(2);
        p6=complex(p6real,p6imag);
       
        z=zeros(1,7);
        z(t1)=p1;
        z(t2)=p2;
        z(t3)=p3;
        z(t4)=p4;
        z(t5)=p5;
        z(t6)=p6;
             
        z_fft=(1/sqrt(Nfft))*fft(z,Nfft);
        tx_sym_freqchannel=conv(z,partoser_sym);
        d=tx_sym_freqchannel;
        
        % adding noise %
        tx_sym_freqchannel = tx_sym_freqchannel(1:Nfft+Ncp);
        N0 = 1/(10.^(i/10));
        n = sqrt(N0)*(randn(1,length(tx_IFFT_CP))+(1i*randn(1,length(tx_IFFT_CP))))/sqrt(2);
        tx_sym_freqchannelAWGN = tx_sym_freqchannel + n;
        
        %Receiver
        
         tx_sym_freqchannelAWGN =  tx_sym_freqchannelAWGN(Ncp + 1:end);
         tx_symfreqchanneleq = ifft(fft( tx_sym_freqchannelAWGN)./z_fft);
         
         sertopar_sym_rx = reshape(tx_symfreqchanneleq,Nfft,1);
         rx_fft = (1/sqrt(Nfft))*fft(sertopar_sym_rx,Nfft);
         partoser_sym_rx = reshape(rx_fft,1,Nfft);
         
         rx_Nleft =  partoser_sym_rx(Nleft+1:end);
         diff = Nfft-Nright-Nleft;
         rx_Nright = rx_Nleft(1:diff);  %removing Nleft and Nright
         
         y = rx_Nright(1:(Nused/2));
         z = rx_Nright((Nused/2)+2:end);
         rx = [y z];
         
         rx_pilot = rx(sym);
         rx_pilot = rx_pilot*sqrt(10);
         demod=[];
         op = qamdemod(rx_pilot,M_QAM);
         
         demod = [demod op];
         Count_error = Count_error + biterr(op,ip);
    end
    error(l) = Count_error/(Ndata*4000);
    l = l+1;
    
end


figure(1);
semilogy(SNR,error,'-g^');
grid on;
hold on;
legend('BPSK','QPSK','16 QAM');
title('BER Plot for Vehicular Multipath Channel Profile - Channel A');
xlabel('SNR (dB)');
ylabel('Bit Error Rate');
