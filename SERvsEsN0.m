%Theoretical & Simulated SER for 16 QAM modulation

Eb = 1;                                 %Bit energy per symbol
M = 16;                                 %No of signal points
k = log2(M);                            %No of bits in each signal point

bit_input = 2*10^5;                       %No of symbols

%The signal constellation for QAM (2 4-PAM Signal)
%Real_Const = [-(2*sqrt(M)/2-1):2:-1 1:2:2*sqrt(M)/2-1];
%Imag_Const = [-(2*sqrt(M)/2-1):2:-1 1:2:2*sqrt(M)/2-1];
Real_Const = [-3 -1 1 3];
Esim_N_dB = [4.62: 20];
Eb_N_dB = [0:20];                       %Sampled Eb/No Values from 0 - 20 dB
Es_N_dB = Eb_N_dB + 10*log10(k);        %The Symbol Energy Noise ratio in dB

%THEORETICAL & Approximate COMPUTATION OF SYMBOL ERROR RATE (SER) USING Q FUNCTION

for i=1:length(Eb_N_dB)
    Eb_N = 10.^(0.1*Eb_N_dB(i));                                                %The Decimal value of bit SNR
    Es_N = 10.^(0.1*Es_N_dB(i));                                                %The Decimal value of Symbol SNR
    N0 = Eb/Eb_N;                                                               %The AWGN Noise 
    SER(i) = (3 * Qfunc(sqrt(0.8*Eb_N))) - (2.25*(Qfunc(sqrt(0.8*Eb_N)))^2);    %Theoretical SER
    SER_Approx(i) = (3 * Qfunc(sqrt(0.8*Eb_N)));                                %Approximate SER
    BER (i) = SER(i)/k;                                                         %The BER 
    SER_Ubound(i) = (M-1) * Qfunc(sqrt(0.2*Es_N));                              %The SER Union Bound
    
end


%SIMULATED COMPUTATION OF SYMBOL ERROR RATE (SER) USING RANDOMIZED INPUT
%BITS
ipHat = zeros(1,bit_input);
for ii = 1:length(Eb_N_dB)
    ip = randsrc(1,bit_input,Real_Const) + j*randsrc(1,bit_input,Real_Const);
    s = (1/sqrt(10))*ip; % Energy normalized to 1
    n = 1/sqrt(2)*[randn(1,bit_input) + j*randn(1,bit_input)]; % white guassian noise, 0dB variance

    y = s + 10^(-Eb_N_dB(ii)/20)*n; % additive white gaussian noise

    % demodulation
    y_re = real(y); % real part of the random noise generated
    y_im = imag(y); % imaginary part of the random noise generated

    ipHat_re(find(y_re< -2/sqrt(10)))           = -3; %Corresponding error points as relate to the received signal of the real part in respect to signal constellation
    ipHat_re(find(y_re > 2/sqrt(10)))           =  3;
    ipHat_re(find(y_re>-2/sqrt(10) & y_re<=0))  = -1;
    ipHat_re(find(y_re>0 & y_re<=2/sqrt(10)))   =  1;

    ipHat_im(find(y_im< -2/sqrt(10)))           = -3; %Corresponding error points as relate to the received signal of the Imaginary part respect to signal constellation
    ipHat_im(find(y_im > 2/sqrt(10)))           =  3;
    ipHat_im(find(y_im>-2/sqrt(10) & y_im<=0))  = -1;
    ipHat_im(find(y_im>0 & y_im<=2/sqrt(10)))   =  1;
    ipHat = ipHat_re + j*ipHat_im;
    nErr(ii) = size(find([ip- ipHat]),2);                           %The number of errors are counted
end

simBer = (nErr/bit_input);


figure(1)
semilogy(Es_N_dB, SER, '-rs','linewidth', 2); hold on           %Theoretical plot of SER
semilogy(Eb_N_dB, simBer, 'mx-','linewidth', 2); hold on        %Simulated plot of SER
semilogy(Es_N_dB, SER_Ubound, '-b*','linewidth', 2); hold on    %Union Bound plot of SER
semilogy(Es_N_dB, SER_Approx, '-g','linewidth', 2); hold on     %Approximate plot of SER

hold off
axis([6 27 10^-20 10^1]);       
grid on
hold off
title('Symbol Error rate for 16 Level QAM Modulation');

xlabel('Es/N0 in dB');
ylabel('Symbol Error Rate');
legend('SER','SER-Sim','SER Union Bound', 'SER Approximate');    %Legend used in differentiating plots
grid
    



