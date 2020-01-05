%Theoretical & Simulated SER for 16 QAM modulation

                                 
k = log2(M); 
Eb = 1;                                 
M = 16;

inputbits = 2*10^5;                       


Const = [-3 -1 1 3];
Eb_N_dB = [0:20];                      
Es_N_dB = Eb_N_dB + 10*log10(k);        


for i=1:length(Eb_N_dB)
    Eb_N = 10.^(0.1*Eb_N_dB(i));
    Es_N = 10.^(0.1*Es_N_dB(i));
    N0 = Eb/Eb_N;
    SER(i) = (3 * Qfunc(sqrt(0.8*Eb_N))) - (2.25*(Qfunc(sqrt(0.8*Eb_N)))^2);
    SER_App(i) = (3 * Qfunc(sqrt(0.8*Eb_N)));
    BER (i) = SER(i)/k;
    SER_bound(i) = (M-1) * Qfunc(sqrt(0.2*Es_N)); 
    
end



Hatip = zeros(1,inputbits);     %randomized bits
for ii = 1:length(Eb_N_dB)
    ipvalue = randsrc(1,inputbits,Const) + j*randsrc(1,inputbits,Const);
    s = (1/sqrt(10))*ipvalue; 
    noise = 1/sqrt(2)*[randn(1,inputbits) + j*randn(1,inputbits)]; 

    y = s + 10^(-Eb_N_dB(ii)/20)*noise; 

   
    y_re = real(y); % real part
    y_im = imag(y); % imaginary part

    Hatip_re(find(y_re< -2/sqrt(10)))           = -3;
    Hatip_re(find(y_re > 2/sqrt(10)))           =  3;
    Hatip_re(find(y_re>-2/sqrt(10) & y_re<=0))  = -1;
    Hatip_re(find(y_re>0 & y_re<=2/sqrt(10)))   =  1;

    Hatip_im(find(y_im< -2/sqrt(10)))           = -3;
    Hatip_im(find(y_im > 2/sqrt(10)))           =  3;
    Hatip_im(find(y_im>-2/sqrt(10) & y_im<=0))  = -1;
    Hatip_im(find(y_im>0 & y_im<=2/sqrt(10)))   =  1;
    Hatip = Hatip_re + j*Hatip_im;
    nErr(ii) = size(find([ipvalue- Hatip]),2); 
end

BerSimulation = (nErr/inputbits);


figure(1)
semilogy(Es_N_dB, SER, '-rs','linewidth', 2); hold on
semilogy(Eb_N_dB, BerSimulation, 'mx-','linewidth', 2); hold on
semilogy(Es_N_dB, SER_bound, '-b*','linewidth', 2); hold on
semilogy(Es_N_dB, SER_App, '-g','linewidth', 2); hold on

hold off      
grid on
hold off
title('Symbol Error rate');
grid
xlabel('Es/N0 in dB');
ylabel('Symbol Error Rate');
legend('SER','Simulation SER','SER Bound', 'SER Approx');




