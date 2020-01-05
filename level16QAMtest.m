clc;

close all;

input_bits = 10^5; % number of symbols

nConst = 16; % no. of constellation points

k = log2(nConst); % symbol bits

 

% real and imaginary constellation for 16-QAM

ConstRe = [-(2*sqrt(nConst)/2-1):2:-1 1:2:2*sqrt(nConst)/2-1];

ConstIm = [-(2*sqrt(nConst)/2-1):2:-1 1:2:2*sqrt(nConst)/2-1];

c = 1/sqrt(10);

 

Eb_N_dB = [0:15]; % generation of Es/N values

Es_N_dB = Eb_N_dB + 10*log10(k);

 

% Symbol Mapping for code conversion

base = [0:k-1];
SymMap = bitxor(base,floor(base/2));
[tt ind] = sort(SymMap);   

 

for i = 1:length(Eb_N_dB)

  

% generation of symbols

% ------------------

inBits = rand(1,input_bits*k,1)>0.5; % genarating 1's and 0's randomly

inBitsReshape = reshape(inBits,k,input_bits).';

bin2DecMat = ones(input_bits,1)*(2.^[(k/2-1):-1:0]) ; % binary to decimal conversion

% generating real part

inBitsRe = inBitsReshape(:,[1:k/2]);

inDecRe = sum(inBitsRe.*bin2DecMat,2);

inGrayDecRe = bitxor(inDecRe,floor(inDecRe/2));

% genaration of imaginary part

inBitsIm = inBitsReshape(:,[k/2+1:k]);

inDecIm = sum(inBitsIm.*bin2DecMat,2);

inGrayDecIm = bitxor(inDecIm,floor(inDecIm/2));

% mapping of symbols using defined constellation

Re = ConstRe(inGrayDecRe+1);

Im = ConstIm(inGrayDecIm+1);

% complex constellation

mod = Re + j*Im;

np = c*mod; % transmit power normalization means equated to unity

  

% noise

% -----

n = 1/sqrt(2)*[randn(1,input_bits) + j*randn(1,input_bits)]; % genarating white guassian noise with 0dB variance
an = np + 10^(-Es_N_dB(i)/20)*n; % genaration of additive white gaussian noise

 

% demodulation

% ------------

an_re = real(an)/c; % Seperating real part

an_im = imag(an)/c; % Seperating imaginary part

 

% rounding to the result to the nearest alphabet

inHRe = 2*floor(an_re/2)+1;

inHRe(find(inHRe>max(ConstRe))) = max(ConstRe);

inHRe(find(inHRe<min(ConstRe))) = min(ConstRe);

inHIm = 2*floor(an_im/2)+1;

inHIm(find(inHIm>max(ConstIm))) = max(ConstIm);

inHIm(find(inHIm<min(ConstIm))) = min(ConstIm);

 

% Constellation mapping conversion to decimal

inDecHRe = ind(floor((inHRe+4)/2+1))-1;

inDecHIm = ind(floor((inHIm+4)/2+1))-1;

 

% converting mapping conversion to binary string

inBinHRe = dec2bin(inDecHRe,k/2);

inBinHIm = dec2bin(inDecHIm,k/2);

 

% converting binary string to number

inBinHRe = inBinHRe.';

inBinHRe = inBinHRe(1:end).';

inBinHRe = reshape(str2num(inBinHRe).',k/2,input_bits).' ;

  

inBinHIm = inBinHIm.';

inBinHIm = inBinHIm(1:end).';

inBinHIm = reshape(str2num(inBinHIm).',k/2,input_bits).' ;

 

% error calculation for real and imaginary

Err(i) = size(find([inBitsRe- inBinHRe]),1) + size(find([inBitsIm - inBinHIm]),1) ;

 

end

BERSim = Err/(input_bits*k);
%Pe_th(i)=Qfunc(sqrt(2*EbN0));

BERTh = 3*Qfunc(sqrt(0.8*(10.^(Eb_N_dB/10))));

 

close all; figure

semilogy(Eb_N_dB,BERTh,'mx-','LineWidth',2);

hold on

semilogy(Eb_N_dB,BERSim,'bs-','LineWidth',2);

axis([0 15 10^-5 1])

grid on

legend('theory', 'simulation');

xlabel('Eb/N, dB')

ylabel('Bit Error Rate')

title('Bit Error Performance of 16QAM System')