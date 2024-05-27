TxPwtot = 20*10^3;                   % Total Transmit Power [mwatts]
TxPwtotdBm = 10*log10(TxPwtot);
TxPwPerUsdBm = 33;                     % Transmit power for each link [dBm]
TxPwPerUs = 10^(TxPwPerUsdBm/10);

TxPwtot = 10^2;                   % Total Transmit Power [mwatts]
TxPwtotdBm = 10*log10(TxPwtot);
TxPwPerUsdBm = 21;                     % Transmit power for each link [dBm]
TxPwPerUs = 10^(TxPwPerUsdBm/10);
H = [];
itti =30;
%%%%%%%%%%%%%%%%%%%%%%%%%%%

RatePerUs = [] ;

fmRatePerUs = [] ;


[fmRatePerUs]= calculproprefemtocell(TxPwtot, TxPwPerUsdBm);
[RatePerUs] = calculpropreHSDPA(TxPwtot, TxPwPerUsdBm);

for iitti= 1: itti 
    
if RatePerUs(iitti)<= fmRatePerUs(iitti)
    H(iitti)=0 ;
else
    H(iitti)=1 ;
end
    
end 
H
figure(1)
hold on
plot(RatePerUs,'m')
plot(fmRatePerUs,'g')
figure (2)
hold on 

plot (H,'b +')