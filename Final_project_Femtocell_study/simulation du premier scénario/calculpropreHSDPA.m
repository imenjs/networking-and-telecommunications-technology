                                                                     
                                                                     
                                                                 





function [RatePerUs] = calculproprefemtocell(TxPwtot,TPerUsdBmxPw)
% Simulation  parameters setting

%--------- Base station Model
TxPwtot = 20*10^3;                   % Total Transmit Power [mwatts]
TxPwtotdBm = 10*log10(TxPwtot);
TxPwPerUsdBm= 33;                     % Transmit power for each link [dBm]
TxPwPerUs = 10^(TxPwPerUsdBm/10);
BW = 5*10^6;                       % BandWidth [Hz]
TTI = 2*10^-3;                          % [s]
SF = 16;                           % Spreading Factor
TotCodes = 15;                      % Number of total available codes
OrthoFact = 0.5;                  % Orthogonality factor
PwComCH = 2e3;                          % Power allocated to commun channels (CPICH, SCH,...) [Watt]
PwComCHdB = 10*log10(PwComCH);          % [dBm]

%--------- Propagation conditions
N0dBmPerH = -174 ;                  % [dBm/Hz]
NoisedBm = N0dBmPerH + 10*log10(BW);
NoisePw = 10^(NoisedBm/10);
StdShdB = 8;                        % Standard deviation of the Shadowing [dB]
StdSh = 10^(StdShdB/10);

%-------- System layout and parameters
%== Us, Bs
NoSpMax = 1;           % Number of snapshots
NoSimulStpMax = 1;   % Oservation period
itti = 30;              % Number of Mobiles (Users)

NoUsVct = [itti];
Radius = 1.5 ;         % Radius of the circular cell [m]
NoBs = 6;              % Number of Base Stations
MinReqRatePerUs = 3;       % Minimum offered rate per user [Mbps]
MaxReqRatePerUs = 5;       % Minimum offered rate per user [Mbps]
%
UsActVct = ones(1,itti);    %users' activity (here, all users are active)
ActFactMean = 1;
%
%==Position of the 19 Bases Stations
XBs=[ 2/3,3/2*Radius, -3/2*Radius, 1*Radius, 1/2*Radius, -1/2*Radius ];
YBs=[ sqrt(1/4),sqrt(3)/2*Radius, sqrt(3)*Radius, -sqrt(3)/2*Radius, -sqrt(3)*Radius, -sqrt(3)/2*Radius ];
Bs(:,1)=XBs'; Bs(:,2)=YBs';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               %
%       Main Program            %
%                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for iSnap = 1: NoSpMax

    XMs = [];
    YMs = [];
    dBsiUsj = [];
    Ms= [];
    SINRperUs = [];
    SINRperUsdB = [];
  
    % Initilisation of the position of the Mobiles Ms(XMs,YMs)
    
   anglMs = rand(1,itti);
    rMs = Radius.*rand(1,itti);
    XMs = 0.0001.*cos(anglMs);
    YMs = 0.0003.*sin(anglMs);
    Ms(:,1)=XMs'; Ms(:,2)=YMs'; 
    Ms(:,1)=XMs'; Ms(:,2)=YMs'; % useful for MsPosTest
    %plot(XMs,YMs,'m+')


    %___ Compute the SINR for each PRIMARY user

    % Compute the link gain for all users in the cell of reference
    for iBs = 1: NoBs
        dBsiUsj(iBs,:) = sqrt((XMs-XBs(iBs)).^2+(YMs-YBs(iBs)).^2);
    end
    xPthLssVctdB = 137.4 + 35.2.*log10(dBsiUsj);
    %ShadowdB = StdShdB*randn(NoBs,itti+1);
    LinkBsiMsidB = xPthLssVctdB; % + ShadowdB;
    LinkBsiMsi = 10.^(LinkBsiMsidB/10);

    % Compute_SINR for each primary user
    RxPwdBm_Bs1 = TxPwPerUsdBm - LinkBsiMsidB(1,:);
    RxPw_Bs1= 10.^(RxPwdBm_Bs1/10);

    TotalTxPw_Bs1 = sum(TxPwPerUs*UsActVct)+PwComCH;
    TotalTxPw_Bs1dBm = 10*log10(TotalTxPw_Bs1);

    IntraInterf = (1-OrthoFact)*TotalTxPw_Bs1./LinkBsiMsi(1,:);

    RxPwdBm_Bsk = TxPwtotdBm*ActFactMean - LinkBsiMsidB(2:NoBs,1:itti); %The neiboring cells are "not" loaded.
    RxPw_Bsk= 10.^(RxPwdBm_Bsk/10);

    for iUs = 1: itti
        InterInterf(iUs)= sum(RxPw_Bsk(:,iUs));
    end

    %--- SINR Computation
    SINRperUs= RxPw_Bs1./(IntraInterf + InterInterf + NoisePw);
    SINRperUsdB = 10.*log10(SINRperUs);

    %____________ Applying Modulation/Coding according to the SINR value

    for iUs = 1:itti
        if SINRperUsdB(iUs) <= -16
            CQIPerUs(iUs) = 0;
        elseif SINRperUsdB(iUs) < 14
            CQIPerUs(iUs) = floor(SINRperUsdB(iUs)/1.02+16.62);
        else
            CQIPerUs(iUs) = 30;
        end
    end

    %____________ Computing the total required  rate for each PRIMARY User

    MaxRatePerCode = [461 931 1742 2583]/TTI/10^6; 

    for iUs=1:itti
        if CQIPerUs(iUs)==0; %
            NoCodePerUs(iUs) = 0;
            RatePerUs(iUs) = 0;
        elseif CQIPerUs(iUs)==1; %
            NoCodePerUs(iUs) = 1;
            RatePerUs(iUs) = 377/TTI/10^6;
        elseif CQIPerUs(iUs)==2; %
            NoCodePerUs(iUs) = 1;
            RatePerUs(iUs) = 461/TTI/10^6;
        elseif CQIPerUs(iUs)==3; %
            NoCodePerUs(iUs) = 1;
            RatePerUs(iUs) = 650/TTI/10^6;
        elseif CQIPerUs(iUs)==4; %
            NoCodePerUs(iUs) = 1;
            RatePerUs(iUs) = 792/TTI/10^6;
        elseif CQIPerUs(iUs)==5; %
            NoCodePerUs(iUs) = 1;
            RatePerUs(iUs) = 931/TTI/10^6;
        elseif CQIPerUs(iUs)==6; %
            NoCodePerUs(iUs) = 1;
            RatePerUs(iUs) = 1262/TTI/10^6;
        elseif CQIPerUs(iUs)==7; %
            NoCodePerUs(iUs) = 2;
            RatePerUs(iUs) = 1483/TTI/10^6;
        elseif CQIPerUs(iUs)==8; %
            NoCodePerUs(iUs) = 2;
            RatePerUs(iUs) = 1742/TTI/10^6;
        elseif CQIPerUs(iUs)==9; %
            NoCodePerUs(iUs) = 2;
            RatePerUs(iUs) = 2279/TTI/10^6;
        elseif CQIPerUs(iUs)==10; %
            NoCodePerUs(iUs) = 3;
            RatePerUs(iUs) = 2583/TTI/10^6;
        elseif CQIPerUs(iUs)==11; %
            NoCodePerUs(iUs) = 3;
            RatePerUs(iUs) = 3319/TTI/10^6;
        elseif CQIPerUs(iUs)==12; %
            NoCodePerUs(iUs) = 3;
            RatePerUs(iUs) = 3565/TTI/10^6;
        elseif CQIPerUs(iUs)==13; %
            NoCodePerUs(iUs) = 4;
            RatePerUs(iUs) = 4189/TTI/10^6;
        elseif CQIPerUs(iUs)==14; %
            NoCodePerUs(iUs) = 4;
            RatePerUs(iUs) = 4664/TTI/10^6;
        elseif CQIPerUs(iUs)==15; %
            NoCodePerUs(iUs) = 5;
            RatePerUs(iUs) = 5287/TTI/10^6;
        elseif CQIPerUs(iUs)==16; %
            NoCodePerUs(iUs) = 5;
            RatePerUs(iUs) = 6554/TTI/10^6;
        elseif CQIPerUs(iUs)==17; %
            NoCodePerUs(iUs) = 5;
            RatePerUs(iUs) = 7168/TTI/10^6;
        elseif CQIPerUs(iUs)==18; %
            NoCodePerUs(iUs) = 5;
            RatePerUs(iUs) = 7168/TTI/10^6;
        elseif CQIPerUs(iUs)==19; %
            NoCodePerUs(iUs) = 5;
            RatePerUs(iUs) =7168/TTI/10^6;
        elseif CQIPerUs(iUs)==20; %
            NoCodePerUs(iUs) = 5;
            RatePerUs(iUs) = 7168/TTI/10^6;
        elseif CQIPerUs(iUs)==21; %
            NoCodePerUs(iUs) = 5;
            RatePerUs(iUs) = 7168/TTI/10^6;
        elseif CQIPerUs(iUs)==22; %
            NoCodePerUs(iUs) = 5;
            RatePerUs(iUs) = 7168/TTI/10^6;
        elseif CQIPerUs(iUs)==23; %
            NoCodePerUs(iUs) = 5;
            RatePerUs(iUs) = 7168/TTI/10^6;
        elseif CQIPerUs(iUs)==24; %
            NoCodePerUs(iUs) = 5;
            RatePerUs(iUs) = 7168/TTI/10^6;
        elseif CQIPerUs(iUs)==25; %
            NoCodePerUs(iUs) = 5;
            RatePerUs(iUs) = 7168/TTI/10^6;
        elseif CQIPerUs(iUs)==26; %
            NoCodePerUs(iUs) = 5;
            RatePerUs(iUs) = 7168/TTI/10^6;
        elseif CQIPerUs(iUs)==27; %
            NoCodePerUs(iUs) = 5;
            RatePerUs(iUs) = 7818/TTI/10^6;
        elseif CQIPerUs(iUs)==28; %
            NoCodePerUs(iUs) = 5;
            RatePerUs(iUs) = 7818/TTI/10^6;
        elseif CQIPerUs(iUs)==29; %
            NoCodePerUs(iUs) = 5;
            RatePerUs(iUs) = 7818/TTI/10^6;
        elseif CQIPerUs(iUs)==30; %
            NoCodePerUs(iUs) = 5;
            RatePerUs(iUs) = 7818/TTI/10^6;
        end
    end

end % iSnap


RatePerUs