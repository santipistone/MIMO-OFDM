    clear
    close all 
    clc
%% Sweep parameters
    SNR_dBV     = 3:1:15;            % vector of SNR values in dB
    SNR_dBVL    = length(SNR_dBV);   % length of SNR vector
    nMonteCarlo = 1e1;               % number of Monte Carlo to find the value of each point in the figure
    ofdmIn.Nt   = 2;                 % number of transmit antennas
    ofdmIn.Nr   = 3;                 % number of recieve antennas
    %%ofdmIn.ifDisplayResults    = 0;  % turn off the display
    % other parameters of ofdm can also be set. see help of MIMO_OFDM_LSE_CHAN_EST
%% Outputs
    MSE_CHAN_SIM = zeros(nMonteCarlo,SNR_dBVL);     % MSE of LSE channel estimation in simulation
    MSE_CHAN_THR = zeros(nMonteCarlo,SNR_dBVL);     % MSE of LSE channel estimation in theory
    MSE_CHAN_BER = zeros(nMonteCarlo,SNR_dBVL);     % BER of the MIMO OFDM with LSE channel estimation
%% Loop
cnt3 = 0;
for cnt1 = 1 : SNR_dBVL
    chanIn.SNR_dB              = SNR_dBV(cnt1); % load the SNR value    
    for cnt2 = 1 : nMonteCarlo
        [ofdm chan] = MIMO_OFDM_LSE_CHAN_EST(ofdmIn,chanIn);
        MSE_CHAN_SIM(cnt2,cnt1) = chan.MSE_Simulation;
        MSE_CHAN_THR(cnt2,cnt1) = chan.MSE_Theory;
        MSE_CHAN_BER(cnt2,cnt1) = ofdm.BER;
        
        % update the loop counter
        cnt3 = cnt3 + 1;
        disp([num2str(round(cnt3/(SNR_dBVL*nMonteCarlo)*1000)/10),' is done...'])        
    end   
end
%% Figures
    figure(1)
    clf
    semilogy(SNR_dBV,mean(MSE_CHAN_THR),'-r.',...
             SNR_dBV,mean(MSE_CHAN_SIM),'bx ')
    xlabel('SNR [dB]')
    ylabel('MSE of LSE channel estimation')
    grid on
    legend('Simulation','Theory')
    title(['Nt = ',num2str(ofdm.Nt),', Nr = ',num2str(ofdm.Nr)])
    figure(2)
    clf
    semilogy(SNR_dBV,mean(MSE_CHAN_BER),'b')
    xlabel('SNR [dB]')
    ylabel('Bit error rate (BER)')
    grid on