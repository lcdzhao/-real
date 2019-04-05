clc
clear
% Initialize constants, settings =========================================
settings = initSettings();


% ����α�����,��cacode.m
w_code=cacode(settings.PRN);
%��CA����в���
samplecacodes = makeCaTable(settings.PRN,settings.codeLength,settings.codeFreqBasis ,settings.samplingFreq);
% ��Ƶ��Ӧ�õ����ɢ��������
spread_code= zeros(0,0);
little_spread_code = [ samplecacodes samplecacodes samplecacodes samplecacodes samplecacodes];
for i = 1:201
    spread_code = [spread_code little_spread_code];
end
%figure(3);
%plot(spread_code(1:500));%���ע��ֻ��ȡ��5000������ʵ������38192*2000��(2000ms������)
%title('��Ƶ�������')

%����
t = (0:(length(spread_code) - 1))/settings.samplingFreq;
sendeddataL1=spread_code.*cos(2*pi*settings.IF1.*t);     %L1,����α��
sendeddataL2=cos(2*pi*settings.IF2.*t);                  %L2,������α��
sendeddata = sendeddataL1 + sendeddataL2;
% ������
data= awgn(sendeddata, -10); %��-20db�ֱ��İ�����

%�����źŲ���L1��L2�ֿ�
[receivedL1,receivedL2] = separateSignal(data,settings.samplingFreq);
%����L1����ȡα����ʼ���0.5
acqResult = acquisition(data,settings);


%channel = preRun(acqResults,settings);
%showChannelStatus(channel,settings);
trackResult = tracking(1,acqResult,settings,data);
%[subFrameStart, activeChnList] = findPreambles(trackResults, settings);
 


finalDistance = calculatePseudoranges(...
            trackResults, ...
            acqResults.codePhase(PRN), ...
           settings);
disp(finalDistance);
%codeError= test(PRN,data)


