clc
clear
% Initialize constants, settings =========================================
settings = initSettings();


% 产生伪随机码,看cacode.m
w_code=cacode(settings.PRN);
%对CA码进行采样
samplecacodes = makeCaTable(settings.PRN,settings.codeLength,settings.codeFreqBasis ,settings.samplingFreq);
% 扩频，应该点乘离散的数据码
spread_code= zeros(0,0);
little_spread_code = [ samplecacodes samplecacodes samplecacodes samplecacodes samplecacodes];
for i = 1:201
    spread_code = [spread_code little_spread_code];
end
%figure(3);
%plot(spread_code(1:500));%这块注意只是取了5000个数据实际上有38192*2000个(2000ms的数据)
%title('扩频后的数据')

%调制
t = (0:(length(spread_code) - 1))/settings.samplingFreq;
sendeddataL1=spread_code.*cos(2*pi*settings.IF1.*t);     %L1,搭载伪码
sendeddataL2=cos(2*pi*settings.IF2.*t);                  %L2,不搭载伪码
sendeddata = sendeddataL1 + sendeddataL2;
% 加噪声
data= awgn(sendeddata, -10); %加-20db分贝的白噪声

%接收信号并把L1和L2分开
[receivedL1,receivedL2] = separateSignal(data,settings.samplingFreq);
%捕获L1并获取伪码起始点和0.5
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


