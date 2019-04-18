clc
clear
% Initialize constants, settings =========================================
settings = initSettings();
distenses = input('Please enter distences:');
delay_times = distenses/settings.c;
%delay_points = round((delay_times/settings.CA_Period) * settings.samplesPerCode);
% ����α�����,��cacode.m
w_code=cacode(settings.PRN,settings);


% ��Ƶ��Ӧ�õ����ɢ��������
for delay_point_index = 1:length(delay_times)
    %��CA����в���
    samplecacodes = makeCaTable(delay_times(delay_point_index),...
        settings.PRN,settings.codeLength,settings.codeFreqBasis ,settings.samplingFreq,settings);
    spread_code= zeros(0,0);            
    little_spread_code = [ samplecacodes...
        samplecacodes samplecacodes samplecacodes samplecacodes];
    for i = 1:(settings.msToProcess * settings.codeSplitSpace/5)+1
        spread_code = [spread_code little_spread_code];
    end
    %figure(3);
    %plot(spread_code(1:500));%���ע��ֻ��ȡ��5000������ʵ������38192*2000��(2000ms������)
    %title('��Ƶ�������')

    %����
    t = (0:(length(spread_code) - 1))/settings.samplingFreq + delay_times(delay_point_index);
    sendeddataL1=spread_code.*cos(2*pi*settings.IF1.*t);     %L1,����α�� 
    sendeddataL2=cos(2*pi*settings.IF2.*t);                  %L2,������α��
    sendeddata = sendeddataL1 + sendeddataL2;
    % ������
    data= awgn(sendeddata, -10); 

    acqResult = acquisition(data,settings);

    trackResult1 = tracking(1,acqResult,settings,data);
    trackResult2 = tracking2(1,acqResult,settings,data);
    finalDistances = calculatePseudoranges(...
                trackResult1, ...
                trackResult2,...
               settings);
    fprintf("��ʵ���� %f , α���þ��� %f�����Ϊ %f�� \n",...
        distenses(delay_point_index),finalDistances.pseudorange1 , ...
        finalDistances.pseudorange1 - distenses(delay_point_index));
    fprintf("��ʵ���� %f , ˫Ƶα���þ��� %f�����Ϊ %f�� \n",...
        distenses(delay_point_index),finalDistances.pseudorange2 , ...
        finalDistances.pseudorange2 - distenses(delay_point_index));
end