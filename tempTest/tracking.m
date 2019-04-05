function trackResult = tracking(fid, acqResult, settings,data)
%% test
loop_para = loop_canshu_calculate(settings);%计算环路滤波器参数

code_table = cacode(settings.PRN);          %调用函数，产生伪随机码

fll_nco_adder = 0;                   %fll  NCO加法器，应该在滤波器中使用
carrier_nco_sum = 0;                 %给这个值再乘2*pi就是相位了
pll_nco_adder = 0;                   %pll  NCO加法器，应该在滤波器中使用
loop_count = 0;                      %循环次数
code_nco_sum = 0;                    %大概是一个和相位很像的东西
code_nco_adder = 0;                  %ddll  NCO加法器，应该在滤波器中使用
n_IQ = 2;                            
n = 3;
output_fll(2:3) = 0;                 %这里面大概是fll鉴频器的输出
output_filter_fll(1:3) = 0;          %fll 环路滤波器的输出
output_filter_pll(1:3) = 0;          %pll 环路滤波器的输出
output_pll(2:3) = 0;                 %pll 鉴频器的输出
output_filter_ddll(1:3) = 0;         %ddll 环路滤波器的输出
pll_after_filter = 0;

Tcoh = settings.Tcoh;          %积分清零时间
global early_code_nco;      %是接收端的，因为早码的话，即时码和晚码都可以从其中而来
early_code_nco = ((1 - (acqResult.codePhase-3)/settings.samplesPerCode)...
    *settings.codeLength) * 2^settings.nco_Length;
early_code_nco = mod(early_code_nco,2^settings.nco_Length*1023);
local_early_code_last = local_earlycode_initial(settings,code_table); %产生本地超前码，接收端使用，因为早码的话，即时码和晚码都可以从其中而来

trackResult.carrier_nco_phase = zeros(1,settings.msToProcess);       %每次积分清零前载波的nco相位
trackResult.code_nco_phase = zeros(1,settings.msToProcess);          %每次积分清零前B1码的nco相位
trackResult.carrierFreq = zeros(1,settings.msToProcess);
blksize = settings.samplesPerCode;

for loop_num = 1 : settings.msToProcess
    carrier_nco_phase = mod(carrier_nco_sum/2^settings.nco_Length,1) * 2 * pi;  
    if carrier_nco_phase > pi
        trackResult.carrier_nco_phase(loop_num) = ...
            2*pi - carrier_nco_phase;
    else
        trackResult.carrier_nco_phase(loop_num) = carrier_nco_phase;
    end
    trackResult.code_nco_phase(loop_num) = ...
        ((((early_code_nco/2^settings.nco_Length)/settings.codeLength)*settings.samplesPerCode...
         -3) /settings.samplesPerCode)*settings.codeLength;       
    trackResult.carrierFreq(loop_num) = ...
         (settings.middle_freq_nco1 + fll_nco_adder + pll_nco_adder)/settings.transfer_coef;
    trackResult.flag(loop_num) = settings.PLL_flag;         %标识该次循环有没有进行PLL锁定
    %fd_plot(loop_num) = settings.dup_freq;      %记录该次的多普勒频移的真实值
    


    % Read in the appropriate number of samples to process this
    % interation 
    receive_signal = data(fid:fid + blksize - 1);
    fid = fid + blksize ;

          
    %code_nco_sum = code_nco_adder + settings.code_word + fll_nco_adder*(1/763);  %本地再生码环NCO     
  
    
    
    
    
    
    
    
    
    %[signal_modulate_code,settings.signal_phase] = signalcode(settings,code_table);  %产生采样后的CB1码，初始相位为settings.modulate_code_bias_phsae



    
%     receive_signal =  original_signal;  %扩频后的信号
    %产生本地再生载波
    for demond_num = 1:settings.Ncoh 
        local_cos(demond_num) = cos(2*pi*carrier_nco_sum/2^settings.nco_Length);
        local_sin(demond_num) = -sin(2*pi*carrier_nco_sum/2^settings.nco_Length);
        carrier_nco_sum = carrier_nco_sum + settings.middle_freq_nco1 + fll_nco_adder + pll_nco_adder ;%本地再生载波NCO
        
        %carrier_nco_sum大概是一个和相位有关系的东西
    end
    
    code_nco_sum = code_nco_adder + settings.code_word ...              %本地再生码环NCO,和上面的carrier_nco_sum作用一样,2048为基准
        + fll_nco_adder*settings.cofe_FLL_auxi_DDLL;                
  
      %产生本地超前，即时，滞后码
    [local_early_code,local_prompt_code,local_late_code,settings.local_phase]=localcode_generate(local_early_code_last,code_nco_sum,code_table,settings);
    local_early_code_last = local_early_code;
    %载波解调    
    I_demon_carrier = local_cos.*receive_signal;
    Q_demon_carrier = local_sin.*receive_signal;
%     save_I_demon_carrier = [save_I_demon_carrier I_demon_carrier];
%     save_Q_demon_carrier = [save_Q_demon_carrier Q_demon_carrier];
    
    %信号解扩并积分清除
    I_E_final = sum(I_demon_carrier.*local_early_code);
    Q_E_final = sum(Q_demon_carrier.*local_early_code);
    I_P_final(n_IQ) = sum(I_demon_carrier.*local_prompt_code);
    Q_P_final(n_IQ) = sum(Q_demon_carrier.*local_prompt_code);
    I_L_final = sum(I_demon_carrier.*local_late_code);
    Q_L_final = sum(Q_demon_carrier.*local_late_code);
    
    
%     I_P_final(n_IQ) = sum(I_demon_carrier);
%     Q_P_final(n_IQ) = sum(Q_demon_carrier);
    
    
    if  1 == loop_num
        I_P_final(n_IQ - 1) = I_P_final(n_IQ);
        Q_P_final(n_IQ - 1) = Q_P_final(n_IQ);
    else
% %         四象限反正切鉴频器
        dot_fll = I_P_final(n_IQ - 1) * I_P_final(n_IQ) + Q_P_final(n_IQ - 1) * Q_P_final(n_IQ);
        cross_fll = I_P_final(n_IQ - 1) * Q_P_final(n_IQ) - I_P_final(n_IQ) * Q_P_final(n_IQ - 1);
        output_fll(n) = atan2(cross_fll,dot_fll)/(Tcoh*2*pi); 
        result_discriminator_Fll(loop_num) = output_fll(n);
        
        output_filter_fll(n) = (loop_para.cofeone_FLL * output_fll(n)) + (loop_para.cofetwo_FLL * output_fll(n - 1)) + (2 * output_filter_fll(n - 1)) - output_filter_fll(n - 2);
        fll_after_filter(loop_num) = output_filter_fll(n);
        
        fll_nco_adder = output_filter_fll(n) * settings.transfer_coef ;  %频率字转换      
        output_fll(n - 1)=output_fll(n);
        output_filter_fll(n - 2)=output_filter_fll(n - 1);
        output_filter_fll(n - 1)=output_filter_fll(n);
        
         if settings.PLL_flag == 1
            %锁相环鉴相器
            output_pll(n) = atan2(Q_P_final(n_IQ),I_P_final(n_IQ)); 
            output_filter_pll(n) = loop_para.cofeone_PLL*output_pll(n) + loop_para.cofetwo_PLL*output_pll(n-1)+loop_para.cofethree_PLL*output_pll(n-2)+2*output_filter_pll(n-1)-output_filter_pll(n-2);
            result_discriminator_Pll(loop_num) = output_pll(n);
            pll_after_filter(loop_num) = output_filter_pll(n);
            pll_nco_adder = (output_filter_pll(n)/(2*pi)) * settings.transfer_coef;  %频率字转换
            
%             output_pll(1:2) = output_pll(2:3);
%             output_filter_pll(1:2) = output_filter_pll(2:3);
            output_pll(n-2) = output_pll(n-1);
            output_pll(n-1) = output_pll(n);
            output_filter_pll(n-2) = output_filter_pll(n-1);
            output_filter_pll(n-1) = output_filter_pll(n);
         end
        
        I_P_final(n_IQ - 1) = I_P_final(n_IQ);
        Q_P_final(n_IQ - 1) = Q_P_final(n_IQ);
       if 0 == settings.PLL_flag  && abs(output_fll(n))<10  %锁频环工作状态下，信号与本地频差小于10时
            loop_count = loop_count + 1;
            if  loop_count>200            
                   settings.PLL_flag = 1;
            end
       elseif  1 == settings.PLL_flag && abs(output_fll(n))>30      %在锁相环工作状态下，锁频环所鉴出的信号与本地频差大于30时
            loop_count = loop_count-1;
            if  0 == loop_count
                settings.PLL_flag = 0;
            end
       end
    end
 %码环鉴别器
    output_ddll(n) = ((I_E_final - I_L_final)*I_P_final(n_IQ) + (Q_E_final - Q_L_final)*Q_P_final(n_IQ) )/((I_P_final(n_IQ)^2 + Q_P_final(n_IQ)^2)*2);  % DDLL_discri_1
       
    result_ddll(loop_num) = output_ddll(n);
    %码环滤波器（二阶）
    output_filter_ddll(n) = output_filter_ddll(n -1) + (loop_para.cofeone_DDLL*output_ddll(n)) + loop_para.cofetwo_DDLL*output_ddll(n - 1);
    result_DDLL_filter(loop_num) = output_filter_ddll(n);
    % 转换成频率控制字
    code_nco_adder = output_filter_ddll(n) * settings.transfer_coef ; %频率字转换
%     Code_NCO=0;
%     C(loop_num)=code_nco_adder;
    %替换
    output_ddll(n - 1)=output_ddll(n);
    output_filter_ddll(n - 1) = output_filter_ddll(n);
    code_phase_discrim(loop_num) = settings.signal_phase - settings.local_phase ;
end

