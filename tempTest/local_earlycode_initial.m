function local_earlycode = local_earlycode_initial(settings,code_table) %%产生初始本地超前码
    global early_code_nco;          %这个是本地早码的相位
    code_word = settings.code_word;
    early_code_temp=[];
    Ncoh = settings.Ncoh;
    for n=1:Ncoh
        early_code_nco = early_code_nco+ code_word ;
        early_code_nco = mod(early_code_nco,2^32*1023);
        index=1+fix(early_code_nco/2^32);
        c=code_table(index);
        early_code_temp=[early_code_temp,c];
    end
    local_earlycode=early_code_temp;
        