function [out_concated_sum_mean] = block_product_adc_rpt_(G0pn,G1pn,Y_bin,signY,division,block_b_y,ratio_reram,volt,block_size,rpt,adc_bit)
    % [m,n,b] x [n,1,b] = [m,1,b]


    
    %  G0pn, G1pn are the reram values with read and write noise
    %  block_size is the divided block size


        
    
    
    [m0,nx0,bx0,grpt0] = size(G0pn);
    [m1,nx1,bx1,grpt1] = size(G1pn);
    [ny,k,by] = size(Y_bin);
    assert(nx0==nx1 && nx0==ny);
    assert(bx0==bx1 && m0==m1);
    assert(grpt0==rpt)
    assert(rpt>=1)
    if block_size == 0
        block_size = nx0;
    end
    batch_number = ceil(nx0/block_size);
    Y_bin_signed = signY.*Y_bin;

    adc_bit_tmp = zeros(bx0,1);

    for i = 1:bx0%loop x
        adc_bit_tmp(i) = ceil(log2(block_size))+division(i);% here use round or use ceil instead
    end
    if ~exist('adc_bit', 'var') || adc_bit==0
        adc_bit = adc_bit_tmp;
    end
    for i = 1:bx0%loop x
        assert(adc_bit_tmp(i)>=adc_bit(i))
        assert(adc_bit(i)>1)
    end

    out_concated_sum = zeros(m0,k);
    for r = 1:rpt
        out_concated = zeros(m0,k);
    for i = 1:bx0%loop x
        weight_bit = sum(division(i+1:end));
        l_x = 2.^weight_bit;
        Imax_adc = (2.^(adc_bit_tmp(i))-1)*ratio_reram(i)*volt;
        Imin_adc = 0;

        for j = 1:by%loop y
            l_y = (2.^block_b_y).^(by-j);
            for t = 1:batch_number
                if t == batch_number
                    temp_out = (G0pn(:,(1+(t-1)*block_size):(end),i,r) * (Y_bin_signed((1+(t-1)*block_size):(end),:,j)))-(G1pn(:,(1+(t-1)*block_size):(end),i,r) * (Y_bin_signed((1+(t-1)*block_size):(end),:,j)));
                else
                    temp_out = (G0pn(:,(1+(t-1)*block_size):(t*block_size),i,r) * (Y_bin_signed((1+(t-1)*block_size):(t*block_size),:,j)))-(G1pn(:,(1+(t-1)*block_size):(t*block_size),i,r) * (Y_bin_signed((1+(t-1)*block_size):(t*block_size),:,j)));
                end
                [temp_out_adc,I_res,res_bit] = ADC(temp_out,Imin_adc,Imax_adc,adc_bit(i));

%                 max_error = max(temp_out-temp_out_adc/(2.^res_bit)*I_res);

                out_concated = out_concated +  temp_out_adc*l_x * l_y;
            end
        end
    end
    out_concated_sum = out_concated_sum + out_concated;
    end
    out_concated_sum_mean = out_concated_sum/rpt;
end