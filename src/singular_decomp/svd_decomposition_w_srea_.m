function [V,S,U] = svd_decomposition_w_srea_(dataX,division,Gmin,Gmax,Wr_E,miu_Rd_E,sig_Rd_E,voltage,block_size,vrpt_time,iter)
    if ~exist('iter', 'var')
        iter = 8;
    end
    if vrpt_time<=0
        vrpt_time = 1;
    end


    [V,S,U] = power_iter_quant(dataX,division,Gmin,Gmax,Wr_E,miu_Rd_E,sig_Rd_E,voltage,block_size,vrpt_time,iter);


    end



function [Vout,Sout,Uout] = power_iter_quant(dataX,division,Gmin,Gmax,Wr_E,miu_Rd_E,sig_Rd_E,voltage,block_size,vrpt_time,iter)
    if ~exist('iter', 'var')
        iter = 20;
    end
    if block_size <= 0
        block_size = 0;
    end

     
    bit = sum(division);
    sizeX = size(dataX);
    Vout = zeros(sizeX(2),sizeX(2));
    Uout = zeros(sizeX(1),sizeX(2));
    Dout = zeros(1,sizeX(2));
    reram_blocks = max(size(division));

    
    [X_quant,X_bin_multi_block_dec,signX,X,ratio_quantX] = quantX_(dataX,bit,division);

    [G0_n,G1_n,G0,G1,ratio_project] = memristor_gen_(X_bin_multi_block_dec,signX,division,Gmin,Gmax,Wr_E);

    
    lambda0 = 0;
    lambdaq0 = 0;
    lambda = lambda0;
    lambdaq = lambdaq0;
    block_n = ndims(division);
    ratio_projectV = zeros(ones(1,block_n));

    for i = 1:sizeX(2)

        x0_r = rand(sizeX(2),1);
        norm_x0 = norm(x0_r);
        x0 = x0_r/norm_x0;
        u0_r = rand(sizeX(1),1);
        norm_u0 = norm(u0_r);
        u0 = u0_r/norm_u0;


        if(i>1)
            [Vec_quant,Vec_bin_multi_block_dec,signVec,Vec,float_bit_Vec] = quantXV_(V(:,i-1)',bit,bit,division,bit);
            [Vq_quant,Vq_bin_multi_block_dec,signVq,Vecq,float_bit_Vq] = quantXV_(Vq(:,i-1)',bit,bit,division,bit);
            if (i-1)==1
                [GV0_n,GV1_n,GV0,GV1,ratio_projectV] = memristor_gen_(Vec_bin_multi_block_dec,signVec,division,Gmin,Gmax,Wr_E);
            else
                assert(ratio_projectV(1)~=0)
                [GV0_n,GV1_n,GV0,GV1] = IR_trans_(Vec_bin_multi_block_dec,signVec,division,Gmin,Gmax,Wr_E,ratio_projectV);
            end
            if reram_blocks ==1
                G_vec_n0((i-1),:,1) = GV0_n;
                G_vec_n1((i-1),:,1) = GV1_n;

            else
                G_vec_n0((i-1),:,:) = GV0_n;
                G_vec_n1((i-1),:,:) = GV1_n;

            end

            for vrpt = 1:vrpt_time
                if (i-1)==1
                    [GVq0_n,GVq1_n,GVq0,GVq1,ratio_projectVq] = memristor_gen_(Vq_bin_multi_block_dec,signVq,division,Gmin,Gmax,Wr_E);
                else
                    assert(ratio_projectVq(1)~=0)
                    [GVq0_n,GVq1_n,GVq0,GVq1] = IR_trans_(Vq_bin_multi_block_dec,signVq,division,Gmin,Gmax,Wr_E,ratio_projectVq);
                end

                if reram_blocks == 1
                    G_Vq_n0((i-1),:,1,vrpt) = GVq0_n;
                    G_Vq_n1((i-1),:,1,vrpt) = GVq1_n;

                else
                    G_Vq_n0((i-1),:,:,vrpt) = GVq0_n;
                    G_Vq_n1((i-1),:,:,vrpt) = GVq1_n;

                end
            end
        end
        v = x0;
        v2_quant_norm = x0;
        for k = 1:iter


            [G0_n_rd] = add_rd_noise_(G0_n,miu_Rd_E,sig_Rd_E);
            [G1_n_rd] = add_rd_noise_(G1_n,miu_Rd_E,sig_Rd_E);

            block_n_y = 8;
            block_b_y = 1;

            adc_bit = 0;
            expected_bit_1 = block_n_y + ceil(log2(sizeX(2))) + sum(division,"all");
            minused_bits = expected_bit_1-8-1;
            [v_quant,v_bin_multi_block_dec,signv,v2_quant_norm,float_bit_v] = quantV_(v2_quant_norm,block_n_y,block_n_y,block_n_y,block_b_y,8);
            v_bin = v_bin_multi_block_dec*voltage;


            v1_quant = block_product_adc_(G0_n_rd,G1_n_rd,v_bin,signv,division,block_b_y,ratio_project,voltage,block_size,adc_bit);
            v1_quant_proj = v1_quant/2.^(minused_bits);
            uq = v1_quant/norm(v1_quant);
            [v1q_out,v1q_bin_multi_block_dec,signv1q,v1q,float_bit_v1q] = quantV_(v1_quant_proj,block_n_y,block_n_y,block_n_y,block_b_y,0);
            v1q_bin = v1q_bin_multi_block_dec * voltage;

            G0_n_rd_T = transpose_block(G0_n_rd);
            G1_n_rd_T = transpose_block(G1_n_rd);
            
            v2_quant_tmp = block_product_adc_(G0_n_rd_T,G1_n_rd_T,v1q_bin,signv1q,division,block_b_y,ratio_project,voltage,block_size,adc_bit);


            if(i>1)

                [G_Vq_n0_rd] = add_rd_noise_(G_Vq_n0,miu_Rd_E,sig_Rd_E);
                [G_Vq_n1_rd] = add_rd_noise_(G_Vq_n1,miu_Rd_E,sig_Rd_E);

                vq1 = block_product_adc_rpt_(G_Vq_n0_rd,G_Vq_n1_rd,v_bin,signv,division,block_b_y,ratio_projectVq,voltage,block_size,vrpt_time,adc_bit);
                vq1_proj = vq1/2.^(minused_bits);
                vq1_proj_d = -D_quant'.* vq1_proj;

                extra_bit = 6;
                block_n_y_v = block_n_y + extra_bit;
                [vq1_out,vq1_bin_multi_block_dec,signvq1,vq1_proj_d_q,float_bit_vq1] = quantV_(vq1_proj_d,block_n_y_v,block_n_y_v,block_n_y_v,block_b_y,0); 
                vq1_bin = vq1_bin_multi_block_dec * voltage;

                [G_Vq_n0_rd] = add_rd_noise_(G_Vq_n0,miu_Rd_E,sig_Rd_E);
                [G_Vq_n1_rd] = add_rd_noise_(G_Vq_n1,miu_Rd_E,sig_Rd_E);
                G_Vq_n0_rd_T = transpose_block(G_Vq_n0_rd);
                G_Vq_n1_rd_T = transpose_block(G_Vq_n1_rd);

                vq2 = block_product_adc_rpt_(G_Vq_n0_rd_T,G_Vq_n1_rd_T,vq1_bin,signvq1,division,block_b_y,ratio_projectVq,voltage,block_size,vrpt_time,adc_bit);
                vq2_proj = vq2/((2.^float_bit_Vq).^2)*(ratio_quantX.^2);
                v2_quant = v2_quant_tmp + vq2_proj;


            else
                v2_quant = v2_quant_tmp;
            end
            v2_quant_norm = v2_quant/norm(v2_quant);

            lambdaq0 = lambdaq;
            [v1_quant_l,v1_bin_multi_block_dec_l,signv1_l,v1_l,float_bit_v1_l] = quantV_(v2_quant_norm,block_n_y,block_n_y,block_n_y,block_b_y,8);
            v1_bin_l = v1_bin_multi_block_dec_l*voltage;
            v1_quant_l = block_product_adc_(G0_n_rd,G1_n_rd,v1_bin_l,signv1_l,division,block_b_y,ratio_project,voltage,block_size,adc_bit);
            v1_quant_proj_l = v1_quant_l/2.^(minused_bits);
            [v1q_out_l,v1q_bin_multi_block_dec_l,signv1q_l,v1q_l,float_bit_v1q_l] = quantV_(v1_quant_proj_l,block_n_y,block_n_y,block_n_y,block_b_y,0);
            v1q_bin_l = v1q_bin_multi_block_dec_l * voltage;


            G0_n_rd_T = transpose_block(G0_n_rd);
            G1_n_rd_T = transpose_block(G1_n_rd);
            v2_quant_tmp_l = block_product_adc_(G0_n_rd_T,G1_n_rd_T,v1q_bin_l,signv1q_l,division,block_b_y,ratio_project,voltage,block_size,adc_bit);



            if(i>1)

                [G_Vq_n0_rd] = add_rd_noise_(G_Vq_n0,miu_Rd_E,sig_Rd_E);
                [G_Vq_n1_rd] = add_rd_noise_(G_Vq_n1,miu_Rd_E,sig_Rd_E);

                vq1_l = block_product_adc_rpt_(G_Vq_n0_rd,G_Vq_n1_rd,v1_bin_l,signv1_l,division,block_b_y,ratio_projectVq,voltage,block_size,vrpt_time,adc_bit);
                vq1_proj_l = vq1_l/2.^(minused_bits);
                vq1_proj_d_l = -D_quant'.* vq1_proj_l;

                extra_bit = 6;
                block_n_y_v = block_n_y + extra_bit;

                [vq1_out_l,vq1_bin_multi_block_dec_l,signvq1_l,vq1_proj_d_q_l,float_bit_vq1_l] = quantV_(vq1_proj_d_l,block_n_y_v,block_n_y_v,block_n_y_v,block_b_y,0);
                vq1_bin_l = vq1_bin_multi_block_dec_l * voltage;


                [G_Vq_n0_rd] = add_rd_noise_(G_Vq_n0,miu_Rd_E,sig_Rd_E);
                [G_Vq_n1_rd] = add_rd_noise_(G_Vq_n1,miu_Rd_E,sig_Rd_E);
                G_Vq_n0_rd_T = transpose_block(G_Vq_n0_rd);
                G_Vq_n1_rd_T = transpose_block(G_Vq_n1_rd);

                vq2_l = block_product_adc_rpt_(G_Vq_n0_rd_T,G_Vq_n1_rd_T,vq1_bin_l,signvq1_l,division,block_b_y,ratio_projectVq,voltage,block_size,vrpt_time,adc_bit);
                
                         
                vq2_proj_l = vq2_l/((2.^float_bit_Vq).^2)*(ratio_quantX.^2);
                v2_quant_l = v2_quant_tmp_l + vq2_proj_l;
                v2_quant_recover_l = v2_quant_l/(ratio_quantX.^2)*2.^(-float_bit_v1_l+minused_bits);
                lambdaq = dot(v2_quant_recover_l,v2_quant_norm);

            else
                v2_quant_tmp_recover_l = v2_quant_tmp_l/(ratio_quantX.^2)*2.^(-float_bit_v1_l+minused_bits);
                lambdaq = dot(v2_quant_tmp_recover_l,v2_quant_norm);
            end

            err = abs(lambdaq-lambdaq0);

            if err<1e-1
                break;
            end
        end

    vector = v;
    vector_quant = v2_quant_norm;
    value_quant = lambdaq;
    V(:,i) = vector;
    Vq(:,i) = vector_quant;
    D_quant(i) = value_quant;
    if value_quant>0
        sq = sqrt(value_quant);
    else
        sq = 0;
    end
    Sq(i) = sq;
    U(:,i) = uq;
    end
    sizeVq = size(Vq);
    sizeUq = size(U);
    sizeSq = size(Sq);

    Vout(1:sizeVq(1),1:sizeVq(2)) = Vq;
    Uout(1:sizeUq(1),1:sizeUq(2)) = U;
    Sout(1:sizeSq(2)) = Sq;
    
end


