function [G0_n,G1_n,G0,G1] = IR_trans_(I_bin,signI,division,Gmin,Gmax,Wr_E,ratio)

    block_n = max(size(division));
    ndimIbin = ndims(I_bin);
    indexingCell = repmat({':'}, 1, ndimIbin-1);
    I = zeros(size(I_bin));
    for i=1 : block_n
        if block_n ==1
            I = I_bin*ratio(i);
        else
            I(indexingCell{:},i) = I_bin(indexingCell{:},i)*ratio(i);
        end
    end

    G0 = zeros(size(I));
    G1 = zeros(size(I));
    G0(signI<=0)=Gmin;
    G1(signI>=0)=Gmin;
    G0(signI>=0) = I(signI>=0)+ G1(signI>=0);
    G1(signI<=0) = G0(signI<=0)+ I(signI<=0);
    Wr_EG0 = 2*Wr_E*(rand(size(I))-0.5);
    Wr_EG1 = 2*Wr_E*(rand(size(I))-0.5);
    G0(G0>Gmax) = Gmax;
    G1(G1>Gmax) = Gmax;
    G0_n = G0+Wr_EG0;
    G1_n = G1+Wr_EG1;
    G0_n(G0_n<Gmin) = Gmin;
    G1_n(G1_n<Gmin) = Gmin;

end