function [G0_n,G1_n,G0,G1,ratio] = memristor_gen_(X_bin,signX,division,Gmin,Gmax,Wr_E)

    block_n = max(size(division));
    ndimXbin = ndims(X_bin);
    ratio = zeros(ones(1,block_n));
    indexingCell = repmat({':'}, 1, ndimXbin-1);
    X_bin_proj = zeros(size(X_bin));
    for i=1 : block_n
        ratio(i) = (Gmax-Gmin)/(2.^(division(i))-1);
        if division(i)==0
            ratio(i) = 1;
        end
        if block_n==1
            X_bin_proj = X_bin*ratio(i);
        else
            X_bin_proj(indexingCell{:},i) = X_bin(indexingCell{:},i)*ratio(i);
        end
    end


    if(~exist('Wr_E','var'))
        Wr_E = 0;
    end
    
    G0 = zeros(size(X_bin));
    G1 = zeros(size(X_bin));


    G0(signX<=0)=Gmin;
    G1(signX>=0)=Gmin;
    G0(signX>=0) = X_bin_proj(signX>=0)+ G1(signX>=0);
    G1(signX<=0) = G0(signX<=0)+ X_bin_proj(signX<=0);

    Wr_EG0 = 2*Wr_E*(rand(size(X_bin))-0.5);
    Wr_EG1 = 2*Wr_E*(rand(size(X_bin))-0.5);


    G0_n = G0+Wr_EG0;
    G1_n = G1+Wr_EG1;

    G0_n(G0_n<Gmin) = Gmin;
    G1_n(G1_n<Gmin) = Gmin;

end