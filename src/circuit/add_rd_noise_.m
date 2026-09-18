function [Rout] = add_rd_noise_(Rin,miu_Rd_E,sig_Rd_E,loop)
if ~exist('loop', 'var')
    loop = 0;
end

sizeRin =  size(Rin);

Rd_ER_O = normrnd(miu_Rd_E,sig_Rd_E,sizeRin);


Rout = Rin+Rd_ER_O;

Rout(Rout<0)=0;

end