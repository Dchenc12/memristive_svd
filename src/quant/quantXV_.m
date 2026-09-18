function [V_quant,V_bin_multi_block_dec,signV,V,float_bit] = quantXV_(V,quant_bit,cut_bit,division,float_bit_q)
assert(quant_bit>2)
if ~exist('float_bit_q', 'var')
    float_bit = floor(quant_bit/2)+4;
else
    float_bit = float_bit_q;
end
assert(float_bit<=quant_bit)
nDimsV = ndims(V);
% block_n = ndims(division);
block_n = max(size(division));
signV = sign(V);
signV = repmat(signV,[ones(1,ndims(V)),block_n]);

V_abs = abs(V);
% T = numerictype('WordLength', quant_bit, 'FractionLength', float_bit, 'Signed', false);
% V_bin_obj = fi(V_abs, T);
% V_bin_str = bin(V_bin_obj);
% if no toolbox Fixed-point Designer###########################################
scale_factor = 2^float_bit;

V_scaled = round(V_abs * scale_factor);
max_val = 2^quant_bit - 1;
V_scaled = min(max(V_scaled, 0), max_val);
V_int = uint32(V_scaled);
V_bin_str = dec2bin(V_int, quant_bit);
% if numel(V_abs) > 1
%     V_bin_str = strings(size(V_abs));
%     for i = 1:numel(V_abs)
%         V_scaled = round(V_abs(i) * scale_factor);
%         V_scaled = min(max(V_scaled, 0), max_val);
%         V_bin_str(i) = dec2bin(uint32(V_scaled), quant_bit);
%     end
% end

%############################################
% split_str = strsplit(binary_str);
split_str_matrix = cellfun(@strsplit, cellstr(V_bin_str), 'UniformOutput', false);

V_bin = string(vertcat(split_str_matrix{:}));
remain_bit = cut_bit + 1;
assert(remain_bit<=quant_bit+1)
V_bin_cut = extractBefore(V_bin, remain_bit);

V_bin_reshape = reshape(V_bin_cut, [size(V)]);


V_bin_multi_block = repmat(V_bin_reshape,[ones(1,nDimsV),block_n]);

indexingCell = repmat({':'}, 1, nDimsV);

for i = 1:block_n
    V_bin_multi_block(indexingCell{:},i) = extractBetween(V_bin_reshape, sum(division(1:i-1))+1, sum(division(1:i)));
end


V_bin_multi_block_dec = bin2dec(V_bin_multi_block);




% V_quant = V_bin_obj;
% V_quant = dec(V_quant);
% V_quant = strsplit(V_quant);
% V_quant = str2double(V_quant); 

V_quant = bin2dec(V_bin_str);


end