function [V_quant,V_bin_multi_block_dec,signV,V,float_bit] = quantV_(V,quant_bit,cut_bit,block_n,block_b,float_bit_q)
assert(quant_bit>2)
if ~exist('float_bit_q', 'var')
    float_bit = floor(quant_bit/2)+4;
else
    float_bit = float_bit_q;
end
assert(float_bit<=quant_bit)
maxV = max(abs(V),[],"all");
% assert(maxV<=2.^(quant_bit-float_bit))
quant_bit_out = quant_bit;
overflow = false;
if maxV>2^(quant_bit-float_bit) && block_b ==1
    quant_bit = ceil(log2(maxV));
    block_n = quant_bit;
    overflow = true;
%     disp('overflow')
end
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
if overflow
    remain_bit = quant_bit+1;
    V_bin_cut = extractBefore(V_bin, remain_bit);
    mask = zeros(size(signV));
    mask(:,:,1:cut_bit) = 1;
end

split_into_parts = @(str) arrayfun(@(k) str((k-1)*block_b+1:k*block_b), 1:numel(str)/block_b, 'UniformOutput', false);

V_bin_splited = cellfun(split_into_parts, cellstr(V_bin_cut), 'UniformOutput', false);

dims = size(V);


V_bin_splited = string(vertcat(V_bin_splited{:}));
V_bin_multi_block = reshape(V_bin_splited, [dims, block_n]);


V_bin_multi_block_dec = bin2dec(V_bin_multi_block);
if overflow
    V_bin_multi_block_dec = V_bin_multi_block_dec.*mask;
end
% V_quant = V_bin_obj;
% V_quant = bin(V_quant);

V_quant = V_bin_str;

V_quant = bin2dec(V_quant);

end