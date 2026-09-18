function [X_quant,X_bin_multi_block_dec,signX,X,ratio] = quantX_(X,bit,division)

% % division shows how to divide the block. 
% % division = [3, 5] means the first block has 3 bits, the second block has 5 bits.
% % n-bit for presentation, if include sign bit actually is n+1 bit. (bit=8, represents [-256, 255])
% % this part of quantization should not be out of range.

assert(sum(division)==bit)

assert(bit>2)

upper_bound = 2.^(bit)-1;
lower_bound = -2.^(bit);



block_n = max(size(division));
nDimsX = ndims(X);
assert(block_n>=1)

x_abs_max = max(abs(X),[],'all');

ratio = upper_bound/x_abs_max;
X_proj = X*ratio;
X_proj_abs = abs(X_proj);
signX = sign(X);
signX(signX==0)=1;
signX = repmat(signX,[ones(1,ndims(X_proj)),block_n]);

% T = numerictype('WordLength', bit, 'FractionLength', 0, 'Signed', false);
% X_bin_obj = fi(X_proj_abs, T);
% X_bin_str = bin(X_bin_obj);
% if no toolbox Fixed-point Designer###########################################
max_val = 2^bit - 1;  
X_bin = min(max(round(X_proj_abs), 0), max_val); 
X_bin = uint32(X_bin); 
X_bin_str = dec2bin(X_bin, bit); 
% if numel(X_bin) > 1
%     X_bin_str = strings(size(X_bin));
%     for i = 1:numel(X_bin)
%         X_bin_str(i) = dec2bin(X_bin(i), bit);
%     end
% end
%############################################

split_str_matrix = cellfun(@strsplit, cellstr(X_bin_str), 'UniformOutput', false);

X_bin = string(vertcat(split_str_matrix{:}));

x_bin_reshape = reshape(X_bin, [size(X)]);

X_bin_multi_block = repmat(x_bin_reshape,[ones(1,nDimsX),block_n]);

indexingCell = repmat({':'}, 1, nDimsX);

for i = 1:block_n
    X_bin_multi_block(indexingCell{:},i) = extractBetween(x_bin_reshape, sum(division(1:i-1))+1, sum(division(1:i)));
end

X_bin_multi_block_dec = bin2dec(X_bin_multi_block);
X_quant = floor(X_proj);




X_quant(X_quant<lower_bound)=lower_bound;
X_quant(X_quant>upper_bound)=upper_bound;


end