function [output] = transpose_block(InTensor)

% sizeIn = size(InTensor);

order = [2,1,3:ndims(InTensor)];
output = permute(InTensor,order);
end

