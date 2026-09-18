function [Vin] = get_direction_block(Vec,Vin)
% singular vectors can be solved in two opposite directions
% this function is designed to choose the direction
    sizeVec = size(Vec);
    sizeVin = size(Vin);
    assert(sizeVec(1)==sizeVin(1));
    for vc = 1:sizeVin(2)
        if Vin(:,vc)'*Vec<0
            Vin(:,vc) = -Vin(:,vc);
        end
    end
    