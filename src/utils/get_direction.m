function [Vin] = get_direction(Vo,Vin)
% singular vectors can be solved in two opposite directions
% this function is designed to choose the direction
    sizeVo = size(Vo);
    sizeVin = size(Vin);
    assert(sizeVo(1)==sizeVin(1) && sizeVo(2)==sizeVin(2));
    for vc = 1:sizeVo(2)
        if Vin(:,vc)'*Vo(:,vc)<0
            Vin(:,vc) = -Vin(:,vc);
        end
    end
    