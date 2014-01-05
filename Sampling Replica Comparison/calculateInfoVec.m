function [Iout] = calculateInfoVec(fNoise,q0Vec,L,N, dblErr) 

    if(min(size(q0Vec))>1)
        error('wrong dims for q0vec');
    end
    
    Iout = zeros(1,length(q0Vec));
    
    for cnt=1:length(q0Vec)
        q0 = q0Vec(cnt);
        Iout(cnt) = calculateInfo(fNoise,q0,L,N, dblErr);
    end
end

