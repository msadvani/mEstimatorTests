function [Iout] = calculateInfoMat(fNoise,q0Mat,L,N, dblErr) 

    Iout = zeros(size(q0Mat));
    
    for cnt = 1:numel(q0Mat)
        q0 = q0Mat(cnt);
        Iout(cnt) = calculateInfo(fNoise,q0,L,N, dblErr);
    end
end

