function [Iout] = calculateInfo(fNoise,q0,L,N, dblErr) 
%Inputs: function describing noise,q0 - variance of gaussian convolution, L - half length of integral interval, N - number of points in integral, error allowed w/ doubling 


    I11 = InfoApprox(fNoise,q0,L,N);
    I12 = InfoApprox(fNoise,q0,L,2*N);
   

    errN= abs(I12 - I11);
    %if(errN>dblErr || I12==0)
    if(errN>dblErr)
        Iout = calculateInfo(fNoise,q0,L,2*N,dblErr);
    else
        %In case you need to resize the interval
        I22 = InfoApprox(fNoise,q0,2*L,2*N);
        errL = abs(I22-I11);

        if(errL>dblErr)
            Iout = calculateInfo(fNoise,q0,2*L,2*N,dblErr);
        else    
            Iout = I22;
        end  
    end
    
    if (Iout==0)
        error('Iout likely inaccurate (Info==0)')
    end
end

