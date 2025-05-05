function s = rad2pi(angle_rad, tol)
%Converte un angolo espresso in radianti in stringa k*pi, ossia in pigreco

    if nargin < 2 
        tol = 1e-6;
    end
    
    [N, D] = rat(angle_rad/pi, tol);
    N = double(N);
    disp(N)
    if D==1 
        s = sprintf('%d*pi', N);
    else
        s = sprintf('%d/%d*pi', N, D);
    end
end