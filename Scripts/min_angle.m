function diff=min_angle(th_d, th)
    nd = [cos(th_d), sin(th_d), 0];
    n = [cos(th), sin(th), 0];
    diff_abs = acos(n*nd');
    diff_sign = cross(n,nd);
    if diff_sign(3)>0
        diff=diff_abs;
    else
        diff=-diff_abs;
    end
end

