function ang = wrapPi(ang)

    %%normalizza gli angoli nell'intervallo (-pi,pi] ossia quello usato in
    %%robotica.

    ang = mod(ang+pi, 2*pi)-pi;

    mask = ang <= -pi;
    ang(mask) = ang(mask)+2*pi;
end

