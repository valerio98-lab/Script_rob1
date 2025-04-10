Rin=calcRotationMatrix(pi/2, pi/4, -pi/4, 'ZXY');
disp('Matrice di rotazione calcolata:');
disp(Rin);

check_rotation_matrix(Rin)

r = [0, -sqrt(2)/2, sqrt(2)/2];
R_rod = RodriguezMatrix(r, pi/6);
disp('Matrice di rotazione calcolata con il metodo di Rodriguez:');
disp(R_rod);

check_rotation_matrix(R_rod, true, pi/6)

Rif = Rin'*R_rod;
disp('Matrice di rotazione tra init frame e final frame:');
disp(Rif)

%%Calcolo la sequenza di rotazioni dalle quali estrarre gli angoli RPY. 

syms phi chi

Ry = [cos(phi), 0 ,sin(phi);
            0, 1, 0; 
            -sin(phi), 0 cos(phi)]; 

Rx = [0, 0, 1;
            0, cos(chi), sin(chi);
            0 -sin(chi), cos(chi)];

R_YXY = Ry*Rx*Ry;
disp(R_YXY)


