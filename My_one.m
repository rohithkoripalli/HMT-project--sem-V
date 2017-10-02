L = 20E-3;             % total length of domain, m
d = 50E-3;            % total thickness of domain, m
rho = 8000;             % density of material, kg/m3
cp = 316;               % heat capacity, J/kg-K
k = 15;                 % thermal conductivity, W/m-K

T_inf_o = 25;             % free stream temperature, degC
T_init = T_inf_o;         % initial temperature, degC
T_inf_I = 200;
ho = 2500;                % convective heat transfer coefficient, W/m2-K
hi = 500;

nx = 11;               % number of points in x-direction
ny = 26;                 % number of points in y-direction

t_final = 750;         % total simulation time, s
nTimes = 1001;            % total number of time steps

nTemps = nx * ny;       % number of temperatures at each time step
dt = t_final/(nTimes-1);% time step, s

dx = L/(nx-1);          % discretization length in x-direction, m
dy = d/(ny-1);          % discretization length in y-direction, m

x = 0:dx:L;             % vector of x locations (for plotting), m
y = 0:dy:d;             % vector of y-locations (for plotting), m
t = 0:dt:t_final;       % vector of times (for plotting), m

T = zeros(nTemps, nTimes);
A = zeros(nTemps, nTemps);
B = A;
C = zeros(nTemps, 1);

T(:, 1) = T_init;

for z=0:1:9
    for q=0:1:4
        T(117+q+z*nx,:)= T_inf_I;
    end
end


for iPoint = 1:nTemps
    if (iPoint==10*nx-5)
        A(iPoint, iPoint-1) = k * dy/ (dx);
        A(iPoint, iPoint+1) = k * dy/ (2*dx);
        A(iPoint, iPoint - nx) = k * dx/ (dy);
        A(iPoint, iPoint+ nx) = k * dx/ (2*dy);
        A(iPoint, iPoint) = - k * dy/ (dx) - k * dy/ (2*dx) - k * dx/ (dy) - ...
                    k * dx/ (2*dy) - 3*rho * cp * dx * dy/(4*dt);
        B(iPoint, iPoint) = - 3*rho * cp * dx * dy/(4*dt); 
        C(iPoint) = -hi * T_inf_I*(dx/2 + dy/2);

    elseif ((10*nx-5<iPoint)&&(iPoint<10*nx))
        A(iPoint, iPoint-nx) = k * dx/ (dy);
        A(iPoint, iPoint+1) = k * dy/ (2*dx);
        A(iPoint, iPoint-1) = k * dy/ (2*dx);
        A(iPoint, iPoint) = - k * dx/ (dy) - k * dy/ (2*dx) - k * dy/ (2*dx)- ...
                    rho * cp * (dx/2) * dy/dt;
        B(iPoint, iPoint) = - rho * cp * (dx/2) * dy/dt;
        C(iPoint) = -hi * T_inf_I*dx;
    elseif (iPoint==10*nx)
        A(iPoint,iPoint-nx) = k * dx/ (2*dy);
        A(iPoint, iPoint-1) = k * dy/ (2*dx);
        A(iPoint, iPoint) = - k * dx/ (2*dy) - k * dy/ (2*dx) - ...
                rho * cp * (dx/4) * dy/dt;
        B(iPoint, iPoint) = - rho * cp * (dx/4) * dy/dt;
        C(iPoint) = -hi * T_inf_I*dx/2;
    elseif (iPoint==21*nx-5)
        A(iPoint, iPoint-1) = k * dy/ (dx);
        A(iPoint, iPoint+1) = k * dy/ (2*dx);
        A(iPoint, iPoint - nx) = k * dx/ (2*dy);
        A(iPoint, iPoint+ nx) = k * dx/ dy;
        A(iPoint, iPoint) = - k * dy/ (dx) - k * dy/ (2*dx) - k * dx/ (2*dy) - ...
                k * dx/ dy - 3*rho * cp * dx * dy/(4*dt);
        B(iPoint, iPoint) = - 3*rho * cp * dx * dy/(4*dt); 
        C(iPoint) = -hi * T_inf_I*(dx/2 + dy/2);
    elseif ((21*nx-5<iPoint)&&(iPoint<21*nx))
        A(iPoint, iPoint-1) = k * dy/ (2*dx);
        A(iPoint, iPoint+nx) = k * dx/ (dy);
        A(iPoint, iPoint+1) = k * dy/ (2*dx);
        A(iPoint, iPoint) = - k * dy/ (2*dx) - k * dy/ (2*dx) - k * dx/(dy)-...
                    rho * cp * (dx/2) * dy/dt;
        B(iPoint, iPoint) = - rho * cp * (dx/2) * dy/dt;
        C(iPoint) = -hi * T_inf_I*dx;
    elseif (iPoint==21*nx)    
        A(iPoint, iPoint+nx) = k * dx/ (2*dy);
        A(iPoint, iPoint-1) = k * dy/ (2*dx);
        A(iPoint, iPoint) = - k * dx/ (2*dy) - k * dy/ (2*dx) - ...
                rho * cp * (dx/4) * dy/dt;
        B(iPoint, iPoint) =  rho * cp * (dx/4) * dy/dt;
        C(iPoint) = -hi * T_inf_I*dx/2;
    elseif (iPoint == (11*nx-5)|| iPoint == (12*nx-5)|| iPoint == (13*nx-5)||...
         iPoint == (14*nx-5)|| iPoint == (15*nx-5)|| iPoint == (16*nx-5)|| iPoint == (17*nx-5)||...
         iPoint == (18*nx-5)|| iPoint == (19*nx-5)|| iPoint == (20*nx-5))
        A(iPoint, iPoint-nx) = k * dx/ (2*dy);
        A(iPoint, iPoint-1) = k * dy/ (dx);
        A(iPoint, iPoint+nx) = k * dx/ (2*dy);
        A(iPoint, iPoint) = - k * dx/ (2*dy) - k * dy/ (dx) - k * dx/(2*dy)- ...
            rho * cp * (dx/2) * dy/dt;
        B(iPoint, iPoint) = - rho * cp * (dx/2) * dy/dt;
        C(iPoint) = -hi * T_inf_I*dy;
    elseif ((116<iPoint)&&(iPoint<=121) || (127<iPoint)&&(iPoint<=135) || (138<iPoint)&&(iPoint<=143) ||...
            (149<iPoint)&&(iPoint<=154) || (160<iPoint)&&(iPoint<=165)|| ...
            (171<iPoint)&&(iPoint<=176) || (182<iPoint)&&(iPoint<=187) || (193<iPoint)&&(iPoint<=198) ||...
            (204<iPoint)&&(iPoint<=209) || (215<iPoint)&&(iPoint<=220))
        A(iPoint, iPoint) = 1;
        B(iPoint, iPoint) = 1;
    elseif iPoint == 1
        A(iPoint, iPoint+1) = k * dy/ (2 * dx);
        A(iPoint, iPoint+nx) = k * dx/ (2 * dy);
        A(iPoint,iPoint) = - k * dy/ (2 * dx) - k * dx/ (2 * dy) - ho * dx/2 - ...
            rho * cp * (dx/2) * (dy/2)/ dt;
        B(iPoint, iPoint) = - rho * cp * (dx/2) * (dy/2)/ dt;
        C(iPoint) = - ho * dx/2 * T_inf_o;
    elseif (1<iPoint)&&(iPoint<nx)
        A(iPoint, iPoint-1) = k * dy/ (2*dx);
        A(iPoint, iPoint+1) = k * dy/ (2*dx);
        A(iPoint, iPoint+nx) = k * dx/ dy;
        A(iPoint, iPoint) = - k * dy/ (2*dx) - k * dy/ (2*dx) - k * dx/ dy - ...
            ho * dx - rho * cp * dx * (dy/2)/ dt;
        B(iPoint, iPoint) = - rho * cp * dx * (dy/2)/ dt;
        C(iPoint) = - ho * dx * T_inf_o;
    elseif iPoint == nx
        A(iPoint, iPoint-1) = k * dy/ (2*dx);
        A(iPoint, iPoint+nx) = k * dx/ (2*dy);
        A(iPoint, iPoint) = - k * dy/ (2*dx) - k * dx/ (2*dy) - ho * dx/2- ...
              rho * cp * (dx/2) * (dy/2)/ dt;
        B(iPoint, iPoint) = - rho * cp * (dx/2) * (dy/2)/ dt;
        C(iPoint) = - ho * T_inf_o * (dx/2);
    elseif ( iPoint == (nx+1) || iPoint == (2*nx+1) || iPoint == (3*nx+1)|| iPoint == (4*nx+1)||...
         iPoint == (5*nx+1)|| iPoint == (6*nx+1)|| iPoint == (7*nx+1)|| iPoint == (8*nx+1)|| iPoint == (9*nx+1)||...
         iPoint == (10*nx+1)|| iPoint == (11*nx+1)|| iPoint == (12*nx+1)|| iPoint == (13*nx+1)||...
         iPoint == (14*nx+1)|| iPoint == (15*nx+1)|| iPoint == (16*nx+1)|| iPoint == (17*nx+1)||...
         iPoint == (18*nx+1)|| iPoint == (19*nx+1)|| iPoint == (20*nx+1)|| iPoint == (21*nx+1)||...
         iPoint == (22*nx+1)|| iPoint == (23*nx+1)|| iPoint == (24*nx+1))
      
        A(iPoint, iPoint-nx) = k * dx/ (2 * dy);
        A(iPoint, iPoint+nx) = k * dx/ (2 * dy);
        A(iPoint, iPoint+1) = k * dy/ dx;
        A(iPoint, iPoint) = - k * dx/ (2 * dy) - k * dx/ (2 * dy) - k * dy/ dx-...
             rho * cp * (dx/2) * dy/ dt;
        B(iPoint, iPoint) = - rho * cp * (dx/2) * dy/ dt;    
    
    elseif ( iPoint == (2*nx)||iPoint == (3*nx)|| iPoint == (4*nx)||...
         iPoint == (5*nx)|| iPoint == (6*nx)|| iPoint == (7*nx)|| iPoint == (8*nx)|| iPoint == (9*nx)||...
         iPoint == (22*nx)|| iPoint == (23*nx)|| iPoint == (24*nx)|| iPoint == (25*nx) )
        A(iPoint, iPoint-nx) = k * dx/ (2*dy);
        A(iPoint, iPoint+nx) = k * dx/ (2*dy);
        A(iPoint, iPoint-1) = k * dy/ dx;
        A(iPoint, iPoint) = - k * dx/ (2*dy) - k * dx/ (2*dy) - ...
            k * dy/ dx - rho * cp * (dx/2) * dy/dt;
        B(iPoint, iPoint) = - rho * cp * (dx/2) * dy/dt; 
    
    elseif iPoint == (25*nx+1)
        A(iPoint, iPoint-nx) = k * dx/ (2 * dy);
        A(iPoint, iPoint+1) = k * dy/ (2*dx);
        A(iPoint, iPoint) = - k * dx/ (2 * dy) - k * dy/ (2*dx)-...
             rho * cp * (dx/2) * (dy/2)/ dt;
        B(iPoint, iPoint) = - rho * cp * (dx/2) * (dy/2)/ dt;
    
    elseif ((nx+1 < iPoint) &&  (iPoint < 2*nx)||(2*nx+1< iPoint) &&  (iPoint < 3*nx)||(3*nx+1 < iPoint) &&  (iPoint < 4*nx)||...
            (4*nx+1< iPoint) &&  (iPoint < 5*nx)||(5*nx+1 < iPoint) &&  (iPoint < 6*nx)||(6*nx+1 < iPoint) &&  (iPoint < 7*nx)||...
            (7*nx+1< iPoint) &&  (iPoint < 8*nx)||(8*nx+1 < iPoint) &&  (iPoint < 9*nx)||...
            (9*nx+1 < iPoint) &&  (iPoint < 10*nx-5) ||(10*nx+1 < iPoint) &&  (iPoint < 11*nx-5)||...
            (11*nx+1 < iPoint) &&  (iPoint < 12*nx-5)||(12*nx+1 < iPoint) &&  (iPoint < 13*nx-5)||...
            (13*nx+1 < iPoint) &&  (iPoint < 14*nx-5)||(14*nx+1 < iPoint) &&  (iPoint < 15*nx-5)|| ...
            (15*nx+1 < iPoint) &&  (iPoint < 16*nx-5)||(16*nx+1 < iPoint) &&  (iPoint < 17*nx-5)||...
            (17*nx+1 < iPoint) &&  (iPoint < 18*nx-5)||(18*nx+1 < iPoint) &&  (iPoint < 19*nx-5)||...
            (19*nx+1 < iPoint) &&  (iPoint < 20*nx-5)||(20*nx+1 < iPoint) &&  (iPoint < 21*nx-5)||...
            (21*nx+1< iPoint) &&  (iPoint < 22*nx)||(22*nx+1 < iPoint) &&  (iPoint < 23*nx)||...
            (23*nx+1< iPoint) &&  (iPoint < 24*nx)||(24*nx+1 < iPoint) &&  (iPoint < 25*nx))    
 
   
        A(iPoint, iPoint-1) = k * dy/ (dx);
        A(iPoint, iPoint+1) = k * dy/ (dx);
        A(iPoint, iPoint - nx) = k * dx/ dy;
        A(iPoint, iPoint + nx) = k * dx/ dy;
        A(iPoint, iPoint) = - k * dy/ (dx) - k * dy/ (dx) - k * dx/ dy - ...
            k * dx/ dy - rho * cp * dx * dy/(dt);
        B(iPoint, iPoint) = - rho * cp * dx * dy/(dt);
    elseif ((25*nx+1<iPoint)&&(iPoint<26*nx))
        A(iPoint, iPoint-1) = k * dy/ (2*dx);
        A(iPoint, iPoint+1) = k * dy/ (2*dx);
        A(iPoint, iPoint-nx) = k * dx/ dy;
        A(iPoint, iPoint) = - k * dy/ (2*dx) - k * dy/ (2*dx) - k * dx/ dy - ...
              rho * cp * dx * (dy/2)/ dt;
        B(iPoint, iPoint) = - rho * cp * dx * (dy/2)/ dt;
       
    elseif iPoint==26*nx
        A(iPoint, iPoint-1) = k * dy/ (2*dx);
        A(iPoint, iPoint-nx) = k * dx/ (2*dy);
        A(iPoint, iPoint) = - k * dy/ (2*dx) - k * dx/ (2*dy) - ...
             rho * cp * (dx/2) * (dy/2)/ dt;
        B(iPoint, iPoint) = - rho * cp * (dx/2) * (dy/2)/ dt;
       
    end
end

            
for iTime = 1:1:(nTimes-1)
    T(:,iTime+1)= A\(B*T(:,iTime)+ C);
end

   
surf(T)
 
