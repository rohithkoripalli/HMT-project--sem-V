
L = 220E-3;             % total length of domain, m
d = 1.25E-3;            % total thickness of domain, m
w = 22E-3;              % width of plastic film, m
rho = 7850;             % density of material, kg/m3
cp = 435;               % heat capacity, J/kg-K
k = 60;                 % thermal conductivity, W/m-K

T_inf = 25;             % free stream temperature, degC
T_init = T_inf;         % initial temperature, degC
h = 100;                % convective heat transfer coefficient, W/m2-K

delta_t_ON = 10.0;      % time for which laser is switched on, s
q0_laser = 85E+3;       % laser irradiation, W/m2;
t_final = 30.0;         % total simulation time, s

nx = 111;               % number of points in x-direction
ny = 4;                 % number of points in y-direction

nTimes = 61;            % total number of time steps

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

for iPoint = 1:nTemps
    if iPoint == 1
        for iTime = 1:(nTimes-1)
            T(iPoint,iTime+1) = k * dy *(T(iPoint+1,iTime)-T(iPoint,iTime))/ (2 * dx) + ... 
             k * dx*(T(iPoint+nx,iTime)-T(iPoint,iTime))/ (2 * dy) + ...
             h*(T_inf - T(iPoint,iTime))*dx/2 +  rho * cp * (dx/2) * (dy/2)* T(iPoint,iTime)/ dt  /...
              rho * cp * (dx/2) * (dy/2)/ dt;
        end
        
    elseif iPoint < nx
         for iTime = 1:(nTimes-1)
            T(iPoint,iTime+1) = k*dy*(T(iPoint-1,iTime)- T(iPoint,iTime))/(2*dx)  +...
                k * dy *(T(iPoint+1,iTime)-T(iPoint,iTime))/ (2 * dx) + ... 
             k * dx*(T(iPoint+nx,iTime)-T(iPoint,iTime))/ ( dy) + ...
             h*(T_inf-T(iPoint,iTime))*dx +  rho * cp * (dx) * (dy/2)*T(iPoint,iTime)/ dt  /...
              rho * cp * (dx) * (dy/2)/ dt;
        end
    elseif iPoint == nx
        A(iPoint, iPoint-1) = k * dy/ (2*dx);
        A(iPoint, 2*nx) = k * dx/ (2*dy);
        A(iPoint, iPoint) = - k * dy/ (2*dx) - k * dx/ (2*dy) - h * dx/2 ...
            - h * dy/2 - rho * cp * (dx/2) * (dy/2)/ dt;
        B(iPoint, iPoint) = - rho * cp * (dx/2) * (dy/2)/ dt;
        C(iPoint) = - h * T_inf * (dx/2 + dy/2);
        for iTime = 1:(nTimes-1)
            T(iPoint,iTime+1) = k * dy *(T(iPoint-1,iTime)-T(iPoint,iTime))/ (2 * dx) + ... 
             k * dx*(T(iPoint+nx,iTime)-T(iPoint,iTime))/ (2 * dy) + ...
             h*(T_inf-T(iPoint,iTime))*dx/2 + h*(T_inf-T(iPoint,iTime))*dy/2 + ...
             rho * cp * (dx/2) * (dy/2)*T(iPoint,iTime)/ dt  /...
              rho * cp * (dx/2) * (dy/2)/ dt;
        end
        
    elseif ( iPoint == (nx+1) || iPoint == (2*nx+1) )
        for iTime = 1:(nTimes-1)
            T(iPoint,iTime+1) = k * dy *(T(iPoint-nx,iTime)-T(iPoint,iTime))/ (2 * dx) + ... 
             k * dy*(T(iPoint+1,iTime)-T(iPoint,iTime))/( dx) + ...
             k*dx*(T(iPoint+nx,iTime) - T(iPoint,iTime))/(2*dy) + ...
             rho * cp * (dx) * (dy)*T(iPoint,iTime)/ dt  /...
              rho * cp * (dx) * (dy)/ dt;
        end
    elseif iPoint < (2*nx)
        A(iPoint, iPoint-1) = k * dy/ dx;
        A(iPoint, iPoint+1) = k * dy/ dx;
        A(iPoint, iPoint - nx) = k * dx/ dy;
        A(iPoint, iPoint + nx) = k * dx/ dy;
        A(iPoint, iPoint) = - k * dy/ dx - k * dy/ dx - k * dx/ dy - ...
            k * dx/ dy - rho * cp * dx * dy/dt;
        B(iPoint, iPoint) = - rho * cp * dx * dy/dt; 
        for iTime = 1:(nTimes-1)
            T(iPoint,iTime+1) = k * dx *(T(iPoint-nx,iTime)-T(iPoint,iTime))/(dy) + ... 
             k * dy*(T(iPoint+1,iTime)-T(iPoint,iTime))/( dx) + ...
             k*dx*(T(iPoint+nx,iTime) - T(iPoint,iTime))/(dy) +...  
             k*dy*(T(iPoint-1,iTime)-T(iPoint,iTime))/(dx) + ...
             rho * cp * (dx) * (dy)*T(iPoint,iTime)/ dt  /...
              rho * cp * (dx) * (dy)/ dt;
        end
    elseif ( iPoint == (2*nx) || iPoint == (3*nx) )
        A(iPoint, iPoint-nx) = k * dx/ (2*dy);
        A(iPoint, iPoint+nx) = k * dx/ (2*dy);
        A(iPoint, iPoint-1) = k * dy/ dx;
        A(iPoint, iPoint) = - k * dx/ (2*dy) - k * dx/ (2*dy) - ...
            k * dy/ dx - h * dy - rho * cp * (dx/2) * dy/dt;
        B(iPoint, iPoint) = - rho * cp * (dx/2) * dy/dt;
        C(iPoint) = - h * T_inf * dy;
         for iTime = 1:(nTimes-1)
            T(iPoint,iTime+1) = k * dx *(T(iPoint-nx,iTime)-T(iPoint,iTime))/(2*dy) + ... 
             k * dy*(T(iPoint-1,iTime)-T(iPoint,iTime))/( dx) + ...
             k*dx*(T(iPoint+nx,iTime) - T(iPoint,iTime))/(2*dy) +...  
             h*dy*(T_inf - T(iPoint,iTime)) + ...
             rho * cp * (dx) * (dy)*T(iPoint,iTime)/ 2*dt  /...
              rho * cp * (dx) * (dy)/ 2*dt;
        end
    elseif iPoint < (3*nx)
        A(iPoint, iPoint-1) = k * dy/ dx;
        A(iPoint, iPoint+1) = k * dy/ dx;
        A(iPoint, iPoint - nx) = k * dx/ dy;
        A(iPoint, iPoint + nx) = k * dx/ dy;
        A(iPoint, iPoint) = - k * dy/ dx - k * dy/ dx - k * dx/ dy - ...
            k * dx/ dy - rho * cp * dx * dy/dt;
        B(iPoint, iPoint) = - rho * cp * dx * dy/dt;  
         for iTime = 1:(nTimes-1)
            T(iPoint,iTime+1) = k * dx *(T(iPoint-nx,iTime)-T(iPoint,iTime))/(dy) + ... 
             k * dy*(T(iPoint+1,iTime)-T(iPoint,iTime))/( dx) + ...
             k*dx*(T(iPoint+nx,iTime) - T(iPoint,iTime))/(dy) +...  
             k*dy*(T(iPoint-1,iTime)-T(iPoint,iTime))/(dx) + ...
             rho * cp * (dx) * (dy)*T(iPoint,iTime)/ dt  /...
              rho * cp * (dx) * (dy)/ dt;
        end
    elseif iPoint == (3*nx+1)
        A(iPoint, iPoint-nx) = k * dx/ (2 * dy);
        A(iPoint, iPoint+1) = k * dy/ dx;
        A(iPoint, iPoint) = - k * dx/ (2 * dy) - k * dy/ dx...
            - rho * cp * (dx/2) * dy/ dt - h * dx/2;
        B(iPoint, iPoint) = - rho * cp * (dx/2) * dy/ dt;
        C(iPoint) = - h * T_inf * dx/2;
         for iTime = 1:(nTimes-1)
            T(iPoint,iTime+1) = k * dx *(T(iPoint-nx,iTime)-T(iPoint,iTime))/(2*dy) + ... 
             k * dy*(T(iPoint+1,iTime)-T(iPoint,iTime))/( 2*dx) + ...  
             h*dx*(T_inf - T(iPoint,iTime)/2) + ...
             rho * cp * (dx) * (dy)*T(iPoint,iTime)/ 4*dt  /...
              rho * cp * (dx) * (dy)/ 4*dt;
        end
    elseif iPoint < 4*nx
        A(iPoint, iPoint-1) = k * dy/ (2*dx);
        A(iPoint, iPoint+1) = k * dy/ (2*dx);
        A(iPoint, iPoint-nx) = k * dx/ dy;
        A(iPoint, iPoint) = - k * dy/ (2*dx) - k * dy/ (2*dx) - k * dx/ dy - ...
            h * dx - rho * cp * dx * (dy/2)/ dt;
        B(iPoint, iPoint) = - rho * cp * dx * (dy/2)/ dt;
        C(iPoint) = - h * dx * T_inf;
         for iTime = 1:(nTimes-1)
            T(iPoint,iTime+1) = k * dx *(T(iPoint-nx,iTime)-T(iPoint,iTime))/(dy) + ... 
             k * dy*(T(iPoint+1,iTime)-T(iPoint,iTime))/( 2*dx) + ...  
             k*dy*(T(iPoint-1,iTime) - T(iPoint,iTime)/(2*dx))+...
             h*dx*(T_inf - T(iPoint,iTime)) + ...
             rho * cp * (dx) * (dy)*T(iPoint,iTime)/ 2*dt  /...
              rho * cp * (dx) * (dy)/ 2*dt;
        end
    else
        A(iPoint, iPoint-1) = k * dy/ (2*dx);
        A(iPoint, iPoint-nx) = k * dx/ (2*dy);
        A(iPoint, iPoint) = - k * dy/ (2*dx) - k * dx/ (2*dy) - ...
            h * (dx/2) - h * (dy/2) - rho * cp * (dx/2) * (dy/2)/ dt;
        B(iPoint, iPoint) = - rho * cp * (dx/2) * (dy/2)/ dt;
        C(iPoint) = - h * T_inf * (dx/2 + dy/2);
         for iTime = 1:(nTimes-1)
            T(iPoint,iTime+1) = k * dx *(T(iPoint-nx,iTime)-T(iPoint,iTime))/ (2 * dy) + ... 
             k * dy*(T(iPoint-1,iTime)-T(iPoint,iTime))/ (2 * dx) + ...
             h*(T_inf-T(iPoint,iTime))*dx/2 + h*(T_inf-T(iPoint,iTime))*dy/2 + ...
             rho * cp * (dx/2) * (dy/2)*T(iPoint,iTime)/ dt  /...
              rho * cp * (dx/2) * (dy/2)/ dt;
        end
        
    end
end

for iTime = 1:(nTimes-1)
    C1 = C;
    if t(iTime) <= delta_t_ON
        for iPoint = 1:12
            if ( (iPoint > 1) && (iPoint < 12) )
                C1(iPoint) = C1(iPoint) - q0_laser * dx;
            else
                C1(iPoint) = C1(iPoint) - q0_laser * dx/2;
            end
        end
    end
    
   
end
surf(T)
