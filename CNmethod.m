N = 20;
h = 1/N;
a = 1;
tau = h;
T = 0.2;
M = round(T/tau);
lambda = a*tau/(h^2); 
t_cur = 0;
i = 1:1:(N-1);
u = (sin(pi*i*h))'; % initial condition, change as necessary
error = 0;

for r=1:1:M         % march forward through time      

    % CN Scheme is A_tau*u_n+1=B_tau*u_n=u_int
    u_int = zeros(N-1,1);
    u_int(1) = (1-lambda)*u(1)+(lambda/2)*u(2);
    u_int(N-1) = (1-lambda)*u(N-1)+(lambda/2)*u(N-2);
    for (p = 2:1:(N-2))
        u_int(p) = (1-lambda)*u(p)+(lambda/2)*u(p-1)+(lambda/2)*u(p+1); % this multiplies B_tau and u_n
    end
    % Now to solve the linear system A_tau*u_n+1=u_int
    % tridiagonal system where a_dia is array of diagonal elements and b_off
    % and c_off are arrays of offdiagonal elements
    a_dia = (ones(1,N-1)*(1+lambda))';
    b_dia = (ones(1,N-2)*(-lambda/2))';
    c_dia = (ones(1,N-2)*(-lambda/2))';
    
    % using thomas algorithm
    % creating zeros
    
    for (j = 2:1:(N-1))
        a_dia(j) = a_dia(j) - (b_dia(j-1)/a_dia(j-1))*c_dia(j-1);
        u_int(j) = u_int(j) - (b_dia(j-1)/a_dia(j-1))*u_int(j-1);
    end
    u_int(N-1) = u_int(N-1)/a_dia(N-1);
    
    % backward substitution
    
    for (k = (N-2):(-1):1)
        u_int(k) = (u_int(k)-((c_dia(k))*u_int(k+1)))/(a_dia(k));
    end
    t_cur = t_cur+tau;
    u = u_int; % solution is updated
    u_exact = (sin(pi*i*h)*exp(-pi*pi*t_cur))';
    errmat = (abs(u-u_exact));
    
end
u_exact = (sin(pi*i*h)*exp(-pi*pi*t_cur))';
errmat = (abs(u-u_exact));
error = max(errmat);
