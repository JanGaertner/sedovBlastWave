%% Analytical solution to the Sedov shock wave problem
% Author: Jan Gärtner
% Further details are described in:
% http://www.ita.uni-heidelberg.de/~dullemond/lectures/num_fluid_2011/Chapter_10.pdf
%
%
% License:
% Copyright 2020 Jan Gärtner
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions 
% are met:
% 1. Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright 
%    notice, this list of conditions and the following disclaimer in the 
%    documentation and/or other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
% TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF 
% THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

rho0 = 10;      % Initial density
r    = 0.0011;   % Radius of high energy region
p    = 1E+7;    % Initial pressure
p1   = 2870;    % Ambient pressure
gamma= 1.4;     % Isotropic exponent
t    = 5E-05;   % Time after simulation start
d    = 2;       % 2D or 3D simulation



E0 = p*(d+1)/3*(pi*r.^d)/(gamma-1);

R = (E0*t.^2/rho0).^(1/(d+2));


p2   = p;
rho2 = (gamma+1)/(gamma-1)*rho0;



%% Numerical integration of Sedov-Taylor problem

% these are the start conditions
eta_s = 1;
RES = 1;
iter = 0;
maxIter=1000;
while RES > 1E-6
    deltaEta = (eta_s/100);
    eta = (linspace(0,eta_s,101))';
    A = zeros(numel(eta),1);
    B = zeros(numel(eta),1);
    C = zeros(numel(eta),1);
    D = zeros(numel(eta),1);
    
    A(end) = 1;
    B(end) = 1;
    C(end) = 1;
    D(end) = 1;
    
    
    for k = (flip(1:numel(eta)))
        if (k==1)
            break;
        end
        
        % Coefficient for PDE matrix
        A1 = -eta(k) + (2*eta(k)*C(k))/(gamma+1);
        B1 = 0;
        C1 = (2*eta(k)*A(k)/(gamma+1));
        D1 = 0;
        K1 = 6/(gamma+1)*A(k)*C(k);
        
        % Short hand notation
        beta = 4/(5*(gamma+1));
        
        A2 = 0;
        B2 = (2/5)*((gamma-1)/(gamma+1))*eta(k)/A(k);
        C2 = -(2/5)*eta(k) + beta*C(k)*eta(k);
        D2 = 0;
        K2 = -C(k)+beta*C(k).^2+(2/5)*((gamma-1)/(gamma+1))*(2*B(k))/A(k);
        
        A3 = -(2/5)*eta(k)*D(k) + beta*eta(k)*C(k)*D(k);
        B3 = -(2/5)*eta(k)+beta*eta(k)*gamma*C(k);
        C3 = beta*eta(k)*gamma*B(k) + beta*eta(k)*A(k)*C(k).^2;
        D3 = -(2/5)*eta(k)*A(k) + beta*eta(k)*A(k)*C(k);
        K3 = -2*(B(k)+A(k)*D(k)) + beta*(5*C(k)*(gamma*B(k)+A(k)*C(k).^2));
        
        A4 = 0;
        B4 = 0;
        C4 = -2*C(k);
        D4 = 1;
        K4 = 0;
        
        M = [A1,B1,C1,D1;...
            A2,B2,C2,D2;...
            A3,B3,C3,D3;...
            A4,B4,C4,D4];
        
        b = M*[A(k);B(k);C(k);D(k)]+ [K1;K2;K3;K4]*deltaEta;
        
        x = M\b;
        A(k-1) = x(1);
        B(k-1) = x(2);
        C(k-1) = x(3);
        D(k-1) = x(4);
    end
    
    % Integrate Energy to check solution
    E = 0;
    
    for k = 1:(numel(eta)-1)
        l = (B(k)+A(k)*C(k).^2)*eta(k).^4;
        r = (B(k+1)+A(k+1)*C(k+1).^2)*eta(k+1).^4;
        
        E = E + (l+r)/2*(eta(k+1)-eta(k));
    end
    E = E*(32*pi)/(25*(gamma.^2-1));
    
    if (E < 1)
        eta_s = eta_s*(2-E);
    else
        temp = 1 + (eta_s-1)/E;
        % blend the value
        eta_s = (temp+eta_s)/2;
    end

    % Error
    RES = abs(1-E);
    iter=iter+1;
    if (iter>maxIter)
        fprintf('\n--> Reached Max Iterations\n')
        break;
    end
end

fprintf('eta_s = %f\n',eta_s);
fprintf('With an accuracy of %f\n',1-abs(1-E)*100)


%% Plot Results

r = eta/eta_s.*(E0.*t.^2./rho0).^(1/(d+2));

rho = rho2*A;
P = p2.*(eta./eta_s).^2.*B;

% Extend region for nices plots
r(end+1)=r(end);
r(end+1)=2*r(end);

P(end+1)=p1;
P(end+1)=p1;

rho(end+1)=rho0;
rho(end+1)=rho0;

eta(end+1) = eta(end);
eta(end+1) = r(end)/(E0.*t.^2./rho0).^(1/(d+2))*eta_s;

subplot(2,2,1)
plot(eta,rho)
title('Density')

subplot(2,2,2)
plot(eta,P)
title('Pressure')

subplot(2,2,3)
plot(r,rho);

subplot(2,2,4)
plot(r,P)




%% Write out data:
M = [r,P,rho];
[filename,path] = uiputfile('*.dat','Save results for gnuplot');


dlmwrite(fullfile(path,filename),M,'delimiter',' ');
