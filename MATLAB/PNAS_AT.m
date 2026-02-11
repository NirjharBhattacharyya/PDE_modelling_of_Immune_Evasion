clear;close all; clc;
% --- Physical Parameters ---
theta = 19/20;          
Dm = 8.64E-7;               
chiP = 1E-5; % 
dM = 0.015;              
M0 = 5E-5;               
Kp = 1E-5;               
a_max = 0.2;            

% --- Numerical Params ---
Lx = 1.0; Ly = 1.0; Nx = 60; Ny = 60;
dx = Lx/Nx; dy = Ly/Ny; % For periodic, we divide by N (interval wrapping)
dt = 0.0001; Nt = 600000; % 60 Days

% --- Grid & Initial Conditions ---
[X, Y] = meshgrid(linspace(0, Lx-dx, Nx), linspace(0, Ly-dy, Ny));
M = zeros(Nx, Ny); 
%P = exp(-400*((X-0.7).^2 + (Y-0.7).^2)); 

radius = sqrt(0.05 / pi);
dist_sq = (X - 0.4).^2 + (Y - 0.4).^2;
%P = double(dist_sq <= radius^2);
P = (double(dist_sq <= radius^2))*267E-12;

% Helper function for periodic gradient
p_grad = @(F, d, dim) (circshift(F, -1, dim) - circshift(F, 1, dim)) / (2*d);

% Pre-calculate P gradients (Periodic)
gradPx = p_grad(P, dx, 2);
gradPy = p_grad(P, dy, 1);
alpha_P = a_max * (P ./ (Kp + P)); 

M_history = zeros(1, Nt);

fprintf('Simulation started (Periodic BCs)...\n');
for n = 1:Nt
    %Laplacian (Five-point stencil as in paper)
    L = 0.8 * Dm * ( (circshift(M, [0, -1]) - 2*M + circshift(M, [0, 1]))/dx^2 + ...
                     (circshift(M, [-1, 0]) - 2*M + circshift(M, [1, 0]))/dy^2 );
    
    %Recruitment & Death
    recruitment = alpha_P .* (M0 - M);
    death = dM * M;
    
    %Divergence of Flux (Periodic)
    fluxX = M .* chiP .* gradPx;
    fluxY = M .* chiP .* gradPy;
    
    div_J = p_grad(fluxX, dx, 2) + p_grad(fluxY, dy, 1);
    
    %Update timeStep
    change = (L + recruitment - theta*(div_J + death)) / theta;
    M = M + dt * change;
    
    %preserve positivty
    M = max(0, M); 
    
    %for plotting
    M_history(n) = mean(M(:));
    
    % Visualization
    if mod(n, 2000) == 0
        subplot(1,2,1);
        imagesc(M); colorbar; 
        title(['M (Day ', num2str(n*dt, '%.1f'), ')']);
        axis image;
        
        subplot(1,2,2);
        plot((1:n)*dt, M_history(1:n), 'b');
        xlabel('Days'); ylabel('Macrophage Concentration');
        title('GrowthCurve');
        drawnow;
    end
end