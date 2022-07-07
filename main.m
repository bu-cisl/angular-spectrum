% physical parameters
pix_xyz = [0.1, 0.1, 0.1]; % pixel size, um
lambda0 = 0.6; % wavelength in vaccum, um
n0 = 1.0; % background refractive index
res = pix_xyz./lambda0.*n0; % note: only 'res' is used in the following code

% initialization
shape = [512 512]; % xy FOV, 51.2um x 51.2um
mu = multipliers(shape, res);

% angular frequency of z
kz = mu.gamma * 2 * pi * res(3);
% Remove evanescent wave
eva = exp(min((mu.gamma - 0.2) * 5, 0));


%% Example 1: a tilted gaussian beam
% You can see the gaussian spot moving along the wave vector, and all the
% energy is absorbed by the boundary finally.
u = mu.tilt([sin(pi/6), 0]) .* mu.gaussian([0.5 0.5], [50 50]);
figure;
imagesc(abs(u).^2); caxis([0, 1]); colorbar; axis image; drawnow;

steps = 100;
dz = 6; % px
p_mat = exp(1i * dz * kz) .* eva;

for i = 1:steps
    fprintf("z = %.1f um\n", i*dz*pix_xyz(3));
    a = fft2(u);
    a = a .* p_mat;
    u = ifft2(a);
    u = u .* mu.soft_crop(0.2, 1, steps); % change strength smaller and you can see energy leak from the other side
    % pause(0.05);
    imagesc(abs(u).^2); caxis([0, 1]); colorbar; axis image; drawnow;
end
pause(1)


%% Example 2: negative circular diffraction, on-axis illumination
% You can first see the famous Arago spot (Poisson spot) and then a
% Fraunhofer diffraction pattern. Boundary absorption ensures there is no
% fake interference pattern of the energy that goes out of the FOV.
u = mu.tilt([0, 0]);
current_phase = 1;
% lazy way to make a circular shadow, r=30pix=3um
u = u .* (mu.gaussian([0.4, 0.4], [30, 30]) < 1/exp(1));
figure;
imagesc(abs(u).^2); caxis([0, 2]); colorbar; axis image; drawnow;

steps = 100;
dz = 6; % px
p_mat = exp(1i * dz * kz) .* eva;

for i = 1:steps
    fprintf("z = %.1f um\n", i*dz*pix_xyz(3));
    a = fft2(u);
    a = a .* p_mat;
    u = ifft2(a);
    current_phase = current_phase * p_mat(1,1);
    u = (u-current_phase) .* mu.soft_crop(0.2, 1, steps) + current_phase;
    % pause(0.05);
    imagesc(abs(u).^2); caxis([0, 2]); colorbar; axis image; drawnow;
end
pause(1)


%% Example 3: negative rectangle diffraction, off-axis illumination
% You can see the combined result of diffraction, pattern moving because of
% off-axis illumination, and boundary absorption.
u = mu.tilt([sin(pi/6), sin(pi/6)]);
current_bg = u;

u(250:290, 150:300) = 0;
figure;
imagesc(abs(u).^2); caxis([0, 2]); colorbar; axis image; drawnow;

steps = 100;
dz = 6; % px
p_mat = exp(1i * dz * kz) .* eva;

for i = 1:steps
    fprintf("z = %.1f um\n", i*dz*pix_xyz(3));
    a = fft2(u);
    a = a .* p_mat;
    u = ifft2(a);
    current_bg = ifft2(fft2(current_bg) .* p_mat);
    u = (u-current_bg) .* mu.soft_crop(0.2, 1, steps) + current_bg;
    % pause(0.05);
    imagesc(abs(u).^2); caxis([0, 2]); colorbar; axis image; drawnow;
end
