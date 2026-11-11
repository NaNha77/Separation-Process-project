% Stage 7: top-tray design following lecture 7, Example 21.6.
% DWSIM 28-position case, R_D = 0.99264, Stage 2 properties.
distillateKgS = 3.5560084416;
R_D = 0.99264;
rhoV = 1.15483;       % kg/m^3
rhoL = 754.931;       % kg/m^3
sigma = 0.0190639;    % N/m

% Example 21.6 sieve-tray geometry and head assumptions.
plateSpacingM = 0.500;
bubbleFraction = 0.70;
plateThicknessIn = 1/8;
holeDiameterIn = 1/4;
holePitchIn = 3/4;    % Triangular pitch.
weirHeightIn = 2;
C0 = 0.73;            % Fig. 21.25, with the example's hole geometry.
beta = 0.60;
phi_d = 0.50;
h_f_L = 10;           % mm of liquid
g = 9.81;             % m/s^2

%% (a) Flooding limit, vapor flow, and diameter
plateSpacingIn = plateSpacingM / 0.0254;
plateSpacingMm = plateSpacingIn * 25.4;

% A total condenser splits the condensed vapor into reflux and distillate.
w_L = distillateKgS * R_D;
w_V = distillateKgS * (1 + R_D);
F_LV = (w_L / w_V) * sqrt(rhoV / rhoL);

% Lygeros-Magoulas fit in DWSIM's guide: spacing in mm, K1 in m/s.
% https://dwsim.org/help/main/24-unit-operations.html#entrainment-flooding
fairCapacity = @(f) 0.0105 + 8.127e-4 * plateSpacingMm^0.755 ...
    .* exp(-1.463 * f.^0.842);
K1 = fairCapacity(F_LV);
K_v = K1 / 0.3048;    % ft/s, for the lecture equation.
sigmaDynCm = 1000 * sigma;
u_c_ft_s = K_v * sqrt((rhoL - rhoV) / rhoV) * (sigmaDynCm / 20)^0.2;
u_c = u_c_ft_s * 0.3048;

% Use the bubbling-area convention and limiting speed of Example 21.6.
Q_V = w_V / rhoV;
A_b = Q_V / u_c;
A_c = A_b / bubbleFraction;
D_c = sqrt(4 * A_c / pi);

%% (b) Pressure drop across the top tray
holeFraction = (pi / (2 * sqrt(3))) * (holeDiameterIn / holePitchIn)^2;
A_h = holeFraction * A_b;
u_0 = Q_V / A_h;
h_d = 51.0 * (u_0 / C0)^2 * (rhoV / rhoL);  % mm of liquid

h_w = weirHeightIn * 25.4;
L_w = 1.62 * D_c / 2;  % Example 21.6 chord factor for 15% area on each side.
q_L = 60 * w_L / rhoL; % m^3/min for the Francis equation.
h_ow = 43.4 * (q_L / L_w)^(2/3);
h_l = beta * (h_w + h_ow);
h_t = h_d + h_l;
deltaP_kPa = rhoL * g * (h_t / 1000) / 1000;

%% (c) Downcomer froth height and tray spacing
Z_c = h_l + h_t + h_f_L;
Z = Z_c / phi_d;
frothMarginMm = plateSpacingMm - Z;
passesFrothCheck = Z < plateSpacingMm;

fprintf('Example 21.6 method | top tray (Stage 2)\n');
fprintf('Spacing = %.3f m = %.3f in = %.1f mm\n', plateSpacingM, plateSpacingIn, plateSpacingMm);
fprintf('Distillate = %.2f kg/h | w_L = %.6f kg/s | w_V = %.6f kg/s\n', ...
    3600 * distillateKgS, w_L, w_V);
fprintf('F_LV = %.8f | K1 = %.6f m/s | K_v = %.6f ft/s\n', F_LV, K1, K_v);
fprintf('u_c = %.4f ft/s = %.4f m/s | Q_V = %.4f m^3/s\n', u_c_ft_s, u_c, Q_V);
fprintf('A_b = %.4f m^2 | A_c = %.4f m^2 | D_c = %.4f m\n', A_b, A_c, D_c);
fprintf('Tray thickness / hole diameter = %.2f | A_h / A_b = %.4f %%\n', ...
    plateThicknessIn / holeDiameterIn, 100 * holeFraction);
fprintf('A_h = %.4f m^2 | u_0 = %.4f m/s | h_d = %.2f mm\n', A_h, u_0, h_d);
fprintf('L_w = %.4f m | q_L = %.5f m^3/min | h_ow = %.2f mm\n', L_w, q_L, h_ow);
fprintf('h_w = %.2f mm | h_l = %.2f mm | h_t = %.2f mm | dP = %.4f kPa/tray\n', ...
    h_w, h_l, h_t, deltaP_kPa);
fprintf('Z_c = %.2f mm | Z = %.2f mm | spacing - Z = %.2f mm\n', Z_c, Z, frothMarginMm);
fprintf('Downcomer froth criterion (1 = passes): %d\n', passesFrothCheck);

%% Fair curve at the actual tray spacing
curveFlow = logspace(log10(min(0.012, F_LV)), log10(max(0.100, F_LV)), 250);
fig = figure('Name', 'Fair flooding capacity', 'NumberTitle', 'off', 'Color', 'white');
theme(fig, 'light');
ax = axes(fig);
curve = loglog(ax, curveFlow, fairCapacity(curveFlow) / 0.3048, ...
    'Color', [0.10 0.28 0.43], 'LineWidth', 2);
hold(ax, 'on');
point = scatter(ax, F_LV, K_v, 75, [0.69 0.29 0.09], 'filled');
grid(ax, 'on');
xlabel(ax, 'Flow parameter, F_{LV}');
ylabel(ax, 'Flooding capacity factor, K_v (ft/s)');
title(ax, sprintf('Fair capacity factor | %.3f m (%.2f in) tray spacing', ...
    plateSpacingM, plateSpacingIn));
subtitle(ax, 'Lygeros-Magoulas fit documented by DWSIM');
legend(ax, [curve point], {'Fair equation', sprintf('Top tray: K_v = %.3f ft/s', K_v)}, ...
    'Location', 'southwest');

plotPath = fullfile(fileparts(mfilename('fullpath')), 'fair_0p5m_example_21_6.png');
temporaryPlotPath = [tempname '.png'];
exportgraphics(fig, temporaryPlotPath, 'Resolution', 200, 'Padding', 18);
[plotSaved, plotSaveMessage] = movefile(temporaryPlotPath, plotPath, 'f');
assert(plotSaved, 'Could not save the Fair chart: %s', plotSaveMessage);
