% Stage 8: Murphree effective-equilibrium stepping for methanol-water.
% The embedded DWSIM UNIFAC bubble-point data are at 101.3 kPa.
x_F = 0.500;
x_D = 0.9983;
x_B = 0.0010;
R_D = 0.99264;
eta_M = 0.453;           % First estimate from the O'Connell overall efficiency.
efficiencyShift = 0.05;  % Absolute efficiency change for the sensitivity cases.

assert(0 < x_B && x_B < x_F && x_F < x_D && x_D < 1 && R_D > 0, ...
    'Expected 0 < x_B < x_F < x_D < 1 and positive reflux.');
efficiencies = [eta_M; eta_M-efficiencyShift; eta_M+efficiencyShift; 1];
assert(all(efficiencies > 0 & efficiencies <= 1), ...
    'All Murphree efficiencies must be in (0, 1].');

vle = unifacVle();
yEq = griddedInterpolant(vle(:,1), vle(:,2), 'pchip', 'none');

% Saturated liquid feed (q = 1). The operating lines meet at x_F.
rectSlope = R_D/(R_D+1);
yAtFeed = rectSlope*x_F + x_D/(R_D+1);
stripSlope = (yAtFeed-x_B)/(x_F-x_B);
yOp = @(x) yAtFeed + stripSlope.*(x-x_F) ...
    + (rectSlope-stripSlope).*max(x-x_F, 0);
gridX = linspace(x_B, x_D, 2001);
assert(all(yEq(gridX) > yOp(gridX)), ...
    'An operating line touches or crosses equilibrium in the separation range.');

cases = ["Base"; "Lower efficiency"; "Higher efficiency"; "Ideal trays"];
for k = numel(efficiencies):-1:1
    results(k) = countTrays(efficiencies(k), x_F, x_D, x_B, yEq, yOp);
end
trayCount = [results.trayCount].';
feedTray = [results.feedTray].';
summary = table(cases, efficiencies, trayCount, feedTray, trayCount+1, ...
    trayCount+2, feedTray+1, [results.xBottom].', 'VariableNames', ...
    {'Case','Eta_M','Trays','Feed_tray','Trays_plus_reboiler', ...
    'DWSIM_positions','DWSIM_feed_position','xB_estimated'});
disp(summary);

base = results(1);
section = repmat("Rectifying",base.trayCount,1);
section(base.feedTray) = "Feed";
section(base.feedTray+1:end) = "Stripping";
stageResult = table((1:base.trayCount).', base.xTray, base.yTray, ...
    base.yOut, section, 'VariableNames', ...
    {'Tray','x_liquid_leaving','y_vapor_leaving','y_vapor_entering','Section'});
disp(stageResult);
fprintf('Reboiler: x_B = %.8f, y = %.8f; target x_B <= %.6f\n', ...
    base.xBottom, base.yReboiler, x_B);

% Plot the entire construction and a magnified view of the top trays.
fig = figure('Name','Murphree McCabe-Thiele', 'NumberTitle','off', ...
    'Color','white','Position',[80 100 1260 640]);
theme(fig,'light');
layout = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
title(layout, sprintf('Methanol-water at 101.3 kPa | R_D = %.5f | eta_M = %.3f', R_D, eta_M));
ax = nexttile(layout);
handles = drawDiagram(ax, base, eta_M, x_F, x_D, x_B, yEq, yOp);
axis(ax,'square');
xlim(ax,[0 1]); ylim(ax,[0 1]);
title(ax,sprintf('%d trays + reboiler | feed tray %d',base.trayCount,base.feedTray));
legend(ax,handles,{'UNIFAC equilibrium','Effective equilibrium','Operating lines', ...
    'y = x','Stage steps','Feed tray'},'Location','southeast','FontSize',8);

detail = nexttile(layout);
drawDiagram(detail, base, eta_M, x_F, x_D, x_B, yEq, yOp);
axis(detail,'square');
xlim(detail,[0.96 1]); ylim(detail,[0.96 1]);
title(detail,'Top-tray detail');

plotPath = fullfile(fileparts(mfilename('fullpath')),'murphree_mccabe_thiele.png');
temporaryPlot = [tempname '.png'];
exportgraphics(fig,temporaryPlot,'Resolution',200,'Padding',18);
[saved,saveMessage] = movefile(temporaryPlot,plotPath,'f');
assert(saved,'Could not save the diagram: %s',saveMessage);

function result = countTrays(eta, x_F, x_D, x_B, yEq, yOp)
    yEffective = @(x) yOp(x) + eta*(yEq(x)-yOp(x));
    yAtBottomTarget = yEq(x_B);
    trayLimit = 1000;
    xTray = zeros(trayLimit,1);
    yTray = zeros(trayLimit,1);
    yOut = zeros(trayLimit,1);
    y = x_D;
    n = 0;
    % Stop when one true-equilibrium reboiler step reaches the bottoms target.
    while y > yAtBottomTarget && n < trayLimit
        n = n+1;
        xTray(n) = fzero(@(x) yEffective(x)-y, [x_B x_D]);
        yTray(n) = y;
        yOut(n) = yOp(xTray(n));
        assert(yOut(n) < y, 'Stage stepping failed to advance.');
        y = yOut(n);
    end
    assert(y <= yAtBottomTarget, 'Could not reach the bottoms target in 1000 trays.');
    xTray = xTray(1:n); yTray = yTray(1:n); yOut = yOut(1:n);
    feedTray = find(xTray <= x_F,1);
    assert(~isempty(feedTray), 'No tray crossed the saturated-liquid feed line.');
    xBottom = fzero(@(x) yEq(x)-y,[0 x_B]);
    result = struct('trayCount',n,'feedTray',feedTray,'xTray',xTray, ...
        'yTray',yTray,'yOut',yOut,'xBottom',xBottom,'yReboiler',y);
end

function handles = drawDiagram(ax, result, eta, x_F, x_D, x_B, yEq, yOp)
    hold(ax,'on');
    x = linspace(0,1,2401);
    xEffective = linspace(x_B,x_D,2401);
    yEffective = yOp(xEffective) + eta*(yEq(xEffective)-yOp(xEffective));
    handles(1) = plot(ax,x,yEq(x),'Color',[0.12 0.35 0.61],'LineWidth',1.9);
    handles(2) = plot(ax,xEffective,yEffective,'--','Color',[0.69 0.25 0.15],'LineWidth',1.8);
    handles(3) = plot(ax,[x_B x_F x_D],yOp([x_B x_F x_D]), ...
        'Color',[0.17 0.50 0.29],'LineWidth',1.5);
    handles(4) = plot(ax,[0 1],[0 1],':','Color',[0.4 0.4 0.4]);
    plot(ax,[x_F x_F],[x_F yOp(x_F)],'-.','Color',[0.40 0.40 0.40]);
    stairsX = [x_D; reshape([result.xTray,result.xTray].',[],1); result.xBottom; result.xBottom];
    stairsY = [x_D; reshape([result.yTray,result.yOut].',[],1); result.yReboiler; result.xBottom];
    handles(5) = plot(ax,stairsX,stairsY,'Color',[0.12 0.12 0.17],'LineWidth',0.9);
    handles(6) = scatter(ax,result.xTray(result.feedTray), ...
        result.yTray(result.feedTray),42,[0.51 0.23 0.64],'filled');
    xlabel(ax,'Liquid methanol mole fraction, x');
    ylabel(ax,'Vapor methanol mole fraction, y');
    box(ax,'on'); grid(ax,'on');
end

function data = unifacVle()
    % DWSIM UNIFAC at 101300 Pa. Columns: x_e, y_e (methanol mole fractions).
    data = [
        0.00000000 0.000000000000;
        0.00001000 0.000079919121;
        0.00005000 0.000399453016;
        0.00010000 0.000798549838;
        0.00025000 0.001993707785;
        0.00050000 0.003978557060;
        0.00075000 0.005954605361;
        0.00100000 0.007921909749;
        0.00150000 0.011830512648;
        0.00200000 0.015704812826;
        0.00300000 0.023352256633;
        0.00400000 0.030867673121;
        0.00500000 0.038254383325;
        0.01000000 0.073366341529;
        0.01500000 0.105702071782;
        0.02000000 0.135574653220;
        0.02500000 0.163253143010;
        0.03000000 0.188969895585;
        0.03500000 0.212926491774;
        0.04000000 0.235298572362;
        0.04500000 0.256239802179;
        0.05000000 0.275885139583;
        0.05500000 0.294353547574;
        0.06000000 0.311750253376;
        0.06500000 0.328168640822;
        0.07000000 0.343691842532;
        0.07500000 0.358394085385;
        0.08000000 0.372341832263;
        0.08500000 0.385594754799;
        0.09000000 0.398206565284;
        0.09500000 0.410225730721;
        0.10000000 0.421696087863;
        0.10500000 0.432657374725;
        0.11000000 0.443145691376;
        0.11500000 0.453193900665;
        0.12000000 0.462831977706;
        0.12500000 0.472087315564;
        0.13000000 0.480984993331;
        0.13500000 0.489548011859;
        0.14000000 0.497797501543;
        0.14500000 0.505752905942;
        0.15000000 0.513432144400;
        0.15500000 0.520851756410;
        0.16000000 0.528027030037;
        0.16500000 0.534972116408;
        0.17000000 0.541700131981;
        0.17500000 0.548223250083;
        0.18000000 0.554552782982;
        0.18500000 0.560699255633;
        0.19000000 0.566672472033;
        0.19500000 0.572481575041;
        0.20000000 0.578135100404;
        0.20500000 0.583641025605;
        0.21000000 0.589006814125;
        0.21500000 0.594239455593;
        0.22000000 0.599345502261;
        0.22500000 0.604331102201;
        0.23000000 0.609202029536;
        0.23500000 0.613963712041;
        0.24000000 0.618621256334;
        0.24500000 0.623179470939;
        0.25000000 0.627642887398;
        0.25500000 0.632015779625;
        0.26000000 0.636302181684;
        0.26500000 0.640505904119;
        0.27000000 0.644630548988;
        0.27500000 0.648679523695;
        0.28000000 0.652656053764;
        0.28500000 0.656563194615;
        0.29000000 0.660403842444;
        0.29500000 0.664180744300;
        0.30000000 0.667896507392;
        0.30500000 0.671553607732;
        0.31000000 0.675154398140;
        0.31500000 0.678701115675;
        0.32000000 0.682195888546;
        0.32500000 0.685640742530;
        0.33000000 0.689037606951;
        0.33500000 0.692388320244;
        0.34000000 0.695694635143;
        0.34500000 0.698958223525;
        0.35000000 0.702180680925;
        0.35500000 0.705363530758;
        0.36000000 0.708508228265;
        0.36500000 0.711616164211;
        0.37000000 0.714688668333;
        0.37500000 0.717727012587;
        0.38000000 0.720732414182;
        0.38500000 0.723706038428;
        0.39000000 0.726649001411;
        0.39500000 0.729562372506;
        0.40000000 0.732447176737;
        0.40500000 0.735304396999;
        0.41000000 0.738134976146;
        0.41500000 0.740939818959;
        0.42000000 0.743719794003;
        0.42500000 0.746475735370;
        0.43000000 0.749208444330;
        0.43500000 0.751918690885;
        0.44000000 0.754607215238;
        0.44500000 0.757274729184;
        0.45000000 0.759921917414;
        0.45500000 0.762549438760;
        0.46000000 0.765157927370;
        0.46500000 0.767747993812;
        0.47000000 0.770320226131;
        0.47500000 0.772875190845;
        0.48000000 0.775413433887;
        0.48500000 0.777935481502;
        0.49000000 0.780441841097;
        0.49500000 0.782933002047;
        0.50000000 0.785409436463;
        0.50500000 0.787871599916;
        0.51000000 0.790319932131;
        0.51500000 0.792754857644;
        0.52000000 0.795176786427;
        0.52500000 0.797586114482;
        0.53000000 0.799983224411;
        0.53500000 0.802368485949;
        0.54000000 0.804742256484;
        0.54500000 0.807104881544;
        0.55000000 0.809456695262;
        0.55500000 0.811798020825;
        0.56000000 0.814129170895;
        0.56500000 0.816450448016;
        0.57000000 0.818762145000;
        0.57500000 0.821064545298;
        0.58000000 0.823357923352;
        0.58500000 0.825642544931;
        0.59000000 0.827918667454;
        0.59500000 0.830186540301;
        0.60000000 0.832446405103;
        0.60500000 0.834698496029;
        0.61000000 0.836943040052;
        0.61500000 0.839180257210;
        0.62000000 0.841410360854;
        0.62500000 0.843633557883;
        0.63000000 0.845850048973;
        0.63500000 0.848060028793;
        0.64000000 0.850263686216;
        0.64500000 0.852461204518;
        0.65000000 0.854652761573;
        0.65500000 0.856838530032;
        0.66000000 0.859018677505;
        0.66500000 0.861193366728;
        0.67000000 0.863362755728;
        0.67500000 0.865526997978;
        0.68000000 0.867686242549;
        0.68500000 0.869840634256;
        0.69000000 0.871990313793;
        0.69500000 0.874135417870;
        0.70000000 0.876276079343;
        0.70500000 0.878412427333;
        0.71000000 0.880544587350;
        0.71500000 0.882672681405;
        0.72000000 0.884796828120;
        0.72500000 0.886917142837;
        0.73000000 0.889033737716;
        0.73500000 0.891146721839;
        0.74000000 0.893256201300;
        0.74500000 0.895362279300;
        0.75000000 0.897465056233;
        0.75500000 0.899564629773;
        0.76000000 0.901661094954;
        0.76500000 0.903754544249;
        0.77000000 0.905845067651;
        0.77500000 0.907932752741;
        0.78000000 0.910017684759;
        0.78500000 0.912099946680;
        0.79000000 0.914179619270;
        0.79500000 0.916256781158;
        0.80000000 0.918331508893;
        0.80500000 0.920403877008;
        0.81000000 0.922473958073;
        0.81500000 0.924541822753;
        0.82000000 0.926607539863;
        0.82500000 0.928671176417;
        0.83000000 0.930732797682;
        0.83500000 0.932792467224;
        0.84000000 0.934850246955;
        0.84500000 0.936906197181;
        0.85000000 0.938960376644;
        0.85500000 0.941012842564;
        0.86000000 0.943063650685;
        0.86500000 0.945112855306;
        0.87000000 0.947160509330;
        0.87500000 0.949206664294;
        0.88000000 0.951251370409;
        0.88500000 0.953294676592;
        0.89000000 0.955336630505;
        0.89500000 0.957377278583;
        0.90000000 0.959416666068;
        0.90500000 0.961454837040;
        0.91000000 0.963491834448;
        0.91500000 0.965527700138;
        0.92000000 0.967562474880;
        0.92500000 0.969596198397;
        0.93000000 0.971628909393;
        0.93500000 0.973660645575;
        0.94000000 0.975691443681;
        0.94500000 0.977721339503;
        0.95000000 0.979750367910;
        0.95500000 0.981778562874;
        0.96000000 0.983805957488;
        0.96500000 0.985832583989;
        0.97000000 0.987858473780;
        0.97500000 0.989883657450;
        0.98000000 0.991908164792;
        0.98500000 0.993932024825;
        0.99000000 0.995955265809;
        0.99500000 0.997977915267;
        0.99600000 0.998382376547;
        0.99700000 0.998786815450;
        0.99800000 0.999191232188;
        0.99830000 0.999312552918;
        0.99850000 0.999393432310;
        0.99900000 0.999595626969;
        0.99925000 0.999696722258;
        0.99950000 0.999797816190;
        0.99975000 0.999898908770;
        0.99990000 0.999959563670;
        0.99995000 0.999979781862;
        0.99999000 0.999995956377;
        1.00000000 1.000000000000;
    ];
end
