% Compartmentalise field shift maps & ks density histograms|
% KS density estimates - does each method follow a sine wave with respect to orientation?
% Add orthogonal plane cylinder growth to mimic crossing fibres, and a weighting function for fibre populations/densities

clear all; close all; clc;

%% --- Parameters -----------------------------------------------
% Environment
voxelSize         = 200; % in µm
axonTargetFill    = 0.15; % percentage target volume
sphereTargetFill  = axonTargetFill * 0.005; % target volume for spheres (modulated by axon volume)

% Geometry
gRatioMean = 0.65;        % G-ratio mean 
gRatioStd  = 0.05;        % G-ratio standard deviation
rMean      = 4;           % Mean axon radius (in µm)
rStd       = 1;         % Standard deviation of axon radius
undAmp     = 1.5;         % Amplitude of sinusoidal undulations in the axon's radial axis
undFreq    = 0.01;        % Frequency of undulations
numRings   = 200;         % Number of cross-sectional 'rings' along each axon
circleRes  = 64;          % Number of points used to define each ring
sphereRadiusRange = [0.4, 2];  % Range of radii for the randomly placed spheres (e.g., cells, obstacles)

segmentLength = round(numRings / 35);   % Number of rings per axon segment (n segments total per axon)
demyelinationProb = 0.2;                % Probability that a segment along an axon is demyelinated
rng(1);  % Set random seed for reproducibility (same seed for each iteration)
%rng('shuffle');  % Random seed


%% --- Initialization -----------------------------------------------------
axonMeshes       = {};
boundingSpheres  = [];
totalVol         = voxelSize^3;
axonFilledVol    = 0;
sphereFilledVol  = 0;
axonCount        = 0;
z                = linspace(0, voxelSize, numRings);

%% --- Compute cylinder face connectivity -------------------
faces = [];
for i = 1:(numRings - 1)
    for j = 1:circleRes
        i1 = (i-1)*circleRes + j;
        i2 = (i-1)*circleRes + mod(j,circleRes)+1;
        i3 = i1 + circleRes;
        i4 = i2 + circleRes;
        faces = [faces; i1 i3 i4; i1 i4 i2];
    end
end

%% --- Grow cylinders ------------------------------------------------
h = waitbar(0,'Growing axons...');
while (axonFilledVol/totalVol < axonTargetFill)

    rAxon   = max(0.2, rMean + randn * rStd);
    gRatio = min(max(normrnd(gRatioMean, gRatioStd), 0.3), 0.95);
    rMyelin = rAxon / gRatio;
    pathLen = voxelSize; % Length is equal to the size of the simulated voxel (i.e 'infinite' cylinders)

    % Random origin and phase (start orientation)
    x0 = rand * voxelSize;
    y0 = rand * voxelSize;
    phaseX = 2*pi*rand;
    phaseY = 2*pi*rand;

    % Generate radial axis
    x = x0 + undAmp * sin(2*pi*undFreq*z + phaseX);
    y = y0 + undAmp * sin(2*pi*undFreq*z + phaseY);
    centerline = [x; y; z];

    % Bounding sphere
    cx = mean(x); cy = mean(y); cz = mean(z);
    radiusBound = rMyelin + undAmp;

    if ~isempty(boundingSpheres)
        dists = sqrt((boundingSpheres(:,1)-cx).^2 + (boundingSpheres(:,2)-cy).^2 + (boundingSpheres(:,3)-cz).^2);
        overlaps = dists < (boundingSpheres(:,4) + radiusBound);
        if any(overlaps)
            continue
        end
    end

    % Build mesh along radial axis
    theta = linspace(0,2*pi,circleRes+1); theta(end)=[];
    circle = [cos(theta); sin(theta)];
    vertsMyelin = zeros(numRings*circleRes,3);
    vertsAxon   = zeros(numRings*circleRes,3);
    tangents    = zeros(numRings,3);
    ringCenters = zeros(numRings,3);

    % --- Segment-based demyelination
    demyelinMap = true(1, numRings);
    startIdxs = 1:segmentLength:numRings;
    for si = startIdxs
        if rand < demyelinationProb
            idxEnd = min(numRings, si + segmentLength - 1);
            demyelinMap(si:idxEnd) = false;
        end
    end
    demyelinatedRings = find(~demyelinMap);

    for i = 1:numRings
        p = centerline(:,i);
        if i==1; t = centerline(:,i+1)-p;
        elseif i==numRings; t = p-centerline(:,i-1);
        else; t = centerline(:,i+1)-centerline(:,i-1);
        end
        t = t/norm(t);
        tangents(i,:) = t'; ringCenters(i,:) = p';
        ref = [0;0;1]; if abs(dot(t,ref))>0.9, ref=[1;0;0]; end
        n1 = cross(t,ref); n1 = n1/norm(n1);
        n2 = cross(t,n1);

        % Myelin or demyelinated at this ring
        if demyelinMap(i)
            rSegment = rMyelin;
        else
            rSegment = rAxon;  % No myelin — same as axon
        end

        ringMyelin = rSegment*(n1*circle(1,:)+n2*circle(2,:))+p;
        ringAxon   = rAxon*(n1*circle(1,:)+n2*circle(2,:))+p;

        vertsMyelin((i-1)*circleRes+(1:circleRes),:) = ringMyelin';
        vertsAxon((i-1)*circleRes+(1:circleRes),:)   = ringAxon';
    end

    axonMeshes{end+1} = struct('vertices',vertsMyelin,'faces',faces, 'innerVerts',vertsAxon,'rAxon',rAxon,'rMyelin',rMyelin, 'tangents',tangents,'ringCenters',ringCenters, 'demyelinatedRings', demyelinatedRings);

    boundingSpheres = [boundingSpheres; cx cy cz radiusBound];
    volMyelin = pi*rMyelin^2*pathLen;
    axonFilledVol = axonFilledVol + volMyelin;
    axonCount = axonCount + 1;

    waitbar(axonFilledVol/totalVol, h, sprintf('Axons: %d | Fill: %.1f%%', axonCount,100*axonFilledVol/totalVol));
end
close(h);

%% --- Sphere placement ---------------------------------------------------
sphereCenters = [];
sphereRadii   = [];
attempts = 0;
h = waitbar(0, 'Placing spheres...');

while (sphereFilledVol/totalVol < sphereTargetFill) && (attempts < 10000)
    attempts = attempts + 1;

    % Random focal center and radius
    r = rand*(sphereRadiusRange(2) - sphereRadiusRange(1)) + sphereRadiusRange(1);
    cx = rand*voxelSize;
    cy = rand*voxelSize;
    cz = rand*voxelSize;
    center = [cx, cy, cz];

    % --- Check against all existing bounding spheres (axon + myelin)
    if ~isempty(boundingSpheres)
        delta = boundingSpheres(:,1:3) - center;              % [N x 3] distances
        dist2 = sum(delta.^2, 2);                             % squared distances
        radiiSum = boundingSpheres(:,4) + r;                  % [N x 1]
        if any(dist2 < (radiiSum.^2))                         % overlap?
            continue
        end
    end

    % --- Check against spheres
    if ~isempty(sphereCenters)
        delta = sphereCenters - center;
        dist2 = sum(delta.^2, 2);
        radiiSum = sphereRadii + r;
        if any(dist2 < (radiiSum.^2))
            continue
        end
    end

    % --- Check against cylinders
    if checkSphereAxonCollision(center, r, axonMeshes)
        continue
    end

    % --- Accept this sphere
    sphereCenters(end+1,:) = center;
    sphereRadii(end+1,1)   = r;
    sphereFilledVol = sphereFilledVol + (4/3)*pi*r^3;

    waitbar(sphereFilledVol/totalVol, h, sprintf('Spheres: %d | Fill: %.2f%%', size(sphereCenters,1), 100*sphereFilledVol/totalVol));
end
close(h);

rAxons   = cellfun(@(a) a.rAxon, axonMeshes);
rMyelins = cellfun(@(a) a.rMyelin, axonMeshes);
gRatios  = rAxons ./ rMyelins;

%% --- Geometry Histograms ---------------------------------------------

figure('Name','Distributions of Axon, Myelin, g-Ratio, and Spheres','NumberTitle','off', 'Color', 'w');
subplot(2,2,1);
histogram(rAxons, 100, 'Normalization', 'pdf', 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'none');
hold on;
[f1, xi1] = ksdensity(rAxons);
plot(xi1, f1, 'r-', 'LineWidth', 2);
hMean1 = xline(mean(rAxons), '--k');
hMedian1 = xline(median(rAxons), ':k');
xlabel('Axon Radius (µm)');
ylabel('Probability Density');
title('Axon Radius');
legend([hMean1 hMedian1], {'Mean', 'Median'}, 'Location', 'northeast');
set(gca, 'FontSize', 12);
subplot(2,2,2);
histogram(rMyelins, 100, 'Normalization', 'pdf', 'FaceColor', [0.3 0.7 0.4], 'EdgeColor', 'none');
hold on;
[f2, xi2] = ksdensity(rMyelins);
plot(xi2, f2, 'r-', 'LineWidth', 2);
xline(mean(rMyelins), '--k');
xline(median(rMyelins), ':k');
xlabel('Myelin Radius (µm)');
ylabel('Probability Density');
title('Myelin Radius');
set(gca, 'FontSize', 12);
subplot(2,2,3);
histogram(gRatios, 100, 'Normalization', 'pdf', 'FaceColor', [0.7 0.5 0.8], 'EdgeColor', 'none');
hold on;
[f3, xi3] = ksdensity(gRatios);
f3 = f3 / max(f3);  % Normalize to max of 1
plot(xi3, f3, 'r-', 'LineWidth', 2);
xline(mean(gRatios), '--k');
xline(median(gRatios), ':k');
xlabel('g-Ratio');
ylabel('Normalized Density');
title('g-Ratio Distribution');
xlim([0 1]);
ylim([0 1.05]); % Optional: give slight headroom above 1
set(gca, 'FontSize', 12);
subplot(2,2,4);
histogram(sphereRadii, 50, 'Normalization', 'pdf', 'FaceColor', [0.8 0.4 0.2], 'EdgeColor', 'none');
hold on;
[f4, xi4] = ksdensity(sphereRadii);
plot(xi4, f4, 'r-', 'LineWidth', 2);
xline(mean(sphereRadii), '--k');
xline(median(sphereRadii), ':k');
xlabel('Sphere Radius (µm)');
ylabel('Probability Density');
title('Sphere Radius Distribution');
set(gca, 'FontSize', 12);
sgtitle('Distributions of Axon Radius, Myelin Radius, g-Ratio, and Sphere Radius', 'FontSize', 14);

%% --- Parameter Summary -------------------------------------------------------
axonTotal = 0;
myelinTotal = 0;

for i = 1:length(axonMeshes)
    rAx = axonMeshes{i}.rAxon;
    rMy = axonMeshes{i}.rMyelin;
    h = voxelSize;

    volAxon = pi * rAx^2 * h;
    volMyelin = pi * (rMy^2 - rAx^2) * h;

    axonTotal = axonTotal + volAxon;
    myelinTotal = myelinTotal + volMyelin;
end

paramText = sprintf([
    '=== Parameter Summary ===\n' ...
    'Voxel size (µm)                : %g\n' ...
    'Target axon fill fraction      : %.3f\n' ...
    'Target sphere fill fraction    : %.6f\n' ...
    'G-ratio                       : %.3f ± %.3f (mean ± std)\n' ...
    'Axon radius (µm)               : %.3f ± %.3f (mean ± std)\n' ...
    'Myelin radius (µm)             : %.3f ± %.3f (mean ± std)\n' ...
    'g-Ratio (axon/myelin)          : %.3f ± %.3f (mean ± std)\n' ...
    'Axon radius range (µm)         : [%.3f, %.3f] (min, max)\n' ...
    'Myelin radius range (µm)       : [%.3f, %.3f] (min, max)\n' ...
    'Undulation amplitude (µm)      : %.3f (mean)\n' ...
    'Undulation frequency (1/µm)    : %.5f (mean)\n' ...
    'Number of rings per axon        : %d\n' ...
    'Circle resolution (points/ring): %d\n' ...
    'Segment length (rings)          : %d\n' ...
    'Demyelination probability      : %.3f\n' ...
    'Total axons generated          : %d\n' ...
    'Total axon volume              : %.3e µm^3 (%.2f%%)\n' ...
    'Total myelin volume            : %.3e µm^3 (%.2f%%)\n' ...
    'Total spheres placed           : %d\n' ...
    'Sphere radius (µm)             : %.3f ± %.3f (mean ± std)\n' ...
    'Sphere radius range (µm)       : [%.3f, %.3f] (min, max)\n' ...
    'Total sphere volume            : %.3e µm^3 (%.2f%%)\n' ...
    'Total volume fill fraction     : %.2f%%\n'], ...
    voxelSize, axonTargetFill, sphereTargetFill, ...
    gRatioMean, gRatioStd, ...
    mean(rAxons), std(rAxons), ...
    mean(rMyelins), std(rMyelins), ...
    mean(rAxons./rMyelins), std(rAxons./rMyelins), ...
    min(rAxons), max(rAxons), ...
    min(rMyelins), max(rMyelins), ...
    undAmp, undFreq, ...
    numRings, circleRes, segmentLength, demyelinationProb, ...
    axonCount, ...
    axonTotal, 100*axonTotal/totalVol, ...
    myelinTotal, 100*myelinTotal/totalVol, ...
    length(sphereRadii), ...
    mean(sphereRadii), std(sphereRadii), ...
    min(sphereRadii), max(sphereRadii), ...
    sphereFilledVol, 100*sphereFilledVol/totalVol, ...
    100*(axonTotal + myelinTotal + sphereFilledVol)/totalVol);

figure('Name','Parameter Summary','NumberTitle','off','Color','w');
annotation('textbox',[0 0 1 1],'Interpreter','none', ...
           'FontName','Courier New','FontSize',10, ...
           'EdgeColor','none','HorizontalAlignment','left', ...
           'VerticalAlignment','top','String',paramText);

%% --- Render -------------------------------------------------------
figure('Color', 'w');
hold on;
axonHandle = patch(NaN(2), NaN(2), [0.2 0.2 0.8], 'FaceAlpha', 1, 'EdgeColor', 'none');
myelinHandle = patch(NaN(2), NaN(2), [0.85 0.85 0.85], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
ironHandle = surf(nan(2), nan(2), nan(2), ...
    'FaceColor', [0.9 0.3 0.3], 'EdgeColor', 'none', 'FaceAlpha', 1);
% Cylinders
for i = 1:axonCount
    patch('Vertices', axonMeshes{i}.innerVerts, 'Faces', axonMeshes{i}.faces, ...
        'FaceColor', [0.2 0.2 0.8], 'EdgeColor', 'none', 'FaceAlpha', 1);
    patch('Vertices', axonMeshes{i}.vertices, 'Faces', axonMeshes{i}.faces, ...
        'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end

% Spheres
[sX, sY, sZ] = sphere(16);
for i = 1:size(sphereCenters, 1)
    r = sphereRadii(i); c = sphereCenters(i,:);
    surf(r*sX + c(1), r*sY + c(2), r*sZ + c(3), ...
        'FaceColor', [0.9 0.3 0.3], 'EdgeColor', 'none', 'FaceAlpha', 1);
end
axis equal;
axis vis3d;
xlabel('X (µm)', 'FontSize', 14);
ylabel('Y (µm)', 'FontSize', 14);
zlabel('Z (µm)', 'FontSize', 14);
xlim([-10 voxelSize+10]);
ylim([-10 voxelSize+10]);
zlim([-10 voxelSize+10]);
grid on;
view(3);
lighting flat;
camlight headlight;
set(gca, 'FontSize', 10);

%% Bounding sphere render
figure;
hold on;
axis equal;
xlabel('X (µm)');
ylabel('Y (µm)');
zlabel('Z (µm)');
grid on;
view(3);
camlight headlight;
lighting gouraud;
set(gca, 'FontSize', 14);
sphereHandles = gobjects(size(boundingSpheres,1),1);
lineHandles = gobjects(length(axonMeshes),1);
for i = 1:size(boundingSpheres,1)
    [sx, sy, sz] = sphere(16);
    cx = boundingSpheres(i,1);
    cy = boundingSpheres(i,2);
    cz = boundingSpheres(i,3);
    r = boundingSpheres(i,4);
    sphereHandles(i) = surf(r*sx + cx, r*sy + cy, r*sz + cz, ...
        'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', [0 0.7 0.7]);
end
for i = 1:length(axonMeshes)
    centers = axonMeshes{i}.ringCenters;
    lineHandles(i) = plot3(centers(:,1), centers(:,2), centers(:,3), 'k--', 'LineWidth', 0.5);
end
legend([sphereHandles(1), lineHandles(1)], {'Bounding spheres', 'Cylinder Axes'}, 'Location', 'northeast', 'Fontsize', 14);
