%% Parameters
numRings    = 80;
circleRes   = 24;
radiusAxon  = 1;
gRatio      = 0.65;
radiusMyelin= radiusAxon / gRatio;
undAmp      = 1.5;
undFreq     = 0.1;
pathLen     = 20;

z = linspace(0, pathLen, numRings);
phaseX = 2*pi*rand;
phaseY = 2*pi*rand;
x = undAmp * sin(2*pi*undFreq*z + phaseX);
y = undAmp * cos(2*pi*undFreq*z + phaseY);
centerline = [x; y; z];

cx = mean(x);
cy = mean(y);
cz = mean(z);
radiusBound = radiusMyelin + undAmp;

theta = linspace(0, 2*pi, circleRes+1);
theta(end) = [];
circle = [cos(theta); sin(theta)];

vertsMyelin = zeros(numRings*circleRes, 3);

for i = 1:numRings
    p = centerline(:,i);
    if i == 1
        t = centerline(:,i+1) - p;
    elseif i == numRings
        t = p - centerline(:,i-1);
    else
        t = centerline(:,i+1) - centerline(:,i-1);
    end
    t = t / norm(t);

    ref = [0; 0; 1];
    if abs(dot(t, ref)) > 0.9
        ref = [1; 0; 0];
    end
    n1 = cross(t, ref);
    n1 = n1 / norm(n1);
    n2 = cross(t, n1);

    ringMyelin = p + radiusMyelin * (n1 * circle(1,:) + n2 * circle(2,:));
    vertsMyelin((i-1)*circleRes + (1:circleRes), :) = ringMyelin';
end

faces = [];
for i = 1:(numRings - 1)
    for j = 1:circleRes
        i1 = (i-1)*circleRes + j;
        i2 = (i-1)*circleRes + mod(j, circleRes) + 1;
        i3 = i1 + circleRes;
        i4 = i2 + circleRes;
        faces = [faces; i1 i3 i4; i1 i4 i2];
    end
end

%% Render
figure('Color', 'w'); hold on; axis equal; grid on;
xlabel('X (µm)');
ylabel('Y (µm)');
zlabel('Z (µm)');
view(3);
patch('Vertices', vertsMyelin, 'Faces', faces, 'FaceColor', [0 0 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
plot3(centerline(1,:), centerline(2,:), centerline(3,:), 'k--', 'LineWidth', 1.5);
[Xs, Ys, Zs] = sphere(50);
Xs = radiusBound * Xs + cx;
Ys = radiusBound * Ys + cy;
Zs = radiusBound * Zs + cz;
surf(Xs, Ys, Zs, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'FaceColor', [1 0 0]);
legend({'Myelin sheath', 'Centerline', 'Bounding sphere'}, 'Location', 'best');
