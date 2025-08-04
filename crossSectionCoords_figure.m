% Parameters
numRings   = 8;
circleRes  = 12;
radius     = 1;
height     = 3;
undAmp     = 1.5;
undFreq    = 0.1;

z = linspace(0, height, numRings);
x = undAmp * sin(2*pi*undFreq*z);
y = undAmp * cos(2*pi*undFreq*z);
centerline = [x; y; z];

vertices = zeros(numRings * circleRes, 3);
tangents = zeros(numRings, 3);
ringCenters = centerline';

theta = linspace(0, 2*pi, circleRes+1); theta(end) = [];
circle = [cos(theta); sin(theta)];

for i = 1:numRings
    p = centerline(:, i);
    
    % Tangent vector
    if i == 1
        t = centerline(:, i+1) - p;
    elseif i == numRings
        t = p - centerline(:, i-1);
    else
        t = centerline(:, i+1) - centerline(:, i-1);
    end
    t = t / norm(t);
    tangents(i, :) = t';

    % Orthonormal basis (n1, n2)
    ref = [0; 0; 1];
    if abs(dot(t, ref)) > 0.9
        ref = [1; 0; 0];
    end
    n1 = cross(t, ref); n1 = n1 / norm(n1);
    n2 = cross(t, n1);

    ring = p + radius * (n1 * circle(1,:) + n2 * circle(2,:));
    vertices((i-1)*circleRes + (1:circleRes), :) = ring';
end

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

%% Render
figure('Color','w');
hold on;
patch('Vertices',vertices, 'Faces',faces, 'FaceColor',[0.8 0.9 1], 'EdgeColor','k', 'FaceAlpha', 0.6);
plot3(x, y, z, 'k--', 'LineWidth', 2);
scatter3(ringCenters(:,1), ringCenters(:,2), ringCenters(:,3), 60, 'k', 'filled');
scale = 1.2 * radius;
for i = 1:numRings
    p = ringCenters(i,:)';
    t = tangents(i,:)';
    ref = [0;0;1]; if abs(dot(t,ref)) > 0.9, ref = [1;0;0]; end
    n1 = cross(t, ref); n1 = n1 / norm(n1);
    n2 = cross(t, n1);

    quiver3(p(1), p(2), p(3), t(1), t(2), t(3), scale, 'r', 'LineWidth', 1.5);
    quiver3(p(1), p(2), p(3), n1(1), n1(2), n1(3), scale, 'g', 'LineWidth', 1.5);
    quiver3(p(1), p(2), p(3), n2(1), n2(2), n2(3), scale, 'b', 'LineWidth', 1.5);
end
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
legend({'Mesh Surface', 'Centerline', 'Centers (c_j)', 't_j (red)', 'n_1 (green)', 'n_2 (blue)'});
view(3);
grid on;
hold off;
