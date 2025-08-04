% Parameters
numRings = 4;
circleRes = 8;
radius = 1;
height = 3;

z = linspace(0, height, numRings);
theta = linspace(0, 2*pi, circleRes+1); theta(end) = [];

vertices = zeros(numRings * circleRes, 3);
for i = 1:numRings
    for j = 1:circleRes
        idx = (i-1)*circleRes + j;
        x = radius * cos(theta(j));
        y = radius * sin(theta(j));
        vertices(idx,:) = [x, y, z(i)];
    end
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

% Render
figure('Color','w');
hold on;
patch('Vertices',vertices, 'Faces',faces, 'FaceColor',[0.8 0.9 1], 'EdgeColor','k', 'FaceAlpha', 0.7);
scatter3(vertices(:,1), vertices(:,2), vertices(:,3), 40, 'filled', 'k');
for i = 1:size(vertices,1)
    text(vertices(i,1), vertices(i,2), vertices(i,3), sprintf('%d', i), 'FontSize', 10, 'Color', 'r', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
end
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
view(3);
grid on;
