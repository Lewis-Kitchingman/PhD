%% Parameters
numRings    = 80;
circleRes   = 24;
radiusAxon  = 1;
gRatio      = 0.65;
radiusMyelin= radiusAxon / gRatio;
undAmp      = 1.5;
undFreq     = 0.1;
segmentLength = 10;
demyelinationProb = 0.2;

z = linspace(0, 10, numRings);

x = undAmp * sin(2*pi*undFreq*z);
y = undAmp * cos(2*pi*undFreq*z);
centerline = [x; y; z];

theta = linspace(0, 2*pi, circleRes+1);
theta(end) = [];
circle = [cos(theta); sin(theta)];

demyelinMap = true(1,numRings);
startIdxs = 1:segmentLength:numRings;
for si = startIdxs
    if rand < demyelinationProb
        idxEnd = min(numRings, si + segmentLength - 1);
        demyelinMap(si:idxEnd) = false;
    end
end

vertsMyelin = zeros(numRings*circleRes, 3);
vertsAxon = zeros(numRings*circleRes, 3);
ringCenters = zeros(numRings, 3);
tangents = zeros(numRings, 3);
ringColors = zeros(numRings, 3);

for i = 1:numRings
    p = centerline(:,i);
    ringCenters(i,:) = p';

    if i == 1
        t = centerline(:,i+1) - p;
    elseif i == numRings
        t = p - centerline(:,i-1);
    else
        t = centerline(:,i+1) - centerline(:,i-1);
    end
    t = t / norm(t);
    tangents(i,:) = t';

    ref = [0;0;1];
    if abs(dot(t,ref)) > 0.9
        ref = [1;0;0];
    end
    n1 = cross(t,ref); n1 = n1/norm(n1);
    n2 = cross(t,n1);

    if demyelinMap(i)
        rOuter = radiusMyelin;
        ringColors(i,:) = [0 0 0.8]; % myelinated
    else
        rOuter = radiusAxon;
        ringColors(i,:) = [0.7 0.7 0.7]; % demyelinated
    end

    ringMyelin = p + rOuter * (n1*circle(1,:) + n2*circle(2,:));
    ringAxon   = p + radiusAxon * (n1*circle(1,:) + n2*circle(2,:));

    vertsMyelin((i-1)*circleRes + (1:circleRes), :) = ringMyelin';
    vertsAxon((i-1)*circleRes + (1:circleRes), :) = ringAxon';
end

faces = [];
for i = 1:(numRings - 1)
    for j = 1:circleRes
        i1 = (i-1)*circleRes + j;
        i2 = (i-1)*circleRes + mod(j,circleRes) + 1;
        i3 = i1 + circleRes;
        i4 = i2 + circleRes;
        faces = [faces; i1 i3 i4; i1 i4 i2];
    end
end

%% Render
figure('Color','w');
hold on;
for i = 1:numRings-1
    vertsIdx = [(i-1)*circleRes + (1:circleRes), i*circleRes + (1:circleRes)];
    faceIdx = [];
    for j = 1:circleRes
        i1 = (j-1) + 1;
        i2 = mod(j,circleRes) + 1;
        faceIdx = [faceIdx; i1 i1+circleRes i2+circleRes; i1 i2+circleRes i2];
    end
    vertsSegment = vertsMyelin(vertsIdx,:);
    c = (ringColors(i,:) + ringColors(i+1,:))/2;
    patch('Vertices', vertsSegment, 'Faces', faceIdx, 'FaceColor', c, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
end
plot3(centerline(1,:), centerline(2,:), centerline(3,:), 'k--', 'LineWidth', 1.5);
hBlue = patch(NaN,NaN,[0 0 0.8],'FaceAlpha',0.4);
hRed = patch(NaN,NaN,[0.7 0.7 0.7],'FaceAlpha',0.4);
legend([hBlue hRed], {'Myelinated segment','Demyelinated segment'}, 'Location', 'northeast');
axis equal;
grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');
view(3);
hold off;
