function overlap = checkAxonCollision(centerlineNew, radiusNew, existingMeshes)
% Returns true if the proposed axon collides with any previously placed axons.
% centerlineNew : [3 x N] path of new axon
% radiusNew     : scalar (outer radius including undulation)
% existingMeshes: cell array of structs with .ringCenters and .rMyelin

overlap = false;

for i = 1:length(existingMeshes)
    ringCenters = existingMeshes{i}.ringCenters; % [N x 3]
    radiusExisting = existingMeshes{i}.rMyelin;  % scalar

    % Fast bounding check
    distCenters = pdist2(centerlineNew', ringCenters);
    minDist = min(distCenters(:));

    if minDist < (radiusNew + radiusExisting)
        overlap = true;
        return;
    end
end
end
