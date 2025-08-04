function collided = checkSphereAxonCollision(center, r, axonMeshes)
collided = false;
% Loop through axons
for i = 1:numel(axonMeshes)
    ax = axonMeshes{i};
    % Loop through segments
    for j = 1:size(ax.ringCenters,1)
        dist = norm(center - ax.ringCenters(j,:));
        % If distance is less than current radius + myelin a collision occurs
        if dist < (r + ax.rMyelin)
            collided = true;
            return;
        end
    end
end
end