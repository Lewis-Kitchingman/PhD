function IronB0_total = compute_iron_B0_3D(xg, yg, zg, B0, sphereCenters, sphereRadii, DeltaChiFe)

    IronB0_total = zeros(size(xg));
    numIron = size(sphereCenters, 1);

    h = waitbar(0, 'Computing ΔB0 from spheres...');
    for i = 1:numIron
        center = sphereCenters(i, :);
        radius = sphereRadii(i);

        dx = xg - center(1);
        dy = yg - center(2);
        dz = zg - center(3);

        r = sqrt(dx.^2 + dy.^2 + dz.^2);

        inside = r <= radius;
        outside = r > radius;

        % --- Inside field (uniform)
        dB_inside = (1/3) * B0 * DeltaChiFe;

        % --- Outside field (dipole)
        r_out = r(outside);
        dz_out = dz(outside);

        dB_outside = (B0 * DeltaChiFe * radius^3) ./ (r_out.^3) .* (1 - 3*(dz_out./r_out).^2);

        IronB0_total(inside) = IronB0_total(inside) + dB_inside;
        IronB0_total(outside) = IronB0_total(outside) + dB_outside;

        waitbar(i/numIron, h, sprintf('Sphere %d of %d', i, numIron));
    end
    close(h);
end
