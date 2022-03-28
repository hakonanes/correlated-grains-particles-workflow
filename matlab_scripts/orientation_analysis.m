clear variables, close all

% MTEX configuration
plotx2east
plotzIntoPlane

% export_fig configuration
res = '-r200';

% Crystal and specimen symmetry
cs = {'notIndexed',...
    crystalSymmetry('m-3m', [4.04 4.04 4.04], 'mineral', 'al')};
cs_al = cs{2};
ssO = specimenSymmetry('orthorhombic');

% Directory and file names
dir_data = '/home/hakon/phd/data/p/prover/300c/3';
dir_kp = fullfile(dir_data, 'kp');
dir_mtex = fullfile(dir_data, 'mtex');
fname_ori = 'xmap_refori2.ang';

% Read orientation data
ebsd = EBSD.load(fullfile(dir_kp, fname_ori), cs, 'columnNames', ...
    {'phi1', 'Phi', 'phi2', 'x', 'y', 'ci', 'iq', 'Phase', 'fit',...
    'detector_signal', 'n_particles', 'r_shift', 'c_shift', 'bse'},...
    'radians');
% Align kikuchipy's to MTEX' crystal reference frame
rot_tsl2mtex = rotation.byAxisAngle(xvector - yvector, 180 * degree);
ebsd = rotate(ebsd, rot_tsl2mtex, 'keepXY');

% Remove unnecessary properties
ebsd.prop = rmfield(ebsd.prop, {'c_shift', 'detector_signal',...
    'r_shift', 'ci', 'iq', 'fit'});

% Step size
dx = ebsd.gridify.dx;

% Fix not getting non-indexed points from .ang file properly
ebsd.phaseMap(1) = -1;

% Set phase of particles to 'notIndexed'
ebsd(ebsd.n_particles > 0).phase = -1;

% Misorientation angle threshold (mat) for grain reconstruction
mat = 1 * degree;

% How high is a high angle grain boundary (HAB)?
hab = 15; % Degrees

% Particle classification threshold in microns
dispersoid_threshold = 0.24;
constituent_particle_threshold = 1;

% Orientation color keys
om_al = ipfHSVKey(cs_al);
om_al.CS2 = ssO;

%% Add GND density (computed elsewhere) to properties
% Read GND density
gnd_struct = load(fullfile(dir_mtex, 'gnd.mat'));
gnd = gnd_struct.gnd;

% Read GND density from denoised data
gnd_denoise_struct = load(fullfile(dir_mtex, 'gnd_denoise.mat'));
gnd_denoise = gnd_denoise_struct.gnd;

%% Add property (must be done before grain reconstruction!)
ebsdg = ebsd.gridify;
ebsdg.prop.gnd = gnd;
ebsdg.prop.gnd_denoise = gnd_denoise;
ebsd = ebsdg(~isnan(ebsdg.oldId));

%% Grain reconstruction
ebsd2 = ebsd;

[grains, ebsd2.grainId, ebsd2.mis2mean] = calcGrains(ebsd2, 'angle',...
    mat, 'boundary', 'tight', 'unitCell');

% Assign small Al grains to surrounding grains
ebsd3 = ebsd2(grains(grains.grainSize < 5));
ebsd3 = ebsd3('al');
ebsd2(ismember(ebsd2.id, ebsd3.id)) = [];
[grains2, ebsd2.grainId, ebsd2.mis2mean] = calcGrains(ebsd2, 'angle',...
    mat, 'boundary', 'tight', 'unitCell');

% Add equivalent circular diameter (ECD)
grains2.prop.ecd = 0.816 * 2 * grains2.equivalentRadius;

% Smooth grain boundaries
grains2 = smooth(grains2, 5);

% Remove unnecessary variables
clear grains ebsd3

%% Boundaries and misorientation angles
gb2 = grains2.boundary;
mori = gb2.misorientation.angle ./ degree;
bin_edges = [mat hab max(mori)];
[~, ~, gb2Id] = histcounts(mori, 'NumBins', 2, 'BinEdges', bin_edges);

%% Plot GNDs, boundaries and particles
figure
plot(ebsd, ebsd.gnd, 'micronBar', 'off')
mtexColorMap LaboTeX
caxis([1e12 1e15])
hold on
plot(ebsd('notIndexed'), 'facecolor', 'k')
plot(gb2(gb2Id == 1), 'linecolor', [0.7 0.7 0.7], 'linewidth', 1)
plot(gb2(gb2Id == 2), 'linecolor', [0 0 0], 'linewidth', 1);
legend('hide')
hold off
export_fig(fullfile(dir_mtex, 'maps_gnd_gb.png'), res)

%% Plot GNDs (denoised), boundaries and particles
figure
plot(ebsd, ebsd.gnd_denoise, 'micronBar', 'off')
caxis([1e12 1e15])
mtexColorMap LaboTeX
hold on
plot(ebsd('notIndexed'), 'facecolor', 'k')
plot(gb2(gb2Id == 1), 'linecolor', [0.7 0.7 0.7], 'linewidth', 1)
plot(gb2(gb2Id == 2), 'linecolor', [0 0 0], 'linewidth', 1);
legend('hide')
hold off
export_fig(fullfile(dir_mtex, 'maps_gnd_denoise_gb.png'), res)

%% Orientation maps
directions = {xvector, yvector, zvector};
titles = {'nd', 'rd', 'td'};
for i=1:length(directions)
    om_al.inversePoleFigureDirection = directions{i};
    figure
    plot(ebsd('al'), om_al.orientation2color(ebsd('al').orientations),...
        'micronBar', 'off')
    hold on
    plot(ebsd('notIndexed'), 'facecolor', 'k')
    legend('hide')
    hold off
    export_fig(fullfile(dir_mtex, ['maps_om_ipf_' titles{i} '.png']), res)
end

%% RD orientation map with grain boundaries overlayed
om_al.inversePoleFigureDirection = yvector;
figure
plot(ebsd('al'), om_al.orientation2color(ebsd('al').orientations),...
    'micronBar', 'off')
hold on
plot(ebsd('notIndexed'), 'facecolor', 'k')
plot(gb2(gb2Id == 1), 'linecolor', [0.7 0.7 0.7], 'linewidth', 1)
plot(gb2(gb2Id == 2), 'linecolor', [0 0 0], 'linewidth', 1);
legend('hide')
hold off
export_fig(fullfile(dir_mtex, 'maps_om_ipf_rd_gb.png'), res)

%% Ideal orientations
% Ideal grain texture components, rotated to match the global reference
% frame: X (east) = ND, Y (north) = RD, Z (out of plane) = TD
rot_scan2global = rotation.byMatrix([0 0 1; 1 0 0; 0 1 0]);

br = rot_scan2global * orientation.byMiller([0 1 1], [2 -1 1], cs_al, ssO);
cu = rot_scan2global * orientation.byMiller([1 1 2], [1 1 -1], cs_al, ssO);
s = rot_scan2global * orientation.byMiller([1 2 3], [6 3 4], cs_al, ssO);
cube = rot_scan2global * orientation.byEuler(0, 0, 0, cs_al, ssO);
cubend = rot_scan2global * orientation.byMiller([0 0 1], [3 1 0], cs_al, ssO);
p = rot_scan2global * orientation.byMiller([0 1 1], [-5 -6 6], cs_al, ssO);
goss = rot_scan2global * orientation.byEuler(0, 45*degree, 0, cs_al, ssO);

ideal_oris = {br, cu, s, cube, cubend, p, goss};
ideal_colors = {'m', 'b', 'g', 'r', [1 0.55 0], 'c', 'y'};
ideal_markers = {'d', '^', 'p', 's', 's', '>', 'o'};
ideal_oris_labels = {'br', 'cu',  's', 'cube', 'cubend', 'p', 'goss'};
n_ideal = length(ideal_oris);

%% Inverse pole figure key with components annotated
om_al.inversePoleFigureDirection = yvector; % RD

figure
plot(om_al)
hold on
ms = 20;
for i=1:length(ideal_oris)
    annotate(ideal_oris{i}, 'marker', ideal_markers{i}, 'markersize',...
        ms, 'markerfacecolor', ideal_colors{i});
end
export_fig(fullfile(dir_mtex, 'ipf_annotated.png'), '-r600')

%% Assign a texture component to each grain (without overlap)
mori_threshold_deg = 15;

% Mean orientations of grains
grain_oris = grains2.meanOrientation;
grain_oris.SS = ssO;
n_grains = length(grains2);

% Get misorientation angle between each grain and component
grain_ideal_mangle = zeros(n_ideal, n_grains);
for i=1:n_ideal
    grain_ideal_mangle(i, :) = angle(grain_oris, ideal_oris{i}) / degree;
end

% Get index of minimum misorientation angle
[~, idx] = nanmin(grain_ideal_mangle);

% Assign ideal orientation and ideal orientation ID to each grain
default_ori = orientation.byEuler(1 * degree, 0, 0, cs_al, ssO);
ideal_ori = repmat(default_ori, n_grains, 1);
ideal_ori_id = zeros(n_grains, 1);
for i=1:length(grains2)
    idx_i = idx(i);
    mangle = grain_ideal_mangle(idx_i, i);
    if ~isnan(mangle) && mangle <= mori_threshold_deg
        ideal_ori(i) = ideal_oris{idx_i};
        ideal_ori_id(i) = idx_i;
    end
end
grains2.prop.ideal_ori = ideal_ori;
grains2.prop.ideal_ori_id = ideal_ori_id;

%% Plot of grains with grain boundaries per component
figure
for i=1:n_ideal
    grains_i = grains2(ismember(grains2.ideal_ori_id, i));
    if ~isempty(grains_i)
        plot(grains_i, 'facecolor', ideal_colors{i}, 'micronBar', 'off')
    end
    hold on
end
plot(ebsd('notIndexed'), 'facecolor', 'k')
plot(gb2(gb2Id == 1), 'linecolor', [0.7 0.7 0.7], 'linewidth', 1)
plot(gb2(gb2Id == 2), 'linecolor', [0 0 0], 'linewidth', 1)
legend('hide')
export_fig(fullfile(dir_mtex, 'maps_grains_ideal_particles.png'), res)

%% ---------------------------------- DISPERSOIDS CLOSE TO GRAIN BOUNDARIES
distance_threshold = dx; % um
distance_considered = 1; % um. Sufficiently large, typically ECD dependent
pad = round(distance_considered / dx);

% Grain boundaries of Al-Al not on the map edges
gb2_al = gb2('al', 'al');
gb2_al = gb2_al(~any(gb2_al.grainId == 0, 2));
gb2_idx = 1:size(gb2_al);

% Extract dispersoids to loop over
dispersoid_condition = grains2.phase == -1 & grains2.ecd <=...
    dispersoid_threshold;
grains_particles = grains2(dispersoid_condition);

% Number of dispersoids
n = length(grains_particles);

% Indices of boundaries each particle is within the distance threshold to.
% 300 deemed a sufficiently large number.
boundary_idx = repmat(zeros(1, 300, 'int64'), [n, 1]);

% Keep track of minimum distance to grain boundary for each particle
grains2.prop.min_distance_to_gb = -ones(size(grains2));

h = waitbar(0, 'Finding boundary segments particles are close to');
for i=1:n
    waitbar(i / n)

    particle = grains_particles(i);

    % Extract coordinates of all particle boundary segment
    gb_particle = particle.boundary;
    gb_particle_midpoint = gb_particle.midPoint;
    gb_particle_midpoint_x = gb_particle_midpoint(:, 1);
    gb_particle_midpoint_y = gb_particle_midpoint(:, 2);

    % Calculate extent around particle
    x_min = min(gb_particle_midpoint_x);
    x_max = max(gb_particle_midpoint_x);
    x_extent = round((x_max - x_min) / dx);
    y_min = min(gb_particle_midpoint_y);
    y_max = max(gb_particle_midpoint_y);
    y_extent = round((y_max - y_min) / dx);
    extent = [x_min - distance_considered, y_min - distance_considered,...
        (x_extent + 2 * pad - 1) * dx, (y_extent + 2 * pad - 1) * dx];

    % Extract coordinates of grain boundary segments of interest
    % surrounding particle
    roi = inpolygon(gb2_al, extent);
    gb_of_interest = gb2_al(logical(roi));
    gb_of_interest_midpoint = gb_of_interest.midPoint;
    gb_of_interest_midpoint_x = gb_of_interest_midpoint(:, 1);
    gb_of_interest_midpoint_y = gb_of_interest_midpoint(:, 2);

    % Calculate distance to particle boundary coordinates for each grain
    % boundary coordinate
    distances = sqrt(...
        (gb_particle_midpoint_x - gb_of_interest_midpoint_x').^2 +...
        (gb_particle_midpoint_y - gb_of_interest_midpoint_y').^2);

    % Extract minimum distances for all grain boundary segments
    min_distance = squeeze(min(distances, [], 1))';

    % Extract indices of all boundaries within ROI within distance
    % threshold
    mask = logical(min_distance <= distance_threshold);
    idx_of_interest = gb2_idx(roi);
    boundary_idx(i, 1:sum(mask)) = idx_of_interest(mask);

    % Assign minimum distance to a grain boundary to the particle
    if ~isempty(min_distance)
        grains2(particle.id).min_distance_to_gb = min(min_distance);
    end

end
close(h)

% Trim array of grain boundary indices to non-zero elements by finding
% the index of max non-zero boundary index across all particles
[~, c] = find(boundary_idx);
max_c = max(c);
boundary_idx2 = boundary_idx(:, 1:max_c);

% Fill sub-array for dispersoids into full array for all grains
boundary_idx_all = repmat(zeros(1, max_c), [length(grains2), 1]);
boundary_idx_all(dispersoid_condition, :) = boundary_idx2;

% Assign boundary indices to all grains
grains2.prop.boundary_idx = boundary_idx_all;

% Assign number of dispersoid particles per boundary segment
boundary_idx = grains2.boundary_idx;
boundary_idx_nz = nonzeros(boundary_idx);
[unique_boundary_idx, ~, ic] = unique(boundary_idx_nz);
particle_counts_present = accumarray(ic, 1);
particle_counts_all = zeros(size(gb2));
particle_counts_all(ismember(gb2_idx, unique_boundary_idx)) =...
    particle_counts_present;
gb2.prop.n_particles_close = particle_counts_all;

%% Assign size of dispersoid particles per boundary segment
particles_close_size = zeros(size(gb2_al)); % Al GBs
particle_sizes = grains_particles.ecd; % Dispersoid sizes

for i=1:n
    gb_idx_i = nonzeros(boundary_idx2(i, :));
    if ~isempty(gb_idx_i)
        for j=1:size(gb_idx_i)
            k = gb_idx_i(j);
            particles_close_size(k) = particles_close_size(k) +...
                particle_sizes(i);
        end
    end
end

particles_close_size_all = zeros(size(gb2)); % All GBs
particles_close_size_all(gb2.isIndexed) = particles_close_size;
gb2.prop.particles_close_size = particles_close_size_all;

%% Sanity check of particle locations and boundaries of interest
gb2_al = gb2('al', 'al');
gb3 = gb2_al(ismember(gb2_idx, grains2.boundary_idx));

% Check that grain boundary indices are correct by correlating particle
% locations with boundary locations
figure
plot(grains_particles, 'facecolor', 'r')
hold on
plot(gb3('al', 'al'), 'linewidth', 2)

% Check that number of particles per boundary is correct
figure
plot(grains_particles, 'facecolor', 'r')
hold on
plot(gb2_al, gb2_al.n_particles_close)

% Check assigned minimum distance to grain boundary by ensuring that
% particles on boundaries and particles within boundaries have correct
% hues
figure
plot(grains2('notIndexed'), grains2('notIndexed').min_distance_to_gb)
hold on
plot(gb2_al)

%% Some grain boundary characteristics
% Number of dispersoids per boundary of interest
gb_component1 = zeros([1, length(gb2)]);
gb_component2 = zeros([1, length(gb2)]);
gb2_ids = gb2.grainId;
gb2_ids1 = gb2_ids(:, 1);
gb2_ids2 = gb2_ids(:, 2);
for i=1:n_ideal
    comp_ids = grains2(grains2.ideal_ori == ideal_oris{i}).id;
    gb_component1(ismember(gb2_ids1, comp_ids)) = i;
    gb_component2(ismember(gb2_ids2, comp_ids)) = i;
end
gb2.prop.component1 = gb_component1;
gb2.prop.component2 = gb_component2;

% Special boundary misorientation
gb2.prop.is_csl3 = angle(CSL(3, cs_al), gb2.misorientation) <= 15 * degree;
gb2.prop.is_csl7 = angle(CSL(7, cs_al), gb2.misorientation) <= 15 * degree;

%% Store misorientation of all boundaries and boundaries with particles
gb2_al = gb2('al', 'al');
export(gb2_al.misorientation, fullfile(dir_mtex, 'mori_gb_all.txt'), 'quaternion', 'radians');
mori_particles = gb2_al(gb2_al.n_particles_close > 0).misorientation;
export(mori_particles, fullfile(dir_mtex, 'mori_gb_with_particles.txt'), 'quaternion', 'radians');

%% -------------------------- SPECIAL ORIENTATION RELATIONSHIP: 40 deg<111>
figure
for i=1:n_ideal
    grains_i = grains2(ismember(grains2.ideal_ori_id, i));
    if ~isempty(grains_i)
        plot(grains_i, 'facecolor', ideal_colors{i}, 'facealpha', 0.25)
    end
    hold on
end
plot(ebsd2('notIndexed'), 'facecolor', 'k', 'micronBar', 'off')
plot(gb2, 'linecolor', [0.5, 0.5, 0.5], 'linewidth', 1)
plot(gb2(gb2.is_csl7), 'linecolor', 'r', 'linewidth', 2)
legend('hide')
export_fig(fullfile(dir_mtex, 'grains_ideal_csl7.png'), res)

figure
for i=1:n_ideal
    grains_i = grains2(ismember(grains2.ideal_ori_id, i));
    if ~isempty(grains_i)
        plot(grains_i, 'facecolor', ideal_colors{i}, 'facealpha', 0.25)
    end
    hold on
end
plot(ebsd2('notIndexed'), 'facecolor', 'k', 'micronBar', 'off')
plot(gb2, 'linecolor', [0.5, 0.5, 0.5], 'linewidth', 1)
plot(gb2(gb2.is_csl3), 'linecolor', 'r', 'linewidth', 2)
legend('hide')
export_fig(fullfile(dir_mtex, 'grains_ideal_csl3.png'), res)

%% ----------------- MISORIENTATIONS OF GRAINS AROUND CONSTITUENT PARTICLES
grains_particle = grains2('notIndexed');

% Whether Al grains neighbor a particle
pairs1 = grains_particle.neighbors('full');
pairs1 = unique(pairs1);
grains2.prop.by_particle = (ismember(grains2.id, pairs1)) &...
    (grains2.phase == 1);

% Whether Al grains neighbor a constituent particle
grains_particle_constituent = grains_particle(grains_particle.ecd >=...
    constituent_particle_threshold);
pairs1 = grains_particle_constituent.neighbors('full');
pairs1 = unique(pairs1);
grains2.prop.by_constituent_particle = (ismember(grains2.id, pairs1)) &...
    (grains2.phase == 1);

% Whether grain boundary segments neighbor a constituent particle
gb_is_by_constituent_particle = all(ismember(gb2_ids,...
    grains2(grains2.by_constituent_particle).id), 2);
gb2.prop.by_constituent_particle = gb_is_by_constituent_particle;

%% Misorientation angles of Al grains around constituent particles
gb_constituent = grains2(grains2.by_constituent_particle).boundary('al', 'al');
gb_other = grains2(~grains2.by_constituent_particle).boundary('al', 'al');
ids_between = grains2(grains2.by_constituent_particle).id;
mask2d = ismember(gb_constituent.grainId, ids_between);
mask1d = all(mask2d, 2);
gb_between = gb_constituent(mask1d);
gb_outwards = gb_constituent(~mask1d);

%% Misorientation angles of grains by constituent particles
fid = fopen(fullfile(dir_mtex, 'gb_by_constituent_mori_angles.txt'), 'w+');
fprintf(fid, '%.5f\n', gb_constituent.misorientation.angle');
fclose(fid);

% Misorientation angles between grains by constituent particles
fid = fopen(fullfile(dir_mtex, 'gb_by_constituent_between_mori_angles.txt'), 'w+');
fprintf(fid, '%.5f\n', gb_between.misorientation.angle');
fclose(fid);

% Misorientation angles outward of grains by constituent particles
fid = fopen(fullfile(dir_mtex, 'gb_by_constituent_outward_mori_angles.txt'), 'w+');
fprintf(fid, '%.5f\n', gb_outwards.misorientation.angle');
fclose(fid);

% Misorientation angles of grains NOT by constituent particles
fid = fopen(fullfile(dir_mtex, 'gb_not_by_constituent_mori_angles.txt'), 'w+');
fprintf(fid, '%.5f\n', gb_other.misorientation.angle');
fclose(fid);

%% Write all relevant data to file
fid = fopen(fullfile(dir_mtex, 'grains.txt'), 'w+');
fprintf(fid, ['#id,phase,size,ideal,gos,gam,gnd,gnd_denoise,'...
    'by_particle,by_constituent_particle,dist_to_gb\n']);
dataMat = [...
    grains2.id,...
    grains2.phase,...
    grains2.grainSize,...
    grains2.ideal_ori_id,...
    grains2.GOS,...
    ebsd2.grainMean(ebsd2.KAM),...
    ebsd2.grainMean(ebsd2.gnd),...
    ebsd2.grainMean(ebsd2.gnd_denoise),...
    grains2.by_particle,...
    grains2.by_constituent_particle,...
    grains2.min_distance_to_gb,...
];
fprintf(fid, '%i,%i,%i,%i,%.5f,%.5f,%.5f,%.5f,%i,%i,%.5f\n', dataMat');
fclose(fid);

% Grain neighbors
fid = fopen(fullfile(dir_mtex, 'grain_neighbors.txt'), 'w+');
fprintf(fid, '%i,%i\n', grains2.neighbors');
fclose(fid);

% Grain boundaries
al_ids = grains2('al').id;
gb2.prop.is_al = ismember(gb2.grainId(:, 1), al_ids) .*...
    ismember(gb2.grainId(:, 2), al_ids);
fid = fopen(fullfile(dir_mtex, 'grain_boundaries.txt'), 'w+');
fprintf(fid, ['id1,id2,angle,is_csl3,is_csl7,n_dispersoids_close,'...
    'dispersoids_close_size,by_constituent_particle,component1,'...
    'component2,is_al,length\n']);
dataMat = [...
    gb2.grainId(:, 1),...
    gb2.grainId(:, 2),...
    gb2.misorientation.angle,...
    gb2.is_csl3,...
    gb2.is_csl7,...
    gb2.n_particles_close,...
    gb2.particles_close_size,...
    gb2.by_constituent_particle,...
    gb2.component1',...
    gb2.component2',...
    gb2.is_al,...
    gb2.segLength,...
];
fprintf(fid, '%i,%i,%.10f,%i,%i,%i,%.5f,%i,%i,%i,%i,%10f\n', dataMat');
fclose(fid);

%% --------------------- ANALYSIS OF COMBINED DATA FROM ALL THREE DATA SETS
% Misorientations of boundaries

% All
abcd = [];
for i=1:3
    file_i = fullfile(dir_data(1:end-1), num2str(i), 'mtex/mori_gb_all.txt');
    fid_i = fopen(file_i, 'r');
    abcd_i = textscan(fid_i, '%f', 'headerlines', 1);
    fclose(fid_i);
    abcd = [abcd; abcd_i{1}];
end
mori_all = orientation(reshape(abcd, [4, length(abcd) / 4])', cs_al, cs_al);
mori_all.antipodal = 1;

% With particles
abcd = [];
for i=1:3
    file_i = fullfile(dir_data(1:end-1), num2str(i), 'mtex/mori_gb_with_particles.txt');
    fid_i = fopen(file_i, 'r');
    abcd_i = textscan(fid_i, '%f', 'headerlines', 1);
    fclose(fid_i);
    abcd = [abcd; abcd_i{1}];
end
mori_with_particles = orientation(reshape(abcd, [4, length(abcd) / 4])', cs_al, cs_al);
mori_with_particles.antipodal = 1;

%% Misorientations of all boundaries and boundaries with particles in
% axis-angle space

% Axis-angle space of all boundaries
fz = fundamentalRegion(cs_al, cs_al);
figure
plot(fz)
hold on
plot(mori_all, 'markersize', 1, 'markercolor', 'k', 'points', length(mori_with_particles))

% Axis-angle space of boundaries with particles
figure
plot(fz)
hold on
plot(mori_with_particles, 'markersize', 1, 'markercolor', 'k', 'all')

%% Misorientation axis of all boundaries and boundaries with particles
figure
plot(mori_all.axis, 'fundamentalRegion', 'contourf')
mtexColorbar
caxis([0 2])
nextAxis
plot(mori_with_particles.axis, 'fundamentalRegion', 'contourf')
caxis([0 2])

%% Boundary axes and angles
angle_all = mori_all.angle / degree;
mori_all0_15 = mori_all(angle_all <= 15);
mori_all15_30 = mori_all((angle_all > 15) & (angle_all <= 30));
mori_all30_45 = mori_all((angle_all > 30) & (angle_all <= 45));
mori_all45_60 = mori_all(angle_all > 45);

angle_part = mori_with_particles.angle / degree;
mori_part0_15 = mori_with_particles(angle_part <= 15);
mori_part15_30 = mori_with_particles((angle_part > 15) & (angle_part <= 30));
mori_part30_45 = mori_with_particles((angle_part > 30) & (angle_part <= 45));
mori_part45_60 = mori_with_particles(angle_part > 45);

% Scatter
figure
newMtexFigure('layout', [2, 4]);
plot(mori_all0_15.axis, 'fundamentalRegion', 'markersize', 1, 'markercolor', 'k')
mtexTitle('0--15$^{\circ}$', 'interpreter', 'latex')
nextAxis(1, 2)
plot(mori_all15_30.axis, 'fundamentalRegion', 'markersize', 1, 'markercolor', 'k')
mtexTitle('15--30$^{\circ}$', 'interpreter', 'latex')
nextAxis(1, 3)
plot(mori_all30_45.axis, 'fundamentalRegion', 'markersize', 1, 'markercolor', 'k')
mtexTitle('30--45$^{\circ}$', 'interpreter', 'latex')
nextAxis(1, 4)
plot(mori_all45_60.axis, 'fundamentalRegion', 'markersize', 1, 'markercolor', 'k')
mtexTitle('$>$ 45$^{\circ}$', 'interpreter', 'latex')
nextAxis(2, 1)
plot(mori_part0_15.axis, 'fundamentalRegion', 'markersize', 1, 'markercolor', 'k')
nextAxis(2, 2)
plot(mori_part15_30.axis, 'fundamentalRegion', 'markersize', 1, 'markercolor', 'k')
nextAxis(2, 3)
plot(mori_part30_45.axis, 'fundamentalRegion', 'markersize', 1, 'markercolor', 'k')
nextAxis(2, 4)
plot(mori_part45_60.axis, 'fundamentalRegion', 'markersize', 1, 'markercolor', 'k')
export_fig(fullfile(dir_data(1:end-1), 'misorientation_boundary_axes_angles_scatter.png'), res)

%% Contour
clim = [0 3];
figure
newMtexFigure('layout', [2, 4]);
plot(mori_all0_15.axis, 'fundamentalRegion', 'contourf')
mtexTitle('0--15$^{\circ}$', 'interpreter', 'latex')
caxis(clim)
nextAxis(1, 2)
plot(mori_all15_30.axis, 'fundamentalRegion', 'contourf')
mtexTitle('15--30$^{\circ}$', 'interpreter', 'latex')
caxis(clim)
nextAxis(1, 3)
plot(mori_all30_45.axis, 'fundamentalRegion', 'contourf')
mtexTitle('30--45$^{\circ}$', 'interpreter', 'latex')
caxis(clim)
nextAxis(1, 4)
plot(mori_all45_60.axis, 'fundamentalRegion', 'contourf')
mtexTitle('$>$ 45$^{\circ}$', 'interpreter', 'latex')
caxis(clim)
nextAxis(2, 1)
plot(mori_part0_15.axis, 'fundamentalRegion', 'contourf')
caxis(clim)
nextAxis(2, 2)
plot(mori_part15_30.axis, 'fundamentalRegion', 'contourf')
caxis(clim)
nextAxis(2, 3)
plot(mori_part30_45.axis, 'fundamentalRegion', 'contourf')
caxis(clim)
nextAxis(2, 4)
plot(mori_part45_60.axis, 'fundamentalRegion', 'contourf')
caxis(clim)
mtexColorbar
export_fig(fullfile(dir_data(1:end-1), 'misorientation_boundary_axes_angles_contour.png'), res)