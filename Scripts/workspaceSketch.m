function [points, hullPoints, radii] = workspaceGeneric(DH, jointTypes, qLimits, resolution, useParfor)
% WORKSPACEGENERIC  Compute the primary workspace of a planar serial-link robot
%
%   [points, hullPoints, radii] = workspaceGeneric(DH, jointTypes, qLimits, resolution, useParfor)
%
%   Inputs:
%     DH           - Nx4 matrix of DH parameters [a, alpha, d, offset]
%     jointTypes   - cell array of length N: 'R' for revolute, 'P' for prismatic
%     qLimits      - Nx2 matrix of joint limits [min max] (rad for R, units for P)
%     resolution   - scalar or Nx1 vector: number of samples per joint (default = 50)
%     useParfor    - (optional) true to enable parallel loop for 3 joints (requires Parallel Toolbox)
%
%   Outputs:
%     points       - all reachable end-effector [x y] positions (M x 2)
%     hullPoints   - vertices of the convex hull outlining the workspace (K x 2)
%     radii        - unique radii of circles approximating hull boundary (L x 1)

% Default parameters
if nargin < 5 || isempty(useParfor)
    useParfor = false;
end
if nargin < 4 || isempty(resolution)
    resolution = 50;
end

% Number of joints
nJ = size(DH,1);

% Ensure resolution vector
if isscalar(resolution)
    res = repmat(resolution, nJ, 1);
else
    res = resolution(:);
end

% Build robot links
L(nJ) = Link();
for i = 1:nJ
    sigma = strcmpi(jointTypes{i}, 'P'); % 1 for prismatic, 0 for revolute
    % Link parameters: [sigma, theta_offset, d_offset, a, alpha]
    L(i) = Link([sigma, DH(i,4), DH(i,3), DH(i,1), DH(i,2)]);
    L(i).qlim = qLimits(i,:);
end
robot = SerialLink(L, 'name', 'generic');

% Prepare joint samples
tSamples = cell(nJ,1);
sz = zeros(nJ,1);
for i = 1:nJ
    tSamples{i} = linspace(qLimits(i,1), qLimits(i,2), res(i));
    sz(i) = numel(tSamples{i});
end

% Preallocate points matrix
total = prod(sz);
points = zeros(total, 2);

% Enumerate configurations and compute forward kinematics
if nJ == 2
    % Fast nested loops for 2 joints
    t1 = tSamples{1}; t2 = tSamples{2};
    idx = 1;
    for i1 = 1:sz(1)
        for i2 = 1:sz(2)
            T = robot.fkine([t1(i1), t2(i2)]);
            p = transl(T);
            points(idx,:) = p(1:2)';
            idx = idx + 1;
        end
    end
elseif nJ == 3 && useParfor
    % Parallel loop for exactly 3 joints
    parfor k = 1:total
        [i1, i2, i3] = ind2sub(sz, k);
        q = [tSamples{1}(i1), tSamples{2}(i2), tSamples{3}(i3)];
        T = robot.fkine(q);
        p = transl(T);
        points(k,:) = p(1:2)';
    end
elseif nJ == 3
    % Standard nested loops for 3 joints
    t1 = tSamples{1}; t2 = tSamples{2}; t3 = tSamples{3};
    idx = 1;
    for i1 = 1:sz(1)
        for i2 = 1:sz(2)
            for i3 = 1:sz(3)
                T = robot.fkine([t1(i1), t2(i2), t3(i3)]);
                p = transl(T);
                points(idx,:) = p(1:2)';
                idx = idx + 1;
            end
        end
    end
else
    % Generic enumeration for any number of joints
    idx = 1;
    for k = 1:total
        subs = cell(1,nJ);
        [subs{:}] = ind2sub(sz, k);
        q = zeros(1, nJ);
        for j = 1:nJ
            q(j) = tSamples{j}(subs{j});
        end
        T = robot.fkine(q);
        p = transl(T);
        points(idx,:) = p(1:2)';
        idx = idx + 1;
    end
end

% Compute convex hull of reachable points
kIdx = convhull(points(:,1), points(:,2));
hullPoints = points(kIdx, :);

% Compute unique radii from base
dists = hypot(hullPoints(:,1), hullPoints(:,2));
radii = unique(round(dists, 6));
end
