function d2map_zxy = edt_3d_regFixed(ref_cell, mov_cell, mov_shift)
% distance transform for two 3d cells. One is reference cell and the other
% is moving cell.
% ref_cell: 3d binary bounding box of the reference cell(fg cell)
% mov_cell: 3d binary bounding box of the moving cell(voxels where 
% distances needed)
% mov_shift: nx3, yxz. the coordinate difference between ref_cell and the 
% mov_cell in the global coordinate system. (e.g. the coordinate of 
% mov_cell(1,1,1) globally euqals to the coordinate of ref_cell(1,1,1) plus
% mov_shift)

% ccwang@vt.edu, 9/23/2020

[sy1,sx1,sz1] = size(ref_cell);
[sy2,sx2,sz2] = size(mov_cell);

% We choose z->x->y direction. To be more precise, one can select from the 
% minimum values from the follolwing 6 choices:
% zyx: sy1*sx1*sz2 + sy2*sx1*sz2
% zxy: sy1*sx1*sz2 + sy1*sx2*sz2 (our choice)
% xyz: sy1*sx2*sz1 + sy2*sx2*sz1
% xzy: sy1*sx2*sz1 + sy1*sx2*sz2
% yxz: sy2*sx1*sz1 + sy2*sx2*sz1
% yzx: sy2*sx1*sz1 + sy2*sx1*sz2


% transform along z-direction
d2map_z = nan(sy1, sx1, sz2,'single');
target_idx = single((1:sz2) + mov_shift(3));
parabolas_mu_candidates = uint16(1:sz1); % should be dense
for y = 1:sy1
    for x = 1:sx1
        parabolas_mu = parabolas_mu_candidates(ref_cell(y,x,:));% z'
        if ~isempty(parabolas_mu)
            d = edt_1d(parabolas_mu, zeros(size(parabolas_mu),'single'), target_idx);
            %d = edt_1dMex(parabolas_mu, zeros(size(parabolas_mu),'single'), target_idx);
            d2map_z(y, x, :) = d;
        end
    end
end
% transform along x-direction
d2map_zx = nan(sy1, sx2, sz2,'single');
target_idx = single((1:sx2) + mov_shift(2));
parabolas_mu_candidates = uint16(1:sx1); % should be dense
for y = 1:sy1
    for z = 1:sz2
        %parabolas_mu = find(~isnan(d2map_z(y, :, z)));
        parabolas_mu = parabolas_mu_candidates(~isnan(d2map_z(y, :, z))); 
        % parabolas_mu should never be empty, since the bounding boxes are
        % tight; there is only one chance that parabolas_mu is empty, which
        % is that the cell contains two disjoint connected component
        if ~isempty(parabolas_mu)
            d = edt_1d(parabolas_mu, d2map_z(y, parabolas_mu, z), target_idx);
            %d = edt_1dMex(parabolas_mu, d2map_z(y, parabolas_mu, z), target_idx);
            d2map_zx(y, :, z) = d';
        end
    end
end

% transform along y-direction
d2map_zxy = nan(sy2, sx2, sz2,'single');
target_idx = single((1:sy2) + mov_shift(1));
parabolas_mu = uint16(1:sy1); % should be dense with no nan values
valid_ids = find(isnan(d2map_zx(:, 1, 1)));
if ~isempty(valid_ids)
    parabolas_mu(valid_ids) = [];
    d2map_zx(valid_ids, :, :) = [];
end
for x = 1:sx2
    for z = 1:sz2
        d = edt_1d(parabolas_mu, d2map_zx(:, x, z), target_idx);
        %d = edt_1dMex(parabolas_mu, d2map_zx(:, x, z), target_idx);
        d2map_zxy(:, x, z) = d';
    end
end
%% output
% if narginout > 1
%     distance_map = zeros(size(binary_data));
%     distance_map(mov_reg_idx) = distances;
% end
end