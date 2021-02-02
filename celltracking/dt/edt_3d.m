function d2map = edt_3d(binary_data)%, ref_reg_idx, mov_reg_idx
% euclidean distance transform for 3d data
% ref_reg_idx: is the region as foreground: ref_reg_idx = find(binary_data);
% mov_reg_idx: is the region to get distances

[h, w, zSlice] = size(binary_data);
d2map = nan(h, w, zSlice);
d2map(binary_data>0) = 0;
% transform along z-direction
for x = 1:w
    for y = 1:h
        parabolas_mu = find(binary_data(y,x,:));% z'
        if ~isempty(parabolas_mu)
            target_idx = find(~binary_data(y,x,:));
            d = edt_1d(parabolas_mu, zeros(size(parabolas_mu)), target_idx);

            d2map(y,x, target_idx) = d;
        end
    end
end
% transform along x-direction
target_idx = 1:w;
for y = 1:h
%     if y==8
%         keyboard;
%     end
    for z = 1:zSlice
        parabolas_mu = find(~isnan(d2map(y, :, z)));
        if ~isempty(parabolas_mu)
            d = edt_1d(parabolas_mu, d2map(y, parabolas_mu, z), target_idx);

            d2map(y, :, z) = d';
        end
    end
end

% transform along y-direction
target_idx = 1:h;
for x = 1:w
    for z = 1:zSlice
        parabolas_mu = find(~isnan(d2map(:, x, z))); % impossible to be empty
        d = edt_1d(parabolas_mu, d2map(parabolas_mu, x, z), target_idx);
        
        d2map(:, x, z) = d';
    end
end
%% output
% if narginout > 1
%     distance_map = zeros(size(binary_data));
%     distance_map(mov_reg_idx) = distances;
% end
end