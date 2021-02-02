function d2map = edt_2d(binary_data)%, ref_reg_idx, mov_reg_idx
% euclidean distance transform for 3d data
% ref_reg_idx: is the region as foreground: ref_reg_idx = find(binary_data);
% mov_reg_idx: is the region to get distances


[h, w] = size(binary_data);
d2map = nan(h, w);
d2map(binary_data>0) = 0;
% transform along columns
for x = 1:w
    parabolas_mu = find(binary_data(:,x));% y'
    if ~isempty(parabolas_mu)
        target_idx = find(~binary_data(:,x));
        d = edt_1d(parabolas_mu, zeros(size(parabolas_mu)), target_idx);

        d2map(target_idx, x) = d;
    end
end
% transform along rows
%d2map = nan(h, w);
target_idx = 1:w;
for y=1:h
    parabolas_mu = find(~isnan(d2map(y,:)));
    d = edt_1d(parabolas_mu, d2map(y, parabolas_mu), target_idx);
    
    d2map(y,:) = d';
end


%% output
% if narginout > 1
%     distance_map = zeros(size(binary_data));
%     distance_map(mov_reg_idx) = distances;
% end
end