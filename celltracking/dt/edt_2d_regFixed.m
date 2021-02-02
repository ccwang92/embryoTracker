function d2map = edt_2d_regFixed(binary_data, ref_yx, mov_yx)%, ref_reg_idx
% euclidean distance transform for 3d data
% mov_reg_idx: is the region to get distances
[~, od] = sort(ref_yx(:,2),'ascend'); % sort by x in ac
ref_yx = ref_yx(od,:);

[~, od] = sort(mov_yx(:,1),'ascend'); % sort by y in ac
mov_yx = mov_yx(od,:);

[h, w] = size(binary_data);

d2map = nan(h, w);
d2map(sub2ind([h,w], ref_yx(:,1), ref_yx(:,2))) = 0;
% transform along columns
locs = find(ref_yx(2:end,2) - ref_yx(1:end-1,2));
locs_st = [1; locs+1];
locs_end = [locs; size(ref_yx,1)];

target_idx = mov_yx(1,1) : mov_yx(end,1);
for i = 1:length(locs_st)
    x = ref_yx(locs_st(i),2);
    parabolas_mu = sort(ref_yx(locs_st(i):locs_end(i),1));% y'
    
    f = zeros(size(parabolas_mu));
    d = edt_1d(parabolas_mu, f, target_idx);
    
    d2map(target_idx, x) = d;
end
% transform along rows
%d2map = nan(h, w);
target_idx = min(mov_yx(:,2)) : max(mov_yx(:,2));
parabolas_mu = ref_yx(1,2):ref_yx(end,2);
for y=mov_yx(1,1):mov_yx(end,1)
    d = edt_1d(parabolas_mu, d2map(y, parabolas_mu), target_idx);
    d2map(y,target_idx) = d';
end


%% output
% if narginout > 1
%     distance_map = zeros(size(binary_data));
%     distance_map(mov_reg_idx) = distances;
% end
end