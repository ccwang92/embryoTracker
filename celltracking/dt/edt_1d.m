function d2 = edt_1d(parabolas_mu, parabolas_f, target_idx)
% euclidean distance transform for 1d data
% parabolas_mu: the locations of foreground sorted in ascending order
% parabolas_f: vector of size n, as a function of location from 1 to n
% target_idx: interested places sorted in ascending order; if we are given 
% a specific index set, test only the specific area.

% d2: squared distances
% ref: Felzenszwalb and Huttenlocher, Distance transforms of sampled 
% functions, 2012.

n = length(parabolas_mu);
v = zeros(n,1); % location of parabolas (used to locate a parabolas)
z = zeros(n+1,1); % [z(i),z(i+1)] is the range dominated by the v(i) parabolas

k = 1;
v(1) = 1;
z(1) = -inf;
z(2) = inf;
parabolas_mu = single(parabolas_mu);
for q = 2:n
    cur_mu = parabolas_mu(q);
    s  = ((parabolas_f(q)+cur_mu^2) - ...
        (parabolas_f(v(k))+parabolas_mu(v(k))^2)) / (2*cur_mu-2*parabolas_mu(v(k)));
    while (s <= z(k))
        k = k - 1;
        s  = ((parabolas_f(q)+cur_mu^2) - ...
        (parabolas_f(v(k))+parabolas_mu(v(k))^2)) / (2*cur_mu-2*parabolas_mu(v(k)));
    end
    k = k + 1;
    v(k) = q;
    z(k) = s;
    z(k+1) = inf;
end

d2 = nan(length(target_idx),1); % distances for output
k = 1;
for q = 1:length(target_idx)
    while z(k+1) < target_idx(q)
        k = k + 1;
    end
    d2(q) = (target_idx(q)-parabolas_mu(v(k)))^2 + parabolas_f(v(k));
end

end