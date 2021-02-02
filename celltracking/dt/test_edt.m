profile on;
for i=1:1000
    a = rand(1,65534)>0.5;
    kk = find(a);
    d1 = edt_1dMex(uint16(kk), zeros(size(kk),'single'), single(1:length(a)));
    d1 = sqrt(d1);
    
    d2 = bwdist(a);
    
    if max(abs(d1-d2))>0.01
        disp(i);
    end
end
profile viewer;

profile off;


a(:,:,4) = 0;

c1 = edt_3d_regFixed(b, a, [-1 0.5 0.3]);
c2 = edt_3dMex(b, a, single([-1 0.5 0.3]));
disp(max(abs(c1(:)-c2')));

c1 = edt_3d_regFixed(a, b, [-1 0.5 0.3]);
c2 = edt_3dMex(a, b, single([-1 0.5 0.3]));
disp(max(abs(c1(:)-c2')));


c1 = edt_3dMex(ref_cell, mov_cell, mov_shift);
c2 = edt_3d_regFixed(ref_cell, mov_cell, mov_shift);
c1 = reshape(c1, size(c2));
disp(max(abs(c1(:)-c2(:))));

c1 = edt_3dMex(mov_cell, ref_cell, -mov_shift);
c2 = edt_3d_regFixed(mov_cell, ref_cell, -mov_shift);
c1 = reshape(c1, size(c2));
disp(max(abs(c1(:)-c2(:))));


c1 = edt_2dMex(mov_cell, ref_cell, single([0 0]));
c2 = edt_3d_regFixed(mov_cell, ref_cell, [0 0 0]);
c1 = reshape(c1, size(c2));
disp(max(abs(c1(:)-c2(:))));
