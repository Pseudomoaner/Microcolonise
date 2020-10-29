function g2 = im_hessstrflt3(init_im, scale)

%Create normalized kernals - x^2, y^2 and xy (horizontal, vertical and diagonal edge detectors)
[x,y,z]=meshgrid(-3:6/scale:3,-3:6/scale:3,-3:6/scale:3);
g2xx=(2*(x.^2)-1).*exp(-(x.^2+y.^2+z.^2));
g2xx=g2xx/sum(abs(g2xx(:)));

g2xy=2*x.*y.*exp(-(x.^2+y.^2+z.^2));
g2xy=g2xy/sum(abs(g2xy(:)));

g2yy=(2*(y.^2)-1).*exp(-(x.^2+y.^2+z.^2));
g2yy=g2yy/sum(abs(g2yy(:)));

g2xz=2*x.*z.*exp(-(x.^2+y.^2+z.^2));
g2xz=g2xz/sum(abs(g2xz(:)));

g2yz=2*y.*z.*exp(-(x.^2+y.^2+z.^2));
g2yz=g2yz/sum(abs(g2yz(:)));

g2zz=(2*(z.^2)-1).*exp(-(x.^2+y.^2+z.^2));
g2zz=g2zz/sum(abs(g2zz(:)));

%Filter image with edge detectors
g2xx_rst=imfilter(init_im, g2xx, 'symmetric', 'same');
g2(:,:,:,1)=g2xx_rst;
g2yy_rst=imfilter(init_im, g2yy, 'symmetric', 'same');
g2(:,:,:,2)=g2yy_rst;
g2zz_rst=imfilter(init_im, g2zz, 'symmetric', 'same');
g2(:,:,:,3)=g2zz_rst;
g2xy_rst=imfilter(init_im, g2xy, 'symmetric', 'same');
g2(:,:,:,4)=g2xy_rst;
g2xz_rst=imfilter(init_im, g2xz, 'symmetric', 'same');
g2(:,:,:,5)=g2xz_rst;
g2yz_rst=imfilter(init_im, g2yz, 'symmetric', 'same');
g2(:,:,:,6)=g2yz_rst;