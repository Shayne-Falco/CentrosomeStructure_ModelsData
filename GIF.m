az = -37.5;
el = 30;
degStep = 2.5;
detlaT = 0.1;
fCount = 71;
f = getframe(gcf);
[im,map] = rgb2ind(f.cdata,256,'nodither');
im(1,1,1,fCount) = 0;
k = 1;
% spin 45°
for i = -37.5:-degStep:-90
  az = i;
  view([az,el])
  f = getframe(gcf);
  im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
  k = k + 1;
end
% tilt down
% for i = 30:-degStep:-15
%   el = i;
%   view([az,el])
%   f = getframe(gcf);
%   im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
%   k = k + 1;
% end
% spin left
% for i = az:-degStep:-90
%   az = i;
%   view([az,el])
%   f = getframe(gcf);
%   im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
%   k = k + 1;
% end
% spin right
for i = -90:degStep:0
  az = i;
  view([az,el])
  f = getframe(gcf);
  im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
  k = k + 1;
end
for i = 0:-degStep:-37.5
  az = i;
  view([az,el])
  f = getframe(gcf);
  im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
  k = k + 1;
end
% tilt up to original
% for i = el:degStep:90
%   el = i;
%   view([az,el])
%   f = getframe(gcf);
%   im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
%   k = k + 1;
% end
imwrite(im,map,'AnimationParams.gif','DelayTime',detlaT,'LoopCount',inf)