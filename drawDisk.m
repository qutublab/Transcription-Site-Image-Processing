
% In the current figure window draws a circular disk at center (cx, cy) with
% radius r and color rgb.

function drawDisk(cx, cy, r, rgb)
x = cx - r;
y = cy - r;
d = r * 2;
rectangle('Position', [x, y, d, d], 'Curvature', [1, 1], 'FaceColor', rgb);

end
