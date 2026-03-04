function loc = location3d(posx, posy, rfmt, xres, yres)

% rfmt: [row, col, start_slice, step]

if length(rfmt) == 4
    step = rfmt(4);
else
    step = 1;
end
zloc = rfmt(3) + (floor(posy/yres)*rfmt(2) + floor(posx/xres))*step;
xloc = mod(posx, xres);
yloc = mod(posy, yres);

loc = [zloc yloc xloc];