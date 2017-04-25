clc; clear all; close all;

surfpack_load('lift','build_points_lift.dat',2,1);

mach = linspace(1.1,8.0,100);
aoa = linspace(0.0,30.0,100);

[MACH,AOA] = meshgrid(mach,aoa);
LIFT = zeros(size(MACH));


size = size(MACH);
for k = 1:size(1)
    for l = 1:size(2)
       LIFT(k,l) = surfpack_eval('lift',[MACH(k,l),AOA(k,l)]);
    end
end

surf(MACH,AOA,LIFT)
