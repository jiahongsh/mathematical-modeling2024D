function [d,I]=x5
sigma = 120; % 标准差
L = 100; % 潜艇长度
R = 20; % 杀伤半径
W = 20; % 潜艇宽度
H = 25; % 高度
sigma_z = 40; % Z的标准差
l1 = 120;
h0 = 150;

f = @(x, y) (1 / (2 * pi * sigma^2)) * exp(-(x.^2 + y.^2) / (2 * sigma^2));
Phi = @(x) normcdf(x, 0, 1);
dm=1/(1 - Phi((l1 - h0) / sigma_z));
% 定义函数 g(z)
g_z = @(z) (1/sigma_z)*dm * (1 / sqrt(2 * pi) ) * exp(-((z - h0).^2) / (2 * sigma_z^2));
%test=integral(@(z) g_z(z),120,200);
fun = @(x,y,z) f(x,y).*g_z(z);

d = 87.5:1:100;

%I1 = arrayfun(@(d) integral3(@(x, y, z) f(x, y) .* g_z(z), -L/2, L/2, -W/2, W/2, l1, d-R-H/2), d);
%I2=[];
%I3=[];
%I4=[];
% I5=arrayfun(@(d) integral(@(z) g_z(z), l1,d+H/2), d);
% I5=0.083734*I5;
I6=[];
I7=[];
I8=[];

for i=1:length(d)
    dx=0.5;
    dy=0.5;
    dz=0.5;
    %%%以下计算I6
    % 初始化黎曼和
    sum=0;
    % 计算黎曼和
    zmin = l1;
    zmax = d(i) +R+ 0.5 * H;
    xmin = @(z) -L/2-sqrt(R^2 - (d(i) - z + H/2).^2);
    xmax = @(z) -L/2;
    ymin = @(x,z) -W/2-sqrt(R^2 - (d(i) - z + H/2).^2-(x+L/2).^2);
    ymax = @(x,z) W/2+sqrt(R^2 - (d(i) - z + H/2).^2-(x+L/2).^2);
    for z = zmin:dz:zmax
        xl=xmin(z);
        xu=xmax(z);
        for x =xl:dx:xu
            yl=ymin(x,z);
            yu=ymax(x,z);
            for y =yl:dy:yu
                sum = sum + fun(x,y,z) * dx * dy * dz;
            end
        end
    end
    I6=[I6 sum];

    %%%以下计算I7
    % 初始化黎曼和
    sum=0;
    % 计算黎曼和
    zmin = l1;
    zmax = d(i) +R+ 0.5 * H;
    xmin = @(z) -L/2;
    xmax = @(z) L/2;
    ymin = @(x,z) -W/2-sqrt(R^2 - (d(i) - z + H/2).^2);
    ymax = @(x,z) W/2+sqrt(R^2 - (d(i) - z + H/2).^2);
    for z = zmin:dz:zmax
        xl=xmin(z);
        xu=xmax(z);
        for x =xl:dx:xu
            yl=ymin(x,z);
            yu=ymax(x,z);
            for y =yl:dy:yu
                sum = sum + fun(x,y,z) * dx * dy * dz;
            end
        end
    end
    I7=[I7 sum];

    %%%以下计算I8
    % 初始化黎曼和
    sum=0;
    % 计算黎曼和
    zmin = l1;
    zmax = d(i) +R+ 0.5 * H;
    xmin = @(z) L/2;
    xmax = @(z) L/2+R;
    ymin = @(x,z) -W/2-sqrt(R^2 - (d(i) - z + H/2).^2-(x-L/2).^2);
    ymax = @(x,z) W/2+sqrt(R^2 - (d(i) - z + H/2).^2-(x-L/2).^2);
    for z = zmin:dz:zmax
        xl=xmin(z);
        xu=xmax(z);
        for x =xl:dx:xu
            yl=ymin(x,z);
            yu=ymax(x,z);
            for y =yl:dy:yu
                sum = sum + fun(x,y,z) * dx * dy * dz;
            end
        end
    end
    I8=[I8 sum];
end
I=I6+I7+I8;
[peak,i]=max(I);
peak_d=d(i);
