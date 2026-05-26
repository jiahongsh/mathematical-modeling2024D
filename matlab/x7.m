function [d,I]=x7
sigma = 120; % 标准差
L = 100; % 潜艇长度
R = 20; % 杀伤半径
W = 20; % 潜艇宽度
H = 25; % 高度
sigma_z = 40; % Z的标准差
l1 = 120;
h0 = 150;
a = L + 2 * R;
b = W + 2 * R;

d = 152.5:1:165;

% 计算四个代表性位置的深弹命中概率
I_center = compute_single_depth_charge(0, 0, d, sigma, L, R, W, H, sigma_z, l1, h0);
I_corner = compute_single_depth_charge(a, b, d, sigma, L, R, W, H, sigma_z, l1, h0);
I_updown = compute_single_depth_charge(0, b, d, sigma, L, R, W, H, sigma_z, l1, h0);
I_leftright = compute_single_depth_charge(a, 0, d, sigma, L, R, W, H, sigma_z, l1, h0);

% 联合命中率是各部分之和（由于互斥非重叠性）
I = I_center + 4 * I_corner + 2 * I_updown + 2 * I_leftright;
[peak,i]=max(I);
peak_d=d(i);
end

function I_out = compute_single_depth_charge(x_offset, y_offset, d, sigma, L, R, W, H, sigma_z, l1, h0)
% 基于坐标平移，单深弹在局部坐标系中的积分计算器

f = @(x, y) (1 / (2 * pi * sigma^2)) * exp(-((x + x_offset).^2 + (y + y_offset).^2) / (2 * sigma^2));
Phi = @(x) normcdf(x, 0, 1);
dm=1/(1 - Phi((l1 - h0) / sigma_z));
% 定义函数 g(z)
g_z = @(z) (1/sigma_z)*dm * (1 / sqrt(2 * pi) ) * exp(-((z - h0).^2) / (2 * sigma_z^2));
fun = @(x,y,z) f(x,y).*g_z(z);

I1 = arrayfun(@(di) integral3(@(x, y, z) f(x, y) .* g_z(z), -L/2, L/2, -W/2, W/2, l1, di-R-H/2), d);
I2=[];
I3=[];
I4=[];
I5=arrayfun(@(di) integral(@(z) g_z(z), di-H/2,di+H/2), d);

% 单枚深弹在局部坐标系中的水平覆盖积分概率常数（对应I5）
% 它可以直接利用 integral2 快速算出来！
C_k = integral2(f, -R-L/2, -L/2, @(x)-W/2-sqrt(max(0, R^2 - (x + L/2).^2)), @(x)W/2+sqrt(max(0, R^2 - (x + L/2).^2)));
C_k = C_k + integral2(f, -L/2, L/2, -W/2-R, W/2+R);
C_k = C_k + integral2(f, L/2, R+L/2, @(x)-W/2-sqrt(max(0, R^2 - (x - L/2).^2)), @(x)W/2+sqrt(max(0, R^2 - (x - L/2).^2)));

I5 = C_k * I5;
I6=[];
I7=[];
I8=[];

for i=1:length(d)
    dx=0.5;
    dy=0.5;
    dz=0.5;
    di = d(i);
    
    %%%以下计算I2
    sum2=0;
    zmin = di - R - 0.5 * H;
    zmax = di - 0.5 * H;
    xmin = @(z) -L/2-sqrt(R^2 - (di - z - H/2).^2);
    xmax = @(z) -L/2;
    ymin = @(x,z) -W/2-sqrt(R^2 - (di - z - H/2).^2-(x+L/2).^2);
    ymax = @(x,z) W/2+sqrt(R^2 - (di - z - H/2).^2-(x+L/2).^2);
    for z = zmin:dz:zmax
        xl=xmin(z);
        xu=xmax(z);
        for x =xl:dx:xu
            yl=ymin(x,z);
            yu=ymax(x,z);
            for y =yl:dy:yu
                sum2 = sum2 + fun(x,y,z) * dx * dy * dz;
            end
        end
    end
    I2=[I2 sum2];

    %%%以下计算I3
    sum3=0;
    xmin = -0.5*L;
    xmax = 0.5*L;
    ymin = @(z) -W/2-sqrt(R^2 - (di - z - H/2).^2);
    ymax = @(z) W/2+sqrt(R^2 - (di - z - H/2).^2);
    for z = zmin:dz:zmax
        xl=xmin;
        xu=xmax;
        for x =xl:dx:xu
            yl=ymin(z);
            yu=ymax(z);
            for y =yl:dy:yu
                sum3 = sum3 + fun(x,y,z) * dx * dy * dz;
            end
        end
    end
    I3=[I3 sum3];

    %%%以下计算I4
    sum4=0;
    xmin = @(z) L/2;
    xmax = @(z) L/2+sqrt(R^2 - (di - z - H/2).^2);
    ymin = @(x,z) -W/2-sqrt(R^2 - (di - z - H/2).^2-(x-L/2).^2);
    ymax = @(x,z) W/2+sqrt(R^2 - (di - z - H/2).^2-(x-L/2).^2);
    for z = zmin:dz:zmax
        xl=xmin(z);
        xu=xmax(z);
        for x =xl:dx:xu
            yl=ymin(x,z);
            yu=ymax(x,z);
            for y =yl:dy:yu
                sum4 = sum4 + fun(x,y,z) * dx * dy * dz;
            end
        end
    end
    I4=[I4 sum4];

    %%%以下计算I6
    sum6=0;
    zmin = di + 0.5 * H;
    zmax = di +R+ 0.5 * H;
    xmin = @(z) -L/2-sqrt(R^2 - (di - z + H/2).^2);
    xmax = @(z) -L/2;
    ymin = @(x,z) -W/2-sqrt(R^2 - (di - z + H/2).^2-(x+L/2).^2);
    ymax = @(x,z) W/2+sqrt(R^2 - (di - z + H/2).^2-(x+L/2).^2);
    for z = zmin:dz:zmax
        xl=xmin(z);
        xu=xmax(z);
        for x =xl:dx:xu
            yl=ymin(x,z);
            yu=ymax(x,z);
            for y =yl:dy:yu
                sum6 = sum6 + fun(x,y,z) * dx * dy * dz;
            end
        end
    end
    I6=[I6 sum6];

    %%%以下计算I7
    sum7=0;
    xmin = @(z) -L/2;
    xmax = @(z) L/2;
    ymin = @(x,z) -W/2-sqrt(R^2 - (di - z + H/2).^2);
    ymax = @(x,z) W/2+sqrt(R^2 - (di - z + H/2).^2);
    for z = zmin:dz:zmax
        xl=xmin(z);
        xu=xmax(z);
        for x =xl:dx:xu
            yl=ymin(x,z);
            yu=ymax(x,z);
            for y =yl:dy:yu
                sum7 = sum7 + fun(x,y,z) * dx * dy * dz;
            end
        end
    end
    I7=[I7 sum7];

    %%%以下计算I8
    sum8=0;
    xmin = @(z) L/2;
    xmax = @(z) L/2+R;
    ymin = @(x,z) -W/2-sqrt(R^2 - (di - z + H/2).^2-(x-L/2).^2);
    ymax = @(x,z) W/2+sqrt(R^2 - (di - z + H/2).^2-(x-L/2).^2);
    for z = zmin:dz:zmax
        xl=xmin(z);
        xu=xmax(z);
        for x =xl:dx:xu
            yl=ymin(x,z);
            yu=ymax(x,z);
            for y =yl:dy:yu
                sum8 = sum8 + fun(x,y,z) * dx * dy * dz;
            end
        end
    end
    I8=[I8 sum8];
end

I_out = I1 + I2 + I3 + I4 + I5 + I6 + I7 + I8;
end
