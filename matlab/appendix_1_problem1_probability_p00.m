sigma = 120; %标准差，单位：米
L = 100; %潜艇长度，单位：米
R = 20; %杀伤半径，单位：米
W = 20; %潜艇宽度，单位：米

%定义被积函数
f = @(x, y) exp(-(x.^2 + y.^2) / (2 * sigma^2));

%计算第一个积分项（左侧部分）
I1 = integral2(f, -R-L/2, -L/2, @(x)-W/2-sqrt(R^2 - (x + L/2).^2), @(x)W/2+sqrt(R^2 - (x + L/2).^2));

%计算第二个积分项（中间部分）
I2 = integral2(f, -L/2, L/2, -W/2-R, W/2+R);

%计算第三个积分项（右侧部分）
I3 = integral2(f, L/2, R+L/2, @(x)-W/2-sqrt(R^2 - (x - L/2).^2), @(x)W/2+sqrt(R^2 - (x - L/2).^2));

%计算总积分
poe = (1 / (2 * pi * sigma^2)) * (I1 + I2 + I3);

%显示结果
disp(['The probability p(0,0) is: ', num2str(poe)]);

%保存结果，供Python测试对齐
save('matlab/appendix_1_result.mat', 'poe');
