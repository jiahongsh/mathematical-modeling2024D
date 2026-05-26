function [d,I,peak,peak_d]=appendix_2_1_problem2_depth_152_5_to_180
% 附录2.1：第二问，d = 152.5:1:180 的分段积分计算。
[d,I]=x2;
[peak,i]=max(I);
peak_d=d(i);
