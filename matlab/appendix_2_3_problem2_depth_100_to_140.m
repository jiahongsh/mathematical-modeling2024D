function [d,I,peak,peak_d]=appendix_2_3_problem2_depth_100_to_140
% 附录2.3：第二问，d = 100:1:140 的分段积分计算。
[d,I]=x4;
[peak,i]=max(I);
peak_d=d(i);
