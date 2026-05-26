function [d,I,peak,peak_d]=appendix_2_5_problem2_fine_depth_156_to_160
% 附录2.5：第二问精搜，d = 156:0.01:160 的分段积分计算。
[d,I]=x6;
[peak,i]=max(I);
peak_d=d(i);
