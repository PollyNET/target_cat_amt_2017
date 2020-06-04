function [data] = smooth_vert(C,smooth_points_vert)
A=C;
B=ones(smooth_points_vert,1);
N=floor(size(A,2)/smooth_points_vert);
for e=1:N
  data(:,e)=(A(:,(e-1)*smooth_points_vert+1:(e-1)*smooth_points_vert+smooth_points_vert)*B);
end
data=data;
end