function data = smooth_temp(C,smooth_points_temp)
%smooth_points_temp=3
%A=[1,2,3,4,5,11,14,6,7,8,9,10,12,15]

%A=size(C');
A=C';
B=ones(smooth_points_temp,1);

N=floor(size(A,2)/smooth_points_temp);
for e=1:N
  data(:,e)=(A(:,(e-1)*smooth_points_temp+1:(e-1)*smooth_points_temp+smooth_points_temp)*B);
  %if e<2
  %  f=d ;
  %else
  %  f=[f,d];
  %end
end

data=data';
end
