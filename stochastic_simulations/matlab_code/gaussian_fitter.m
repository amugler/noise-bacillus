%{ This code creates a function which is used to create a Gaussian fit of the periodograms generated in the other codes
%}

function z = gaussian_fitter(x)
global xvals yvals


c=zeros(1,numel(yvals));
for j =1:numel(xvals)
% z = z + (yvals(j) - x(1)/sqrt(2*pi*x(3))*exp(-(xvals(j)-x(2))^2/2/x(3)))^2;
c(j)=(yvals(j) - x(1)/sqrt(2*pi*x(3))*exp(-(xvals(j)-x(2))^2/2/x(3)))^2;
end

z=sum(c);
end
