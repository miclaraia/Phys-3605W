l=3;
w=2;

dial_error=6
voltage_error=.01

A=csvread('calibration.csv',1)

x=A(:,1)
yerr=dial_error*ones(23,1)
y=A(:,2)

X= [ones(23,1) x]
%%
[a,aerr,b,berr]=my_fit(x, y, yerr)
B=[ a b ]';

xerr=sqrt(yerr.^2/b^2+aerr^2/b^2+berr^2*y.^2/b^4)

subplot(l,w,1);
p=plot(x, y, 'b.');
p.MarkerSize=10;
hold on;
line=errorbar(x,y, yerr);
line.LineStyle='none';
line.Color='b';
line=errorbar(x,y, xerr, 'horizontal');
line.Color='b';
line.LineStyle='none';
%errorbar(x, y, yerr, 'b', 'horizontal')
fit=linspace(min(x),max(x),10)';
plot(fit, [ones(size(fit)) fit]*B, '-r');
ylabel('Dial Setting');
xlabel('Wavelength');
hold off;
Chi=(y-(a+b.*x))./yerr
Chi2=Chi.^2

ax=subplot(l,w,2)
hist(ax, Chi)
%hold on
% plot(x, Chi2)
Chi2=sum(Chi2)
%%

 % Curve 1, .5 mm
data=csvread('curve_1.csv',1)
x=data(:,1);
y=data(:,2);
d=size(x)
xerr=dial_error*ones(d);
yerr=voltage_error*ones(d);

[a,aerr,b,berr]=my_fit(x, y, yerr)
B1=[ a b ]';

subplot(l,w,3);
plot(x, y, '.');
hold on
%l=errorbar(x,y, yerr);
%l.LineStyle='none';
%l.Color='b';
%l=errorbar(x,y, xerr, 'horizontal');
%l.LineStyle='none';
%l.Color='b';


data=csvread('curve_2.csv',1)
x=data(:,1);
y=data(:,2);
d=size(x)
xerr=dial_error*ones(d);
yerr=voltage_error*ones(d);

[a,aerr,b,berr]=my_fit(x, y, yerr)
B2=[ a b ]';
plot(x, y, 'm.');

fit=linspace(min(x),max(x),10)';
plot(fit, [ones(size(fit)) fit]*B1, '-r');
plot(fit, [ones(size(fit)) fit]*B2, '-c');
hold off

subplot(l,w,5)
Chi=(y-(a+b.*x))./yerr;
[x xerr y yerr Chi]
hist(Chi)
%hold on
% plot(x, Chi2)
Chi2=sum(Chi.^2)
 %%%%%%%%%%%%%%%%%%%%%%%


data=csvread('curve_3.csv',1)
x=data(:,1);
y=data(:,2);
d=size(x)
xerr=dial_error*ones(d);
yerr=voltage_error*ones(d);

[a,aerr,b,berr]=my_fit(x, y, yerr)
B3=[ a b ]';

subplot(l,w,4)
plot(x, y, 'm.');
hold on

fit=linspace(min(x),max(x),10)';
plot(fit, [ones(size(fit)) fit]*B3, '-c');
hold off

subplot(l,w,6)
Chi=(y-(a+b.*x))./yerr;
[x xerr y yerr Chi]
hist(Chi)
%hold on
% plot(x, Chi2)
Chi2=sum(Chi.^2)


print('second', '-dpng')
%%


function [a,aerr,b,berr]=my_fit(x, y, yerr)
    delta = sum((y.^2)./(yerr.^2))*sum(1./(yerr.^2))-(sum(y./yerr.^2))^2;
    M=[sum(1./yerr.^2) sum(x./yerr.^2);
       sum(x./yerr.^2) sum(x.^2./yerr.^2)];
    Y=[sum(y./yerr.^2); sum(x.*y./yerr.^2)];

    a = det([Y M(:,2)])/det(M);
    b = det([M(:,1) Y])/det(M);
    B = [a b]';

    aerr=sqrt(sum(x.^2./yerr.^2)/det(M));
    berr=sqrt(sum(1./yerr.^2)/det(M));
end
