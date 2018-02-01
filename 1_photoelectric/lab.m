more off
A=[579.05 369; 579.05 369; 579.05 369.5;
    576.96 377.5; 576.96 376.5; 576.96 376;
    546.06 486.5; 546.06 487.5; 546.06 488; 546.06 488;
    491.6 681; 491.6 686.5; 491.6 685;
    435.84 884; 435.84 886; 435.84 885;
    434.75 888.5;
    433.92 891.5;
    407.78 992; 407.78 985;
    404.66 999.5; 404.66 995; 404.66 995.5;
    ]
A=csvread('calibration.csv',1)
dial_error=6

x=A(:,1)
yerr=dial_error*ones(23,1)
y=A(:,2)

X= [ones(23,1) x]
%%
[a,aerr,b,berr]=my_fit(x, y, yerr)
B=[ a b ]';

xerr=sqrt(yerr.^2/b^2+aerr^2/b^2+berr^2*y.^2/b^4)

subplot(2,2,1);
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

ax=subplot(2,2,2)
hist(ax, Chi)
%hold on
% plot(x, Chi2)
Chi2=sum(Chi2)
%%
data=csvread('curve_1.csv',1)
voltage_error=.01
x=data(:,1);
y=data(:,2);
d=size(x)
xerr=dial_error*ones(d);
yerr=voltage_error*ones(d);

[a,aerr,b,berr]=my_fit(x, y, yerr)
B=[ a b ]';

subplot(2,2,3);
plot(x, y, '.');
hold on
l=errorbar(x,y, yerr);
l.LineStyle='none';
l.Color='b';
l=errorbar(x,y, xerr, 'horizontal');
l.LineStyle='none';
l.Color='b';

fit=linspace(min(x),max(x),10)';
plot(fit, [ones(size(fit)) fit]*B, '-r');
hold off

subplot(2,2,4)
Chi=(y-(a+b.*x))./yerr;
[x xerr y yerr Chi]
hist(Chi)
%hold on
% plot(x, Chi2)
Chi2=sum(Chi.^2)
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