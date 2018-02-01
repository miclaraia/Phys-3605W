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

x=A(:,1)
yerr=6*ones(23,1)
y=A(:,2)

X= [ones(23,1) x]
x2 = linspace(400, 600, 10)'
X2 = [ones(10,1) x2]
%%
[a,aerr,b,berr]=my_fit(x, y, yerr)

xerr=sqrt(yerr.^2/b^2+aerr^2/b^2+berr^2*y.^2/b^4)

subplot(2,1,1);
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
plot(x2, X2*B, '-r');
ylabel('Dial Setting');
xlabel('Wavelength');
hold off;
Chi=(y-(a+b.*x))./yerr
Chi2=Chi.^2

ax=subplot(2,1,2)
hist(ax, Chi)
%hold on
% plot(x, Chi2)
Chi2=sum(Chi2)
%%

function [a,aerr,b,berr]=my_fit(x, y, yerr)
    delta = sum((y.^2)./(yerr.^2))*sum(1./(yerr.^2))-(sum(y./yerr.^2))^2
    M=[sum(1./yerr.^2) sum(x./yerr.^2);
       sum(x./yerr.^2) sum(x.^2./yerr.^2)]
    Y=[sum(y./yerr.^2); sum(x.*y./yerr.^2)]

    a = det([Y M(:,2)])/det(M)
    b = det([M(:,1) Y])/det(M)
    B = [a b]'

    aerr=sqrt(sum(x.^2./yerr.^2)/det(M))
    berr=sqrt(sum(1./yerr.^2)/det(M))
end