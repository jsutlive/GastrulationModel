interval = 32;

file2 = 'Distribution.txt';
B = dlmread(file2);
n = size(B,1);
xy = zeros(n,2);

for i = n:-1:1
    xy(i,:) = B(i,:)-B(1,:); % use coordinate relative to the first point
end
length = xy(n,1) - xy(1,1);

maxval = 0;
for i = 1:1:n
    xy(i,1) = xy(i,1)/length;
    if xy(i,2)<maxval
        maxval = xy(i,2);
    end
end
for i = 1:1:n
    xy(i,2) = xy(i,2)/maxval;
end
plot(xy(:,1),xy(:,2),'r','LineWidth',2);
xlabel('spatial distribution');ylabel('experimental');
set(gca,'FontSize',30,'FontWeight','bold','FontName','Times New Roman');


x = 0 : 1/(interval-1) : 1
xq = interp1(xy(:,1),xy(:,2),x);
figure(2)
plot(x,xq,'r','LineWidth',2);
xlabel('spatial distribution');ylabel('experimental');
set(gca,'FontSize',30,'FontWeight','bold','FontName','Times New Roman');

xy1 = zeros(interval/2,1);
for i = interval/2+1:1:interval
    xy1(i-interval/2) = (xq(i)+xq(interval+1-i))/2;
end

figure(3)
plot(xy1,'r','LineWidth',2);
xlabel('spatial distribution');ylabel('experimental');
set(gca,'FontSize',30,'FontWeight','bold','FontName','Times New Roman');

s = join(['distribution-', int2str(interval), '.txt']);

fileID = fopen(s,'w');
for i = 1:1:size(xy1,1)
    fprintf(fileID,'%8d\n',xy1(i));
end
fclose(fileID);

figure(4)
gau = exp(-0.5*(((x-0.5)*5).^2));
plot(x,gau,'r','LineWidth',2);
xlabel('spatial distribution');ylabel('Gaussian');
set(gca,'FontSize',30,'FontWeight','bold','FontName','Times New Roman');


figure(5)
amp = 1.0+(10-1.0)*x.^3.*(10.0-15.0*x+6.0*x.^2);
%amp = 10.0 - 9*(x-1).*(x-1);
plot(x,amp,'r','LineWidth',2);
ylim([0,10]);
xlabel('time');ylabel('polynomial');
set(gca,'FontSize',30,'FontWeight','bold','FontName','Times New Roman');

figure(6)
amp = 10 - 9*(1-x).^2
plot(x,amp,'r','LineWidth',2);
xlabel('time');ylabel('parabolic');
set(gca,'FontSize',30,'FontWeight','bold','FontName','Times New Roman');

figure(7)
amp = exp(2.3*x);
plot(x,amp,'r','LineWidth',2);
xlabel('time');ylabel('exponential');
set(gca,'FontSize',30,'FontWeight','bold','FontName','Times New Roman');

