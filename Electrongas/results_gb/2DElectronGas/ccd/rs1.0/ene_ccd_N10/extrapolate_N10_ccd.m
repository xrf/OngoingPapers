clear all;
close all;
format long;

fileID = fopen('N10_ccd_alpha0.5.txt','r');
formatSpec = '%d %d %f %f';
sizeM = [4 10];
M = fscanf(fileID,formatSpec,sizeM);
M = M'
fclose(fileID);

x = [1./M(:,2)]
y = [M(:,4)]
n = 2;

[p,S] = polyfit(x,y,n)

hold on;
plot(x,y,'Color','g','LineStyle','-','Marker','o')

x2 = [0:0.00001:0.0005];
%y2 = 0.936311464260502*x2 - 0.198842056852821*x2.^0.0;
y2 = 7.778835165525041*x2.^2 + 0.774991224012795*x2 - 0.198778813173882*x2.^0.0;

plot(x2,y2,'Color','r','LineStyle','--')