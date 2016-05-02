clear all;
close all;
clc;


% for i = 1:10
% x = 1:.1:i;
% y = i*x.^2;
% figure(1)
% hold on;
% plot(x,y)
% end
% 
% z = 0:.1:15;
% w = 4.*z;
% 
% figure(1)
% plot(z,w);

s = {' HI','BYE','TOO'};
x = 1:.1:10;
y = x.^2;
t  = 5;
for i=1:10
    if i<5
        t = 6;
    end
hold on;
figure(t)
plot(x,i*x.^2);
end

