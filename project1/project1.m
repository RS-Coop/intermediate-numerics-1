%APPM 4650 Summer 2020 Project 1
clear all;
close all;
format long;

start = 0;
finish = 10;
h = 0.1;
%Pr_numbers = [0.2,0.5,2,5,10];
Pr_numbers = [0.1,0.2,0.5,1,2,5,10];
Pr_numbers_legend = [];

y3init = findy3init(@y1prime,@y2prime,@y3prime,0,0,0,1,0,10,0.1,1e-15);
fprintf('Using the garden hose method, y3(0) = %f \n', y3init);

nu_thermals = [];
figure(1);
hold on;
for Pr = Pr_numbers
    y5init = findy5init(@y1prime,@y2prime,@y3prime,@y4prime,@y5prime,0,0,y3init,1,-1,0.5,Pr,start,finish,h,1e-15);
    fprintf('Using the garden hose method, y5(0) = %f for Pr = %f \n',y5init,Pr);
      
    [y1,y2,y3,y4,y5] = rk4five(@y1prime,@y2prime,@y3prime,@y4prime,@y5prime,0,0,y3init,1,y5init,Pr,start,finish,h);
    plot([start:h:finish],y4);
    
    fprintf('etaM = %f for Pr = %f \n',find(y2 > 0.95,1),Pr);
    nu_thermals = [nu_thermals, find(y4 < 0.01,1)];
    Pr_numbers_legend = [Pr_numbers_legend,compose('Pr = %f',Pr)];
end
title('G(x) vs. x for various Pr#');
xlabel('x');
ylabel('G(x)');
legend(Pr_numbers_legend);
hold off;

figure(2);
hold on;
scatter(Pr_numbers,nu_thermals);
title('Nu Thermal vs. Pr number');
xlabel('Pr#');
ylabel('Nu Thermal');
axis([0 max(Pr_numbers) 0 max(nu_thermals)]);
hold off;

y5init = findy5init(@y1prime,@y2prime,@y3prime,@y4prime,@y5prime,0,0,y3init,1,-1,1,5,start,finish,h,1e-15);
[y1,y2,y3,y4,y5] = rk4five(@y1prime,@y2prime,@y3prime,@y4prime,@y5prime,0,0,y3init,1,y5init,Pr,start,finish,h);
figure(3);
hold on;
plot(start:h:finish,y2);
plot(start:h:finish,(1/2).*((start:h:finish).*y2 - y1));
legend('Fprime(x)','(1/2)*(x*Fprime - F)');
title('Velocity Profiles');
xlabel('x');
hold off;

etaT = find(y2 > 0.95,1);



function y3init = findy3init(y1prime,y2prime,y3prime,y1init,y2init,y3min,y3max,start,finish,h,tol)
    diff = 1; %distance between y2(last) and 1
    increment = (y3max - y3min)/10;
    while diff > tol
        [~,y2atmin,~] = rk4three(y1prime,y2prime,y3prime,y1init,y2init,y3min,start,finish,h);
        y2assmin = y2atmin(end);
        [~,y2atmax,~] = rk4three(y1prime,y2prime,y3prime,y1init,y2init,y3max,start,finish,h);
        y2assmax = y2atmax(end);
        prod = (1-y2assmin)*(1-y2assmax);
        assert(prod < 0);
        
        while (1-y2assmin)*(1-y2assmax) < 0 && y3min < y3max
            [~,y2atmin,~] = rk4three(y1prime,y2prime,y3prime,y1init,y2init,y3min,start,finish,h);
            y2assmin = y2atmin(end);
            y3min = y3min + increment;
        end
        y3min = y3min - 2*increment;
        [~,y2atmin,~] = rk4three(y1prime,y2prime,y3prime,y1init,y2init,y3min,start,finish,h);
        y2assmin = y2atmin(end);
        while (1-y2assmin)*(1-y2assmax) < 0 && y3min < y3max
            [~,y2atmax,~] = rk4three(y1prime,y2prime,y3prime,y1init,y2init,y3max,start,finish,h);
            y2assmax = y2atmax(end);
            y3max = y3max - increment;
        end
        y3max = y3max + 2*increment;
        increment = increment / 10;
        
        diff = abs(1-((y2assmax + y2assmin)/2));
    end
    y3init = (y3max + y3min)/2;
end

function [y1,y2,y3] = rk4three(y1prime,y2prime,y3prime,y1init,y2init,y3init,start,finish,h)
    y1 = y1init;
    y2 = y2init;
    y3 = y3init;
    
    for x = start:h:finish - h
        k11 = h*y1prime(x,y1(end),y2(end),y3(end));
        k12 = h*y2prime(x,y1(end),y2(end),y3(end));
        k13 = h*y3prime(x,y1(end),y2(end),y3(end));
        
        k21 = h*y1prime(x+(h/2),y1(end)+(k11/2),y2(end)+(k12/2),y3(end)+(k13/2));
        k22 = h*y2prime(x+(h/2),y1(end)+(k11/2),y2(end)+(k12/2),y3(end)+(k13/2));
        k23 = h*y3prime(x+(h/2),y1(end)+(k11/2),y2(end)+(k12/2),y3(end)+(k13/2));
        
        k31 = h*y1prime(x+(h/2),y1(end)+(k21/2),y2(end)+(k22/2),y3(end)+(k23/2));
        k32 = h*y2prime(x+(h/2),y1(end)+(k21/2),y2(end)+(k22/2),y3(end)+(k23/2));
        k33 = h*y3prime(x+(h/2),y1(end)+(k21/2),y2(end)+(k22/2),y3(end)+(k23/2));
        
        k41 = h*y1prime(x+h,y1(end)+k31,y2(end)+k32,y3(end)+k33);
        k42 = h*y2prime(x+h,y1(end)+k31,y2(end)+k32,y3(end)+k33);
        k43 = h*y3prime(x+h,y1(end)+k31,y2(end)+k32,y3(end)+k33);
        
        y1 = [y1, y1(end)+(k11+2*k21+2*k31+k41)/6];
        y2 = [y2, y2(end)+(k12+2*k22+2*k32+k42)/6];
        y3 = [y3, y3(end)+(k13+2*k23+2*k33+k43)/6];
    end
end

function y5init = findy5init(y1prime,y2prime,y3prime,y4prime,y5prime,y1init,y2init,y3init,y4init,y5min,y5max,Pr,start,finish,h,tol)
    diff = 1; %distance between y4(last) and 0
    increment = (y5max - y5min)/10;
    while diff > tol
        [~,~,~,y4atmin,~] = rk4five(y1prime,y2prime,y3prime,y4prime,y5prime,y1init,y2init,y3init,y4init,y5min,Pr,start,finish,h);
        y4assmin = y4atmin(end);
        [~,~,~,y4atmax,~] = rk4five(y1prime,y2prime,y3prime,y4prime,y5prime,y1init,y2init,y3init,y4init,y5max,Pr,start,finish,h);
        y4assmax = y4atmax(end);
        prod = y4assmin*y4assmax;
        assert(prod < 0);
        
        while y4assmin*y4assmax < 0 && y5min < y5max
            [~,~,~,y4atmin,~] = rk4five(y1prime,y2prime,y3prime,y4prime,y5prime,y1init,y2init,y3init,y4init,y5min,Pr,start,finish,h);
            y4assmin = y4atmin(end);
            y5min = y5min + increment;
        end
        y5min = y5min - 2*increment;
        [~,~,~,y4atmin,~] = rk4five(y1prime,y2prime,y3prime,y4prime,y5prime,y1init,y2init,y3init,y4init,y5min,Pr,start,finish,h);
        y4assmin = y4atmin(end);
        while y4assmin*y4assmax < 0 && y5min < y5max
            [~,~,~,y4atmax,~] = rk4five(y1prime,y2prime,y3prime,y4prime,y5prime,y1init,y2init,y3init,y4init,y5max,Pr,start,finish,h);
            y4assmax = y4atmax(end);
            y5max = y5max - increment;
        end
        y5max = y5max + 2*increment;
        increment = increment / 10;
        
        diff = (abs(y4assmax) + abs(y4assmin))/2;
    end
    y5init = (y5max + y5min)/2;
end

function [y1,y2,y3,y4,y5] = rk4five(y1prime,y2prime,y3prime,y4prime,y5prime,y1init,y2init,y3init,y4init,y5init,Pr,start,finish,h)
    y1 = y1init;
    y2 = y2init;
    y3 = y3init;
    y4 = y4init;
    y5 = y5init;
    
    for x = start:h:finish - h
        k11 = h*y1prime(x,y1(end),y2(end),y3(end));
        k12 = h*y2prime(x,y1(end),y2(end),y3(end));
        k13 = h*y3prime(x,y1(end),y2(end),y3(end));
        k14 = h*y4prime(x,y1(end),y2(end),y3(end),y4(end),y5(end));
        k15 = h*y5prime(x,y1(end),y2(end),y3(end),y4(end),y5(end),Pr);
        
        k21 = h*y1prime(x+(h/2),y1(end)+(k11/2),y2(end)+(k12/2),y3(end)+(k13/2));
        k22 = h*y2prime(x+(h/2),y1(end)+(k11/2),y2(end)+(k12/2),y3(end)+(k13/2));
        k23 = h*y3prime(x+(h/2),y1(end)+(k11/2),y2(end)+(k12/2),y3(end)+(k13/2));
        k24 = h*y4prime(x+(h/2),y1(end)+(k11/2),y2(end)+(k12/2),y3(end)+(k13/2),y4(end)+(k14/2),y5(end)+(k15/2));
        k25 = h*y5prime(x+(h/2),y1(end)+(k11/2),y2(end)+(k12/2),y3(end)+(k13/2),y4(end)+(k14/2),y5(end)+(k15/2),Pr);
        
        k31 = h*y1prime(x+(h/2),y1(end)+(k21/2),y2(end)+(k22/2),y3(end)+(k23/2));
        k32 = h*y2prime(x+(h/2),y1(end)+(k21/2),y2(end)+(k22/2),y3(end)+(k23/2));
        k33 = h*y3prime(x+(h/2),y1(end)+(k21/2),y2(end)+(k22/2),y3(end)+(k23/2));
        k34 = h*y4prime(x+(h/2),y1(end)+(k21/2),y2(end)+(k22/2),y3(end)+(k23/2),y4(end)+(k24/2),y5(end)+(k25/2));
        k35 = h*y5prime(x+(h/2),y1(end)+(k21/2),y2(end)+(k22/2),y3(end)+(k23/2),y4(end)+(k24/2),y5(end)+(k25/2),Pr);
        
        k41 = h*y1prime(x+h,y1(end)+k31,y2(end)+k32,y3(end)+k33);
        k42 = h*y2prime(x+h,y1(end)+k31,y2(end)+k32,y3(end)+k33);
        k43 = h*y3prime(x+h,y1(end)+k31,y2(end)+k32,y3(end)+k33);
        k44 = h*y4prime(x+h,y1(end)+k31,y2(end)+k32,y3(end)+k33,y4(end)+k34,y5(end)+k35);
        k45 = h*y5prime(x+h,y1(end)+k31,y2(end)+k32,y3(end)+k33,y4(end)+k34,y5(end)+k35,Pr);
        
        y1 = [y1, y1(end)+(k11+2*k21+2*k31+k41)/6];
        y2 = [y2, y2(end)+(k12+2*k22+2*k32+k42)/6];
        y3 = [y3, y3(end)+(k13+2*k23+2*k33+k43)/6];
        y4 = [y4, y4(end)+(k14+2*k24+2*k34+k44)/6];
        y5 = [y5, y5(end)+(k15+2*k25+2*k35+k45)/6];
    end
end

function f = y1prime(x,y1,y2,y3)
    f = y2;
end

function f = y2prime(x,y1,y2,y3)
    f = y3;
end

function f = y3prime(x,y1,y2,y3)
    f = -(y1*y3)/2;
end

function f = y4prime(x,y1,y2,y3,y4,y5)
    f = y5;
end

function f = y5prime(x,y1,y2,y3,y4,y5,Pr)
    f = -(y1*y5*Pr)/2;
end