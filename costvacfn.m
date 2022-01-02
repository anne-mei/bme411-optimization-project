function f = costvacfn(u1,u2,u3,h_bar);

T = length(u1);
costu1 = 0;
costu2 = 0;
costu3 = 0;
for t=1:T
    if u1 (t) > 0 %add fixed cost for each time
        costu1= 3+costu1+ (10*u1(t))^2;
    end
    if u2 (t) > 0 %add fixed cost for each time updated
        costu2= 2+costu2+ (5*u2(t))^2;
    end
    if u3 (t) > 0 %add fixed cost for each time updated
        costu3= 8+costu2+ (20*u3(t))^2;
    end
    f=costu1+costu2+costu3+2000*h_bar;
end

