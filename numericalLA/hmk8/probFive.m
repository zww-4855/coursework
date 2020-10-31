% cosh(x)=(e^x + e^-x)/2
%sinh(x)=(e^x - e^-x)/2


%uses the fact that e^x+e^-x/(1+e^x-e^-x) == (1+e^-2x)/(1+[1/e^x]-e^-2x)
t=top(1000)
b=bottom(1000)

fprintf('numerator and denominator are: %d! and %d\n', t,b);
fprintf('Ratio is: %d \n', t/b);

function add = top(x)
    add=0.0;
    x
    for i=0:x%x:-1:0
        temp=1+exp(-2*i)
        add=add+temp;
    end
    return
end
    
function addD = bottom(x)
    addD=0.0;
    for i=0:x %x:-1:1
        temp=1-exp(-2*i) +(1.0/exp(i))
        addD=addD+temp;
    end
    return
end


