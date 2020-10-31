% uses the fact that the prduct can be reduced to: 
%  [(n*(n-2)*(n-4)*...) * ((n-1)*(n-3)*...)^-1] = (n)/(n-1) * (n-2)/(n-3) * (n-4)/(n-5) * ....
%
fprintf('Total sum for R(n=%d) is:    %d \n', 100,top(100));
fprintf('Total sum for R(n=%d) is:    %d \n', 400,top(400));
fprintf('Total sum for R(n=%d) [4 million]is:    %d \n', 4000000 ,top(4000000));

function add = top(x)
    add=1.000000000;
    for i=x:-2:2
    tempT=i/(i-1);
    add=tempT*add;
    end
end

%total=1.000000;
%for i=4000000:-2:2
%    tempT=i/(i-1);
%    total=tempT*total;
%end

