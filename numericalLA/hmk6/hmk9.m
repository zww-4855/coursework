%part a
x=(-128:128)'/128;
A=[x.^0 x.^1 x.^2 x.^3]
[Q,R]=qr(A,0);
scale=Q(257,:);
Q=Q*diag(1./scale);
%plot(Q)
%xlabel('Grid Point') 
%ylabel('Approximated P_j (x)')
%legend('P_0','P_1','P_2','P_3')

%part b
m=length(x);
% Definiton of legendre polynominal taken from book
% eq 7.11
realLegPoly=[ones(m,1) x (1.5*x.^2-0.5) (2.5*x.^3-1.5*x)];
error=realLegPoly - Q
%plot(error)
%xlabel('Grid Point') 
%ylabel('Error of Approximated P_j (x) w.r.t Real P_j')
%legend('P_0 Error','P_1 Error','P_2 Error','P_3 Error')
e1=max(error(:,1))
e2=max(error(:,2))
e3=max(error(:,3))
e4=max(error(:,4))
fprintf('Max. error for P_0, P_1, P_2, and P_3, respectively: %d %d %d %d \n',e1,e2,e3,e4)



%part c
v1=4;
v2=25;

error=[];
count=1
pzero=[];
pone=[];
ptwo=[];
pthree=[];
for v = v1:v2
    x=(-2^v:2^v)'/2^v;
    A=[x.^0 x.^1 x.^2 x.^3];
    [Q,R]=qr(A,0);
    scale=Q(2^(v+1)+1,:);
    Q=Q*diag(1./scale);
    
    m=length(x);
    % Definiton of legendre polynominal taken from book
    % eq 7.11
    realLegPoly=[ones(m,1) x (1.5*x.^2-0.5) (2.5*x.^3-1.5*x)];
    error(count)=[max(max(abs(realLegPoly-Q)))];
    
    pzero(count)=max(abs(realLegPoly(:,1)-Q(:,1)));
    pone(count)=max(abs(realLegPoly(:,2)-Q(:,2)));
    ptwo(count)=max(abs(realLegPoly(:,3)-Q(:,3)));
    pthree(count)=max(abs(realLegPoly(:,4)-Q(:,4)));
    
    fprintf('iteration: %i and max. error %d \n',v,error(count))
    count=count+1;
end
final= [pzero' pone' ptwo' pthree']
%scatter(v1:v2,log(error));
plot(log(final));
xlabel('v in Delta(x) = 1/2^v') 
ylabel('Log of Max Error b.t. the Legendre poly Cols.')
legend('P_0 Error','P_1 Error','P_2 Error','P_3 Error')
