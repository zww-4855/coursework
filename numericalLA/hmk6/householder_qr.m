function classic_gs()
% S. Pollock, MAD 6406, fall 2020
% Simple, non-optimized implementation of Classic Gram-Schmidt with
% tests for orthgonality and normalization of the columns of Q.

% Experiment 1: random matrix
% User, set m, n the dimensions of the incoming matrix
%m = 332;
%n = 300;
%A = 256*rand(m,n);
%flip = 0; %% don't transpose in final visualization

% Expriment 2: QR the cat vs. dog picture (uncomment to run it)
tmp = imread('jvsp.jpg');
A = double(tmp(:,:,1))'; % this picture has more columns than rows, so we are
%                         % QR-ing its transpose
flip = 1;                % transpose final visualization
[m,n] = size(A);

% initialize an nxn matrix Rc (c for classic) to contain the upper-triangular R1
% initialize an mxn matrix Qc to contain the vectors q_i
Rc = zeros(n,n);
Qc = zeros(m,n);  
% -- classic Gram-Schmidt
for j = 1:n
  v = A(:,j);
  for i = 1:j-1
    q = Qc(:,i);
    r = q'*A(:,j);
    v = v - q*r;
    Rc(i,j) = r;
  end
  r = norm(v);
  Qc(:,j) = (1/r)*v;
  Rc(j,j) = r;
end

orth_err = zeros(n*(n-1)/2,1);
norm_err = zeros(n,1);
% test for orthogonality
k = 0; %% running index
for i = 1:n-1
  norm_err(i) = abs(Qc(:,i)'*Qc(:,i)-1.0);
  for j = i+1:n
    k = k+1;
    orth_err(k) = abs(Qc(:,i)'*Qc(:,j));
  end
end
norm_err(n) = abs(Qc(:,n)'*Qc(:,n)-1.0);

% -- display some results   
fprintf('\nm = %g, n = %g\n------------------------------------\n',m,n) 
fprintf('maximum orthogonality error: %g \n', max(orth_err)) 
fprintf('average orthogonality error: %g \n', sum(orth_err)*2/(n*(n-1)) ) 
fprintf('maximum normalization error: %g \n', max(norm_err)) 
fprintf('average normalization error: %g \n', sum(norm_err)/n ) 

% Adding Householder QR factorization to Gram-Schmidt - based on Algo 10.1
Rc = zeros(n,n);
Qc = zeros(m,n);  
Qold=eye(m);
% -- classic Gram-Schmidt
for k = 1:n 
  x = A(k:m,k);
  size(x,1)
  unitvec=zeros(size(x,1),1);
  unitvec(1)=1
  minus=-norm(x)*unitvec -x;
  plus=norm(x)*unitvec +x;
  if (abs(minus)>abs(plus))
      v=minus;
  else
      v=plus;
  end
  v=v/norm(v);
  A(k:m,k:n)=A(k:m,k:n)-2*v*(v'*A(k:m,k:n));
  %Build eq. 10.2 matrix for Q
  identityMat=eye(k-1);
  F=zeros(m-k+1);
  FI=eye(m-k+1);
  F=FI-2*(v*v')/(v'*v);
  zeroTop=zeros(k-1,m-(k-1));
  zeroBottom=zeros(m-(k-1),k-1);
  BuiltQ= [identityMat, zeroTop ; zeroBottom, F];
  Qnew=Qold*BuiltQ;
  Qold=Qnew;
end
Rc=A;
Qc=Qold;

orth_err = zeros(n*(n-1)/2,1);
norm_err = zeros(n,1);
% test for orthogonality
k = 0; %% running index
for i = 1:n-1
  norm_err(i) = abs(Qc(:,i)'*Qc(:,i)-1.0);
  for j = i+1:n
    k = k+1;
    orth_err(k) = abs(Qc(:,i)'*Qc(:,j));
  end
end
norm_err(n) = abs(Qc(:,n)'*Qc(:,n)-1.0);

% -- display some results   
fprintf('\nm = %g, n = %g\n------------------------------------\n',m,n) 
fprintf('maximum orthogonality error: %g \n', max(orth_err)) 
fprintf('average orthogonality error: %g \n', sum(orth_err)*2/(n*(n-1)) ) 
fprintf('maximum normalization error: %g \n', max(norm_err)) 
fprintf('average normalization error: %g \n', sum(norm_err)/n ) 
%


% -- visualize the results (visualize the transpose if we transposed
% -- the original matrix, and marked flip==1 )
if flip==1
  Arecon = (Qc*Rc)';
else
  Arecon = (Qc*Rc);
end
figure(1), imshow(Arecon,[0,255]); title(sprintf('reconstructed array'))

