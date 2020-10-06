function img_svd()
% S. Pollock, math 6406, fall 2019
% simple example of using the SVD for image compression.
% demonstrated with cute animal photos.
% You can use other grayscale images if you prefer.
% You can convert a color image (myimage.jpg) to grayscale by
%tmpCat = imread('myimage.jpg');
%tmpCat = rgb2gray(tmpCat); 

tmpCat = imread('sunnyworking.jpg');  %% input the cat (or other grayscale image)
tmpCat = rgb2gray(tmpCat); 
% try cat.jpg, puppy.jpg, lion.jpg, zebra.jpg (parrot.jpg is big), or your own photo!
Cat = double(tmpCat(:,:,1)); %% tmpCat has 3 layers (rgb), but for grayscale img they are the same.

                             %% convert to double so we can svd
clear tmpCat                 %% we don't need it!

[UC,SC,VC] = svd(Cat);       %% svd the cat
sc = diag(SC);               %% extract diagonal entries of S into a vector

[mcat,ncat] = size(Cat);
ndim = min(mcat,ncat);
kCat = zeros(size(Cat));

done = 0; k = 0;
figure(1), clf
while (k < ndim & done ~= 1)
  k = k+1;
  kCat = kCat + sc(k)*UC(:,k)*VC(:,k)';
  figure(1), imshow(kCat,[0,255]); title(sprintf('k = %g, n = %g',k,ndim))
  figure(2), plot(1:ndim,sc,'b',1,sc(1),'ko',ndim,sc(ndim),'ko',k,sc(k),'ro')
  commandwindow  %% return control to commandwindow
  sdone = input('press 1 to stop, return to continue: ');
  if sdone==1
    done=1;
  end
end
K = k; %% (hold for plot)

reallydone = 0;
while (reallydone ~=1)
  J = input('enter a value of k: ')';
  J = max(min(J,ndim),0);
  JCat = zeros(size(Cat));
  for k = 1:J
    JCat = JCat + sc(k)*UC(:,k)*VC(:,k)';
  end
  figure(3), clf, imshow(JCat,[0,255]); title(sprintf('k = %g, n = %g',k,ndim))
  figure(2), plot(1:ndim,sc,'b',1,sc(1),'bo',ndim,sc(ndim),'bo',K,sc(K),'ro',J,sc(J),'go')

  commandwindow
  yesno = input('Do you want to try another value of K (y/n)? ','s');
  if (yesno=='n')
    reallydone = 1;
  end
end



