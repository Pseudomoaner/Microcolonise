function A = im_hessangle3(im, scale)

%Run ridge detection filtering stages - g2 corresponds to the values of the Hessian matrix for each point in the greyscale intensity matrix.
g2 = im_hessstrflt3(im, scale);

%Construct each Hessian for each pixel in image, and calculate associated
%eigenvalues. Assign largest as A.
A = zeros(size(im));
for i = 1:size(im,1)
    for j = 1:size(im,2)
        for k = 1:size(im,3)
            hessMat = [g2(i,j,k,1),g2(i,j,k,4),g2(i,j,k,5);g2(i,j,k,4),g2(i,j,k,2),g2(i,j,k,6);g2(i,j,k,5),g2(i,j,k,6),g2(i,j,k,3)];
            E = eig(hessMat);
            [~,maxInd] = max(abs(E));
            A(i,j,k) = E(maxInd);
        end
    end
end