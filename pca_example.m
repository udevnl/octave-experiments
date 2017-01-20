%
% Shows how PCA (Principle Component Analysis) is applied 
% to approximate data in a lower dimensional space
%
% This example generates 3D points that are approximately on a plane.
% Next PCA is applied to calculate the eigenvalues which form the n-dimensional
% approximates.
%
% Next the 1D and 2D approximates are calculated and plotted.
%
% The 1D approximate is the line that has the lowest projected distance from
% all the original 3D points.
%
% The 2D approximate is the plane that has the lowest projected distance from
% all the original 3D points. 
%
function pca_example()

  function [mu sigma Xn] = normalize(X)
    mu = mean(X);
    sigma = max(X) - min(X);
    Xn = (X - mu) ./ sigma;
  end
  
  % Calculate the approximate of Xn given the eigenvalues U using
  % a dimensions-dimensional approximation.
  function [Xapprox] = approximate(Xn, U, dimensions)
    Ureduce = U(:,1:dimensions);
    Z = Xn * Ureduce;
    Xapprox = Z * Ureduce';
  end
  
  function X = generateRandomPlane(gridSize, noiceFactor)
    s = (0:(gridSize*gridSize-1))';
    xp = floor(s / gridSize);
    yp = mod(s, gridSize);
    noice = (rand(length(s), 3) - .5) * noiceFactor;
    X = [ (xp) (yp) (xp*.1 + yp*.2) ] + noice; 
  end

  % Random 3D data
  X = generateRandomPlane(30, 2);

  % Normalize X
  [mu sigma Xn] = normalize(X);

  % Plot X
  figure;
  plot3(Xn(:,1), Xn(:,2), Xn(:,3), "ro", "markersize", 5, "linewidth", 1); hold on;
  xlim([-.6 .6]);
  ylim([-.6 .6]);
  zlim([-.6 .6]);
  
  % Calculate PCA
  C = cov(Xn);
  [U S V] = svd(C);
  
  % 1D
  Xapprox = approximate(Xn, U, 1);
  vR = 1 - S(1,1) / sum(S(:));
  legend1D = sprintf("1D approximation (Vr=%1.3f)", vR);
  plot3(Xapprox(:,1), Xapprox(:,2), Xapprox(:,3), "bo", "markersize", 5, "linewidth", 1);
  
  % 2D
  Xapprox = approximate(Xn, U, 2);
  vR = 1 - (S(1,1) + S(2,2)) / sum(S(:));
  legend2D = sprintf("2D approximation (Vr=%1.3f)", vR);
  plot3(Xapprox(:,1), Xapprox(:,2), Xapprox(:,3), "go", "markersize", 5, "linewidth", 1);
  
  % Show legend
  legend("original", legend1D, legend2D);

end
