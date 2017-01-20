%
% Shows how PCA (Principle Component Analysis) is applied 
% to approximate data in a lower dimensional space
%
% This example generates 3D points using 3 different generated data sets.
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
    sigma(sigma == 0) = 1;
    Xn = (X - mu) ./ sigma;
  end
  
  % Calculate the approximate of Xn given the eigenvalues U using
  % a dimensions-dimensional approximation.
  function [Xapprox] = approximate(Xn, U, dimensions)
    Ureduce = U(:,1:dimensions);
    Z = Xn * Ureduce;
    Xapprox = Z * Ureduce';
  end
  
  % Generate plane with noice
  function [X] = generateRandomPlane(gridSize, noiceFactor)
    s = (0:(gridSize*gridSize-1))';
    xp = floor(s / gridSize);
    yp = mod(s, gridSize);
    noice = (rand(length(s), 3) - .5) * noiceFactor;
    X = [ (xp) (yp) (xp*.1 + yp*.2) ] + noice; 
  end
  
  % Generate curved plane
  function [X] = generateRandomCurvedPlane(gridSize, noiceFactor)
    s = (0:(gridSize*gridSize-1))';
    xp = floor(s / gridSize);
    yp = mod(s, gridSize);
    zp = cos(pi * ((xp-(gridSize/2))/gridSize)) + sin(pi * ((yp-(gridSize/2))/gridSize));
    noice = (rand(length(s), 3) - .5) * noiceFactor;
    X = [ (xp) (yp) (zp) ]; %+ noice; 
  end
  
  % Generate sphere
  function [X] = generateBol(count)
    alpha = rand(count, 1) * 2 * pi;
    beta = rand(count, 1) * 2 * pi;
    gamma = rand(count, 1) * 2 * pi;
    
    x = ones(count, 1) * .5 - rand(count,1) * 0;
    y = zeros(count, 1);
    z = zeros(count, 1);
    
    y1 = cos(alpha) .* y - sin(alpha) .* z;
    z1 = sin(alpha) .* y + cos(alpha) .* z;
    x2 = cos(beta) .* x - sin(beta) .* z1;
    z2 = sin(beta) .* x + cos(beta) .* z1;
    x3 = cos(gamma) .* x2 - sin(gamma) .* y1;
    y3 = sin(gamma) .* x2 + cos(gamma) .* y1;
    
    X = [ (x3) (y3) (z2) ];
  end
  
  function runAnalysis(Xn)
    % Plot X
    figure;
    
    plot3(Xn(:,1), Xn(:,2), Xn(:,3), "ro", "markersize", 5, "linewidth", 1); hold on;
    
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
    
    xlim([-.6 .6]);
    ylim([-.6 .6]);
    zlim([-.6 .6]);
    
    % Show legend
    legend("original", legend1D, legend2D);
    
  end
  
  % Run analysis on 3 different data sets
  % ------------------------------------------------------------------------
  
  [mu sigma data] = normalize(generateRandomPlane(30, 2));
  runAnalysis(data);
  title("3D plane with random noise");
  
  [mu sigma data] = normalize(generateRandomCurvedPlane(30, 2));
  runAnalysis(data);
  title("3D curved plane");
  
  data = generateBol(1000);
  runAnalysis(data);
  title("3D sphere edge");

end
