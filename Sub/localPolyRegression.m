function pred = localPolyRegression(...
    Xte, Xtr, Ytr, h, polyOrder, kernelParams, onlyG)
% Performs locally polynomial kernel regression on the data.
% Xte is a num_test_pts x num_dims matrix at which we need to estimate the
% function. 
% Xtr, Ytr: are the data matrix and regressors.
% h: The bandwidth for kernel regression
% polyOrder: is the order of the polynomial. We just take all powers up to
%   polyOrder for each variable. we do not include the cross terms since the
%   number of terms required grows exponentially.
% The implementation uses either Gaussian kernels or kernels constructed using
% Legendre polynomials.

  % Prelims
  num_train_pts = size(Xtr, 1);
  num_test_pts = size(Xte, 1);
  num_dims = size(Xte, 2);

  % Determine the kernel Function here
  if strcmp(kernelParams.kernelType, 'gauss')
    kernelFunc = @(arg1, arg2) GaussKernel(h, arg1, arg2);
  elseif strcmp(kernelParams.kernelType, 'legendre')
    kernelFunc = @(arg1, arg2) legendreKernel(arg1, arg2, h, kernelParams.order);
  elseif strcmp(kernelParams.kernelType, 'epanechnikov')
    kernelFunc = @(arg1, arg2) num_dims/4*(num_dims+2)/pi^(num_dims/2)*gamma(num_dims/2)*max(0,1-dist2(arg1,arg2));
  elseif strcmp(kernelParams.kernelType, 'none')
    kernelFunc = @(arg1, arg2) true;
  end

  % Construct the Feature matrix B for the train and test points
  Btr = ones(num_train_pts, 1 + polyOrder * num_dims);
  Bte = ones(num_test_pts, 1 + polyOrder * num_dims);
  for p_iter = 1:polyOrder
    start_idx = 2 + (p_iter-1)*num_dims;
    end_idx = 1 + p_iter*num_dims;
    Btr(:, start_idx:end_idx) = Xtr.^p_iter;
    Bte(:, start_idx:end_idx) = Xte.^p_iter;
  end

  % Compute the kernels from the test points to the train points
  % K_tetr = kernelFunc(Xte, Xtr);

  warning off; % matlab complains when h is too small
  % Finally compute the predictions
  pred = zeros(num_test_pts, 1);
  for test_iter = 1:num_test_pts
    % kernels from current test point to all of the train pts.
    omega = kernelFunc(Xtr, Xte(test_iter, :));
    % Compute the predictions. 
    % The parantheses in the first line computes b(x)^T (B' Omega B)^-1
    % The second line computes B' Omega Y
    if nargin > 6 && onlyG
        k = (Btr' * bsxfun(@times, Btr, omega)) \ ( Btr' * (omega .* Ytr) );
        pred(test_iter) = k(2);
        continue;
    end
    pred(test_iter) = ...
      ( Bte(test_iter, :) / (Btr' * bsxfun(@times, Btr, omega)) ) * ...
      ( Btr' * (omega .* Ytr) );
  end
  warning on; 

end

