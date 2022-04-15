function [predFunc, opt_h, opt_poly_order] = localPolyRegressionCV( ...
  Xtr, Ytr, h_cands, polyOrder_cands, kernelParams)
% This function implements locally Polynomial Kernel regression and searches for
% the optimal hyper-params (bandwidth and poly order)
% if h_cands is empty, picks 10 values based on the std of X. If polyOrder_cands
% is empty uses poly order 2.
% Picks the hyper parametsr using Kfold CV. Outputs predictions and the optimal
% parameters. For meanings of other variables read localPolyRegression.m

  % prelims
  num_train_pts = size(Xtr, 1);
  num_dims = size(Xtr, 2);
  num_kfoldcv_partitions = min(10, num_train_pts);;

  % Set the kernel
  if ~exist('kernelParams', 'var') | isempty(kernelParams)
    kernelParams = struct;
    kernelParams.kernelType = 'none';
  end

  % specify default values for candidats and poly order if not specified
  if ~exist('h_cands', 'var') | isempty(h_cands)
    silverman_h = 1.06 * median(std(Xtr)) / num_train_pts^( -1/(4 + num_dims) );
    h_cands = logspace(-2, 1, 20)' * silverman_h;
  end
  if ~exist('polyOrder_cands', 'var') | isempty(polyOrder_cands)
    polyOrder_cands = 2; % since Holder(2, L) is a fair assumption ?
  end
  num_po_cands = size(polyOrder_cands, 1);
  num_h_cands = size(h_cands, 1);
  opt_h = h_cands(1);
  opt_poly_order = polyOrder_cands(1);

  % Shuffle the data
  shuffle_order = randperm(num_train_pts);
  Xtr = Xtr(shuffle_order, :);
  Ytr = Ytr(shuffle_order, :);

  % Now iterate through these combinations and obtain the optimal value
  best_cv_error = inf;
  for po_iter = 1:num_po_cands
    for h_iter = 1:num_h_cands
      curr_cv_error = KFoldExperiment(Xtr, Ytr, num_kfoldcv_partitions, ...
                        h_cands(h_iter), polyOrder_cands(po_iter), kernelParams);
      if best_cv_error >= curr_cv_error
        best_cv_error = curr_cv_error;
        opt_h = h_cands(h_iter);
        opt_poly_order = polyOrder_cands(po_iter);
      end
    end
  end

  % Finally use the optimal parameters and all the data to fit a function
  predFunc = @(arg) localPolyRegression(arg, Xtr, Ytr, opt_h, opt_poly_order, ...
                                        kernelParams);

end


function kfold_error = KFoldExperiment(X, y, num_partitions, h, polyOrder, ...
                        kernelParams)
% This function computes the cross validation error for the current candidate
% values for h and polyOrder.

  m = size(X, 1);
  kfold_error = 0;

  for kfold_iter = 1:num_partitions
    test_start_idx = round( (kfold_iter-1)*m/num_partitions + 1 );
    test_end_idx   = round( kfold_iter*m/num_partitions );
    train_indices = [1:test_start_idx-1, test_end_idx+1:m];
    test_indices = [test_start_idx : test_end_idx];
    Xtr = X(train_indices, :);
    ytr = y(train_indices);
    Xte = X(test_indices, :);
    yte = y(test_indices);

    % obtain the predictions
    pred = localPolyRegression(Xte, Xtr, ytr, h, polyOrder, kernelParams);
    % accumulate the errors
    kfold_error = kfold_error + sum( (yte - pred).^2 );
  end

end

