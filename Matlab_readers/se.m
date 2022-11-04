function [kern_result] = se(x1,x2,k_params)

 s = k_params(1); % Signal strength parameter
 h = log(s);
 l = k_params(2); % Scale length parameters

% Scale each dimension using the scale-length l to compute distance between x and xprime
scaled_d = sqrt( sum( ((x1 - x2)/l).^2.0) );

% Now apply the SE formula
kern_result = exp(-.5*scaled_d + h);

end 