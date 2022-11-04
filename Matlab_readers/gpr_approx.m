function [func_approx] = gpr_approx(x,kern_params,gpr_params,arg_dummy)

global gridPointsX gridPointsM nInputs

func_approx = 0.0;

if (arg_dummy >= 0.5) && (arg_dummy < 1.5) % This is the 2D case
    
    for i=1:nInputs
        func_approx = func_approx + gpr_params(i)*se(x,gridPointsX(1:2,i),kern_params);
    end 
    
elseif (arg_dummy < 2.5) && (arg_dummy >= 1.5) % This is the 3D case with x-argument
    for i=1:nInputs
        func_approx = func_approx + gpr_params(i)*se(x,gridPointsX(:,i),kern_params);
    end
else % This is the 3D case with m-argument
    for i=1:nInputs
        func_approx = func_approx + gpr_params(i)*se(x,gridPointsM(:,i),kern_params);
    end
end

end