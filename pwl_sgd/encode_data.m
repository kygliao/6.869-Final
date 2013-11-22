function [encoded_data,seg_data] = encode_data(data,param,seg_data)

if(nargin < 2),
    ndivs = 10;
    fprintf('using %i divisions per dimension\n',ndivs);
end
use_pwl = 1;
if(nargin < 3),
    %compute the seg_data 
    mins = min(data,[],2)';
    maxs = max(data,[],2)';
    nsteps = param.ndivs*ones(1,size(data,1));
    step = (maxs-mins)./nsteps;
    if (sum(step<=0)>0)
      fprintf(' small steps?\n');
    end
    seg_data = single(cat(1,mins,maxs,step,nsteps));
    use_pwl = param.use_pwl;

end

%parameters
use_sqrt = 1;              %normalizes the range of each data
iters    = 0;              %set to zero for just encoding (no solver)

%dummy parameters -- used just by the solver
k              = 0;                       % subset size (PEGASOS style)
class_1_weight = 1;                       % relative weight of class 1
lambda         = 0.001;                   % lambda
restart_iter   = 0;
T              = single(ones(1,size(data,2)));

[w,ii,jj,val] = pwl_sgd( data , seg_data, T, k, iters, ...
                                lambda, use_pwl, restart_iter, use_sqrt,class_1_weight);
                              
encoded_data = sparse(ii,jj,val);                              
end

