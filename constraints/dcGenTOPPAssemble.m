function dc = dcGenTOPPAssemble(states, splines, dcFun, selectionMatrixL, selectionMatrixR)
% Assembles the sparse jacobian matrix for the nonlinear TOPP constraints

P = reshape(states,numel(states)/5,5)';
nPoints = size(P,2);
sp = reshape(splines,numel(splines)/4,4)';
%dc0 = dcFunGen(P(:,1),sp(:,1));
%dc = zeros(length(states),size(dc0,2)*size(P,2));
dcCell = cell(nPoints,1);

for i = 1:nPoints
    dcCell{i} = dcFun(P(:,i),sp(:,i));
end
dc = blkdiag(dcCell{:});

% rearranging rows and columns to align with state representation
dc = selectionMatrixL*dc*selectionMatrixR;