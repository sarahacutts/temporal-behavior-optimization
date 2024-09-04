function MIn = nmi2_opt2(Cx, Cy)

% computes the NMI
%
% Stripped down version of partition_distance
% for binary partitions (2 communities), argument 1 can be array, argument
% 2 must be vector; assumes two arguments (both logical)
%
%   2011-2017, Mika Rubinov, UNSW, Janelia HHMI
%   2022, Olaf Sporns

[n,t] = size(Cx);

Px(2,:) = sum(Cx)./n; Px(1,:) = 1-Px(2,:);
HX = - ((Px(1,:).*log(Px(1,:))) + (Px(2,:).*log(Px(2,:))));

% assumes second argument is a single vector
Py(1) = sum(Cy==0)/n; Py(2) = 1-Py(1);
HY = - sum(Py .* log(Py));  

repCy = repmat(Cy,1,t);
% Pxy(1,:) = sum((Cx==0).*(repCy==0))./n;
% Pxy(2,:) = sum((Cx==1).*(repCy==0))./n;
% Pxy(3,:) = sum((Cx==0).*(repCy==1))./n;
% Pxy(4,:) = sum((Cx==1).*(repCy==1))./n;
Pxy(1,:) = sum(~Cx&~repCy)./n;
Pxy(2,:) = sum(Cx&~repCy)./n;
Pxy(3,:) = sum(~Cx&repCy)./n;
Pxy(4,:) = sum(Cx&repCy)./n;
Hxy = - sum(Pxy.*log(Pxy),1);
MIn = 2.*(HX+HY-Hxy)./(HX+HY);

