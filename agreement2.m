function D = agreement2(ci)
%AGREEMENT      Agreement matrix from clusters
%
%   D = AGREEMENT(CI) takes as input a set of vertex partitions CI of
%   dimensions [vertex x partition]. Each column in CI contains the
%   assignments of each vertex to a class/community/module. This function
%   aggregates the partitions in CI into a square [vertex x vertex]
%   agreement matrix D, whose elements indicate the number of times any two
%   vertices were assigned to the same class.
%
%   This version assumes that every one of the ci vectors has exactly 2
%   communities (i.e. ci's are bisections of the full network).
%
%   Inputs,     CI,     set of (possibly) degenerate partitions
%
%   Outputs:    D,      agreement matrix
%
%   Richard Betzel, Indiana University, 2012
%   Olaf Sporns, IU 2020


%n = size(ci,1);

ind = single([(ci==1) (ci==2)]);

D = ind*ind';

%D = D.*~eye(n);
