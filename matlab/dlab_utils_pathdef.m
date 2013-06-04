function p = dlab_utils_pathdef
% DLAB_UTILS_PATHDEF - adds paths within the repo
%
% In your startup.m (or at the command line), add:
%
% cd ~/git/dlab-utils/matlab % replace with path to folder containing this file
% path(dlab_utils_pathdef,path);
%
% Tom Davidson (tjd@stanford.edu/tjd@alum.mit.edu) 2010-2013
 
path_to_repo = pwd; % since we cd'd here before running this

p = [...
    path_to_repo pathsep...
    ];

