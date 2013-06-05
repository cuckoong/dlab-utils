function [R, B] = target3d(target, armconfig, tilt, rot, dvintersect, quiet)
% TARGET3D: target a location using a tilted, rotated stereotax
%
% Specify the location of your target, and the angle of approach you would
% like to take to that target. This code generates the translation along
% each of the 3 (possibly non-orthogonal) stereotactic axes needed to reach
% the target.
%
% There are many different conventions in use to describe an angle in 3
% dimensions. Here, we specify the approach angle in terms of its tilt away
% from vertical, and its rotation about the same vertical axis. (See 'tilt'
% and 'rot' below for specifics).
%
% Note that many possible approach angles are awkward or infeasible using a
% given stereotaxic apparatus and probe, so it is recommended that you
% first 'rough out' the tilt and rotation that you would like to use to
% reach your target, then use this code to generate precise coordinates.
% (Protip: you can use a second, untilted and unrotated arm to specify your
% target location in space. See:
%  A. J. Greenshaw, J. Neurosci. Methods 78, 169?172 (1997))
%
% N.B. All coordinates are given as [ML AP DV], relative to Bregma, where:
%          ML (mediolateral): animal right is +ve
%          AP (anteroposterior): anterior is +ve
%          DV (dorsoventral): dorsal (up) is +ve
% (Same order and sign as used by Kopf digital display)
%
% INPUTS:
%  target = [ML AP DV] atlas coordinates of target location.
%
%  armconfig = 'LR', or 'FB': is the stereotax arm 'knuckle' configured so
%        that, when a probe is located over bregma, the arm tilts 
%        left-right, or front-back?
%
%  tilt = degrees of tilt away from vertical. Must be positive, and less 
%        than 90 degrees! (direction of tilt is specified by rotation)
%
%  rot = degrees of CW rotation of tilted DV axis away from animal nose.
%        (I.e. 0 means arm tilts towards animals nose, 90 means arm tilts
%        to the right, -90 (or 270) means arm tilts to the left.
%
% OPTIONAL INPUTS:
%  dvintersect = DV heights (in old coordinates) at which to calculate
%        intersection of approach path. Default = []. Useful for 
%        establishing craniotomy locations in unrotated coordinates before 
%        rotating the stereotax arms.
%
%  quiet = true/false: silence text description of outputs, just return
%        coordinates. Default = false.
%
% OUTPUTS:
%  R = coordinates to use to reach the targeted location, using the 
%        rotated/tilted stereotax arm.
%
%  B = location, in *unrotated* coordinates, of the intersection of the 
%        approach path with a DV plane specified by dvintersect. 
%
% EXAMPLES:
%  To target a location at ML 0.25, AP 0.7, DV -6.6, with a stereotax arm 
%  tilted 20 degrees to the right:
%     R = target3d([0.25 0.7 -6.6], 'LR', 20, 90);
%
%     -->R: +2.49 / +0.70 / -6.12 
%
%  To target a location at ML -2, AP -2, DV -3.1, with a stereotax arm 
%  tilted 30 degrees forwards(i.e. leaning towards the animal's nose).
%     R = target3d([-2 -2 -3.1], 'FB', 30, 0);
%
%     -->R: -2.00 / -0.21 / -3.58
%
%  Same target, also calculate the location of the intersection of the
%  approach path with a plane slightly below bregma (e.g. for marking a 
%  craniotomy site before tilting the arm:
%     [R B] = target3d([-2 -2 -3.1], 'FB', 30, 0, -0.3);
%
%     -->R: -2.00 / -0.21 / -3.58
%     -->B: -2.00 / -0.38 / -0.30
%
% Tom Davidson (tjd@stanford.edu) 2010-2013

% TODO: 
% -consider having armconfig as an output, not an input.

%% Check inputs
switch lower(armconfig)
    case 'lr'
        tiltLR = true;
    case 'fb'
        tiltLR = false;
    otherwise
        error('Must provide ''armconfig'': either ''LR'', or ''FB''');
end

if ~exist('dvintersect', 'var')
    dvintersect = [];
end


if ~exist('quiet', 'var')
    quiet = false;
end

% Check tilt inputs
if tilt < 0
    error('''tilt'' parameter must be positive. Specify direction of tilting with ''rot'' parameter');
end

if tilt >90
    error('''tilt'' must be <90 degrees.');
end

if tilt > 40
    warning('''tilt'' greater than 40 degrees is unusual--are you sure?');
end

if mod (tilt * 360,pi) == 0 && tilt ~=0
    warning('You specified ''tilt'' as a fraction of pi. Use degrees, not radians.');
end

% Check rot inputs
if rot <= -360 || rot >=360
    error('Please specify a rotation ''rot'' between (-360,360) degrees');
end

if mod (rot * 360,pi) == 0 && rot ~=0
    warning('You specified ''rot'' as a fraction of pi. Use degrees, not radians.');
end
       
% Confirm that the specified stereotax arm configuration is reasonable for
% the requested approach angle.
if rot < 0,
    rot = rot+360;
end

arm_tilted_frontback = rot >315 || rot <45 || (rot > 135 && rot <225);

if tiltLR,
    if arm_tilted_frontback
        error('Stereotax arm configuration is wrong for requested approach: use ''FB'' (front-back) tilt instead');
    end
else
    if ~arm_tilted_frontback
        error('Stereotax arm configuration is wrong for requested approach: use ''LR'' (left-right) tilt instead');
    end
end
    

%% Convert inputs to polar coordinates

% Convert DV vector specified by 'tilt' and 'rot' to polar coordinates as
% used by Matlab's built-in sph2cart: 
%  phi = elevation of DV axis from x0y0 plane (i.e. from horizontal) 
%  theta = CCW rotation of DV axis in x0y0 plane away from positive 'X' 
%        (i.e. animal right=0, nose=90)
phi = sub_deg2rad(-tilt+90);
theta = sub_deg2rad(-rot+90);


%% Calculate rotation matrix

% Specify x1 (ML) vector in x0y0z0 space

if tiltLR, % Tilt Left-Right  case:
    % Relative to DV, rotation (theta) is same, tilt (phi) is 90deg CW
    [x01(1) x01(2) x01(3)] = sph2cart(theta, phi-pi/2,1);

else % Tilt Front-Back case
    % ML not tilted (phi=0); rotation (theta) 90deg CW relative to DV axis
    [x01(1) x01(2) x01(3)] = sph2cart(theta-pi/2, 0, 1);

end


% Specify y1 (AP) vector in x0y0z0 space 

% (AP stereotax axis can't tilt or rotate)
y01 = [0 1 0];


% Specify z1 (DV) vector in x0y0z0 space 

% (tilted, rotated as per inputs)
[z01(1) z01(2) z01(3)] = sph2cart(theta, phi,1);


% the inverse of this matrix is a rotation matrix that can be used to get 
% coords in x1y1z1 given target in x0y0z0
rotmat = inv([x01' y01' z01']);


%% Generate coordinates in new reference frame
R = rotmat * target';

% output as row vector
R = R';


%% Calculate intersection of approach vector with horizontal planes

% Optionally calculate intersection of approach vector with requested
% planes parallel to x0y0 plane (Bregma plane), in x0y0z0 coordinates.
% (Note: no dependence on rotation matrix, or on rotated coords x1y1z1)

if ~isempty(dvintersect)
  % radial distance from target to Bregma plane
  Br = [-target(3)+dvintersect] / sin(phi);

  % x0y0 distance from target to Bregma plane
  [B(:,1), B(:,2), B(:,3)] = sph2cart(theta,phi,Br);

  % plus x0y0 distance to target
  B = bsxfun(@plus, B, target);
  
else
  B = [];
end
  

%% Pretty-print outputs

if ~quiet
    disp(sprintf('\nAll coords ML(Right+)/AP(Anterior+)/DV(Dorsal+), as on Kopf stereotax digital display.\n'));
    disp(sprintf('R (target location in new coords): %+0.2f / %+0.2f / %+0.2f \n', R));
    if ~isempty(B),
      disp(sprintf('B (intersection with DV plane in old coords): %+0.2f / %+0.2f / %+0.2f \n', B'));
    end
end



function rad = sub_deg2rad(deg)
% DEG2RAD-converts a matrix of angles in degrees into (-pi pi] radians

rad = mod(deg * (pi/180), 2*pi);

radneg = rad > pi;
rad(radneg) = rad(radneg) - 2*pi;
