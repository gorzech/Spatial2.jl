function  f_out = apply_external_forces( parent, Xup, f_in, f_ext )

# apply_external_forces  subtract f_ext from a given cell array of forces
# f_out=apply_external_forces(parent,Xup,f_in,f_ext)  incorporates the
# external forces specified in f_ext into the calculations of a dynamics
# algorithm.  It does this by subtracting the contents of f_ext from an
# array of forces supplied by the calling function (f_in) and returning the
# result.  f_ext has the following format: (1) it must either be an empty
# cell array, indicating that there are no external forces, or else a cell
# array containing NB elements such that f_ext{i} is the external force
# acting on body i; (2) f_ext{i} must either be an empty array, indicating
# that there is no external force acting on body i, or else a spatial or
# planar vector (as appropriate) giving the external force expressed in
# absolute coordinates.  apply_external_forces performs the calculation
# f_out = f_in - transformed f_ext, where f_out and f_in are cell arrays of
# forces expressed in link coordinates; so f_ext has to be transformed to
# link coordinates before use.  The arguments parent and Xup contain the
# parent array and link-to-link coordinate transforms for the model to
# which the forces apply, and are used to work out the coordinate
# transforms.

# Note: the possibility exists to allow various formats for f_ext;
# e.g. 6/3xNB matrix, or structure with shortened cell array f_ext and a
# list of body numbers (f_ext{i} applies to body b(i))

f_out = f_in;

if length(f_ext) > 0
  for i = 1:length(parent)
    if parent(i) == 0
      Xa{i} = Xup{i};
    else
      Xa{i} = Xup{i} * Xa{parent(i)};
    end
    if length(f_ext{i}) > 0
      f_out{i} = f_out{i} - Xa{i}' \ f_ext{i};
    end
  end
end
function  ret = EnerMo( model, q, qd )

# EnerMo  calculate energy, momentum and related quantities
# EnerMo(robot,q,qd)  returns a structure containing the fields KE, PE,
# htot, Itot, mass, cm and vcm.  These fields contain the kinetic and
# potential energies of the whole system, the total spatial momentum, the
# total spatial inertia, total mass, position of centre of mass, and the
# linear velocity of centre of mass, respectively.  Vector quantities are
# expressed in base coordinates.  PE is defined to be zero when cm is
# zero.

for i = 1:model.NB
  [ XJ, S ] = jcalc( model.jtype{i}, q(i) );
  vJ = S*qd(i);
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
  end
  Ic{i} = model.I{i};
  hc{i} = Ic{i} * v{i};
  KE(i) = 0.5 * v{i}' * hc{i};
end

ret.Itot = zeros(size(Ic{1}));
ret.htot = zeros(size(hc{1}));

for i = model.NB:-1:1
  if model.parent(i) ~= 0
    Ic{model.parent(i)} = Ic{model.parent(i)} + Xup{i}'*Ic{i}*Xup{i};
    hc{model.parent(i)} = hc{model.parent(i)} + Xup{i}'*hc{i};
  else
    ret.Itot = ret.Itot + Xup{i}'*Ic{i}*Xup{i};
    ret.htot = ret.htot + Xup{i}'*hc{i};
  end
end

a_grav = get_gravity(model);

if length(a_grav) == 6
  g = a_grav(4:6);			# 3D linear gravitational accn
  h = ret.htot(4:6);			# 3D linear momentum
else
  g = a_grav(2:3);			# 2D gravity
  h = ret.htot(2:3);			# 2D linear momentum
end

[mass, cm] = mcI(ret.Itot);

ret.KE = sum(KE);
ret.PE = - mass * dot(cm,g);
ret.mass = mass;
ret.cm = cm;
ret.vcm = h / mass;
function  Q = fbanim( X, Qr )

# fbanim  Floating Base Inverse Kinematics for Animation
# Q=fbanim(X) and Q=fbanim(X,Qr) calculate matrices of joint position data
# for use in showmotion animations.  Q=fbanim(X) calculates a 6xN matrix Q
# from a 13x1xN or 13xN or 7x1xN or 7xN matrix X, such that each column of
# X contains at least the first 7 elements of a 13-element singularity-free
# state vector used by FDfb and IDfb, and the corresponding column of Q
# contains the position variables of the three prismatic and three revolute
# joints that form the floating base.  Each column of Q is calculated using
# fbkin; but fbanim then adjusts the three revolute joint variables to
# remove the pi and 2*pi jumps that occur when the angles wrap around or
# pass through a kinematic singularity.  As a result, some or all of the
# revolute joint angles may grow without limit.  This is essential for a
# smooth animation.  Q=fbanim(X,Qr) performs the same calculation as just
# described, but then appends Qr to Q.  Qr must be an MxN or Mx1xN matrix
# of joint position data for the real joints in a floating-base mechanism,
# and the resulting (6+M)xN matrix Q contains the complete set of joint
# data required by showmotion.  Note: the algorithm used to remove the
# jumps makes continuity assumptions (specifically, less than pi/2 changes
# from one column to the next, except when passing through a singularity)
# that are not guaranteed to be true.  Therefore, visible glitches in an
# animation are still possible.

if length(size(X)) == 3			# collapse 3D -> 2D array
  tmp(:,:) = X(:,1,:);
  X = tmp;
end

# apply kinematic transform using fbkin

for i = 1 : size(X,2)
  Q(:,i) = fbkin( X(1:7,i) );
end

# This code removes wrap-arounds and step-changes on passing through a
# singularity.  Whenever q6 or q4 wrap, they jump by 2*pi.  However, when
# q5 hits pi/2 (or -pi/2), q4 and q6 both jump by pi, and q5 turns around.
# To undo this, the code looks to see whether n is an even or odd number,
# and replaces q5 with pi-q5 whenever n is odd.  q4 is calculated via the
# sum or difference of q4 and q6 on the grounds that the sum well defined
# at the singularity at q5=pi/2 and the difference is well defined at
# q5=-pi/2.

for i = 2 : size(X,2)
  n = round( (Q(6,i-1) - Q(6,i)) / pi );
  q6 = Q(6,i) + n*pi;
  if Q(5,i) >= 0
    q46 = Q(4,i) + Q(6,i);
    q46 = q46 + 2*pi * round( (Q(4,i-1)+Q(6,i-1) - q46) / (2*pi) );
    Q(4,i) = q46 - q6;
  else
    q46 = Q(4,i) - Q(6,i);
    q46 = q46 + 2*pi * round( (Q(4,i-1)-Q(6,i-1) - q46) / (2*pi) );
    Q(4,i) = q46 + q6;
  end
  Q(6,i) = q6;
  if mod(n,2) == 0
    q5 = Q(5,i);
  else
    q5 = pi - Q(5,i);
  end
  Q(5,i) = q5 + 2*pi * round( (Q(5,i-1) - q5) / (2*pi) );
end

# add the rest of the joint data, if argument Qr has been supplied

if nargin == 2
  if length(size(Qr)) == 3		# collapse 3D -> 2D array
    tmp = zeros(size(Qr,1),size(Qr,3));
    tmp(:,:) = Qr(:,1,:);
    Qr = tmp;
  end
  Q = [ Q; Qr ];
end
function  [o1, o2, o3] = fbkin( i1, i2, i3 )

# fbkin  Forward and Inverse Kinematics of Floating Base
# [x,xd]=fbkin(q,qd,qdd) calculates the forward kinematics, and
# [q,qd,qdd]=fbkin(x,xd) calculates the inverse kinematics, of a spatial
# floating base constructed by floatbase.  The vectors q, qd and qdd are
# the joint position, velocity and acceleration variables for the six
# joints forming the floating base; x is the 13-element singularity-free
# state vector used by FDfb and IDfb; and xd is its derivative.  The
# component parts of x are: a unit quaternion specifying the orientation of
# the floating-base frame relative to the fixed base, a 3D vector giving
# the position of the origin of the floating-base frame, and the spatial
# velocity of the floating base expressed in fixed-base coordinates.  If
# the acceleration calculation is not needed then use the simpler calls
# x=fbkin(q,qd) and [q,qd]=fbkin(x); and if the velocity calculation is
# also not needed then use p=fbkin(q) and q=fbkin(p) (or q=fbkin(x)), where
# p contains the first 7 elements of x.  If the length of the first
# argument is 7 or 13 then it is assumed to be p or x, otherwise it is
# assumed to be q.  All vectors are expected to be column vectors.  The
# returned value of q(5) is normalized to the range -pi/2 to pi/2; and the
# returned values of q(4) and q(6) are normalized to the range -pi to pi.
# The calculation of qd and qdd fails if the floating base is in one of its
# two kinematic singularities, which occur when q(5)==+/-pi/2.  However,
# the calculation of q is robust, and remains accurate both in and close to
# the singularities.

if length(i1) == 13 || length(i1) == 7
  if nargout == 3
    [o1, o2, o3] = invkin( i1, i2 );
  elseif nargout == 2
    [o1, o2] = invkin( i1 );
  else
    o1 = invkin( i1 );
  end
else
  if nargout == 2
    [o1, o2] = fwdkin( i1, i2, i3 );
  elseif nargin == 2
    o1 = fwdkin( i1, i2 );
  else
    o1 = fwdkin( i1 );
  end
end



function  [x, xd] = fwdkin( q, qd, qdd )

c4 = cos(q(4));  s4 = sin(q(4));
c5 = cos(q(5));  s5 = sin(q(5));
c6 = cos(q(6));  s6 = sin(q(6));

E = [  c5*c6,  c4*s6+s4*s5*c6,  s4*s6-c4*s5*c6;
      -c5*s6,  c4*c6-s4*s5*s6,  s4*c6+c4*s5*s6;
        s5,       -s4*c5,           c4*c5 ];

qn = rq(E);				# unit quaternion fixed-->floating
r = q(1:3);				# position of floating-base origin

x(1:4,1) = qn;
x(5:7) = r;

if nargin > 1				# do velocity calculation

  S = [ 1  0    s5;
	0  c4  -s4*c5;
	0  s4   c4*c5 ];

  omega = S*qd(4:6);
  rd = qd(1:3);				# lin vel of floating-base origin

  v = [ omega; rd+cross(r,omega) ];	# spatial vel in fixed-base coords

  x(8:13) = v;
end

if nargout == 2				# calculate xd

  c4d = -s4*qd(4);  s4d = c4*qd(4);
  c5d = -s5*qd(5);  s5d = c5*qd(5);

  Sd = [ 0  0     s5d;
	 0  c4d  -s4d*c5-s4*c5d;
	 0  s4d   c4d*c5+c4*c5d ];

  omegad = S*qdd(4:6) + Sd*qd(4:6);
  rdd = qdd(1:3);
					# spatial accn in fixed-base coords
  a = [ omegad; rdd+cross(rd,omega)+cross(r,omegad) ];

  xd(1:4,1) = rqd( omega, rq(E) );
  xd(5:7) = rd;
  xd(8:13) = a;
end



function  [q, qd, qdd] = invkin( x, xd )

E = rq(x(1:4));				# coord xfm fixed-->floating
r = x(5:7);				# position of floating-base origin

q(1:3,1) = r;

q(5) = atan2( E(3,1), sqrt(E(1,1)*E(1,1)+E(2,1)*E(2,1)) );
q(6) = atan2( -E(2,1), E(1,1) );
if E(3,1) > 0
  q(4) = atan2( E(2,3)+E(1,2), E(2,2)-E(1,3) ) - q(6);
else
  q(4) = atan2( E(2,3)-E(1,2), E(2,2)+E(1,3) ) + q(6);
end
if q(4) > pi
  q(4) = q(4) - 2*pi;
elseif q(4) < -pi
  q(4) = q(4) + 2*pi;
end

if nargout > 1				# calculate qd

  c4 = cos(q(4));  s4 = sin(q(4));
  c5 = cos(q(5));  s5 = sin(q(5));

  S = [ 1  0    s5;
	0  c4  -s4*c5;
	0  s4   c4*c5 ];

  omega = x(8:10);
  rd = x(11:13) - cross(r,omega);	# lin vel of floating-base origin

  qd(1:3,1) = rd;
  qd(4:6) = S \ omega;			# this will fail at a singularity
end

if nargout == 3				# calculate qdd

  c4d = -s4*qd(4);  s4d = c4*qd(4);
  c5d = -s5*qd(5);  s5d = c5*qd(5);

  Sd = [ 0  0     s5d;
	 0  c4d  -s4d*c5-s4*c5d;
	 0  s4d   c4d*c5+c4*c5d ];

  omegad = xd(8:10);
  rdd = xd(11:13) - cross(rd,omega) - cross(r,omegad);

  qdd(1:3,1) = rdd;
  qdd(4:6) = S \ (omegad - Sd*qd(4:6));	# this will fail at a singularity
end
function  qdd = FDab( model, q, qd, tau, f_ext )

# FDab  Forward Dynamics via Articulated-Body Algorithm
# FDab(model,q,qd,tau,f_ext,grav_accn)  calculates the forward dynamics of
# a kinematic tree via the articulated-body algorithm.  q, qd and tau are
# vectors of joint position, velocity and force variables; and the return
# value is a vector of joint acceleration variables.  f_ext is an optional
# argument specifying the external forces acting on the bodies.  It can be
# omitted if there are no external forces.  The format of f_ext is
# explained in the source code of apply_external_forces.

a_grav = get_gravity(model);

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
  vJ = S{i}*qd(i);
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
    c{i} = zeros(size(a_grav));		# spatial or planar zero vector
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    c{i} = crm(v{i}) * vJ;
  end
  IA{i} = model.I{i};
  pA{i} = crf(v{i}) * model.I{i} * v{i};
end

if nargin == 5
  pA = apply_external_forces( model.parent, Xup, pA, f_ext );
end

for i = model.NB:-1:1
  U{i} = IA{i} * S{i};
  d{i} = S{i}' * U{i};
  u{i} = tau(i) - S{i}'*pA{i};
  if model.parent(i) ~= 0
    Ia = IA{i} - U{i}/d{i}*U{i}';
    pa = pA{i} + Ia*c{i} + U{i} * u{i}/d{i};
    IA{model.parent(i)} = IA{model.parent(i)} + Xup{i}' * Ia * Xup{i};
    pA{model.parent(i)} = pA{model.parent(i)} + Xup{i}' * pa;
  end
end

for i = 1:model.NB
  if model.parent(i) == 0
    a{i} = Xup{i} * -a_grav + c{i};
  else
    a{i} = Xup{i} * a{model.parent(i)} + c{i};
  end
  qdd(i,1) = (u{i} - U{i}'*a{i})/d{i};
  a{i} = a{i} + S{i}*qdd(i);
end
function  qdd = FDcrb( model, q, qd, tau, f_ext )

# FDcrb  Forward Dynamics via Composite-Rigid-Body Algorithm
# FDcrb(model,q,qd,tau,f_ext)  calculates the forward dynamics of a
# kinematic tree via the composite-rigid-body algorithm.  q, qd and tau are
# vectors of joint position, velocity and force variables; and the return
# value is a vector of joint acceleration variables.  f_ext is an optional
# argument specifying the external forces acting on the bodies.  It can be
# omitted if there are no external forces.  The format of f_ext is
# explained in the source code of apply_external_forces.

if nargin == 4
  [H,C] = HandC( model, q, qd );
else
  [H,C] = HandC( model, q, qd, f_ext );
end

qdd = H \ (tau - C);
function  [xdfb,qdd] = FDfb( model, xfb, q, qd, tau, f_ext )

# FDfb  Floating-Base Forward Dynamics via Articulated-Body Algorithm
# [xdfb,qdd]=FDfb(model,xfb,q,qd,tau,f_ext) calculates the forward dynamics
# of a floating-base kinematic tree via the articulated-body algorithm.
# This function avoids the kinematic singularity in the six-joint chain
# created by floatbase to mimic a true 6-DoF joint.  xfb is a 13-element
# column vector containing: a unit quaternion specifying the orientation of
# the floating base (=body 6)'s coordinate frame relative to the fixed
# base; a 3D vector specifying the position of the origin of the floating
# base's coordinate frame in fixed-base coordinates; and a spatial vector
# giving the velocity of the floating base in fixed-base coordinates.  The
# return value xdfb is the time-derivative of this vector.  The arguments
# q, qd and tau contain the position, velocity and force variables for the
# real joints in the system (i.e., joints 7 onwards in the system model);
# so q(i), qd(i) and tau(i) all apply to joint i+6.  The return value qdd
# is the time-derivative of qd.  f_ext is an optional argument specifying
# the external forces acting on the bodies.  It can be omitted if there are
# no external forces.  If supplied, it must be a cell array of length
# model.NB, of which the first 5 elements are ignored, and f_ext{6} onward
# specify the forces acting on the floating base (body 6) onward.  The
# format of f_ext is explained in the source code of apply_external_forces.

a_grav = get_gravity(model);

qn = xfb(1:4);				# unit quaternion fixed-->f.b.
r = xfb(5:7);				# position of f.b. origin
Xup{6} = plux( rq(qn), r );		# xform fixed --> f.b. coords

vfb = xfb(8:end);
v{6} = Xup{6} * vfb;			# f.b. vel in f.b. coords

IA{6} = model.I{6};
pA{6} = crf(v{6}) * model.I{6} * v{6};

for i = 7:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i-6) );
  vJ = S{i}*qd(i-6);
  Xup{i} = XJ * model.Xtree{i};
  v{i} = Xup{i}*v{model.parent(i)} + vJ;
  c{i} = crm(v{i}) * vJ;
  IA{i} = model.I{i};
  pA{i} = crf(v{i}) * model.I{i} * v{i};
end

if nargin == 6 && length(f_ext) > 0
  prnt = model.parent(6:end) - 5;
  pA(6:end) = apply_external_forces(prnt, Xup(6:end), pA(6:end), f_ext(6:end));
end

for i = model.NB:-1:7
  U{i} = IA{i} * S{i};
  d{i} = S{i}' * U{i};
  u{i} = tau(i-6) - S{i}'*pA{i};
  Ia = IA{i} - U{i}/d{i}*U{i}';
  pa = pA{i} + Ia*c{i} + U{i} * u{i}/d{i};
  IA{model.parent(i)} = IA{model.parent(i)} + Xup{i}' * Ia * Xup{i};
  pA{model.parent(i)} = pA{model.parent(i)} + Xup{i}' * pa;
end

a{6} = - IA{6} \ pA{6};			# floating base accn without gravity

qdd = zeros(0,1);			# avoids a matlab warning when NB==6

for i = 7:model.NB
  a{i} = Xup{i} * a{model.parent(i)} + c{i};
  qdd(i-6,1) = (u{i} - U{i}'*a{i})/d{i};
  a{i} = a{i} + S{i}*qdd(i-6);
end

qnd = rqd( vfb(1:3), qn );		# derivative of qn
rd = Vpt( vfb, r );			# lin vel of flt base origin
afb = Xup{6} \ a{6} + a_grav;		# true f.b. accn in fixed-base coords

xdfb = [ qnd; rd; afb ];
function  qdd = FDgq( model, q, qd, tau, f_ext )

# FDgq  Forward Dynamics via CRBA + constraint function gamma_q
# FDgq(model,q,qd,tau,f_ext)  calculates the forward dynamics of a
# kinematic tree, subject to the kinematic constraints embodied in the
# function model.gamma_q, via the composite-rigid-body algorithm.  q, qd
# and tau are vectors of joint position, velocity and force variables; and
# the return value is a vector of joint acceleration variables.  q and qd do
# not have to satisfy the constraints exactly, but are expected to be
# close.  qdd will typically contain a constraint-stabilization component
# that tends to reduce constraint violation over time.  f_ext is an optional
# argument specifying the external forces acting on the bodies.  It can be
# omitted if there are no external forces.  The format of f_ext is
# explained in the source code of apply_external_forces.  A detailed
# description of gamma_q appears at the end of this source code.

[q,qd,G,g] = model.gamma_q( model, q, qd );

if nargin == 4
  [H,C] = HandC( model, q, qd );
else
  [H,C] = HandC( model, q, qd, f_ext );
end

qdd = G * ((G'*H*G) \ (G'*(tau-C-H*g))) + g;       # cf. eq 3.20 in RBDA




#{
How to Create Your Own gamma_q
------------------------------

[See also Section 3.2 of "Rigid Body Dynamics Algorithms", and maybe also
Section 8.3.]

The purpose of this function is to define a set of algebraic constraints among
the elements of q, and therefore also among the elements of qd and qdd.  These
constraints could be due to kinematic loops, gears, pulleys, and so on.  To
create a function gamma_q, proceed as follows:

1.  Identify a set of independent variables, y.  Typically, y will be a subset
    of the variables in q.

2.  Define a function, gamma, that maps y to q; i.e., q=gamma(y).  The
    calculated value of q must satisfy the constraints exactly.

3.  Define an inverse function, gamma^-1, that maps q to y.  This function
    must satisfy gamma^-1(gamma(y))==y for all y, and gamma(gamma^-1(q))==q
    for all q that satisfy the constraints exactly.

4.  Optional:  For some kinds of constraints, gamma may be ambiguous (i.e.,
    multi-valued).  In these cases, gamma can be modified to take a second
    argument, q0, which is a vector satisfying gamma^-1(q0)==y.  This vector
    plays the role of a disambiguator---it identifies which one of multiple
    possible values of gamma(y) is the correct one in the present context.
    For differentiation purposes, this extra argument is regarded as a
    constant.

5.  Having defined gamma, G is the matrix partial d/dy gamma (i.e., it is
    the Jacobian of gamma).

6.  G provides the following relationship between qd and yd: qd = G * yd,
    which can be used to calculate qd from yd.  Alternatively, qd could be
    calculated directly from d/dt gamma(y).

7.  To obtain yd from qd, use yd = d/dt gamma^-1(q).  In general, the
    right-hand side will be a function of both q and qd.  However, if y is a
    subset of q then yd is the same subset of qd.

8.  Given that qd = G * yd, it follows that qdd = G * ydd + dG/dt * yd.
    However, it is advisable to add a stabilization term to this formula, so
    that qdd = G * ydd + g, where g = dG/dt * yd + gs, and gs is a
    stabilization term defined as follows: gs = 2/Ts * qderr + 1/Ts^2 * qerr.
    Ts is a stabilization time constant, and qerr and qderr are measures of
    the degree to which the (given) prevailing values of q and qd fail to
    satisfy the constraints.  If q0 and qd0 are the prevailing values then
    qerr = gamma(gamma^-1(q0)) - q0 and qderr = G*yd - qd0 (where yd is
    calculated from qd0).

The function gamma_q can now be defined as follows:

function call:  [q,qd,G,g] = gamma_q( model, q0, qd0 )

where q0 and qd0 are the current values of the joint position and velocity
variables, q and qd are new values that exactly satisfy the constraints, and G
and g are as described above.

function body:

y = gamma^-1(q0);
q = gamma(y);			# or gamma(y,q0)

G = Jacobian of gamma;

yd = d/dt gamma^-1(q0);		# a function of q0 and qd0
qd = G * yd;			# or d/dt gamma(y)

Ts = some suitable value, such as 0.1 or 0.01;

gs = 2/Ts * (qd - qd0) + 1/Ts^2 * (q - q0);

g = dG/dt * yd + gs;

Tip: gamma_q can be used at the beginning of a simulation run to initialize q
and qd to values that satisfy the constraints exactly.
#}
function  a_grav = get_gravity( model )

# get_gravity  spatial/planar gravitational accn vector for given model
# get_gravity(model) returns the gravitational acceleration vector to be
# used in dynamics calculations for the given model.  The return value is
# either a spatial or a planar vector, according to the type of model.  It
# is computed from the field model.gravity, which is a 2D or 3D (row or
# column) vector specifying the linear acceleration due to gravity.  If
# this field is not present then get_gravity uses the following defaults:
# [0,0,-9.81] for spatial models and [0,0] for planar.

if isfield( model, 'gravity' )
  g = model.gravity;
else
  g = [0;0;-9.81];
end

if size(model.Xtree{1},1) == 3		# is model planar?
  a_grav = [0;g(1);g(2)];
else
  a_grav = [0;0;0;g(1);g(2);g(3)];
end
function  [H,C] = HandC( model, q, qd, f_ext )

# HandC  Calculate coefficients of equation of motion
# [H,C]=HandC(model,q,qd,f_ext)  calculates the coefficients of the
# joint-space equation of motion, tau=H(q)qdd+C(d,qd,f_ext), where q, qd
# and qdd are the joint position, velocity and acceleration vectors, H is
# the joint-space inertia matrix, C is the vector of gravity,
# external-force and velocity-product terms, and tau is the joint force
# vector.  Algorithm: recursive Newton-Euler for C, and
# Composite-Rigid-Body for H.  f_ext is an optional argument specifying the
# external forces acting on the bodies.  It can be omitted if there are no
# external forces.  The format of f_ext is explained in the source code of
# apply_external_forces.

a_grav = get_gravity(model);

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
  vJ = S{i}*qd(i);
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
    avp{i} = Xup{i} * -a_grav;
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    avp{i} = Xup{i}*avp{model.parent(i)} + crm(v{i})*vJ;
  end
  fvp{i} = model.I{i}*avp{i} + crf(v{i})*model.I{i}*v{i};
end

if nargin == 4
  fvp = apply_external_forces( model.parent, Xup, fvp, f_ext );
end

for i = model.NB:-1:1
  C(i,1) = S{i}' * fvp{i};
  if model.parent(i) ~= 0
    fvp{model.parent(i)} = fvp{model.parent(i)} + Xup{i}'*fvp{i};
  end
end

IC = model.I;				# composite inertia calculation

for i = model.NB:-1:1
  if model.parent(i) ~= 0
    IC{model.parent(i)} = IC{model.parent(i)} + Xup{i}'*IC{i}*Xup{i};
  end
end

H = zeros(model.NB);

for i = 1:model.NB
  fh = IC{i} * S{i};
  H(i,i) = S{i}' * fh;
  j = i;
  while model.parent(j) > 0
    fh = Xup{j}' * fh;
    j = model.parent(j);
    H(i,j) = S{j}' * fh;
    H(j,i) = H(i,j);
  end
end
function  [qdd_out,tau_out] = HD( model, fd, q, qd, qdd, tau, f_ext )

# HD  Articulated-Body Hybrid Dynamics Algorithm
# [qdd_out,tau_out]=HD(model,fd,q,qd,qdd,tau,f_ext) calculates the hybrid
# dynamics of a kinematic tree using the articulated-body algorithm.  fd is
# an array of boolean values such that fd(i)==1 if joint i is a
# forward-dynamics joint, and fd(i)==0 otherwise.  If fd(i)==1 then tau(i)
# contains the given force at joint i, and the value of qdd(i) is ignored;
# and if fd(i)==0 then qdd(i) contains the given acceleration at joint i,
# and the value of tau(i) is ignored.  Likewise, if fd(i)==1 then
# qdd_out(i) contains the calculated acceleration at joint i, and
# tau_out(i) contains the given force copied from tau(i); and if fd(i)==0
# then tau_out(i) contains the calculated force and qdd_out(i) the given
# acceleration copied from qdd(i).  Thus, the two output vectors are always
# fully instantiated.  f_ext is an optional argument specifying the
# external forces acting on the bodies.  It can be omitted if there are no
# external forces.  The format of f_ext is explained in the source code of
# apply_external_forces.

a_grav = get_gravity(model);

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
  vJ = S{i}*qd(i);
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
    c{i} = zeros(size(a_grav));		# spatial or planar zero vector
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    c{i} = crm(v{i}) * vJ;
  end
  if fd(i) == 0
    c{i} = c{i} + S{i} * qdd(i);
  end
  IA{i} = model.I{i};
  pA{i} = crf(v{i}) * model.I{i} * v{i};
end

if nargin == 7
  pA = apply_external_forces( model.parent, Xup, pA, f_ext );
end

for i = model.NB:-1:1
  if fd(i) == 0
    if model.parent(i) ~= 0
      Ia = IA{i};
      pa = pA{i} + IA{i}*c{i};
      IA{model.parent(i)} = IA{model.parent(i)} + Xup{i}' * Ia * Xup{i};
      pA{model.parent(i)} = pA{model.parent(i)} + Xup{i}' * pa;
    end
  else
    U{i} = IA{i} * S{i};
    d{i} = S{i}' * U{i};
    u{i} = tau(i) - S{i}'*pA{i};
    if model.parent(i) ~= 0
      Ia = IA{i} - U{i}/d{i}*U{i}';
      pa = pA{i} + Ia*c{i} + U{i} * u{i}/d{i};
      IA{model.parent(i)} = IA{model.parent(i)} + Xup{i}' * Ia * Xup{i};
      pA{model.parent(i)} = pA{model.parent(i)} + Xup{i}' * pa;
    end
  end
end

for i = 1:model.NB
  if model.parent(i) == 0
    a{i} = Xup{i} * -a_grav + c{i};
  else
    a{i} = Xup{i} * a{model.parent(i)} + c{i};
  end
  if fd(i) == 0
    qdd_out(i,1) = qdd(i);
    tau_out(i,1) = S{i}'*(IA{i}*a{i} + pA{i});
  else
    qdd_out(i,1) = (u{i} - U{i}'*a{i})/d{i};
    tau_out(i,1) = tau(i);
    a{i} = a{i} + S{i}*qdd_out(i);
  end
end
function  tau = ID( model, q, qd, qdd, f_ext )

# ID  Inverse Dynamics via Recursive Newton-Euler Algorithm
# ID(model,q,qd,qdd,f_ext) calculates the inverse dynamics of a kinematic
# tree via the recursive Newton-Euler algorithm.  q, qd and qdd are vectors
# of joint position, velocity and acceleration variables; and the return
# value is a vector of joint force variables.  f_ext is an optional
# argument specifying the external forces acting on the bodies.  It can be
# omitted if there are no external forces.  The format of f_ext is
# explained in the source code of apply_external_forces.

a_grav = get_gravity(model);

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
  vJ = S{i}*qd(i);
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
    a{i} = Xup{i}*(-a_grav) + S{i}*qdd(i);
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    a{i} = Xup{i}*a{model.parent(i)} + S{i}*qdd(i) + crm(v{i})*vJ;
  end
  f{i} = model.I{i}*a{i} + crf(v{i})*model.I{i}*v{i};
end

if nargin == 5
  f = apply_external_forces( model.parent, Xup, f, f_ext );
end

for i = model.NB:-1:1
  tau(i,1) = S{i}' * f{i};
  if model.parent(i) ~= 0
    f{model.parent(i)} = f{model.parent(i)} + Xup{i}'*f{i};
  end
end
function  [xdfb,tau] = IDfb( model, xfb, q, qd, qdd, f_ext )

# IDfb  Floating-Base Inverse Dynamics (=Hybrid Dynamics)
# [xdfb,tau]=IDfb(model,xfb,q,qd,qdd,f_ext) calculates the inverse dynamics
# of a floating-base kinematic tree via the algorithm in Table 9.6 of RBDA
# (which is really a special case of hybrid dynamics), using the same
# singularity-free representation of the motion of the floating base as
# used by FDfb.  xfb is a 13-element column vector containing: a unit
# quaternion specifying the orientation of the floating base (=body 6)'s
# coordinate frame relative to the fixed base; a 3D vector specifying the
# position of the origin of the floating base's coordinate frame in
# fixed-base coordinates; and a spatial vector giving the velocity of the
# floating base in fixed-base coordinates.  The return value xdfb is the
# time-derivative of xfb.  The arguments q, qd and qdd contain the
# position, velocity and acceleration variables for the real joints in the
# system (i.e., joints 7 onwards in the system model); so q(i), qd(i) and
# qdd(i) all apply to joint i+6.  The return value tau is the vector of
# force variables required to produce the given acceleration qdd.  f_ext is
# an optional argument specifying the external forces acting on the bodies.
# It can be omitted if there are no external forces.  If supplied, it must
# be a cell array of length model.NB, of which the first 5 elements are
# ignored, and f_ext{6} onward specify the forces acting on the floating
# base (body 6) onward.  The format of f_ext is explained in the source
# code of apply_external_forces.

a_grav = get_gravity(model);

qn = xfb(1:4);				# unit quaternion fixed-->f.b.
r = xfb(5:7);				# position of f.b. origin
Xup{6} = plux( rq(qn), r );		# xform fixed --> f.b. coords

vfb = xfb(8:end);
v{6} = Xup{6} * vfb;			# f.b. vel in f.b. coords

a{6} = zeros(6,1);

IC{6} = model.I{6};
pC{6} = model.I{6}*a{6} + crf(v{6})*model.I{6}*v{6};

for i = 7:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i-6) );
  vJ = S{i}*qd(i-6);
  Xup{i} = XJ * model.Xtree{i};
  v{i} = Xup{i}*v{model.parent(i)} + vJ;
  a{i} = Xup{i}*a{model.parent(i)} + S{i}*qdd(i-6) + crm(v{i})*vJ;
  IC{i} = model.I{i};
  pC{i} = IC{i}*a{i} + crf(v{i})*IC{i}*v{i};
end

if nargin == 6 && length(f_ext) > 0
  prnt = model.parent(6:end) - 5;
  pC(6:end) = apply_external_forces(prnt, Xup(6:end), pC(6:end), f_ext(6:end));
end

for i = model.NB:-1:7
  IC{model.parent(i)} = IC{model.parent(i)} + Xup{i}'*IC{i}*Xup{i};
  pC{model.parent(i)} = pC{model.parent(i)} + Xup{i}'*pC{i};
end

a{6} = - IC{6} \ pC{6};			# floating-base acceleration
                                        # without gravity
for i = 7:model.NB
  a{i} = Xup{i} * a{model.parent(i)};
  tau(i-6,1) = S{i}'*(IC{i}*a{i} + pC{i});
end

qnd = rqd( vfb(1:3), qn );		# derivative of qn
rd = Vpt( vfb, r );			# lin vel of flt base origin
afb = Xup{6} \ a{6} + a_grav;		# f.b. accn in fixed-base coords

xdfb = [ qnd; rd; afb ];
function  [Xj,S] = jcalc( jtyp, q )

# jcalc  joint transform and motion subspace matrices.
# [Xj,S]=jcalc(type,q)  returns the joint transform and motion subspace
# matrices for a joint of the given type.  jtyp can be either a string or a
# structure containing a string-valued field called 'code'.  Either way,
# the string contains the joint type code.  For joints that take
# parameters (e.g. helical), jtyp must be a structure, and it must contain
# a field called 'pars', which in turn is a structure containing one or
# more parameters.  (For a helical joint, pars contains a parameter called
# 'pitch'.)  q is the joint's position variable.

if ischar( jtyp )
  code = jtyp;
else
  code = jtyp.code;
end

switch code
  case 'Rx'				# revolute X axis
    Xj = rotx(q);
    S = [1;0;0;0;0;0];
  case 'Ry'				# revolute Y axis
    Xj = roty(q);
    S = [0;1;0;0;0;0];
  case {'R','Rz'}			# revolute Z axis
    Xj = rotz(q);
    S = [0;0;1;0;0;0];
  case 'Px'				# prismatic X axis
    Xj = xlt([q 0 0]);
    S = [0;0;0;1;0;0];
  case 'Py'				# prismatic Y axis
    Xj = xlt([0 q 0]);
    S = [0;0;0;0;1;0];
  case {'P','Pz'}			# prismatic Z axis
    Xj = xlt([0 0 q]);
    S = [0;0;0;0;0;1];
  case 'H'				# helical (Z axis)
    Xj = rotz(q) * xlt([0 0 q*jtyp.pars.pitch]);
    S = [0;0;1;0;0;jtyp.pars.pitch];
  case 'r'				# planar revolute
    Xj = plnr( q, [0 0] );
    S = [1;0;0];
  case 'px'				# planar prismatic X axis
    Xj = plnr( 0, [q 0] );
    S = [0;1;0];
  case 'py'				# planar prismatic Y axis
    Xj = plnr( 0, [0 q] );
    S = [0;0;1];
  otherwise
    error( 'unrecognised joint code ''#s''', code );
end
