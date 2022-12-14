function  vcross = crf( v )

# crf  spatial/planar cross-product operator (force).
# crf(v)  calculates the 6x6 (or 3x3) matrix such that the expression
# crf(v)*f is the cross product of the motion vector v with the force
# vector f.  If length(v)==6 then it is taken to be a spatial vector, and
# the return value is a 6x6 matrix.  Otherwise, v is taken to be a planar
# vector, and the return value is 3x3.

vcross = -crm(v)';
function  vcross = crm( v )

# crm  spatial/planar cross-product operator (motion).
# crm(v)  calculates the 6x6 (or 3x3) matrix such that the expression
# crm(v)*m is the cross product of the motion vectors v and m.  If
# length(v)==6 then it is taken to be a spatial vector, and the return
# value is a 6x6 matrix.  Otherwise, v is taken to be a planar vector, and
# the return value is 3x3.

if length(v) == 6

  vcross = [  0    -v(3)  v(2)   0     0     0    ;
	      v(3)  0    -v(1)   0     0     0    ;
	     -v(2)  v(1)  0      0     0     0    ;
	      0    -v(6)  v(5)   0    -v(3)  v(2) ;
	      v(6)  0    -v(4)   v(3)  0    -v(1) ;
	     -v(5)  v(4)  0     -v(2)  v(1)  0 ];

else

  vcross = [  0     0     0    ;
	      v(3)  0    -v(1) ;
	     -v(2)  v(1)  0 ];
end
function  f = Fpt( fp, p )

# Fpt  forces at points --> spatial/planar forces
# f=Fpt(fp,p)  converts one or more linear forces fp acting at one or more
# points p to their equivalent spatial or planar forces.  In 3D, fp and p
# are 3xn arrays and f is 6xn.  In 2D, fp and p are 2xn and f is 3xn.  In
# both cases, f(:,i) is the equivalent of fp(:,i) acting at p(:,i).

if size(fp,1)==3			# 3D forces at 3D points
  f = [ cross(p,fp,1); fp ];
else					# 2D forces at 2D points
  f = [ p(1,:).*fp(2,:) - p(2,:).*fp(1,:); fp ];
end
function  [o1,o2,o3] = mcI( i1, i2, i3 )

# mcI  rigid-body inertia <--> mass, CoM and rotational inertia.
# rbi=mcI(m,c,I) and [m,c,I]=mcI(rbi) convert between the spatial or planar
# inertia matrix of a rigid body (rbi) and its mass, centre of mass and
# rotational inertia about the centre of mass (m, c and I).  In the spatial
# case, c is 3x1, I is 3x3 and rbi is 6x6.  In the planar case, c is 2x1, I
# is a scalar and rbi is 3x3.  In both cases, m is a scalar.  When c is an
# argument, it can be either a row or a column vector.  If only one
# argument is supplied then it is assumed to be rbi; and if it is a 6x6
# matrix then it is assumed to be spatial.  Otherwise, three arguments must
# be supplied, and if length(c)==3 then mcI calculates a spatial inertia
# matrix.  NOTE: (1) mcI(rbi) requires rbi to have nonzero mass; (2) if |c|
# is much larger than the radius of gyration, or the dimensions of the
# inertia ellipsoid, then extracting I from rbi is numerically
# ill-conditioned.

if nargin == 1
  [o1,o2,o3] = rbi_to_mcI( i1 );
else
  o1 = mcI_to_rbi( i1, i2, i3 );
end


function  rbi = mcI_to_rbi( m, c, I )

if length(c) == 3			# spatial

  C = skew(c);
  rbi = [ I + m*C*C', m*C; m*C', m*eye(3) ];

else					# planar

  rbi = [  I+m*dot(c,c), -m*c(2), m*c(1);
	  -m*c(2),        m,      0;
	   m*c(1),        0,      m  ];
end


function  [m,c,I] = rbi_to_mcI( rbi )

if all(size(rbi)==[6 6])		# spatial

  m = rbi(6,6);
  mC = rbi(1:3,4:6);
  c = skew(mC)/m;
  I = rbi(1:3,1:3) - mC*mC'/m;

else					# planar

  m = rbi(3,3);
  c = [rbi(3,1);-rbi(2,1)]/m;
  I = rbi(1,1) - m*dot(c,c);

end
function  [o1,o2] = plnr( i1, i2 )

# plnr  compose/decompose planar-vector coordinate transform.
# X=plnr(theta,r) and [theta,r]=plnr(X)  compose a planar-vector coordinate
# transform X from its component parts theta and r, and decompose it into
# those parts, respectively.  theta is a scalar and r is a 2D vector.  r is
# returned as a column vector, but it can be supplied as a row or column
# vector.  X is the transform from A to B coordinates, in which frame B is
# located at point r (expressed in A coordinates) and is rotated by an
# angle theta (radians) relative to frame A.  If two arguments are supplied
# then they are assumed to be theta and r, otherwise X.

if nargin == 2				# theta,r --> X

  c = cos(i1);
  s = sin(i1);

  o1 = [        1          0  0 ;
	 s*i2(1)-c*i2(2)   c  s ;
	 c*i2(1)+s*i2(2)  -s  c ];

else					# X --> theta,r

  c = i1(2,2);
  s = i1(2,3);

  o1 = atan2(s,c);
  o2 = [ s*i1(2,1)+c*i1(3,1); s*i1(3,1)-c*i1(2,1) ];

end
function  out = pluho( in )

# pluho  convert Plucker <--> 4x4 homogeneous coordinate transform.
# X=pluho(T) and T=pluho(X) convert between a Plucker coordinate transform
# matrix X and a 4x4 homogeneous coordinate transform matrix T.  If the
# argument is a 6x6 matrix then it is taken to be X, otherwise T.
# NOTE: the 4x4 matrices used in 3D graphics (e.g. OpenGL and Matlab handle
# graphics) are displacement operators, which are the inverses of
# coordinate transforms.  For example, to set the 'matrix' property of a
# Matlab hgtransform graphics object so as to rotate its children by an
# angle theta about the X axis, use inv(pluho(rotx(theta))).

# Formulae:
#   X = [  E   0 ]   T = [ E  -Er ]
#       [ -Erx E ]       [ 0   1  ]

if all(size(in)==[6 6])			# Plucker -> 4x4 homogeneous
  E = in(1:3,1:3);
  mErx = in(4:6,1:3);			# - E r cross
  out = [ E, skew(mErx*E'); 0 0 0 1 ];
else					# 4x4 homogeneous -> Plucker
  E = in(1:3,1:3);
  mEr = in(1:3,4);			# - E r
  out = [ E, zeros(3); skew(mEr)*E, E ];
end
function  [o1,o2] = plux( i1, i2 )

# plux  compose/decompose Plucker coordinate transform.
# X=plux(E,r) and [E,r]=plux(X)  compose a Plucker coordinate transform X
# from its component parts E and r, and decompose it into those parts,
# respectively.  E is a 3x3 rotational coordinate transform and r is a 3D
# vector.  r is returned as a column vector, but it can be supplied as a
# row or column vector.  X is a coordinate transform corresponding to a
# shift of origin by an amount specified by r, followed by a rotation about
# the new origin as specified by E.  For example, plux(rx(1),[2 3 4]) makes
# the same transform as rotx(1)*xlt([2 3 4]).  If two arguments are
# supplied then they are assumed to be E and r, otherwise X.

if nargin == 2				# E,r --> X

  o1 = [ i1, zeros(3); -i1*skew(i2), i1 ];

else					# X --> E,r

  o1 = i1(1:3,1:3);
  o2 = -skew(o1'*i1(4:6,1:3));

end
function  X = rotx( theta )

# rotx  spatial coordinate transform (X-axis rotation).
# rotx(theta)  calculates the coordinate transform matrix from A to B
# coordinates for spatial motion vectors, where coordinate frame B is
# rotated by an angle theta (radians) relative to frame A about their
# common X axis.

c = cos(theta);
s = sin(theta);

X = [ 1  0  0  0  0  0 ;
      0  c  s  0  0  0 ;
      0 -s  c  0  0  0 ;
      0  0  0  1  0  0 ;
      0  0  0  0  c  s ;
      0  0  0  0 -s  c
    ];
function  X = roty( theta )

# roty  spatial coordinate transform (Y-axis rotation).
# roty(theta)  calculates the coordinate transform matrix from A to B
# coordinates for spatial motion vectors, where coordinate frame B is
# rotated by an angle theta (radians) relative to frame A about their
# common Y axis.

c = cos(theta);
s = sin(theta);

X = [ c  0 -s  0  0  0 ;
      0  1  0  0  0  0 ;
      s  0  c  0  0  0 ;
      0  0  0  c  0 -s ;
      0  0  0  0  1  0 ;
      0  0  0  s  0  c
    ];
function  X = rotz( theta )

# rotz  spatial coordinate transform (Z-axis rotation).
# rotz(theta)  calculates the coordinate transform matrix from A to B
# coordinates for spatial motion vectors, where coordinate frame B is
# rotated by an angle theta (radians) relative to frame A about their
# common Z axis.

c = cos(theta);
s = sin(theta);

X = [  c  s  0  0  0  0 ;
      -s  c  0  0  0  0 ;
       0  0  1  0  0  0 ;
       0  0  0  c  s  0 ;
       0  0  0 -s  c  0 ;
       0  0  0  0  0  1
    ];
function  vp = Vpt( v, p )

# Vpt  spatial/planar velocities --> velocities at points
# vp=Vpt(v,p)  calculates the linear velocities vp at one or more points p
# due to one or more spatial/planar velocities v.  In 3D, v is either 6x1
# or 6xn, and p and vp are 3xn.  In 2D, v is either 3x1 or 3xn, and p and
# vp are 2xn.  If v is just a single spatial/planar vector then it applies
# to every point in p; otherwise, vp(:,i) is calculated from v(:,i) and
# p(:,i).

if size(v,2)==1 && size(p,2)>1
  v = repmat(v,1,size(p,2));
end

if size(v,1)==6				# 3D points and velocities
  vp = v(4:6,:) + cross(v(1:3,:),p,1);
else					# 2D points and velocities
  vp = v(2:3,:) + [ -v(1,:).*p(2,:); v(1,:).*p(1,:) ];
end
function  X = xlt( r )

# xlt  spatial coordinate transform (translation of origin).
# xlt(r)  calculates the coordinate transform matrix from A to B
# coordinates for spatial motion vectors, in which frame B is translated by
# an amount r (3D vector) relative to frame A.  r can be a row or column
# vector.

X = [  1     0     0    0  0  0 ;
       0     1     0    0  0  0 ;
       0     0     1    0  0  0 ;
       0     r(3) -r(2) 1  0  0 ;
      -r(3)  0     r(1) 0  1  0 ;
       r(2) -r(1)  0    0  0  1
    ];
function  xp = Xpt( X, p )

# Xpt  apply Plucker/planar coordinate transform to 2D/3D points
# xp=Xpt(X,p)  applies the coordinate transform X to the points in p,
# returning the new coordinates in xp.  If X is a 6x6 matrix then it is
# taken to be a Plucker coordinate transform, and p is expected to be a
# 3xn matrix of 3D points.  Otherwise, X is assumed to be a planar
# coordinate transform and p a 2xn array of 2D points.
  
if all(size(X)==[6 6])			# 3D points
  E = X(1:3,1:3);
  r = -skew(E'*X(4:6,1:3));
else					# 2D points
  E = X(2:3,2:3);
  r = [ X(2,3)*X(2,1)+X(3,3)*X(3,1); X(2,3)*X(3,1)-X(3,3)*X(2,1) ];
end

if size(p,2) > 1
  r = repmat(r,1,size(p,2));
end

xp = E * (p - r);
function v = XtoV( X )

# XtoV  obtain spatial/planar vector from small-angle transform.
# XtoV(X)  interprets X as the coordinate transform from A to B
# coordinates, which implicitly defines the location of frame B relative to
# frame A.  XtoV calculates the velocity of a third frame, C(t), that
# travels at constant velocity from frame A to frame B in one time unit.
# Thus, C(0)=A, C(1)=B and dC/dt is a constant.  The return value, v, is
# the velocity of C, calculated using a small-angle approximation.  It is
# therefore exact only if A and B are parallel.  The return value is a
# spatial vector if X is a 6x6 matrix; otherwise it is a planar vector.
# The return value is an invariant of X (i.e., v=X*v), and can therefore be
# regarded as being expressed in both A and B coordinates.

if all(size(X)==[6 6])			# Plucker xform -> spatial vector

  v = 0.5 * [ X(2,3) - X(3,2);
	      X(3,1) - X(1,3);
	      X(1,2) - X(2,1);
	      X(5,3) - X(6,2);
	      X(6,1) - X(4,3);
	      X(4,2) - X(5,1) ];

else					# planar xform -> planar vector

  v = [  X(2,3);
	 (X(3,1) + X(2,2)*X(3,1) + X(2,3)*X(2,1))/2;
	 (-X(2,1) - X(2,2)*X(2,1) + X(2,3)*X(3,1))/2  ];
end
