function rq(vm)

    # rq  unit quaternion <--> 3x3 coordinate rotation matrix
    # E=rq(q) and q=rq(E)  convert between a unit quaternion q, representing
    # the orientation of a coordinate frame B relative to frame A, and the 3x3
    # coordinate rotation matrix E that transforms from A to B coordinates.
    # For example, if B is rotated relative to A about their common X axis by
    # an angle h, then q=[cos(h/2);sin(h/2);0;0] and rq(q) produces the same
    # matrix as rx(h).  If the argument is a 3x3 matrix then it is assumed to
    # be E, otherwise it is assumed to be q.  rq(E) expects E to be accurately
    # orthonormal, and returns a quaternion in a 4x1 matrix; but rq(q) accepts
    # any nonzero quaternion, contained in either a row or a column vector, and
    # normalizes it before use.  As both q and -q represent the same rotation,
    # rq(E) returns the value that satisfies q[1]>0.  If q[1]==0 then it picks
    # the value such that the largest-magnitude element is positive.  In the
    # event of a tie, the smaller index wins.

    if size(vm) == (3, 3)
        Etoq(vm)
    else
        qtoE(vm)
    end
end

function qtoE(q)

    q = q / norm(q)

    q0s = q[1] * q[1]
    q1s = q[2] * q[2]
    q2s = q[3] * q[3]
    q3s = q[4] * q[4]
    q01 = q[1] * q[2]
    q02 = q[1] * q[3]
    q03 = q[1] * q[4]
    q12 = q[2] * q[3]
    q13 = q[4] * q[2]
    q23 = q[3] * q[4]

    2 * [
        q0s+q1s-0.5 q12+q03 q13-q02
        q12-q03 q0s+q2s-0.5 q23+q01
        q13+q02 q23-q01 q0s+q3s-0.5
    ]
end

function Etoq(E)

    # for sufficiently large q0, this function formulates 2*q0 times the
    # correct return value; otherwise, it formulates 4*|q1| or 4*|q2| or 4*|q3|
    # times the correct value.  The final normalization step yields the correct
    # return value.

    trc = tr(E)# trace is 4*q0^2-1
    v = -skew(E)# v is 2*q0 * [q1;q2;q3]

    if trc > 0
        q = [(trc + 1) / 2; v]
    else
        E = E - (trc - 1) / 2 * I
        E = E + E'
        if E(1, 1) >= E(2, 2) && E(1, 1) >= E(3, 3)
            q = [2 * v(1); E(:, 1)]
        elseif E(2, 2) >= E(3, 3)
            q = [2 * v(2); E(:, 2)]
        else
            q = [2 * v(3); E(:, 3)]
        end
        if q[1] < 0
            q = -q
        end
    end

    q / norm(q)
end

function rqd(in1, in2)
    # rqd  derivative of unit quaternion from angular velocity
    # qd=rqd(wA,q) and qd=rqd(q,wB) calculate the derivative of a unit
    # quaternion, q, representing the orientation of a coordinate frame B
    # relative to frame A, given the angular velocity w of B relative to A.  If
    # w is expressed in A coordinates then use rqd(wA,q); and if w is expressed
    # in B coordinates then use rqd(q,wB).  If the length of the first argument
    # is 4 then it is assumed to be q, otherwise it is assumed to be wA.  The
    # return value is a column vector, but the arguments can be row or column
    # vectors.  It is not necessary for |q| to be exactly 1.  If |q|~=1 then qd
    # contains a magnitude-stabilizing term that will cause |q| to converge
    # towards 1 if q is obtained by numerical integration of qd.

    Kstab = 0.1  # magnitude stabilization constant:
    # value not critical, but K>1 too big

    if length(in1) == 4 # arguments are q and wB
        q = in1
        # w(1:3, 1) = in2
        w = in2
        Q = [q[1] -q[2] -q[3] -q[4]
            q[2] q[1] -q[4] q[3]
            q[3] q[4] q[1] -q[2]
            q[4] -q[3] q[2] q[1]]
    else # arguments are wA and q
        q = in2
        # w(1:3, 1) = in1
        w = in1
        Q = [q[1] -q[2] -q[3] -q[4]
            q[2] q[1] q[4] -q[3]
            q[3] -q[4] q[1] q[2]
            q[4] q[3] -q[2] q[1]]
    end

    0.5 * Q * [Kstab * norm(w) * (1 - norm(q)); w]
end

function rv(vm)

    # rv  3D rotation vector <--> 3x3 coordinate rotation matrix
    # E=rv(v) and v=rv(E) convert between a rotation vector v, whose magnitude
    # and direction describe the angle and axis of rotation of a coordinate
    # frame B relative to frame A, and the 3x3 coordinate rotation matrix E
    # that transforms from A to B coordinates.  For example, if v=[theta;0;0]
    # then rv(v) produces the same matrix as rx(theta).  If the argument is a
    # 3x3 matrix then it is assumed to be E, otherwise it is assumed to be v.
    # rv(E) expects E to be accurately orthonormal, and returns a column vector
    # with a magnitude in the range [0,pi].  If the magnitude is exactly pi
    # then the sign of the return value is unpredictable, since pi*u and -pi*u,
    # where u is any unit vector, both represent the same rotation.  rv(v) will
    # accept a row or column vector of any magnitude.

    if size(vm) == (3, 3)
        Etov(vm)
    else
        vtoE([vm[1], vm[2], vm[3]])
    end
end

function vtoE(v)
    theta = norm(v)
    if theta == 0
        E = Matrix(1.0I, (3, 3))
    else
        s, c = sincos(theta)
        c1 = 2 * sin(theta / 2)^2# 1-cos(h) == 2sin^2(h/2)
        u = v / theta
        E = c * I - s * skew(u) + c1 * u * u'
    end
end

function Etov(E)
    # This function begins by extracting the skew-symmetric component of E,
    # which, as can be seen from the previous function, is -s*skew(v/theta).
    # For angles sufficiently far from pi, v is then calculated directly from
    # this quantity.  However, for angles close to pi, E is almost symmetric,
    # and so extracting the skew-symmetric component becomes numerically
    # ill-conditioned, and provides an increasingly inaccurate value for the
    # direction of v.  Therefore, the component c1*u*u' is extracted as well,
    # and used to get an accurate value for the direction of v.  If the angle
    # is exactly pi then the sign of v will be indeterminate, since +v and -v
    # both represent the same rotation, but the direction will still be
    # accurate.

    w = -skew(E)# w == s/theta * v
    s = norm(w)
    c = (tr(E) - 1) / 2
    theta = atan(s, c)

    if s == 0
        v = [0, 0, 0]
    elseif theta < 0.9 * pi # a somewhat arbitrary threshold
        v = theta / s * w
    else # extract v*v' component from E and
        E = E - c * I # pick biggest column (best chance
        E = E + E'# to get sign right)
        if E[1, 1] >= E[2, 2] && E[1, 1] >= E[3, 3]
            if w[1] >= 0
                v = E[:, 1]
            else
                v = -E[:, 1]
            end
        elseif E[2, 2] >= E[3, 3]
            if w[2] >= 0
                v = E[:, 2]
            else
                v = -E[:, 2]
            end
        else
            if w[3] >= 0
                v = E[:, 3]
            else
                v = -E[:, 3]
            end
        end
        v = theta / norm(v) * v
    end
    return v
end

function rx(theta)
    # rx  3x3 coordinate rotation (X-axis)
    # rx(theta)  calculates the 3x3 rotational coordinate transform matrix from
    # A to B coordinates, where coordinate frame B is rotated by an angle theta
    # (radians) relative to frame A about their common X axis.

    s, c = sincos(theta)

    [1 0 0
        0 c s
        0 -s c]
end

function ry(theta)
    # ry  3x3 coordinate rotation (Y-axis)
    # ry(theta)  calculates the 3x3 rotational coordinate transform matrix from
    # A to B coordinates, where coordinate frame B is rotated by an angle theta
    # (radians) relative to frame A about their common Y axis.

    s, c = sincos(theta)

    [c 0 -s
        0 1 0
        s 0 c]
end

function rz(theta)
    # rz  3x3 coordinate rotation (Z-axis)
    # rz(theta)  calculates the 3x3 rotational coordinate transform matrix from
    # A to B coordinates, where coordinate frame B is rotated by an angle theta
    # (radians) relative to frame A about their common Z axis.

    s, c = sincos(theta)

    [c s 0
        -s c 0
        0 0 1]
end

function skew(vm)
    # skew  convert 3D vector <--> 3x3 skew-symmetric matrix
    # S=skew(v) and v=skew(A) calculate the 3x3 skew-symmetric matrix S
    # corresponding to the given 3D vector v, and the 3D vector corresponding
    # to the skew-symmetric component of the given arbitrary 3x3 matrix A.  For
    # vectors a and b, skew(a)*b is the cross product of a and b.  If the
    # argument is a 3x3 matrix then it is assumed to be A, otherwise it is
    # assumed to be v.  skew(A) produces a column-vector result, but skew(v)
    # will accept a row or column vector argument.

    if size(vm) == (3, 3) # do v = skew(A)
        0.5 * [vm[3, 2] - vm[2, 3], vm[1, 3] - vm[3, 1], vm[2, 1] - vm[1, 2]]
    else # do S = skew(v)
        [0.0 -vm[3] vm[2]
            vm[3] 0 -vm[1]
            -vm[2] vm[1] 0]
    end
end
