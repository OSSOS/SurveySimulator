module elemutils

  use datadec

contains
! \subroutine{coord\_cart}

  subroutine coord_cart (mu, o_m, p, v)

! This routine transforms delaunay variables into cartisian
! variables.
!
! ANGLES ARE GIVEN IN RADIAN !!!!
!
! \subsection{Arguments}
! \subsubsection{Definitions}
! \begin{verse}
! \verb|a, e, inc, node, peri, m| = osculating elements \\
! \verb|x, y, z, vx, vy, vz| = cartesian variables: $X, Y, Z, Px, Py, Pz$
! \end{verse}
!
! Angles are in [radian]
!
! \subsubsection{Declarations}
!
!f2py intent(in) mu
!f2py intent(in) o_m
!f2py intent(out) p
!f2py intent(out) v
    implicit none

    type(t_orb_m) :: o_m
    type(t_v3d) :: p, v
    real (kind=8) :: mu

! \subsection{Variables}
! \subsubsection{Definitions}
! \begin{verse}
! \verb|cos_e, sin_e, cos_i, sin_i| = sinus and cosines of $E$ and
!      $i$ \\
! \verb|delau| = Delaunay variables: $l, \cos(g), \sin(g), \cos(h),
!      \sin(h), L, G, H$ \\
! \verb|e| = eccentric anomaly \\
! \verb|mat| = rotation matrix \\
! \verb|q_vec, qp| = $q$ and $\dot q$ \\
! \verb|tmp| = temporary variable
! \end{verse}
!
! \subsubsection{Declarations}
    integer :: i
    real (kind=8) :: delau(8), cos_e, cos_i, e, mat(3,3), q_vec(2), qp(2), &
         sin_e, sin_i, tmp, signe, f, de, fp, fpp, fppp
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, TwoPi = 2.d0*Pi

! Computation of sinus and cosines of angles.

    signe = 1.d0
    if (o_m%a .lt. 0.d0) signe = -1.d0
    cos_i = dcos(o_m%inc)
    sin_i = dsqrt(1.d0 - cos_i**2)
    delau(2) = dcos(o_m%peri)
    delau(3) = dsin(o_m%peri)
    delau(4) = dcos(o_m%node)
    delau(5) = dsin(o_m%node)
    delau(1) = o_m%m - int(o_m%m/TwoPi)*TwoPi
    delau(6) = signe*dsqrt(mu*o_m%a*signe)
    delau(7) = abs(delau(6))*dsqrt((1.d0 - o_m%e**2)*signe)

! Rotation matrix.
! The rotation matrix is the composition of 3 matrices: $R_{xq} =
!  R_3(-h) \cdot R_1(-i) \cdot R_3(-g)$:
! \begin{displaymath}
!     R_{xq} = \left(\matrix{
!       \cos(h)\cos(g)-\frac{H}{G}\sin(h)\sin(g)&
!       -\cos(h)\sin(g)-\frac{H}{G}\sin(h)\cos(g)&
!       \sqrt{1-\frac{H^2}{G^2}}\sin(h) \cr
!       \sin(h)\cos(g)+\frac{H}{G}\cos(h)\sin(g)&
!       -\sin(h)\sin(g)+\frac{H}{G}\cos(h)\cos(g)&
!       -\sqrt{1-\frac{H^2}{G^2}}\cos(h) \cr
!       \sqrt{1-\frac{H^2}{G^2}}\sin(g)&
!       \sqrt{1-\frac{H^2}{G^2}}\cos(g)&
!       \frac{H}{G} \cr}\right);
! \end{displaymath}

    mat(1,1) = delau(4)*delau(2) - cos_i*delau(5)*delau(3)
    mat(1,2) = -delau(4)*delau(3) - cos_i*delau(5)*delau(2)
    mat(2,1) = delau(5)*delau(2) + cos_i*delau(4)*delau(3)
    mat(2,2) = -delau(5)*delau(3) + cos_i*delau(4)*delau(2)
    mat(3,1) = sin_i*delau(3)
    mat(3,2) = sin_i*delau(2)

! Eccentric anomaly.
! We solve iteratively the equation:
! \begin{displaymath}
!     E - e \sin(E) = l
! \end{displaymath}
! using the accelerated Newton's method (see Danby).

    e = delau(1) + sign(.85d0, dsin(delau(1)))*o_m%e
    i = 0
1000 continue
    sin_e = o_m%e*dsin(e)
    f = e - sin_e - delau(1)
    if (dabs(f) .gt. 1.d-14) then
       cos_e = o_m%e*dcos(e)
       fp = 1.d0 - cos_e
       fpp = sin_e
       fppp = cos_e
       de = -f/fp
       de = -f/(fp + de*fpp/2.d0)
       de = -f/(fp + de*fpp/2.d0 + de*de*fppp/6.d0)
       e = e + de
       i = i + 1
       if (i .lt. 20) goto 1000
       write (6, *) 'COORD_CART: No convergence: ', i, f
       write (6, *) mu, e, de
       write (6, *) o_m%a, o_m%e, o_m%inc
       write (6, *) o_m%node, o_m%peri, o_m%m
       stop
    end if
1100 continue

! Coordinates relative to the orbit.
! The cartisian coordinate are given by $\vec X = R_{xq} \vec q$
! and $\vec P = R_{xq} \dot{\vec q}$, where:
! \begin{eqnarray*}
!     \vec q & = & \left(\frac{L^2}{\mu}(\cos(E) - e),
!                        \frac{GL}{\mu}\sin(E), 0\right), \\
!     \dot{\vec q} & = & \frac{\mu}{L(1 - e\cos(E))}
!                  \left(-\sin(E), \frac{G}{L}\cos(E), 0\right)
! \end{eqnarray*}

    cos_e = dcos(e)
    sin_e = dsin(e)
    q_vec(1) = delau(6)**2*(cos_e - o_m%e)/mu
    q_vec(2) = delau(7)*delau(6)*sin_e/mu
    tmp = mu/(delau(6)*(1.d0 - o_m%e*cos_e))
    qp(1) = -sin_e*tmp
    qp(2) = delau(7)*cos_e*tmp/delau(6)

! Cartisian coordinates

    p%x = mat(1,1)*q_vec(1) + mat(1,2)*q_vec(2)
    p%y = mat(2,1)*q_vec(1) + mat(2,2)*q_vec(2)
    p%z = mat(3,1)*q_vec(1) + mat(3,2)*q_vec(2)
    v%x = mat(1,1)*qp(1) + mat(1,2)*qp(2)
    v%y = mat(2,1)*qp(1) + mat(2,2)*qp(2)
    v%z = mat(3,1)*qp(1) + mat(3,2)*qp(2)

    return
  end subroutine coord_cart

! \subroutine{pos\_cart}

  subroutine pos_cart (o_m, p)

! This routine transforms delaunay variables into cartisian
! variables, positions only.
!
! \subsection{Arguments}
! \subsubsection{Definitions}
! \begin{verse}
! \verb|a, e, inc, node, peri, m| = osculating elements \\
! \verb|x, y, z| = cartesian variables: $X, Y, Z$
! \end{verse}
!
! Angles are in [radian]
!
! \subsubsection{Declarations}
!
!f2py intent(in) o_m
!f2py intent(out) p
    implicit none

    type(t_orb_m) :: o_m
    type(t_v3d) :: p

! \subsection{Variables}
! \subsubsection{Definitions}
! \begin{verse}
! \verb|cos_e, sin_e, cos_i, sin_i| = sinus and cosines of $E$ and
!      $i$ \\
! \verb|delau| = Delaunay variables: $l, \cos(g), \sin(g), \cos(h),
!      \sin(h), L, G, H$ \\
! \verb|e| = eccentric anomaly \\
! \verb|mat| = rotation matrix \\
! \verb|q_vec| = $q$ \\
! \end{verse}
!
! \subsubsection{Declarations}
    integer :: i
    real (kind=8) :: delau(8), cos_e, cos_i, e, mat(3,3), q_vec(2), &
         sin_e, sin_i, signe, f, de, fp, fpp, fppp
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, TwoPi = 2.d0*Pi

! Computation of sinus and cosines of angles.

    signe = 1.d0
    if (o_m%a .lt. 0.d0) signe = -1.d0
    cos_i = dcos(o_m%inc)
    sin_i = dsqrt(1.d0 - cos_i**2)
    delau(2) = dcos(o_m%peri)
    delau(3) = dsin(o_m%peri)
    delau(4) = dcos(o_m%node)
    delau(5) = dsin(o_m%node)
    delau(1) = o_m%m - int(o_m%m/TwoPi)*TwoPi
    delau(6) = signe*dsqrt(o_m%a*signe)
    delau(7) = abs(delau(6))*dsqrt((1.d0 - o_m%e**2)*signe)

! Rotation matrix.
! The rotation matrix is the composition of 3 matrices: $R_{xq} =
!  R_3(-h) \cdot R_1(-i) \cdot R_3(-g)$:
! \begin{displaymath}
!     R_{xq} = \left(\matrix{
!       \cos(h)\cos(g)-\frac{H}{G}\sin(h)\sin(g)&
!       -\cos(h)\sin(g)-\frac{H}{G}\sin(h)\cos(g)&
!       \sqrt{1-\frac{H^2}{G^2}}\sin(h) \cr
!       \sin(h)\cos(g)+\frac{H}{G}\cos(h)\sin(g)&
!       -\sin(h)\sin(g)+\frac{H}{G}\cos(h)\cos(g)&
!       -\sqrt{1-\frac{H^2}{G^2}}\cos(h) \cr
!       \sqrt{1-\frac{H^2}{G^2}}\sin(g)&
!       \sqrt{1-\frac{H^2}{G^2}}\cos(g)&
!       \frac{H}{G} \cr}\right);
! \end{displaymath}

    mat(1,1) = delau(4)*delau(2) - cos_i*delau(5)*delau(3)
    mat(1,2) = -delau(4)*delau(3) - cos_i*delau(5)*delau(2)
    mat(2,1) = delau(5)*delau(2) + cos_i*delau(4)*delau(3)
    mat(2,2) = -delau(5)*delau(3) + cos_i*delau(4)*delau(2)
    mat(3,1) = sin_i*delau(3)
    mat(3,2) = sin_i*delau(2)

! Eccentric anomaly.
! We solve iteratively the equation:
! \begin{displaymath}
!     E - e \sin(E) = l
! \end{displaymath}
! using the accelerated Newton's method (see Danby).

    e = delau(1) + sign(.85d0, dsin(delau(1)))*o_m%e
    i = 0
1000 continue
    sin_e = o_m%e*dsin(e)
    f = e - sin_e - delau(1)
    if (dabs(f) .gt. 1.d-14) then
       cos_e = o_m%e*dcos(e)
       fp = 1.d0 - cos_e
       fpp = sin_e
       fppp = cos_e
       de = -f/fp
       de = -f/(fp + de*fpp/2.d0)
       de = -f/(fp + de*fpp/2.d0 + de*de*fppp/6.d0)
       e = e + de
       i = i + 1
       if (i .lt. 20) goto 1000
       write (6, *) 'POS_CART: No convergence: ', i, f
       write (6, *) e, de
       write (6, *) o_m%a, o_m%e, o_m%inc
       write (6, *) o_m%node, o_m%peri, o_m%m
       stop
    end if
1100 continue

! Coordinates relative to the orbit.
! The cartisian coordinate are given by $\vec X = R_{xq} \vec q$
! and $\vec P = R_{xq} \dot{\vec q}$, where:
! \begin{eqnarray*}
!     \vec q & = & \left(\frac{L^2}{\mu}(\cos(E) - e),
!                        \frac{GL}{\mu}\sin(E), 0\right), \\
!     \dot{\vec q} & = & \frac{\mu}{L(1 - e\cos(E))}
!                  \left(-\sin(E), \frac{G}{L}\cos(E), 0\right)
! \end{eqnarray*}

    cos_e = dcos(e)
    sin_e = dsin(e)
    q_vec(1) = delau(6)**2*(cos_e - o_m%e)
    q_vec(2) = delau(7)*delau(6)*sin_e

! Cartisian coordinates

    p%x = mat(1,1)*q_vec(1) + mat(1,2)*q_vec(2)
    p%y = mat(2,1)*q_vec(1) + mat(2,2)*q_vec(2)
    p%z = mat(3,1)*q_vec(1) + mat(3,2)*q_vec(2)

    return
  end subroutine pos_cart

! \subroutine{coord\_cart}

  subroutine PQ_cart (inc, node, peri, P, Q, R)

! This routine transforms delaunay variables into cartisian
! variables.
!
! ANGLES ARE GIVEN IN RADIAN !!!!
!
! \subsection{Arguments}
! \subsubsection{Definitions}
! \begin{verse}
! \verb|inc, node, peri| = osculating elements \\
! \verb|P, Q, R| = $\vec P$, $\vec Q$ and $\vec P \times \vec Q$ vectors
! \end{verse}
!
! \subsubsection{Declarations}
!
!f2py intent(in) inc
!f2py intent(in) node
!f2py intent(in) peri
!f2py intent(out) P
!f2py intent(out) Q
!f2py intent(out) R
    implicit none

    real (kind=8) :: inc, node, peri, P(3), Q(3), R(3)

! \subsection{Variables}
! \subsubsection{Definitions}
! \begin{verse}
! \verb|cos_i, sin_i| = sinus and cosines of $i$ \\
! \verb|delau| = Delaunay variables: $l, \cos(g), \sin(g), \cos(h),
!      \sin(h), L, G, H$
! \end{verse}
!
! \subsubsection{Declarations}
    integer :: i
    real (kind=8) :: delau(8), cos_i, sin_i

! Computation of sinus and cosines of angles.

    cos_i = dcos(inc)
    sin_i = dsqrt(1.d0 - cos_i**2)
    delau(2) = dcos(peri)
    delau(3) = dsin(peri)
    delau(4) = dcos(node)
    delau(5) = dsin(node)

! Rotation matrix.
! The rotation matrix is the composition of 3 matrices: $R_{xq} =
!  R_3(-h) \cdot R_1(-i) \cdot R_3(-g)$:
! \begin{displaymath}
!     R_{xq} = \left(\matrix{
!       \cos(h)\cos(g)-\frac{H}{G}\sin(h)\sin(g)&
!       -\cos(h)\sin(g)-\frac{H}{G}\sin(h)\cos(g)&
!       \sqrt{1-\frac{H^2}{G^2}}\sin(h) \cr
!       \sin(h)\cos(g)+\frac{H}{G}\cos(h)\sin(g)&
!       -\sin(h)\sin(g)+\frac{H}{G}\cos(h)\cos(g)&
!       -\sqrt{1-\frac{H^2}{G^2}}\cos(h) \cr
!       \sqrt{1-\frac{H^2}{G^2}}\sin(g)&
!       \sqrt{1-\frac{H^2}{G^2}}\cos(g)&
!       \frac{H}{G} \cr}\right);
! \end{displaymath}

    P(1) = delau(4)*delau(2) - cos_i*delau(5)*delau(3)
    Q(1) = -delau(4)*delau(3) - cos_i*delau(5)*delau(2)
    R(1) = sin_i*delau(5)
    P(2) = delau(5)*delau(2) + cos_i*delau(4)*delau(3)
    Q(2) = -delau(5)*delau(3) + cos_i*delau(4)*delau(2)
    R(2) = -sin_i*delau(4)
    P(3) = sin_i*delau(3)
    Q(3) = sin_i*delau(2)
    R(3) = cos_i

    return
  end subroutine PQ_cart

! \subroutine{osc\_el}

  subroutine osc_el (mu, p, v, o_m)

! This routine transforms cartisian variables into delaunay
! variables.
!
! \subsection{Arguments}
! \subsubsection{Definitions}
! \begin{verse}
! \verb|a, e, inc, node, peri, m| = osculating elements \\
! \verb|cart| = cartesian variables: $X, Y, Z, Px, Py, Pz$
! \end{verse}
!
! \subsubsection{Declarations}
!f2py intent(in) mu
!f2py intent(in) p
!f2py intent(in) v
!f2py intent(out) o_m
    implicit none

    type(t_v3d) :: p, v
    type(t_orb_m) :: o_m
    real (kind=8) :: mu

! \subsection{Variables}
! \subsubsection{Definitions}
! \begin{verse}
! \verb|cos_i, sin_i| = sinus and cosines of $i$ \\
! \verb|delau| = Delaunay variables:  $l, \cos(g), \sin(g), \cos(h),
!      \sin(h), L, G, H$ \\
! \verb|e| = eccentric anomaly \\
! \verb|f| = $f$ true anomaly \\
! \verb|g| = $g$ argument of pericenter \\
! \verb|h_vec| = $\vec h = \vec X \times \vec P$ \\
! \verb|p_vec| = $\vec p = -\mu \frac{\vec X}{r}
!      - \vec h \times \vec P$ \\
! \verb|r| = radial distance \\
! \verb|tmp1, tmp2| = temporary variables \\
! \verb|v2| = velocity squared
! \end{verse}
!
! \subsubsection{Declarations}
    real (kind=8) :: delau(8), e, f, h_vec(3), p_vec(3), &
         r, tmp1, tmp2, v2, signe

! Computation of angular momentum and eccentricity vector.
! \begin{eqnarray*}
!     \vec h & = & \vec X \times \vec P, \\
!     \vec p & = & -\mu \frac{\vec X}{|\vec X|}
!                  - \vec h \times \vec P
! \end{eqnarray*}

    h_vec(1) = p%y*v%z - p%z*v%y
    h_vec(2) = p%z*v%x - p%x*v%z
    h_vec(3) = p%x*v%y - p%y*v%x
    r = 1.d0/dsqrt(p%x**2 + p%y**2 + p%z**2)
    p_vec(1) = -mu*p%x*r - h_vec(2)*v%z + h_vec(3)*v%y
    p_vec(2) = -mu*p%y*r - h_vec(3)*v%x + h_vec(1)*v%z
    p_vec(3) = -mu*p%z*r - h_vec(1)*v%y + h_vec(2)*v%x

! Computation of momenta.
! \begin{eqnarray*}
!     L & = & \mu\sqrt{\frac{1}
!                      {\frac{2\mu}{|\vec X|}-|\vec P|^2}}, \\
!     G & = & |\vec h|, \\
!     H & = &  h_z
! \end{eqnarray*}

    v2 = v%x**2 + v%y**2 + v%z**2
    tmp1 = 2.d0*mu*r - v2
    signe = 1.d0
    if (tmp1 .lt. 0.d0) signe = -1.d0
    delau(6) = signe*mu/dsqrt(signe*tmp1)
    delau(7) = dsqrt(h_vec(1)**2 + h_vec(2)**2 + h_vec(3)**2)
    delau(8) = h_vec(3)

    if ((p%z .eq. 0.d0) .and. (v%z .eq. 0.d0)) then
       delau(4) = 1.d0
       delau(5) = 0.d0
       tmp1 = 1.d0/dsqrt(p_vec(1)**2 + p_vec(2)**2)
       delau(2) = p_vec(1)*tmp1
       delau(3) = p_vec(2)*tmp1
    else

! Longitude of node.
! \begin{eqnarray*}
!     \cos(h) & = & -\frac{h_y}{\sqrt{h_x^2 + h_y^2}}, \\
!     \sin(h) & = & \frac{h_x}{\sqrt{h_x^2 + h_y^2}}
! \end{eqnarray*}

       tmp1 = 1.d0/dsqrt(h_vec(1)**2 + h_vec(2)**2)
       delau(4) = -h_vec(2)*tmp1
       delau(5) = h_vec(1)*tmp1

! Argument of pericenter.
! Let us call $\vec N$ the vector derived from $\vec h$, pointing
! to the ascending node:
! \begin{displaymath}
!     \vec N = \left(-h_y, h_x, 0\right)
! \end{displaymath}
! From this, we get:
! \begin{eqnarray*}
!     \cos(g) & = & \frac{\vec N \cdot \vec p}{|\vec N||\vec p|}, \\
!     \sin(g) & = & \frac{\vec N \times \vec p}{|\vec N||\vec p|}
!                   \cdot \frac{\vec h}{|\vec h|}
! \end{eqnarray*}

       tmp2 = 1.d0/dsqrt(p_vec(1)**2 + p_vec(2)**2 + p_vec(3)**2)
       delau(2) = (h_vec(1)*p_vec(2) - h_vec(2)*p_vec(1))*tmp1*tmp2
       delau(3) = ((h_vec(1)**2 + h_vec(2)**2)*p_vec(3) &
            - h_vec(3)*(h_vec(1)*p_vec(1) + h_vec(2)*p_vec(2))) &
            *tmp1*tmp2/delau(7)
    end if

! Mean anomaly
! We define $\vec X_{orb} = R_1(i) \cdot R_3(h) \vec X$. It turns
! out that $\vec X_{orb} = (r \cos(g+f), r \sin(g+f), 0)$. Hence:
! \begin{eqnarray*}
!     \cos(g+f) & = & \cos(h) X + \sin(h) Y, \\
!     \sin(g+f) & = & \cos(i) \left(\cos(h) Y - \sin(h) X\right)
!                   + \sin(i) Z
! \end{eqnarray*}
! Furthermore, we have the relation:
! \begin{displaymath}
!     \tan(\frac{E}{2}) = \sqrt{\frac{1 - e}{1 + e}}
!                         \tan(\frac{f}{2})
! \end{displaymath}
! and finally:
! \begin{displaymath}
!     l = E - e \sin(E)
! \end{displaymath}

    o_m%e = dsqrt(dmax1(1.d0 - signe*(delau(7)/delau(6))**2, 0.d0))
    tmp1 = (p%x*p_vec(1) + p%y*p_vec(2) + p%z*p_vec(3))
    tmp2 = ((p%z*p_vec(2) - p%y*p_vec(3))*h_vec(1) &
         + (p%x*p_vec(3) - p%z*p_vec(1))*h_vec(2) &
         + (p%y*p_vec(1) - p%x*p_vec(2))*h_vec(3)) &
         /delau(7)
    f = datan2(tmp2, tmp1)
    e = 2.d0*datan(dsqrt(signe*(1.d0 - o_m%e)/(1.d0 + o_m%e))*dtan(f/2.d0))
    o_m%m = e - o_m%e*dsin(e)
    o_m%node = datan2(delau(5), delau(4))
    o_m%peri = datan2(delau(3), delau(2))
    o_m%a = signe*delau(6)**2/mu
    o_m%inc = dacos(dmax1(dmin1(delau(8)/delau(7),1.d0),-1.d0))

    return
  end subroutine osc_el

end module elemutils
