
c \subroutine{coord\_cart}

      subroutine coord_cart (mu, a, ecc, inc, capo, smallo, capm,
     $  x, y, z, vx, vy, vz)

c This routine transforms delaunay variables into cartisian
c variables.
c
c ANGLES ARE GIVEN IN RADIAN !!!!
c
c \subsection{Arguments}
c \subsubsection{Definitions}
c \begin{verse}
c \verb|a, e, inc, capo, smallo, capm| = osculating elements \\
c \verb|x, y, z, vx, vy, vz| = cartesian variables: $X, Y, Z, Px, Py, Pz$
c \end{verse}
c
c \subsubsection{Declarations}
c
Cf2py intent(in) mu
Cf2py intent(in) a
Cf2py intent(in) ecc
Cf2py intent(in) inc
Cf2py intent(in) capo
Cf2py intent(in) smallo
Cf2py intent(in) capm
Cf2py intent(out) x
Cf2py intent(out) y
Cf2py intent(out) z
Cf2py intent(out) vx
Cf2py intent(out) vy
Cf2py intent(out) vz

      implicit none

      real*8
     $  a, ecc, inc, capo, smallo, capm, mu, x, y, z, vx, vy, vz

c \subsection{Variables}
c \subsubsection{Definitions}
c \begin{verse}
c \verb|cos_e, sin_e, cos_i, sin_i| = sinus and cosines of $E$ and
c      $i$ \\
c \verb|delau| = Delaunay variables: $l, \cos(g), \sin(g), \cos(h),
c      \sin(h), L, G, H$ \\
c \verb|e| = eccentric anomaly \\
c \verb|mat| = rotation matrix \\
c \verb|q_vec, qp| = $q$ and $\dot q$ \\
c \verb|tmp| = temporary variable
c \end{verse}
c
c \subsubsection{Declarations}

      integer*4
     $  i

      real*8
     $  delau(8), cos_e, cos_i, e, mat(3,3), q_vec(2), qp(2),
     $  sin_e, sin_i, tmp, signe, Pi, TwoPi, f, de, fp, fpp, fppp

      parameter
     $  (Pi = 3.141592653589793238d0, TwoPi = 2.d0*Pi)

c Computation of sinus and cosines of angles.

      signe = 1.d0
      if (a .lt. 0.d0) signe = -1.d0
      cos_i = dcos(inc)
      sin_i = dsqrt(1.d0 - cos_i**2)
      delau(2) = dcos(smallo)
      delau(3) = dsin(smallo)
      delau(4) = dcos(capo)
      delau(5) = dsin(capo)
      delau(1) = capm - int(capm/TwoPi)*TwoPi
      delau(6) = signe*dsqrt(mu*a*signe)
      delau(7) = abs(delau(6))*dsqrt((1.d0 - ecc**2)*signe)

c Rotation matrix.
c The rotation matrix is the composition of 3 matrices: $R_{xq} =
c  R_3(-h) \cdot R_1(-i) \cdot R_3(-g)$:
c \begin{displaymath}
c     R_{xq} = \left(\matrix{
c       \cos(h)\cos(g)-\frac{H}{G}\sin(h)\sin(g)&
c       -\cos(h)\sin(g)-\frac{H}{G}\sin(h)\cos(g)&
c       \sqrt{1-\frac{H^2}{G^2}}\sin(h) \cr
c       \sin(h)\cos(g)+\frac{H}{G}\cos(h)\sin(g)&
c       -\sin(h)\sin(g)+\frac{H}{G}\cos(h)\cos(g)&
c       -\sqrt{1-\frac{H^2}{G^2}}\cos(h) \cr
c       \sqrt{1-\frac{H^2}{G^2}}\sin(g)&
c       \sqrt{1-\frac{H^2}{G^2}}\cos(g)&
c       \frac{H}{G} \cr}\right);
c \end{displaymath}

      mat(1,1) = delau(4)*delau(2) - cos_i*delau(5)*delau(3)
      mat(1,2) = -delau(4)*delau(3) - cos_i*delau(5)*delau(2)
      mat(2,1) = delau(5)*delau(2) + cos_i*delau(4)*delau(3)
      mat(2,2) = -delau(5)*delau(3) + cos_i*delau(4)*delau(2)
      mat(3,1) = sin_i*delau(3)
      mat(3,2) = sin_i*delau(2)

c Eccentric anomaly.
c We solve iteratively the equation:
c \begin{displaymath}
c     E - e \sin(E) = l
c \end{displaymath}
c using the accelerated Newton's method (see Danby).

      e = delau(1) + sign(.85d0, dsin(delau(1)))*ecc
      i = 0
 1000 continue
      sin_e = ecc*dsin(e)
      f = e - sin_e - delau(1)
      if (dabs(f) .gt. 1.d-14) then
         cos_e = ecc*dcos(e)
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
         write (6, *) a, ecc, inc
         write (6, *) capo, smallo, capm
         stop
      end if
 1100 continue

c Coordinates relative to the orbit.
c The cartisian coordinate are given by $\vec X = R_{xq} \vec q$
c and $\vec P = R_{xq} \dot{\vec q}$, where:
c \begin{eqnarray*}
c     \vec q & = & \left(\frac{L^2}{\mu}(\cos(E) - e),
c                        \frac{GL}{\mu}\sin(E), 0\right), \\
c     \dot{\vec q} & = & \frac{\mu}{L(1 - e\cos(E))}
c                  \left(-\sin(E), \frac{G}{L}\cos(E), 0\right)
c \end{eqnarray*}

      cos_e = dcos(e)
      sin_e = dsin(e)
      q_vec(1) = delau(6)**2*(cos_e - ecc)/mu
      q_vec(2) = delau(7)*delau(6)*sin_e/mu
      tmp = mu/(delau(6)*(1.d0 - ecc*cos_e))
      qp(1) = -sin_e*tmp
      qp(2) = delau(7)*cos_e*tmp/delau(6)

c Cartisian coordinates

      x = mat(1,1)*q_vec(1) + mat(1,2)*q_vec(2)
      y = mat(2,1)*q_vec(1) + mat(2,2)*q_vec(2)
      z = mat(3,1)*q_vec(1) + mat(3,2)*q_vec(2)
      vx = mat(1,1)*qp(1) + mat(1,2)*qp(2)
      vy = mat(2,1)*qp(1) + mat(2,2)*qp(2)
      vz = mat(3,1)*qp(1) + mat(3,2)*qp(2)

      return
      end

c \subroutine{pos\_cart}

      subroutine pos_cart (a, ecc, inc, capo, smallo, capm,
     $  x, y, z)

c This routine transforms delaunay variables into cartisian
c variables, positions only.
c
c \subsection{Arguments}
c \subsubsection{Definitions}
c \begin{verse}
c \verb|a, e, inc, capo, smallo, capm| = osculating elements \\
c \verb|x, y, z| = cartesian variables: $X, Y, Z$
c \end{verse}
c
c \subsubsection{Declarations}
c
Cf2py intent(in) a
Cf2py intent(in) ecc
Cf2py intent(in) inc
Cf2py intent(in) capo
Cf2py intent(in) smallo
Cf2py intent(in) capm
Cf2py intent(out) x
Cf2py intent(out) y
Cf2py intent(out) z

      implicit none

      real*8
     $  a, ecc, inc, capo, smallo, capm, x, y, z

c \subsection{Variables}
c \subsubsection{Definitions}
c \begin{verse}
c \verb|cos_e, sin_e, cos_i, sin_i| = sinus and cosines of $E$ and
c      $i$ \\
c \verb|delau| = Delaunay variables: $l, \cos(g), \sin(g), \cos(h),
c      \sin(h), L, G, H$ \\
c \verb|e| = eccentric anomaly \\
c \verb|mat| = rotation matrix \\
c \verb|q_vec| = $q$ \\
c \end{verse}
c
c \subsubsection{Declarations}

      integer*4
     $  i

      real*8
     $  delau(8), cos_e, cos_i, e, mat(3,3), q_vec(2),
     $  sin_e, sin_i, signe, Pi, TwoPi, f, de, fp, fpp, fppp

      parameter
     $  (Pi = 3.141592653589793238d0, TwoPi = 2.d0*Pi)

c Computation of sinus and cosines of angles.

      signe = 1.d0
      if (a .lt. 0.d0) signe = -1.d0
      cos_i = dcos(inc)
      sin_i = dsqrt(1.d0 - cos_i**2)
      delau(2) = dcos(smallo)
      delau(3) = dsin(smallo)
      delau(4) = dcos(capo)
      delau(5) = dsin(capo)
      delau(1) = capm - int(capm/TwoPi)*TwoPi
      delau(6) = signe*dsqrt(a*signe)
      delau(7) = abs(delau(6))*dsqrt((1.d0 - ecc**2)*signe)

c Rotation matrix.
c The rotation matrix is the composition of 3 matrices: $R_{xq} =
c  R_3(-h) \cdot R_1(-i) \cdot R_3(-g)$:
c \begin{displaymath}
c     R_{xq} = \left(\matrix{
c       \cos(h)\cos(g)-\frac{H}{G}\sin(h)\sin(g)&
c       -\cos(h)\sin(g)-\frac{H}{G}\sin(h)\cos(g)&
c       \sqrt{1-\frac{H^2}{G^2}}\sin(h) \cr
c       \sin(h)\cos(g)+\frac{H}{G}\cos(h)\sin(g)&
c       -\sin(h)\sin(g)+\frac{H}{G}\cos(h)\cos(g)&
c       -\sqrt{1-\frac{H^2}{G^2}}\cos(h) \cr
c       \sqrt{1-\frac{H^2}{G^2}}\sin(g)&
c       \sqrt{1-\frac{H^2}{G^2}}\cos(g)&
c       \frac{H}{G} \cr}\right);
c \end{displaymath}

      mat(1,1) = delau(4)*delau(2) - cos_i*delau(5)*delau(3)
      mat(1,2) = -delau(4)*delau(3) - cos_i*delau(5)*delau(2)
      mat(2,1) = delau(5)*delau(2) + cos_i*delau(4)*delau(3)
      mat(2,2) = -delau(5)*delau(3) + cos_i*delau(4)*delau(2)
      mat(3,1) = sin_i*delau(3)
      mat(3,2) = sin_i*delau(2)

c Eccentric anomaly.
c We solve iteratively the equation:
c \begin{displaymath}
c     E - e \sin(E) = l
c \end{displaymath}
c using the accelerated Newton's method (see Danby).

      e = delau(1) + sign(.85d0, dsin(delau(1)))*ecc
      i = 0
 1000 continue
      sin_e = ecc*dsin(e)
      f = e - sin_e - delau(1)
      if (dabs(f) .gt. 1.d-14) then
         cos_e = ecc*dcos(e)
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
         write (6, *) a, ecc, inc
         write (6, *) capo, smallo, capm
         stop
      end if
 1100 continue

c Coordinates relative to the orbit.
c The cartisian coordinate are given by $\vec X = R_{xq} \vec q$
c and $\vec P = R_{xq} \dot{\vec q}$, where:
c \begin{eqnarray*}
c     \vec q & = & \left(\frac{L^2}{\mu}(\cos(E) - e),
c                        \frac{GL}{\mu}\sin(E), 0\right), \\
c     \dot{\vec q} & = & \frac{\mu}{L(1 - e\cos(E))}
c                  \left(-\sin(E), \frac{G}{L}\cos(E), 0\right)
c \end{eqnarray*}

      cos_e = dcos(e)
      sin_e = dsin(e)
      q_vec(1) = delau(6)**2*(cos_e - ecc)
      q_vec(2) = delau(7)*delau(6)*sin_e

c Cartisian coordinates

      x = mat(1,1)*q_vec(1) + mat(1,2)*q_vec(2)
      y = mat(2,1)*q_vec(1) + mat(2,2)*q_vec(2)
      z = mat(3,1)*q_vec(1) + mat(3,2)*q_vec(2)

      return
      end

c \subroutine{osc\_el}

      subroutine osc_el (mu, x, y, z, vx, vy, vz,
     $  a, ecc, inc, capo, smallo, capm)

c This routine transforms cartisian variables into delaunay
c variables.
c
c \subsection{Arguments}
c \subsubsection{Definitions}
c \begin{verse}
c \verb|a, e, inc, capo, smallo, capm| = osculating elements \\
c \verb|cart| = cartesian variables: $X, Y, Z, Px, Py, Pz$
c \end{verse}
c
c \subsubsection{Declarations}

      implicit none

      real*8
     $  a, ecc, inc, capo, smallo, capm, mu, x, y, z, vx, vy, vz

c \subsection{Variables}
c \subsubsection{Definitions}
c \begin{verse}
c \verb|cos_i, sin_i| = sinus and cosines of $i$ \\
c \verb|delau| = Delaunay variables:  $l, \cos(g), \sin(g), \cos(h),
c      \sin(h), L, G, H$ \\
c \verb|e| = eccentric anomaly \\
c \verb|f| = $f$ true anomaly \\
c \verb|g| = $g$ argument of pericenter \\
c \verb|h_vec| = $\vec h = \vec X \times \vec P$ \\
c \verb|p_vec| = $\vec p = -\mu \frac{\vec X}{r}
c      - \vec h \times \vec P$ \\
c \verb|r| = radial distance \\
c \verb|tmp1, tmp2| = temporary variables \\
c \verb|v2| = velocity squared
c \end{verse}
c
c \subsubsection{Declarations}

      real*8
     $  delau(8), e, f, h_vec(3), p_vec(3),
     $  r, tmp1, tmp2, v2, signe, cart(6)

c Computation of angular momentum and eccentricity vector.
c \begin{eqnarray*}
c     \vec h & = & \vec X \times \vec P, \\
c     \vec p & = & -\mu \frac{\vec X}{|\vec X|}
c                  - \vec h \times \vec P
c \end{eqnarray*}

      cart(1) = x
      cart(2) = y
      cart(3) = z
      cart(4) = vx
      cart(5) = vy
      cart(6) = vz
      h_vec(1) = cart(2)*cart(6) - cart(3)*cart(5)
      h_vec(2) = cart(3)*cart(4) - cart(1)*cart(6)
      h_vec(3) = cart(1)*cart(5) - cart(2)*cart(4)
      r = 1.d0/dsqrt(cart(1)**2 + cart(2)**2 + cart(3)**2)
      p_vec(1) = -mu*cart(1)*r - h_vec(2)*cart(6) + h_vec(3)*cart(5)
      p_vec(2) = -mu*cart(2)*r - h_vec(3)*cart(4) + h_vec(1)*cart(6)
      p_vec(3) = -mu*cart(3)*r - h_vec(1)*cart(5) + h_vec(2)
     $  *cart(4)

c Computation of momenta.
c \begin{eqnarray*}
c     L & = & \mu\sqrt{\frac{1}
c                      {\frac{2\mu}{|\vec X|}-|\vec P|^2}}, \\
c     G & = & |\vec h|, \\
c     H & = &  h_z
c \end{eqnarray*}

      v2 = cart(4)**2 + cart(5)**2 + cart(6)**2
      tmp1 = 2.d0*mu*r - v2
      signe = 1.d0
      if (tmp1 .lt. 0.d0) signe = -1.d0
      delau(6) = signe*mu/dsqrt(signe*tmp1)
      delau(7) = dsqrt(h_vec(1)**2 + h_vec(2)**2 + h_vec(3)**2)
      delau(8) = h_vec(3)

      if ((cart(3) .eq. 0.d0) .and. (cart(6) .eq. 0.d0)) then
         delau(4) = 1.d0
         delau(5) = 0.d0
         tmp1 = 1.d0/dsqrt(p_vec(1)**2 + p_vec(2)**2)
         delau(2) = p_vec(1)*tmp1
         delau(3) = p_vec(2)*tmp1
      else

c Longitude of node.
c \begin{eqnarray*}
c     \cos(h) & = & -\frac{h_y}{\sqrt{h_x^2 + h_y^2}}, \\
c     \sin(h) & = & \frac{h_x}{\sqrt{h_x^2 + h_y^2}}
c \end{eqnarray*}

         tmp1 = 1.d0/dsqrt(h_vec(1)**2 + h_vec(2)**2)
         delau(4) = -h_vec(2)*tmp1
         delau(5) = h_vec(1)*tmp1

c Argument of pericenter.
c Let us call $\vec N$ the vector derived from $\vec h$, pointing
c to the ascending node:
c \begin{displaymath}
c     \vec N = \left(-h_y, h_x, 0\right)
c \end{displaymath}
c From this, we get:
c \begin{eqnarray*}
c     \cos(g) & = & \frac{\vec N \cdot \vec p}{|\vec N||\vec p|}, \\
c     \sin(g) & = & \frac{\vec N \times \vec p}{|\vec N||\vec p|}
c                   \cdot \frac{\vec h}{|\vec h|}
c \end{eqnarray*}

         tmp2 = 1.d0/dsqrt(p_vec(1)**2 + p_vec(2)**2 + p_vec(3)**2
     $     )
         delau(2) = (h_vec(1)*p_vec(2) - h_vec(2)*p_vec(1))*tmp1
     $     *tmp2
         delau(3) = ((h_vec(1)**2 + h_vec(2)**2)*p_vec(3)
     $     - h_vec(3)*(h_vec(1)*p_vec(1) + h_vec(2)*p_vec(2)))
     $     *tmp1*tmp2/delau(7)
      end if

c Mean anomaly
c We define $\vec X_{orb} = R_1(i) \cdot R_3(h) \vec X$. It turns
c out that $\vec X_{orb} = (r \cos(g+f), r \sin(g+f), 0)$. Hence:
c \begin{eqnarray*}
c     \cos(g+f) & = & \cos(h) X + \sin(h) Y, \\
c     \sin(g+f) & = & \cos(i) \left(\cos(h) Y - \sin(h) X\right)
c                   + \sin(i) Z
c \end{eqnarray*}
c Furthermore, we have the relation:
c \begin{displaymath}
c     \tan(\frac{E}{2}) = \sqrt{\frac{1 - e}{1 + e}}
c                         \tan(\frac{f}{2})
c \end{displaymath}
c and finally:
c \begin{displaymath}
c     l = E - e \sin(E)
c \end{displaymath}

      ecc = dsqrt(dmax1(1.d0 - signe*(delau(7)/delau(6))**2, 0.d0))
      tmp1 = (cart(1)*p_vec(1) + cart(2)*p_vec(2) + cart(3)*p_vec(3))
      tmp2 = ((cart(3)*p_vec(2) - cart(2)*p_vec(3))*h_vec(1)
     $  + (cart(1)*p_vec(3) - cart(3)*p_vec(1))*h_vec(2)
     $  + (cart(2)*p_vec(1) - cart(1)*p_vec(2))*h_vec(3))
     $  /delau(7)
      f = datan2(tmp2, tmp1)
      e = 2.d0*datan(dsqrt(signe*(1.d0 - ecc)/(1.d0 + ecc))
     $  *dtan(f/2.d0))
      capm = e - ecc*dsin(e)
      capo = datan2(delau(5), delau(4))
      smallo = datan2(delau(3), delau(2))
      a = signe*delau(6)**2/mu
      inc = dacos(dmax1(dmin1(delau(8)/delau(7),1.d0),-1.d0))

      return
      end
