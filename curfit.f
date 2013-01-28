      subroutine curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,
     * wrk,lwrk,iwrk,ier)
c  given the set of data points (x(i),y(i)) and the set of positive
c  numbers w(i),i=1,2,...,m,subroutine curfit determines a smooth spline
c  approximation of degree k on the interval xb <= x <= xe.
c  if iopt=-1 curfit calculates the weighted least-squares spline
c  according to a given set of knots.
c  if iopt>=0 the number of knots of the spline s(x) and the position
c  t(j),j=1,2,...,n is chosen automatically by the routine. the smooth-
c  ness of s(x) is then achieved by minimalizing the discontinuity
c  jumps of the k-th derivative of s(x) at the knots t(j),j=k+2,k+3,...,
c  n-k-1. the amount of smoothness is determined by the condition that
c  f(p)=sum((w(i)*(y(i)-s(x(i))))**2) be <= s, with s a given non-
c  negative constant, called the smoothing factor.
c  the fit s(x) is given in the b-spline representation (b-spline coef-
c  ficients c(j),j=1,2,...,n-k-1) and can be evaluated by means of
c  subroutine splev.
c
c  calling sequence:
c     call curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,
c    * lwrk,iwrk,ier)
c
c  parameters:
c   iopt  : integer flag. on entry iopt must specify whether a weighted
c           least-squares spline (iopt=-1) or a smoothing spline (iopt=
c           0 or 1) must be determined. if iopt=0 the routine will start
c           with an initial set of knots t(i)=xb, t(i+k+1)=xe, i=1,2,...
c           k+1. if iopt=1 the routine will continue with the knots
c           found at the last call of the routine.
c           attention: a call with iopt=1 must always be immediately
c           preceded by another call with iopt=1 or iopt=0.
c           unchanged on exit.
c   m     : integer. on entry m must specify the number of data points.
c           m > k. unchanged on exit.
c   x     : real array of dimension at least (m). before entry, x(i)
c           must be set to the i-th value of the independent variable x,
c           for i=1,2,...,m. these values must be supplied in strictly
c           ascending order. unchanged on exit.
c   y     : real array of dimension at least (m). before entry, y(i)
c           must be set to the i-th value of the dependent variable y,
c           for i=1,2,...,m. unchanged on exit.
c   w     : real array of dimension at least (m). before entry, w(i)
c           must be set to the i-th value in the set of weights. the
c           w(i) must be strictly positive. unchanged on exit.
c           see also further comments.
c   xb,xe : real values. on entry xb and xe must specify the boundaries
c           of the approximation interval. xb<=x(1), xe>=x(m).
c           unchanged on exit.
c   k     : integer. on entry k must specify the degree of the spline.
c           1<=k<=5. it is recommended to use cubic splines (k=3).
c           the user is strongly dissuaded from choosing k even,together
c           with a small s-value. unchanged on exit.
c   s     : real.on entry (in case iopt>=0) s must specify the smoothing
c           factor. s >=0. unchanged on exit.
c           for advice on the choice of s see further comments.
c   nest  : integer. on entry nest must contain an over-estimate of the
c           total number of knots of the spline returned, to indicate
c           the storage space available to the routine. nest >=2*k+2.
c           in most practical situation nest=m/2 will be sufficient.
c           always large enough is  nest=m+k+1, the number of knots
c           needed for interpolation (s=0). unchanged on exit.
c   n     : integer.
c           unless ier =10 (in case iopt >=0), n will contain the
c           total number of knots of the spline approximation returned.
c           if the computation mode iopt=1 is used this value of n
c           should be left unchanged between subsequent calls.
c           in case iopt=-1, the value of n must be specified on entry.
c   t     : real array of dimension at least (nest).
c           on succesful exit, this array will contain the knots of the
c           spline,i.e. the position of the interior knots t(k+2),t(k+3)
c           ...,t(n-k-1) as well as the position of the additional knots
c           t(1)=t(2)=...=t(k+1)=xb and t(n-k)=...=t(n)=xe needed for
c           the b-spline representation.
c           if the computation mode iopt=1 is used, the values of t(1),
c           t(2),...,t(n) should be left unchanged between subsequent
c           calls. if the computation mode iopt=-1 is used, the values
c           t(k+2),...,t(n-k-1) must be supplied by the user, before
c           entry. see also the restrictions (ier=10).
c   c     : real array of dimension at least (nest).
c           on succesful exit, this array will contain the coefficients
c           c(1),c(2),..,c(n-k-1) in the b-spline representation of s(x)
c   fp    : real. unless ier=10, fp contains the weighted sum of
c           squared residuals of the spline approximation returned.
c   wrk   : real array of dimension at least (m*(k+1)+nest*(7+3*k)).
c           used as working space. if the computation mode iopt=1 is
c           used, the values wrk(1),...,wrk(n) should be left unchanged
c           between subsequent calls.
c   lwrk  : integer. on entry,lwrk must specify the actual dimension of
c           the array wrk as declared in the calling (sub)program.lwrk
c           must not be too small (see wrk). unchanged on exit.
c   iwrk  : integer array of dimension at least (nest).
c           used as working space. if the computation mode iopt=1 is
c           used,the values iwrk(1),...,iwrk(n) should be left unchanged
c           between subsequent calls.
c   ier   : integer. unless the routine detects an error, ier contains a
c           not*(7+3*k))-positive value on exit, i.e.
c    ier=0  : normal return. the spline returned has a residual sum of
c             squares fp such that abs(fp-s)/s <= tol with tol a relat-
c             ive tolerance set to 0.001 by the program.
c    ier=-1 : normal return. the spline returned is an interpolating
c             spline (fp=0).
c    ier=-2 : normal return. the spline returned is the weighted least-
c             squares polynomial of degree k. in this extreme case fp
c             gives the upper bound fp0 for the smoothing factor s.
c    ier=1  : error. the required storage space exceeds the available
c             storage space, as specified by the parameter nest.
c             probably causes : nest too small. if nest is already
c             large (say nest > m/2), it may also indicate that s is
c             too small
c             the approximation returned is the weighted least-squares
c             spline according to the knots t(1),t(2),...,t(n). (n=nest)
c             the parameter fp gives the corresponding weighted sum of
c             squared residuals (fp>s).
c    ier=2  : error. a theoretically impossible result was found during
c             the iteration proces for finding a smoothing spline with
c             fp = s. probably causes : s too small.
c             there is an approximation returned but the corresponding
c             weighted sum of squared residuals does not satisfy the
c             condition abs(fp-s)/s < tol.
c    ier=3  : error. the maximal number of iterations maxit (set to 20
c             by the program) allowed for finding a smoothing spline
c             with fp=s has been reached. probably causes : s too small
c             there is an approximation returned but the corresponding
c             weighted sum of squared residuals does not satisfy the
c             condition abs(fp-s)/s < tol.
c    ier=10 : error. on entry, the input data are controlled on validity
c             the following restrictions must be satisfied.
c             -1<=iopt<=1, 1<=k<=5, m>k, nest>2*k+2, w(i)>0,i=1,2,...,m
c             xb<=x(1)<x(2)<...<x(m)<=xe, lwrk>=(k+1)*m+nest*(7+3*k)
c             if iopt=-1: 2*k+2<=n<=min(nest,m+k+1)
c                         xb<t(k+2)<t(k+3)<...<t(n-k-1)<xe
c                       the schoenberg-whitney conditions, i.e. there
c                       must be a subset of data points xx(j) such that
c                         t(j) < xx(j) < t(j+k+1), j=1,2,...,n-k-1
c             if iopt>=0: s>=0
c                         if s=0 : nest >= m+k+1
c             if one of these conditions is found to be violated,control
c             is immediately repassed to the calling program. in that
c             case there is no approximation returned.
c
c  further comments:
c   by means of the parameter s, the user can control the tradeoff
c   between closeness of fit and smoothness of fit of the approximation.
c   if s is too large, the spline will be too smooth and signal will be
c   lost ; if s is too small the spline will pick up too much noise. in
c   the extreme cases the program will return an interpolating spline if
c   s=0 and the weighted least-squares polynomial of degree k if s is
c   very large. between these extremes, a properly chosen s will result
c   in a good compromise between closeness of fit and smoothness of fit.
c   to decide whether an approximation, corresponding to a certain s is
c   satisfactory the user is highly recommended to inspect the fits
c   graphically.
c   recommended values for s depend on the weights w(i). if these are
c   taken as 1/d(i) with d(i) an estimate of the standard deviation of
c   y(i), a good s-value should be found in the range (m-sqrt(2*m),m+
c   sqrt(2*m)). if nothing is known about the statistical error in y(i)
c   each w(i) can be set equal to one and s determined by trial and
c   error, taking account of the comments above. the best is then to
c   start with a very large value of s ( to determine the least-squares
c   polynomial and the corresponding upper bound fp0 for s) and then to
c   progressively decrease the value of s ( say by a factor 10 in the
c   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
c   approximation shows more detail) to obtain closer fits.
c   to economize the search for a good s-value the program provides with
c   different modes of computation. at the first call of the routine, or
c   whenever he wants to restart with the initial set of knots the user
c   must set iopt=0.
c   if iopt=1 the program will continue with the set of knots found at
c   the last call of the routine. this will save a lot of computation
c   time if curfit is called repeatedly for different values of s.
c   the number of knots of the spline returned and their location will
c   depend on the value of s and on the complexity of the shape of the
c   function underlying the data. but, if the computation mode iopt=1
c   is used, the knots returned may also depend on the s-values at
c   previous calls (if these were smaller). therefore, if after a number
c   of trials with different s-values and iopt=1, the user can finally
c   accept a fit as satisfactory, it may be worthwhile for him to call
c   curfit once more with the selected value for s but now with iopt=0.
c   indeed, curfit may then return an approximation of the same quality
c   of fit but with fewer knots and therefore better if data reduction
c   is also an important objective for the user.
c
c  other subroutines required:
c    fpback,fpbspl,fpchec,fpcurf,fpdisc,fpgivs,fpknot,fprati,fprota
c
c  references:
c   dierckx p. : an algorithm for smoothing, differentiation and integ-
c                ration of experimental data using spline functions,
c                j.comp.appl.maths 1 (1975) 165-184.
c   dierckx p. : a fast algorithm for smoothing data on a rectangular
c                grid while using spline functions, siam j.numer.anal.
c                19 (1982) 1286-1304.
c   dierckx p. : an improved algorithm for curve fitting with spline
c                functions, report tw54, dept. computer science,k.u.
c                leuven, 1981.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author:
c    p.dierckx
c    dept. computer science, k.u. leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  creation date : may 1979
c  latest update : march 1987
c
c  ..
c  ..scalar arguments..
      real*8 xb,xe,s,fp
      integer iopt,m,k,nest,n,lwrk,ier
c  ..array arguments..
      real*8 x(m),y(m),w(m),t(nest),c(nest),wrk(lwrk)
      integer iwrk(nest)
c  ..local scalars..
      real*8 tol
      integer i,ia,ib,ifp,ig,iq,iz,j,k1,k2,lwest,maxit,nmin
c  ..
c  we set up the parameters tol and maxit
      maxit = 20
      tol = 0.1e-02
c  before starting computations a data check is made. if the input data
c  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(k.le.0 .or. k.gt.5) go to 50
      k1 = k+1
      k2 = k1+1
      if(iopt.lt.(-1) .or. iopt.gt.1) go to 50
      nmin = 2*k1
      if(m.lt.k1 .or. nest.lt.nmin) go to 50
      lwest = m*k1+nest*(7+3*k)
      if(lwrk.lt.lwest) go to 50
      if(xb.gt.x(1) .or. xe.lt.x(m) .or. w(1).le.0.) go to 50
      do 10 i=2,m
         if(x(i-1).ge.x(i) .or. w(i).le.0.) go to 50
  10  continue
      if(iopt.ge.0) go to 30
      if(n.lt.nmin .or. n.gt.nest) go to 50
      j = n
      do 20 i=1,k1
         t(i) = xb
         t(j) = xe
         j = j-1
  20  continue
      call fpchec(x,m,t,n,k,ier)
      if(ier) 50,40,50
  30  if(s.lt.0.) go to 50
      if(s.eq.0. .and. nest.lt.(m+k1)) go to 50
      ier = 0
c we partition the working space and determine the spline approximation.
  40  ifp = 1
      iz = ifp+nest
      ia = iz+nest
      ib = ia+nest*k1
      ig = ib+nest*k2
      iq = ig+nest*k2
      call fpcurf(iopt,x,y,w,m,xb,xe,k,s,nest,tol,maxit,k1,k2,n,t,c,fp,
     * wrk(ifp),wrk(iz),wrk(ia),wrk(ib),wrk(ig),wrk(iq),iwrk,ier)
  50  return
      end
      subroutine fpback(a,z,n,k,c,nest)
c  subroutine fpback calculates the solution of the system of
c  equations a*c = z with a a n x n upper triangular matrix
c  of bandwidth k.
c  ..
c  ..scalar arguments..
      integer n,k,nest
c  ..array arguments..
      real*8 a(nest,k),z(n),c(n)
c  ..local scalars..
      real*8 store
      integer i,i1,j,k1,l,m
c  ..
      k1 = k-1
      c(n) = z(n)/a(n,1)
      i = n-1
      if(i.eq.0) go to 30
      do 20 j=2,n
        store = z(i)
        i1 = k1
        if(j.le.k1) i1 = j-1
        m = i
        do 10 l=1,i1
          m = m+1
          store = store-c(m)*a(i,l+1)
  10    continue
        c(i) = store/a(i,1)
        i = i-1
  20  continue
  30  return
      end
      subroutine fpbspl(t,n,k,x,l,h)
c  subroutine fpbspl evaluates the (k+1) non-zero b-splines of
c  degree k at t(l) <= x < t(l+1) using the stable recurrence
c  relation of de boor and cox.
c  ..
c  ..scalar arguments..
      real*8 x
      integer n,k,l
c  ..array arguments..
      real*8 t(n),h(6)
c  ..local scalars..
      real*8 f,one
      integer i,j,li,lj
c  ..local arrays..
      real*8 hh(5)
c  ..
      one = 0.1e+01
      h(1) = one
      do 20 j=1,k
        do 10 i=1,j
          hh(i) = h(i)
  10    continue
        h(1) = 0.
        do 20 i=1,j
          li = l+i
          lj = li-j
          f = hh(i)/(t(li)-t(lj))
          h(i) = h(i)+f*(t(li)-x)
          h(i+1) = f*(x-t(lj))
  20  continue
      return
      end
      subroutine fpchec(x,m,t,n,k,ier)
c  subroutine fpchec verifies the number and the position of the knots
c  t(j),j=1,2,...,n of a spline of degree k, in relation to the number
c  and the position of the data points x(i),i=1,2,...,m. if all of the
c  following conditions are fulfilled, the error parameter ier is set
c  to zero. if one of the conditions is violated ier is set to ten.
c      1) k+1 <= n-k-1 <= m
c      2) t(1) <= t(2) <= ... <= t(k+1)
c         t(n-k) <= t(n-k+1) <= ... <= t(n)
c      3) t(k+1) < t(k+2) < ... < t(n-k)
c      4) t(k+1) <= x(i) <= t(n-k)
c      5) the conditions specified by schoenberg and whitney must hold
c         for at least one subset of data points, i.e. there must be a
c         subset of data points y(j) such that
c             t(j) < y(j) < t(j+k+1), j=1,2,...,n-k-1
c  ..
c  ..scalar arguments..
      integer m,n,k,ier
c  ..array arguments..
      real*8 x(m),t(n)
c  ..local scalars..
      integer i,j,k1,k2,l,nk1,nk2,nk3
      real*8 tj,tl
c  ..
      k1 = k+1
      k2 = k1+1
      nk1 = n-k1
      nk2 = nk1+1
      ier = 10
c  check condition no 1
      if(nk1.lt.k1 .or. nk1.gt.m) go to 80
c  check condition no 2
      j = n
      do 20 i=1,k
        if(t(i).gt.t(i+1)) go to 80
        if(t(j).lt.t(j-1)) go to 80
        j = j-1
  20  continue
c  check condition no 3
      do 30 i=k2,nk2
        if(t(i).le.t(i-1)) go to 80
  30  continue
c  check condition no 4
      if(x(1).lt.t(k1) .or. x(m).gt.t(nk2)) go to 80
c  check condition no 5
      if(x(1).ge.t(k2) .or. x(m).le.t(nk1)) go to 80
      i = 1
      l = k2
      nk3 = nk1-1
      if(nk3.lt.2) go to 70
      do 60 j=2,nk3
        tj = t(j)
        l = l+1
        tl = t(l)
  40    i = i+1
        if(i.ge.m) go to 80
        if(x(i).le.tj) go to 40
        if(x(i).ge.tl) go to 80
  60  continue
  70  ier = 0
  80  return
      end
      subroutine fpcurf(iopt,x,y,w,m,xb,xe,k,s,nest,tol,maxit,k1,k2,
     * n,t,c,fp,fpint,z,a,b,g,q,nrdata,ier)
c  ..
c  ..scalar arguments..
      real*8 xb,xe,s,tol,fp
      integer iopt,m,k,nest,maxit,k1,k2,n,ier
c  ..array arguments..
      real*8 x(m),y(m),w(m),t(nest),c(nest),fpint(nest),
     * z(nest),a(nest,k1),b(nest,k2),g(nest,k2),q(m,k1)
      integer nrdata(nest)
c  ..local scalars..
      real*8 acc,con1,con4,con9,cos,half,fpart,fpms,fpold,fp0,f1,f2,f3,
     * one,p,pinv,piv,p1,p2,p3,rn,sin,store,term,wi,xi,yi
      integer i,ich1,ich3,it,iter,i1,i2,i3,j,k3,l,l0,
     * mk1,new,nk1,nmax,nmin,nplus,npl1,nrint,n8
c  ..local arrays..
      real*8 h(7)
c  ..function references
      real*8 abs,fprati
      integer max0,min0
c  ..subroutine references..
c    fpback,fpbspl,fpgivs,fpdisc,fpknot,fprota
c  ..
c  set constants
      one = 0.1e+01
      con1 = 0.1e0
      con9 = 0.9e0
      con4 = 0.4e-01
      half = 0.5e0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  part 1: determination of the number of knots and their position     c
c  **************************************************************      c
c  given a set of knots we compute the least-squares spline sinf(x),   c
c  and the corresponding sum of squared residuals fp=f(p=inf).         c
c  if iopt=-1 sinf(x) is the requested approximation.                  c
c  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
c    if fp <=s we will continue with the current set of knots.         c
c    if fp > s we will increase the number of knots and compute the    c
c       corresponding least-squares spline until finally fp<=s.        c
c    the initial choice of knots depends on the value of s and iopt.   c
c    if s=0 we have spline interpolation; in that case the number of   c
c    knots equals nmax = m+k+1.                                        c
c    if s > 0 and                                                      c
c      iopt=0 we first compute the least-squares polynomial of         c
c      degree k; n = nmin = 2*k+2                                      c
c      iopt=1 we start with the set of knots found at the last         c
c      call of the routine, except for the case that s > fp0; then     c
c      we compute directly the least-squares polynomial of degree k.   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  determine nmin, the number of knots for polynomial approximation.
      nmin = 2*k1
      if(iopt.lt.0) go to 60
c  calculation of acc, the absolute tolerance for the root of f(p)=s.
      acc = tol*s
c  determine nmax, the number of knots for spline interpolation.
      nmax = m+k1
      if(s.gt.0.) go to 45
c  if s=0, s(x) is an interpolating spline.
c  test whether the required storage space exceeds the available one.
      n = nmax
      if(nmax.gt.nest) go to 420
c  find the position of the interior knots in case of interpolation.
  10  mk1 = m-k1
      if(mk1.eq.0) go to 60
      k3 = k/2
      i = k2
      j = k3+2
      if(k3*2.eq.k) go to 30
      do 20 l=1,mk1
        t(i) = x(j)
        i = i+1
        j = j+1
  20  continue
      go to 60
  30  do 40 l=1,mk1
        t(i) = (x(j)+x(j-1))*half
        i = i+1
        j = j+1
  40  continue
      go to 60
c  if s>0 our initial choice of knots depends on the value of iopt.
c  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
c  polynomial of degree k which is a spline without interior knots.
c  if iopt=1 and fp0>s we start computing the least squares spline
c  according to the set of knots found at the last call of the routine.
  45  if(iopt.eq.0) go to 50
      if(n.eq.nmin) go to 50
      fp0 = fpint(n)
      fpold = fpint(n-1)
      nplus = nrdata(n)
      if(fp0.gt.s) go to 60
  50  n = nmin
      fpold = 0.
      nplus = 0
      nrdata(1) = m-2
c  main loop for the different sets of knots. m is a save upper bound
c  for the number of trials.
  60  do 200 iter = 1,m
        if(n.eq.nmin) ier = -2
c  find nrint, tne number of knot intervals.
        nrint = n-nmin+1
c  find the position of the additional knots which are needed for
c  the b-spline representation of s(x).
        nk1 = n-k1
        i = n
        do 70 j=1,k1
          t(j) = xb
          t(i) = xe
          i = i-1
  70    continue
c  compute the b-spline coefficients of the least-squares spline
c  sinf(x). the observation matrix a is built up row by row and
c  reduced to upper triangular form by givens transformations.
c  at the same time fp=f(p=inf) is computed.
        fp = 0.
c  initialize the observation matrix a.
        do 80 i=1,nk1
          z(i) = 0.
          do 80 j=1,k1
            a(i,j) = 0.
  80    continue
        l = k1
        do 130 it=1,m
c  fetch the current data point x(it),y(it).
          xi = x(it)
          wi = w(it)
          yi = y(it)*wi
c  search for knot interval t(l) <= xi < t(l+1).
  85      if(xi.lt.t(l+1) .or. l.eq.nk1) go to 90
          l = l+1
          go to 85
c  evaluate the (k+1) non-zero b-splines at xi and store them in q.
  90      call fpbspl(t,n,k,xi,l,h)
          do 95 i=1,k1
            q(it,i) = h(i)
            h(i) = h(i)*wi
  95      continue
c  rotate the new row of the observation matrix into triangle.
          j = l-k1
          do 110 i=1,k1
            j = j+1
            piv = h(i)
            if(piv.eq.0.) go to 110
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,a(j,1),cos,sin)
c  transformations to right hand side.
            call fprota(cos,sin,yi,z(j))
            if(i.eq.k1) go to 120
            i2 = 1
            i3 = i+1
            do 100 i1 = i3,k1
              i2 = i2+1
c  transformations to left hand side.
              call fprota(cos,sin,h(i1),a(j,i2))
 100        continue
 110      continue
c  add contribution of this row to the sum of squares of residual
c  right hand sides.
 120      fp = fp+yi**2
 130    continue
        if(ier.eq.(-2)) fp0 = fp
        fpint(n) = fp0
        fpint(n-1) = fpold
        nrdata(n) = nplus
c  backward substitution to obtain the b-spline coefficients.
        call fpback(a,z,nk1,k1,c,nest)
c  test whether the approximation sinf(x) is an acceptable solution.
        if(iopt.lt.0) go to 440
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 440
c  if f(p=inf) < s accept the choice of knots.
        if(fpms.lt.0.) go to 250
c  if n = nmax, sinf(x) is an interpolating spline.
        if(n.eq.nmax) go to 430
c  increase the number of knots.
c  if n=nest we cannot increase the number of knots because of
c  the storage capacity limitation.
        if(n.eq.nest) go to 420
c  determine the number of knots nplus we are going to add.
        if(ier.eq.0) go to 140
        nplus = 1
        ier = 0
        go to 150
 140    npl1 = nplus*2
        rn = nplus
        if(fpold-fp.gt.acc) npl1 = rn*fpms/(fpold-fp)
        nplus = min0(nplus*2,max0(npl1,nplus/2,1))
 150    fpold = fp
c  compute the sum((w(i)*(y(i)-s(x(i))))**2) for each knot interval
c  t(j+k) <= x(i) <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
        fpart = 0.
        i = 1
        l = k2
        new = 0
        do 180 it=1,m
          if(x(it).lt.t(l) .or. l.gt.nk1) go to 160
          new = 1
          l = l+1
 160      term = 0.
          l0 = l-k2
          do 170 j=1,k1
            l0 = l0+1
            term = term+c(l0)*q(it,j)
 170      continue
          term = (w(it)*(term-y(it)))**2
          fpart = fpart+term
          if(new.eq.0) go to 180
          store = term*half
          fpint(i) = fpart-store
          i = i+1
          fpart = store
          new = 0
 180    continue
        fpint(nrint) = fpart
        do 190 l=1,nplus
c  add a new knot.
          call fpknot(x,m,t,n,fpint,nrdata,nrint,nest,1)
c  if n=nmax we locate the knots as for interpolation.
          if(n.eq.nmax) go to 10
c  test whether we cannot further increase the number of knots.
          if(n.eq.nest) go to 200
 190    continue
c  restart the computations with the new set of knots.
 200  continue
c  test whether the least-squares kth degree polynomial is a solution
c  of our approximation problem.
 250  if(ier.eq.(-2)) go to 440
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  part 2: determination of the smoothing spline sp(x).                c
c  ***************************************************                 c
c  we have determined the number of knots and their position.          c
c  we now compute the b-spline coefficients of the smoothing spline    c
c  sp(x). the observation matrix a is extended by the rows of matrix   c
c  b expressing that the kth derivative discontinuities of sp(x) at    c
c  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c
c  ponding weights of these additional rows are set to 1/p.            c
c  iteratively we then have to determine the value of p such that      c
c  f(p)=sum((w(i)*(y(i)-sp(x(i))))**2) be = s. we already know that    c
c  the least-squares kth degree polynomial corresponds to p=0, and     c
c  that the least-squares spline corresponds to p=infinity. the        c
c  iteration process which is proposed here, makes use of rational     c
c  interpolation. since f(p) is a convex and strictly decreasing       c
c  function of p, it can be approximated by a rational function        c
c  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  c
c  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      c
c  to calculate the new value of p such that r(p)=s. convergence is    c
c  guaranteed by taking f1>0 and f3<0.                                 c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  evaluate the discontinuity jump of the kth derivative of the
c  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
      call fpdisc(t,n,k2,b,nest)
c  initial value for p.
      p1 = 0.
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      p = 0.
      do 255 i=1,nk1
         p = p+a(i,1)
 255  continue
      rn = nk1
      p = rn/p
      ich1 = 0
      ich3 = 0
      n8 = n-nmin
c  iteration process to find the root of f(p) = s.
      do 360 iter=1,maxit
c  the rows of matrix b with weight 1/p are rotated into the
c  triangularised observation matrix a which is stored in g.
        pinv = one/p
        do 260 i=1,nk1
          c(i) = z(i)
          g(i,k2) = 0.
          do 260 j=1,k1
            g(i,j) = a(i,j)
 260    continue
        do 300 it=1,n8
c  the row of matrix b is rotated into triangle by givens transformation
          do 270 i=1,k2
            h(i) = b(it,i)*pinv
 270      continue
          yi = 0.
          do 290 j=it,nk1
            piv = h(1)
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,g(j,1),cos,sin)
c  transformations to right hand side.
            call fprota(cos,sin,yi,c(j))
            if(j.eq.nk1) go to 300
            i2 = k1
            if(j.gt.n8) i2 = nk1-j
            do 280 i=1,i2
c  transformations to left hand side.
              i1 = i+1
              call fprota(cos,sin,h(i1),g(j,i1))
              h(i) = h(i1)
 280        continue
            h(i2+1) = 0.
 290      continue
 300    continue
c  backward substitution to obtain the b-spline coefficients.
        call fpback(g,c,nk1,k2,c,nest)
c  computation of f(p).
        fp = 0.
        l = k2
        do 330 it=1,m
          if(x(it).lt.t(l) .or. l.gt.nk1) go to 310
          l = l+1
 310      l0 = l-k2
          term = 0.
          do 320 j=1,k1
            l0 = l0+1
            term = term+c(l0)*q(it,j)
 320      continue
          fp = fp+(w(it)*(term-y(it)))**2
 330    continue
c  test whether the approximation sp(x) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 440
c  test whether the maximal number of iterations is reached.
        if(iter.eq.maxit) go to 400
c  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3.ne.0) go to 340
        if((f2-f3).gt.acc) go to 335
c  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p.le.p1) p=p1*con9 + p2*con1
        go to 360
 335    if(f2.lt.0.) ich3=1
 340    if(ich1.ne.0) go to 350
        if((f1-f2).gt.acc) go to 345
c  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3.lt.0.) go to 360
        if(p.ge.p3) p = p2*con1 + p3*con9
        go to 360
 345    if(f2.gt.0.) ich1=1
c  test whether the iteration process proceeds as theoretically
c  expected.
 350    if(f2.ge.f1 .or. f2.le.f3) go to 410
c  find the new value for p.
        p = fprati(p1,f1,p2,f2,p3,f3)
 360  continue
c  error codes and messages.
 400  ier = 3
      go to 440
 410  ier = 2
      go to 440
 420  ier = 1
      go to 440
 430  ier = -1
 440  return
      end
      subroutine fpdisc(t,n,k2,b,nest)
c  subroutine fpdisc calculates the discontinuity jumps of the kth
c  derivative of the b-splines of degree k at the knots t(k+2)..t(n-k-1)
c  ..scalar arguments..
      integer n,k2,nest
c  ..array arguments..
      real*8 t(n),b(nest,k2)
c  ..local scalars..
      real*8 an,fac,prod
      integer i,ik,j,jk,k,k1,l,lj,lk,lmk,lp,nk1,nrint
c  ..local array..
      real*8 h(12)
c  ..
      k1 = k2-1
      k = k1-1
      nk1 = n-k1
      nrint = nk1-k
      an = nrint
      fac = an/(t(nk1+1)-t(k1))
      do 40 l=k2,nk1
        lmk = l-k1
        do 10 j=1,k1
          ik = j+k1
          lj = l+j
          lk = lj-k2
          h(j) = t(l)-t(lk)
          h(ik) = t(l)-t(lj)
  10    continue
        lp = lmk
        do 30 j=1,k2
          jk = j
          prod = h(j)
          do 20 i=1,k
            jk = jk+1
            prod = prod*h(jk)*fac
  20      continue
          lk = lp+k1
          b(lmk,j) = (t(lk)-t(lp))/prod
          lp = lp+1
  30    continue
  40  continue
      return
      end
      subroutine fpgivs(piv,ww,cos,sin)
c  subroutine fpgivs calculates the parameters of a givens
c  transformation .
c  ..
c  ..scalar arguments..
      real*8 piv,ww,cos,sin
c  ..local scalars..
      real*8 dd,one,store
c  ..function references..
      real*8 abs,sqrt
c  ..
      one = 0.1e+01
      store = abs(piv)
      if(store.ge.ww) dd = store*sqrt(one+(ww/piv)**2)
      if(store.lt.ww) dd = ww*sqrt(one+(piv/ww)**2)
      cos = ww/dd
      sin = piv/dd
      ww = dd
      return
      end
      subroutine fpknot(x,m,t,n,fpint,nrdata,nrint,nest,istart)
c  subroutine fpknot locates an additional knot for a spline of degree
c  k and adjusts the corresponding parameters,i.e.
c    t     : the position of the knots.
c    n     : the number of knots.
c    nrint : the number of knotintervals.
c    fpint : the sum of squares of residual right hand sides
c            for each knot interval.
c    nrdata: the number of data points inside each knot interval.
c  istart indicates that the smallest data point at which the new knot
c  may be added is x(istart+1)
c  ..
c  ..scalar arguments..
      integer m,n,nrint,nest,istart
c  ..array arguments..
      real*8 x(m),t(nest),fpint(nest)
      integer nrdata(nest)
c  ..local scalars..
      real*8 an,am,fpmax
      integer ihalf,j,jbegin,jj,jk,jpoint,k,maxbeg,maxpt,
     * next,nrx,number
c  ..
      k = (n-nrint-1)/2
c  search for knot interval t(number+k) <= x <= t(number+k+1) where
c  fpint(number) is maximal on the condition that nrdata(number)
c  not equals zero.
      fpmax = 0.
      jbegin = istart
      do 20 j=1,nrint
        jpoint = nrdata(j)
        if(fpmax.ge.fpint(j) .or. jpoint.eq.0) go to 10
        fpmax = fpint(j)
        number = j
        maxpt = jpoint
        maxbeg = jbegin
  10    jbegin = jbegin+jpoint+1
  20  continue
c  let coincide the new knot t(number+k+1) with a data point x(nrx)
c  inside the old knot interval t(number+k) <= x <= t(number+k+1).
      ihalf = maxpt/2+1
      nrx = maxbeg+ihalf
      next = number+1
      if(next.gt.nrint) go to 40
c  adjust the different parameters.
      do 30 j=next,nrint
        jj = next+nrint-j
        fpint(jj+1) = fpint(jj)
        nrdata(jj+1) = nrdata(jj)
        jk = jj+k
        t(jk+1) = t(jk)
  30  continue
  40  nrdata(number) = ihalf-1
      nrdata(next) = maxpt-ihalf
      am = maxpt
      an = nrdata(number)
      fpint(number) = fpmax*an/am
      an = nrdata(next)
      fpint(next) = fpmax*an/am
      jk = next+k
      t(jk) = x(nrx)
      n = n+1
      nrint = nrint+1
      return
      end
      real*8 function fprati(p1,f1,p2,f2,p3,f3)
c  given three points (p1,f1),(p2,f2) and (p3,f3), function fprati
c  gives the value of p such that the rational interpolating function
c  of the form r(p) = (u*p+v)/(p+w) equals zero at p.
c  ..
c  ..scalar arguments..
      real*8 p1,f1,p2,f2,p3,f3
c  ..local scalars..
      real*8 h1,h2,h3,p
c  ..
      if(p3.gt.0.) go to 10
c  value of p in case p3 = infinity.
      p = (p1*(f1-f3)*f2-p2*(f2-f3)*f1)/((f1-f2)*f3)
      go to 20
c  value of p in case p3 ^= infinity.
  10  h1 = f1*(f2-f3)
      h2 = f2*(f3-f1)
      h3 = f3*(f1-f2)
      p = -(p1*p2*h3+p2*p3*h1+p3*p1*h2)/(p1*h1+p2*h2+p3*h3)
c  adjust the value of p1,f1,p3 and f3 such that f1 > 0 and f3 < 0.
  20  if(f2.lt.0.) go to 30
      p1 = p2
      f1 = f2
      go to 40
  30  p3 = p2
      f3 = f2
  40  fprati = p
      return
      end
      subroutine fprota(cos,sin,a,b)
c  subroutine fprota applies a givens rotation to a and b.
c  ..
c  ..scalar arguments..
      real*8 cos,sin,a,b
c ..local scalars..
      real*8 stor1,stor2
c  ..
      stor1 = a
      stor2 = b
      b = cos*stor2+sin*stor1
      a = cos*stor1-sin*stor2
      return
      end

      subroutine splev(t,n,c,k,x,y,m,ier)
c  subroutine splev evaluates in a number of points x(i),i=1,2,...,m
c  a spline s(x) of degree k, given in its b-spline representation.
c
c  calling sequence:
c     call splev(t,n,c,k,x,y,m,ier)
c
c  input parameters:
c    t    : array,length n, which contains the position of the knots.
c    n    : integer, giving the total number of knots of s(x).
c    c    : array,length n, which contains the b-spline coefficients.
c    k    : integer, giving the degree of s(x).
c    x    : array,length m, which contains the points where s(x) must
c           be evaluated.
c    m    : integer, giving the number of points where s(x) must be
c           evaluated.
c
c  output parameter:
c    y    : array,length m, giving the value of s(x) at the different
c           points.
c    ier  : error flag
c      ier = 0 : normal return
c      ier =10 : invalid input data (see restrictions)
c
c  restrictions:
c    m >= 1
c    t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
c
c  other subroutines required: fpbspl.
c
c  references :
c    de boor c  : on calculating with b-splines, j. approximation theory
c                 6 (1972) 50-62.
c    cox m.g.   : the numerical evaluation of b-splines, j. inst. maths
c                 applics 10 (1972) 134-149.
c    dierckx p. : curve and surface fitting with splines, monographs on
c                 numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1987
c
c  ..scalar arguments..
      integer n,k,m,ier
c  ..array arguments..
      real*8 t(n),c(n),x(m),y(m)
c  ..local scalars..
      integer i,j,k1,l,ll,l1,nk1
      real*8 arg,sp,tb,te
c  ..local array..
      real*8 h(6)
c  ..
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 10
      if(m-1) 100,30,10
  10  do 20 i=2,m
        if(x(i).lt.x(i-1)) go to 100
  20  continue
  30  ier = 0
c  fetch tb and te, the boundaries of the approximation interval.
      k1 = k+1
      nk1 = n-k1
      tb = t(k1)
      te = t(nk1+1)
      l = k1
      l1 = l+1
c  main loop for the different points.
      do 80 i=1,m
c  fetch a new x-value arg.
        arg = x(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
c  search for knot interval t(l) <= arg < t(l+1)
  40    if(arg.lt.t(l1) .or. l.eq.nk1) go to 50
        l = l1
        l1 = l+1
        go to 40
c  evaluate the non-zero b-splines at arg.
  50    call fpbspl(t,n,k,arg,l,h)
c  find the value of s(x) at x=arg.
        sp = 0.
        ll = l-k1
        do 60 j=1,k1
          ll = ll+1
          sp = sp+c(ll)*h(j)
  60    continue
        y(i) = sp
  80  continue
 100  return
      end
