      subroutine gwf(at,zt,ap,zp,in,il,ij,vw,rw,aw,vs,rs,as,em,ex,wf,eo)

      IMPLICIT REAL*8 (A-H,O-Z)

      INTEGER at,zt,ap,zp,in,il,ij
      REAL vw,aw,rw,vs,as,rs,em,ex
      REAL*8 wf(100)

cf2py intent(in) at,zt,ap,zp,in,il,ij,vw,rw,aw,vs,rs,as,em,ex
cf2py intent(out) wf(100), eo

      common/cadif/adif,rr0
      common/ws1/nrp,lp,j2,iat,izt,iap,izp
      common/ws2/vnorm,vz,vzs,vsop,vnormls
      common/cem/emins,emaxs
      common/crad/xkin,rms,vc,xkint,xmst,nt,nb
      common/main2/del,hbsd2m
      common/crmax/rmax_unbound
      common/crmin/emin,emax,vnmin,vnmax


c    set params from input
      iat = at
      izt = zt
      iap = ap
      izp = zp
      nrp = in
      lp = il
      j2 = ij

c444   format('j2=',i10.3)
c      print 444,j2
c      return

cbab meshsize in fm
      del = 0.1
cbab max r for unbound phase shift calculation
      rmax_unbound = 500.

cbab set ws parameters
      vz = vw
      r0 = rw
      adif = aw

      vsop = vs
      r0so = rs
      aso = as

      vzs = -30.0
      vnorm = 1
      vnormls = 1.

cbab back come here is n value is -1
960   continue
      xkint = 0.
      xmst = 0.
      nt = 0
      nb = 0
910   continue
      call resetb

903   continue
      emins = emin
      emaxs = emax

cbab number of steps for scattering energy
      xne = 0.
      xvn = 0.
      if((emax-emin).gt.0.) xne = 100.
      if((vnmax-vnmin).gt.0.) xvn = 100.

800   continue
      if(nrp.eq.-1) go to 960
      if(j2.eq.0) go to 960


cbab emin.le.0. indicates a bound state
555   format(f10.3)
c      print 555,emin
      if(emin.le.0.) go to 810
      go to 811
810   continue
c defualt setup for bound states
      xne = 50.
c      emax = -1.
c      emin = -90.
      emax = ex
      emin = em
c      print 813
c813   format('boundstate')
      call boundstate(emin,emax,xne)
c only doing one set of q numbers so return after one
      return
      go to 800
811   continue

995   continue
      call res(epr,emin,emax,xne,vnmin,vnmax,xvn,g,ire,0)
      if(ire.eq.1) go to 899

cbab loop until the resonance is found
900   continue
      gold = g
      call res(epr,emin,emax,xne,vnmin,vnmax,xvn,g,ire,0)
      if(ire.eq.1) go to 899
      if(abs(g-gold)/g.lt.0.1) go to 902
      go to 900
902   continue

cbab one more call
      if(xne.gt.1.) call res(epr,emin,emax,xne,vnmin,vnmax,xvn,g,ire,1)

cbab this means that vn was varied so go back and read emin,emax for fixed vn
      if(xvn.eq.1.) go to 998
      epr = emin
      emin = epr/2.
      emax = 2.*epr
c      print 7777,epr,emin,emax
c7777  format(3e12.4)
      emins = emin
      emaxs = emax

cbab number of steps for scattering energy
      xne = 0.
      xvn = 0.
      if((emax-emin).gt.0.) xne = 100.
      xvn = 1.
      go to 995
998   continue

c      if(xvn.gt.1.) go to 903
c      print 881,emin,epr,emax
c881   format(3f10.3)
c      open(unit=12,file=filetop)
      t1 = (emax-emin)/6.
      t2 = (emax-emin)/2.
      if(epr.le.0.1) write(12,850) emin-epr,emax-epr,t1,t2,fileplt
850   format('fi',/,e12.4,',',e12.4,',',e12.4,',',e12.4,
     1',E-E_(o) (MeV)',/,'0.,1.1,.1,.2,P(E)',/,'ds line',/,'df ',a13)
      if(epr.gt.0.1) write(12,851) emin-epr,emax-epr,t1,t2,fileplt
851   format('fi',/,f7.3,',',f7.3,',',f7.3,',',f7.3,
     1',E-E_(o) (MeV)',/,'0.,1.1,.1,.2,P(E)',/,'ds line',/,'df ',a13)
      stop

899   continue
      print 898
898   format('no resonance found')

      contains

      subroutine  well(xpot,r0,v,dv)
      IMPLICIT REAL*8 (A-H,O-Z)
      common/cadif/adif,rr0
      dimension cp(300),cn(300),sp(300),sn(300),co(300)
      dimension v(50000),vbab(5,50000)
      common/cwell/vwst
      common/main1/nmax,nq
      common/main2/del,hbsd2m
      common/ws1/nrp,lp,j2,iat,izt,iap,izp
      common/ws2/vnorm,vz,vzs,vsop,vnormls
      common/cifirstp/ifirstp,ifirstv
      at2 = iat
      ztot = izt*izp
      int = iat-izt

c parameters for ws potential
      rr0c = 1.20
      xnz = int-izt
      if(izp.eq.0.) vws = vz*iap-vzs*xnz/at2
      if(izp.gt.0.) vws = vz*iap+vzs*xnz/at2
      if(ifirstp.eq.999) go to 1050
c      print 1051,vws
c1051  format('potential depth = ',f10.3)
      ifirstp = 999
1050  continue
      vwst = vws*vnorm
c      r0=rr0*(iat)**0.3333
c      r0c=rr0c*(iat)**0.3333
c      r0so=rso*(iat)**0.3333
      r0c=2.0
c      print 7777,vws,r0,adif
c7777  format(1x,3f10.3)

      do 10 i=1,nmax
      x=del*i
      v2=lp*(lp+1)*hbsd2m/x**2
      fl=-(lp+1)
      if(j2.gt.2*lp) fl=lp
      if(j2.eq.2*lp) fl=0.0

c woods-saxon potential
      y=(x-r0)/adif
      yso=(x-r0so)/aso  !
      if(x.lt.10.)e=exp(y)
      if(x.ge.10.)e=9999999.
      if(x.lt.10.) f=1./(1.+e)
      if(x.ge.10.) f=0.
      if(x.lt.10.)eso=exp(yso)  !
      if(x.ge.10.)eso=9999999.  !
      if(x.lt.10.) fso=1./(1.+eso)!
      if(x.ge.10.) fso=0.  !
      vcoul1=ztot*1.44/x
      vcoul2 = 0.
      if(x.lt.r0c) vcoul2=ztot*1.44/r0c*(1.5-0.5*(x/r0c)**2)-vcoul1
      vso=vsop*(fso**2)*eso*fl/(x*aso)
c      if(i.eq.10) print 7771,vnormls,vsop,f,e,fl
c7771  format(10f8.3)
      v1=-vws*f

      v(i)=xpot*(vnorm*vws*f+vso+vcoul2)+v2+vcoul1
c      if(i.eq.10) print 666,lp,v(i),vnormls,vsop,f,e,fl,x,adif
c666   format(i3,10f10.3)
      if(ifirstv.eq.999) go to 502
      vbab(1,i) = vnorm*vws*f
      vbab(2,i) = vso
      vbab(3,i) = vcoul2+vcoul1
      vbab(4,i) = v2
      vbab(5,i) = v(i)
502   continue
10    continue

c      if(ifirstv.eq.999) go to 503
c      do 500 ik = 1,5
c      do 500 i = 1,nmax
c      x = del*i
c500   write(40,501) x,vbab(ik,i)
c501   format(f6.1,1x,f8.2)
c503   continue
      ifirstv = 999



      return
      end

      subroutine resetb
      common/reset/ifirstb
      ifirstb = 0
      return
      end

      subroutine boundstate(emin,emax,xne)
      IMPLICIT REAL*8 (A-H,O-Z)
      common/ceb/eb
      common/cwell/vwst
      common/crad/xkin,rms,vc,xkint,xmst,nt,nb
      common/ws1/nrp,lp,j2,iat,izt,iap,izp
      common/reset/ifirstb
      dimension eis(400)
      estep = (emax-emin)/xne
      ei = emin - estep
      i = 0
600   continue
      i = i + 1
      ei = ei + estep
      if(ei.gt.emax) go to 601

      call gb(ei,phi,nodes)
620   continue

c      print 820,i,nodes,nrp,ei
c820   format(3i3,1x,f10.3)
      eis(i) = eb

c eb.gt.0. indicates a failed guess for bound states
      if(eb.eq.0.) go to 600

c skip the output if the same energy is found from a previous try
      if(i.eq.1) go to 622
      do 622 j = 1,i-1
      i1 = 1000*eis(i)
      i2 = 1000*eis(j)
621   if(i1.eq.i2) go to 600
622   continue
      if(nodes.ne.nrp) go to 600

601   continue
      if(eb.gt.0.) go to 830

      if(ifirstb.eq.999) go to 615
      ifirstb = 999
c      xmst = 0.
c      nt = 0
c      nb = 0
c      if(ifirst610.ne.999) print 610
c      if(ifirst610.ne.999) write(10,610)
      if(ifirst610.ne.999) ifirst610 = 999
610   format(1x,'n,l,2j,spe,k,kt,rms,rmst,v0 = ')
615   continue
      nb = nb + 1
      xmst = xmst + (j2+1.)*rms*rms
      nt = nt + (j2+1)
      rmst = sqrt(xmst/nt)
      xkint = xkint + (j2+1.)*xkin
      eo = eb
c     print 614,nrp,lp,j2,eb,xkin,xkint,rms,rmst,vwst
c     write(10,614) nrp,lp,j2,eb,xkin,xkint,rms,rmst,vwst
614   format(1x,3i3,1x,7f10.3,/)

c      if(emin.lt.0.and.nb.eq.0) print 623,lp,j2
c      if(emin.lt.0.and.nb.eq.0) write(10,623) lp,j2
c      if(emin.lt.0.) ire = 1

      return
830   continue
c      write(10,623) nrp,lp,j2
c      print 623,nrp,lp,j2
c623   format(1x,3i3,' no bound state found')

      return
      end




      subroutine res(epr,emin,emax,xne,vnmin,vnmax,xvn,g,ire,iout)
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension eis(400),phis(400),sigs(400),vns(400)
      common/cem/emins,emaxs
      common/ws1/nrp,lp,j2,iat,izt,iap,izp
      common/cwell/vwst
      common/crad/xkin,rms,vc,xkint,xmst,nt,nb
      common/ceb/eb

      common/ws2/vnorm,vz,vzs,vsop,vnormls
      common/reset/ifirstb
      common/cifirst610/ifist610
c      character*2 sym

      ire = 0
      nb = 0
c      if(emax.eq.0) emax = emin
      if(xne.eq.0.) xne = 1.
      estep = (emax-emin)/xne
      ei = emin - estep

      if(xvn.eq.0.) xvn= 1.
      if(xvn.gt.1..and.vnmin.eq.0.) vnmin = 1.
      if(xvn.eq.1..and.vnmin.eq.0.) vnmin = vnorm
      if(vnmax.eq.0.) vnmax = vnmin
      if(xvn.gt.1.) vnstep = (vnmax-vnmin)/xvn
      if(xvn.gt.1.) vnorm = vnmin - vnstep

      i = 0

600   continue
      ei = ei + estep
      if(ei.gt.emax) go to 601
      if(xvn.gt.1.) vnorm = vnorm + vnstep
      if(xvn.gt.10.and.vnorm.gt.vnmax) go to 601
      i = i + 1
c      open(unit=97,file='res.dat')
      call gb(ei,phi,nodes)

      eis(i) = ei
      vns(i) = vnorm
      xk2 = ((2.*938.)/(197.*197.))*ei
      sig = (2.*lp+1)*4.*3.14*sin(phi)**2./xk2
      phis(i) = phi
      sigs(i) = sig
      if(ifirst613.ne.999.and.iout.eq.1) write(10,613)
c      if(ifirst613.ne.999.and.iout.eq.1) print 613
      if(ifirst613.ne.999.and.iout.eq.1) ifirst613 = 999
613   format(/,1x,'n, e, phi, sig=sin(phi)**2/k, vnorm')
c      if(iout.eq.1) print 612,nodes,ei,phi,sig,vnorm
c      if(iout.eq.1) write(97,6612) ei,sig
c6612  format(f9.6,1x,e12.4)
      if(iout.eq.1) write(10,612) nodes,ei,phi,sig,vnorm
612   format(1x,i2,1x,f9.5,1x,f7.3,1x,e12.4,2x,f9.5,1x,e12.4)
      if(iout.eq.1) write(11,614) ei-epr,sig/sigmax
614   format(1x,e12.4,1x,e12.4)
      if(xne.eq.1..and.xvn.eq.1.) go to 601
      go to 600

601   continue
      ifirst613 = 0
      imax = i


c------------------------------------------------------------
      sigmax = 0.
cbab ikeep is the value of i at the maximum
      ikeep = 0
      do 700 i = 2,imax
      if(sigs(i).gt.sigs(i-1).and.sigs(i).gt.sigs(i+1)
     1.and.sigs(i).gt.sigmax) ikeep = i
      if(sigs(i).gt.sigs(i-1).and.sigs(i).gt.sigs(i+1)
     1.and.sigs(i).gt.sigmax) sigmax = sigs(i)
700   continue
c      ire = 1
c      return
c701   continue
      if(ikeep.eq.0) ire = 1
      if(ikeep.eq.0) return
      if(ikeep.eq.imax) ire = 1
      if(ikeep.eq.imax) return
      is = ikeep

      do 200 k = is-1,3,-1
      if(2.*sigs(k).lt.sigs(is)) go to 201
200   continue
201   continue
      g1 = 2.*(eis(is)-eis(k))
      g2 = 2.*(vns(is)-vns(k))

      if(g1.ne.0.) epr = eis(is)
      if(g2.ne.0.) epr = vns(is)

c      if(xne.gt.1.) emin = eis(k-1)
      if(xne.gt.1.) emin = epr-2*g1
      if(emin.lt.emins) emin = emins
      if(xvn.gt.1.) vnmin = vns(k-1)

c      if(xne.gt.1.) emax = epr + (epr-emin)
      if(xne.gt.1.) emax = epr+2*g1
      if(emax.gt.emaxs) emax = emaxs
      if(xvn.gt.1.) vnmax = epr + (epr-vnmin)
      g = g1 + g2
      call resout(epr,g1,g2,emin,emax,vnmin,vnmax)
c      call gb(epr,phi,nodes)
c      close(97)

c      print 750,emin,emax,vnmin,vnmax
c      write(10,750) emin,emax,vnmin,vnmax
c750   format(/,1x,'new emin,emax,vnmin,vnmax = ',4f9.6)

c      if(xvn.gt.1.) print 750,emin,vnmin,vnmax
c      write(10,750) emin,vnmin,vnmax
c750   format(/,1x,'new e,vnmin,vnmax = ',4f9.6)

c      if(xvn.eq.1.) print 751,emin,emax,vnmin
c      write(10,751) emin,emax,vnmin
c751   format(/,1x,'new emin,emax,vn = ',4f9.6)


      if(xvn.gt.1.) vnorm = epr
      return
      end

      subroutine resout(epr,g1,g2,emin,emax,vnmin,vnmax)
      IMPLICIT REAL*8 (A-H,O-Z)
      character*2 sym,sym2
      character*6 sym3
      g = 0.
      if(g1.gt.0.) g = g1
      if(g2.gt.0.) g = g2
      sym3 = 'MeV'
      if(g.lt..001) sym3 = 'keV'
      if(g.lt..001) g = 1000.*g
      if(g.lt..001) sym3 = 'ev'
      if(g.lt..001) g = 1000.*g
      if(g.lt..001) sym3 = '10-3ev'
      if(g.lt..001) g = 1000.*g
      if(g.lt..001) sym3 = '10-6ev'
      if(g.lt..001) g = 1000.*g
      if(g1.ne.0.) sym = 'ei'
      if(g2.ne.0.) sym = 'vn'
      if(g1.ne.0.) sym2 = ' G'
      if(g2.ne.0.) sym2 = 'dn'
      if(g2.ne.0.) sym3 = '   '
c      if(sym.eq.'ei') print 702,sym,epr,sym2,g,sym3,emin,emax
c      if(sym.eq.'ei') write(10,702) sym,epr,sym2,g,sym3,emin,emax
      if(sym.eq.'ei') print 702,sym,epr,sym2,g,sym3
      if(sym.eq.'ei') write(10,702) sym,epr,sym2,g,sym3
c      if(sym.eq.'vn') print 703,sym,epr,vnmin,vnmax
c      if(sym.eq.'vn') write(10,703) sym,epr,vnmin,vnmax
      if(sym.eq.'vn') print 703,sym,epr
      if(sym.eq.'vn') write(10,703) sym,epr
702   format(1x,
     1 'resonance at ',a2,
     1' = ',f7.3,' with ',a2,' = ',f10.5,1x,a6,1x,2e12.4)
703   format(1x,
     1 'resonance at ',a2,
     1' = ',f7.3,29x,2e12.4)
      return
      end


      subroutine gb(ei,phi3,nodes)
      IMPLICIT REAL*8 (A-H,O-Z)
c  version 10-10-85 with phase shift (g. f. bertsch)
c nmax = number of mesh points
c del  = mesh size
c iat  = mass number of target plus projectile
c izt  = proton number of target plus projectile
c ei>0   will find scattering phase shift
c ei<0   will search for bound state
c izp=0  for neutron single-particle state
c izp=1  for proton single-particle state
c lp   = orbital angular momentum for single-particle state
c j2   = twice the j value of the single-particle state
      common/cua/dimension ua(50000),ub(50000),v(50000)
      common/ws1/nrp,lp,j2,iat,izt,iap,izp
      common/ws2/vnorm,vz,vzs,vsop,vnormls
      common/main1/nmax,nq
      common/main2/del,hbsd2m
      common/crmax/rmax_unbound
      ep = ei
      nsk=1
      if(ep.lt.0.) rmax = 10.
      if(ep.gt.0.) rmax = rmax_unbound
c      del = 0.1
      nmax = rmax/del
      fd=nsk*del
c hbsd2m = (hbar**2/2mu), mu = reduced mass
      if(iap.eq.1.and.izp.eq.0) xmu = 1.00866*931.50*(iat)/(iat+1)
      if(iap.eq.1.and.izp.eq.1) xmu = 1.00728*931.50*(iat)/(iat+1)
      if(iap.gt.1) xmu = iap*931.50*(iat)/(iat+iap)
      hbsd2m=197.33**2./(2.*xmu)
      sum2=0.0
      sum1=0.0
5     continue
      nsep=nmax
      babone = 1.
      call well(babone,r0,v,dv)
      if(ep.gt.0..and.v(nmax).gt.ep) print 10
      if(ep.gt.0..and.v(nmax).gt.ep) write(10,10)
10    format(1x,'warning: increase rmax')
      if(ep.lt.0) nsep=r0/del
c      write(6,7777) lp
7777  format(1x,i9)
      call schro(v,ep,nsep,ub,lp,phi)
      nodes = 0
      do 21 i=1,80
      if(ub(i)*ub(i+1).lt.0.) nodes=nodes+1
   21 continue
      call xks(ub)
      call prt
      if(ep.lt.0) go to 99
      babone = 0.
      call well(babone,r0,v,dv)
      call schro(v,ep,nsep,ub,lp,phi2)
      phi3=phi-phi2
      phi3=phi3+int(phi3/3.1416)*3.1416
99    return
      end

      subroutine prt
      IMPLICIT REAL*8 (A-H,O-Z)
      common/cua/dimension ua(50000),ub(50000),v(50000)
      common/main2/del,hbsd2m
c      open(unit=98,file='rad.dat')
      r = 0.
      do 200 i = 1,100
      r = r + del
      wf(i) = ub(i)/r
c      write(98,100) r,ub(i)/r
100   format(f5.2,2f10.3)
200   continue
      close(98)
      return
      end

      subroutine  schro(v,e,nsep,u,l,phi)
      IMPLICIT REAL*8 (A-H,O-Z)
      common/ceb/eb
      dimension v(50000),u(50000),q(50000)
      common/main1/nmax,nq
      common/main2/del,hbsd2m
   12 format(8f10.5)
   14 format(1x,8e13.5)
      k=0
    3 do 1 i=1,nmax
      fi=i
      q(i)=(-e+v(i))*del**2/hbsd2m
    1 continue
      u(1)=0.00001
c     write(6,7777) l,u(1)
7777  format(1x,i7,1x,e12.4)
      u(2)=u(1)*(2.0)**(l+1)
      if(q(nmax).lt.0.0) go to 25
      call dif(q,u,1,nsep,1,up1,sum1,u1)
      u(nmax)=0.0
      u(nmax-1)=0.00001
      call dif(q,u,nmax,nsep,-1,up2,sum2,u2)
      x=u1*u2*(up1*u2-up2*u1)/(del**2*(sum1*u2**2+sum2*u1**2))
      k=k+1
      if(k.gt.20) go to 5
      e=e+x*hbsd2m
      eb = e
      conv=0.1d-4
      if(abs(x).gt.conv) go to 3
    5 t=1.
      nq=nsep-1
      do 101 i=1,nq
  101 u(i)=u(i)*u2/t
      do 102 i=nsep,nmax
  102 u(i)=u(i)*u1/t
      sum=0.0
      do 103 i=1,nmax
  103 sum=sum+u(i)**2
      t=sqrt(sum*del)
  105 do 104 i=1,nmax
  104 u(i)=u(i)/t
   99 return
   25 call dif(q,u,1,nmax,1,up1,sum1,u1)
      pk=sqrt(-q(nmax))
c  the next line normalizes the continuum wave function to
c  sin(k r + delta) / sqrt( k ) at large distances.
      t=sqrt(u1**2+(up1/pk)**2)*sqrt(pk/del)
      phi=atan(u1*pk/up1)
      go to 105
      end

      subroutine dif(q,u,n1,n2,nstep,up1,sum,u1)
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension u(50000),q(50000),s(50000)
      sum=0.0
      y1=u(n1)
      n=n1+nstep
      y2=u(n)
      q1=q(n1)
      q2=q(n)
    1 n=n+nstep
      if((n-n2)/nstep.gt.0) go to 99
      q3=q(n)
      u(n)=(2.*y2-y1+.83333*q2*y2+.08333*q1*y1)/(1.-.08333*q3)
      y1=y2
      y2=u(n)
      q1=q2
      q2=q3
      sum=sum+y2**2
      go to 1
   99 nx=n-3*nstep
      up1=(4.*(y2-y1)-(y2-u(nx)))/(2.*nstep)
      u1=y2
      return
      end

      subroutine xks(ub)
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension ub(50000)
cbab this is the radial me:   <l(i1)| kinetic energy |l(i2)>
cbab expectation value of the kinetic energy in the oscillator basis
      common/crad/xkin,rms,vc,xkint,xmst,nt,nb
      common/ws1/nrp,lp,j2,iat,izt,iap,izp
      common/ws2/vnorm,vz,vzs,vsop,vnormls
      common/main1/nmax,nq
      common/main2/del,hbsd2m
      s1 = 0.
      s2 = 0.
      s3 = 0.
      s4 = 0.
      hbsdm = 197.*197./938.
c      hbsdm = 2.*hbsd2m
      h = del
      xll = lp

      r0c = 1.2*(iat**(1./3.))
      ztot = izt
      do 200 ir = 2,100
      r = h*ir
      x1 = ub(ir)
      x2 = ub(ir)
      ra = r-h
      rb = r
      rc = r+h
      x2a = ub(ir-1)
      x2b = ub(ir)
      x2c = ub(ir+1)
      t1 = (x2a+x2c-2.*x2b)/(h*h)
      t2 = xll*(xll+1)*x2/(r*r)
      x = -x1*(hbsdm/2.)*(t1-t2)
      s1 = s1 + x*h
      s2 = s2 + x1*x2*h
      s3 = s3 + x1*x2*h*r*r

      if(r.ge.r0c) s4 = s4 + x1*x2*h*ztot*(1.44/r)
      if(r.lt.r0c)
     1 s4 = s4 + x1*x2*h*ztot*(1.44/r0c)*(1.5-0.5*(r/r0c)**2)

200   continue
      xkin = s1/s2
      rms = sqrt(s3/s2)
      vc = s4/s2
c      xks = s1
      return
      end
      end subroutine gwf

c      REAL*4 wf(100)

c      call get_wspot_wf(6, 3, 1, 0, 0, 1, 1, 1.0, 0.65, 1.25, wfr,wf)

c676   format('wf = ', f10.3)
c      print 676,wf

c681   format('end code')
c      print 681

      end
