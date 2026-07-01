      subroutine vuhard(
C Read only -
     *     nblock, 
     *     jElem, kIntPt, kLayer, kSecPt, 
     *     lAnneal, stepTime, totalTime, dt, cmname,
     *     nstatev, nfieldv, nprops, 
     *     props, tempOld, tempNew, fieldOld, fieldNew,
     *     stateOld,
     *     eqps, eqpsRate,
C Write only -
     *     yield, dyieldDtemp, dyieldDeqps,
     *     stateNew )
C
      include 'vaba_param.inc'
C
      dimension props(nprops), tempOld(nblock), tempNew(nblock),
     1   fieldOld(nblock,nfieldv), fieldNew(nblock,nfieldv),
     2   stateOld(nblock,nstatev), eqps(nblock), eqpsRate(nblock),
     3   yield(nblock), dyieldDtemp(nblock), dyieldDeqps(nblock,2),
     4   stateNew(nblock,nstatev), jElem(nblock)
C
      character*80 cmname
c      
      real*8 Sy,A,n,eps,arfa,beita,scr,epu,B,Tr,Tref,eta
      real*8 part1,part2,part3t1,part3t2,part3
C
      do 100 i = 1,nblock
          Sy=39
          A=320
          n=0.41
          
          eps=0.001
          arfa=7.5
          beita=0.01
          
          epu=10000
          B=110
          
          Tr=25
          Tref=350
          eta=0.0043
          
          part1=Sy+A*eqps(i)**n
          
          if(eqpsRate(i) .le. eps)then
              part2=0.d0
          elseif ((eqpsRate(i) .gt. eps).and.(eqpsRate(i) .le. epu))then
              part2=(arfa*eqps(i)**n+beita)*log(eqpsRate(i)/eps)
          else
              part2=(arfa*eqps(i)**n+beita)*log(eqpsRate(i)/eps)
     1        +B*log(eqpsRate(i)/epu)
          endif
          
          part3t1=1.d0+exp((Tr-Tref)*eta)
          part3t2=1.d0+exp((tempOld(i)-Tref)*eta)
          part3=part3t1/part3t2
          
          yield(i)=(part1+part2)*part3
          
          if(eqpsRate(i) .le. eps)then
           dyieldDeqps(i,1)=n*A*eqps(i)**(n-1)*part3
           dyieldDeqps(i,2)=0.d0
           dyieldDtemp(i)=-part3/part3t2*eta*exp((tempOld(i)-Tref)*eta)*
     1        part1
          
          elseif ((eqpsRate(i) .gt. eps).and.(eqpsRate(i) .le. epu))then
           dyieldDeqps(i,1)=(n*A*eqps(i)**(n-1)+(n*arfa*eqps(i)**(n-1))*
     1        log(eqpsRate(i)/eps))*part3
           dyieldDeqps(i,2)=(arfa*eqps(i)**n+beita)/eqpsRate(i)*part3
           dyieldDtemp(i)=-part3/part3t2*eta*exp((tempOld(i)-Tref)*eta)*
     1        (part1+part2)
              
          else
           dyieldDeqps(i,1)=(n*A*eqps(i)**(n-1)+(n*arfa*eqps(i)**(n-1))*
     1         log(eqpsRate(i)/eps))*part3         
           dyieldDeqps(i,2)=((arfa*eqps(i)**n+beita)/eqpsRate(i)
     1         +B/eqpsRate(i))*part3
           dyieldDtemp(i)=-part3/part3t2*eta*exp((tempOld(i)-Tref)*eta)*
     1         (part1+part2)              
          endif

  100 continue
C
      return
      end