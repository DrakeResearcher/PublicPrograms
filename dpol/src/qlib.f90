      function qmin1(a,b)
      implicit real*16(a-h,o-z)
      real*16 min1
      qmin1 = min1(a,b)
      return
      end
      function qmax1(a,b)
      implicit real*16(a-h,o-z)
      real*16 max1
      qmax1 = max1(a,b)
      return
      end
      function qlog(x)
      implicit real*16(a-h,o-z)
      real*16 log
      qlog = log(x)
      return
      end
      function qsqrt(x)
      real*16 qsqrt,x
      qsqrt = sqrt(x)
      return
      end
      function qexp(x)
      real*16 qexp,x
      qexp = exp(x) 
      return
      end
      function qabs(x)
      real*16 qabs,x
      qabs = abs(x)
      return
      end
