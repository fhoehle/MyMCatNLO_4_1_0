c The following file collects all source files of the 
c
c   CHAPLIN library -- Stephan Buehler, Claude Duhr, 1106.5739 [hep-ph]
c
c version 1.1. They are collected here in order to limit the number
c of external libraries needed by the user when running MC@NLO.
c The present file is linked only when HPL's are actually needed

c-source-file aux.f
C=============================================================================
C---  Auxiliary functions
C=============================================================================

c ------------------------------------------------
c --- some auxiliary functions:
      integer function s(x)
      implicit none
      double complex x
      integer ris
      ris = 1
      if (dimag(x).lt.0d0) then
         ris = -1
      endif
      s=ris
      return
      end function
c ------------------------------------------------
c$$$c s1(x) = s(x^2)...if the imaginary part of x^2 is zero, it depends
c$$$c on the value of the real part of x if we return -1 or 1 because:
c$$$c (x+I*eps) = x^2 + 2*x*I*eps - eps^2
c$$$      integer function s1(x)
c$$$      implicit none
c$$$      double complex x
c$$$      integer ris
c$$$      ris = 1
c$$$      if (dimag(x**2).lt.0d0) then
c$$$         ris = -1
c$$$      else if (dimag(x**2).eq.0d0.and.dreal(x).lt.0d0) then
c$$$         ris = -1
c$$$      endif
c$$$      s1=ris
c$$$      return
c$$$      end function
c$$$c ------------------------------------------------
c$$$c s2(x) = s(4x/(1+x)^2)...if the imaginary part of the arg is zero, it depends
c$$$c on the value of the real part of x if we return -1 or 1 because:
c$$$c 4*(x+I*eps)/(1+x+I*eps)^2 = 4x/(1+x)^2 + 4*(1-x)*I*eps/(1+x)^3
c$$$      integer function s2(x)
c$$$      implicit none
c$$$      double complex x,arg
c$$$      integer ris
c$$$      arg = 4d0*x/(1d0+x)**2
c$$$      ris = 1
c$$$      if (dimag(arg).lt.0d0) then
c$$$         ris = -1
c$$$      else if (dimag(arg).eq.0d0.and.dreal(1d0-x).lt.0d0) then
c$$$         ris = -1
c$$$      endif
c$$$      s2=ris
c$$$      return
c$$$      end function
c$$$c ------------------------------------------------
c-source-file basisexp.f
C=============================================================================
C---  basis expansions
C=============================================================================

C---- expansion of dilogarithm in y = - log(1-z) with Bernoulli numbers  

      double  complex function bsli2_inside(z)
      implicit none
      integer i, Nmax
      double complex ris, z, zb 
      double precision bern(11)

c     bern(i+1) = BernoulliB(2i)/(2i)!
      data bern /1.D0,0.8333333333333333D-1,-0.1388888888888889D-2
     &,0.3306878306878307D-4,-0.8267195767195767D-6,0.208767569878681D-7
     &,-0.5284190138687493D-9,0.1338253653068468D-10
     &,-0.3389680296322583D-12,0.8586062056277845D-14
     &,-0.2174868698558062D-15/
      parameter (Nmax=11)       ! this is half the order we want (coz odd bernoulli numbers are zero except BernoulliB(1)=-0.5d0)
      
      zb = dcmplx(1d0,0d0)-z
      zb = -log(zb)
      ris = -zb**2/4d0          !accounting for BernoulliB(1) = -0.5d0
      do i=1,Nmax
         ris = ris + zb**(2*i-1)*bern(i)/(2*i-1)
      enddo
      bsli2_inside = ris
      return 
      end

C---- expansion of the dilogarithm in log(z) with Zeta values  
C-------- used for border < |z| < 1
      
      double  complex function bsli2_outside(z)
      implicit none
      integer i, Nmax
      double complex ris, z, zb
      double precision zeta(29),zeta0,zeta2
      
c     zeta(i) = Zeta(2-2i-1)/(2i+1)! i.e. Zeta(-1)/6, Zeta(-3)/120, Zeta(-5)/7!....
      data zeta /-0.01388888888888889d0,0.00006944444444444444d0
     &,-7.873519778281683d-7,1.148221634332745d-8,-1.897886998897100d-10
     &,3.387301370953521d-12,-6.372636443183180d-14,1.246205991295067d-
     &15,-2.510544460899955d-17,5.178258806090624d-19,-1.088735736830085
     &d-20,2.325744114302087d-22,-5.035195213147390d-24,1.10264992943812
     &2d-25,-2.438658550900734d-27,5.440142678856252d-29,-1.222834013121
     &735d-30,2.767263468967951d-32,-6.300090591832014d-34,1.44208683884
     &1848d-35,-3.317093999159543d-37,7.663913557920658d-39,-1.777871473
     &383066d-40,4.139605898234137d-42,-9.671557036081102d-44,2.26671870
     &1676613d-45,-5.327956311328254d-47,1.255724838956433d-48,-2.967000
     &542247094d-50/

      parameter (Nmax=29) ! this is half the order we want (coz even zetaval2 are zero except for 0,2)
      parameter (zeta0 = 1.644934066848226d0)
      parameter (zeta2 = -0.2500000000000000d0)

      zb = log(z)
      ris = dcmplx(zeta0, 0d0) + zb*(1d0 -log(-zb)) 
     &     + zb**2*zeta2
      do i=1,Nmax 
         ris = ris + zb**(2*i+1)*zeta(i)
      enddo
      
      bsli2_outside=ris 
      return 
      end

C---- expansion of trilogarithm in y = - log(1-z) with Bernoulli numbers  
      
      double  complex function bsli3_inside(z)
      implicit none
      integer i, Nmax
      double complex ris, z, zb 
      double precision bern(21)

c     bern(n+1) = Sum(Bern(n-k)*Bern(k)/(k+1)!(n-k)!,k=0..n)
      data bern /1d0,-0.7500000000000000d0,0.2361111111111111d0,-0.0347
     &2222222222222d0,0.0006481481481481481d0,0.0004861111111111111d0,-
     &0.00002393550012597632d0,-0.00001062925170068027d0,7.794784580498
     &866d-7,2.526087595532040d-7,-2.359163915200471d-8,-6.168132746415
     &575d-9,6.824456748981078d-10,1.524285616929085d-10,-1.91690941417
     &4054d-11,-3.791718683693992d-12,5.277408409541286d-13,9.471165533
     &842511d-13,-1.432311114490360d-14,-2.372464515550457d-15,3.846565
     &792753191d-16/
      
      parameter (Nmax=21)

      zb = dcmplx(1d0,0d0)-z
      zb = -log(zb)
      ris = dcmplx(0d0, 0d0)
      do i=1,Nmax
         ris = ris + zb**(i)*bern(i)/(i)
      enddo
      bsli3_inside=ris 
      return 
      end

C---- expansion of the trilogarithm in log(z) with Zeta values  
C-------- used for border < |z| < 1
      
      double  complex function bsli3_outside(z)
      implicit none
      integer i, Nmax
      double complex ris, z, zb
      double precision zeta(29),zeta0,zeta1,zeta3

c     zeta(i) = Zeta(3-2i-2)/(2i+2)! i.e. Zeta(-1)/24, Zeta(-3)/6!, Zeta(-5)/8!,....
      data zeta /-0.003472222222222222d0,0.00001157407407407407d0,-9.84
     &1899722852104d-8,1.148221634332745d-9,-1.581572499080917d-11,2.41
     &9500979252515d-13,-3.982897776989488d-15,6.923366618305929d-17,-1
     &.255272230449977d-18,2.353754002768465d-20,-4.536398903458687d-22
     &,8.945169670392643d-24,-1.798284004695496d-25,3.675499764793738d-
     &27,-7.620807971564795d-29,1.600041964369486d-30,-3.39676114756037
     &6d-32,7.282272286757765d-34,-1.575022647958003d-35,3.433540092480
     &589d-37,-7.538849998089870d-39,1.666068164765360d-40,-3.703898902
     &881387d-42,8.279211796468275d-44,-1.859914814630981d-45,4.1976272
     &25327060d-47,-9.514207698800454d-49,2.165042825786954d-50,-4.9450
     &00903745158d-52/

      parameter (zeta0 = 1.202056903159594d0)
      parameter (zeta1 = 1.644934066848226d0)
      parameter (zeta3 = -0.08333333333333333d0)
      parameter (Nmax=30)       ! again, half the order, coz odd zetas are zero except 1,3
      zb = log(z)
      ris = dcmplx(zeta0, 0d0) + zb*zeta1 + zb**3*zeta3
     &     + zb**2*(1d0 + 0.5d0 -log(-zb))/2d0
      do i=2,Nmax 
         ris = ris + zb**(2*i)*zeta(i-1)
      enddo
      bsli3_outside=ris 
      return 
      end

C---- expansion of tetralogarithm in y = - log(1-z) with Bernoulli numbers  
      
      double  complex function bsli4_inside(z)
      implicit none
      integer i, Nmax
      double complex  ris, z, zb 
      double precision bern(21)

c     bern(n+1) = Sum(Sum(Bern(n-k1)*Bern(k2)*Bern(k1-k2)/(k1+1)(k2+1)!(n-k1)!(k1-k2)!,k2=0..k1),k1=0..n)
      data bern /1d0,-0.8750000000000000d0,0.3495370370370370d0,-0.0792
     &8240740740741d0,0.009639660493827160d0,-0.0001863425925925926d0,-
     &0.0001093680638040048d0,0.000006788098837418565d0,0.0000020618654
     &94287074d0,-2.183261421852692d-7,-4.271107367089217d-8,6.53555052
     &3864399d-9,9.049046773887543d-10,-1.872603276102330d-10,-1.917727
     &902789986d-11,5.216900572839828d-12,4.020087098665104d-13,-1.4261
     &64321965609d-13,-8.256053984896996d-15,3.847254012507184d-15,1.64
     &0607009971150d-16/

      parameter (Nmax=20)
      zb = dcmplx(1d0,0d0)-z
      zb = -log(zb)
      ris = dcmplx(0d0, 0d0)
      do i=0,Nmax
         ris = ris+zb**(i+1)*bern(i+1)/(i+1)
      enddo
      bsli4_inside=ris 
      return 
      end

C---- expansion of tetralogarithm in y = log(z) with Zeta values  
C-------- used for 0.3 < |z| < 1
      
      double  complex function bsli4_outside(z)
      implicit none
      integer i, Nmax
      double complex ris, z, zb
      double precision zeta(28),zeta0,zeta1,zeta2,zeta4

c     zeta(i) = Zeta(4-2i-3)/(2i+3)! i.e. Zeta(-1)/120, Zeta(-3)/7!, Zeta(-5)/9!,....
      data zeta /-0.0006944444444444444d0,0.000001653439153439153d0,-1.
     &093544413650234d-8,1.043837849393405d-10,-1.216594230062244d-12,1
     &.613000652835010d-14,-2.342881045287934d-16,3.643877167529436d-18
     &,-5.977486811666558d-20,1.023371305551507d-21,-1.814559561383475d
     &-23,3.313025803849127d-25,-6.200979326536194d-27,1.18564508541733
     &5d-28,-2.309335748959029d-30,4.571548469627103d-32,-9.18043553394
     &6961d-34,1.867249304296863d-35,-3.841518653556106d-37,7.984976959
     &257184d-39,-1.675299999575527d-40,3.544825882479490d-42,-7.558977
     &352819157d-44,1.623374862052603d-45,-3.509273235152795d-47,7.6320
     &49500594654d-49,-1.669159245403588d-50,3.669564111503313d-52/
      
      parameter (zeta0 = 1.082323233711138d0)
      parameter (zeta1 = 1.202056903159594d0)
      parameter (zeta2 = 0.8224670334241131d0)
      parameter (zeta4 = -0.02083333333333333d0)
      parameter (Nmax=29)! half the order, again
      zb = log(z)
      ris = dcmplx(zeta0, 0d0)+zb*zeta1+zb**2*zeta2+zb**4*zeta4
     &     + zb**3*(1d0 + 0.5d0 + 1d0/3d0 -log(-zb))/6d0
      do i=2,Nmax 
         ris = ris+zb**(2*i+1)*zeta(i-1)
      enddo
      bsli4_outside=ris 
      return 
      end

C---- expansion of H2m2(z) = -Li22(-1,z) in y=-log(1+z)
C---- requires the routine bsh2m2_inside_coeff in li22coeff.F for the coefficients
C---- Nmax is the highest order of the taylor expansion

      double complex function bsh2m2_inside(z)
      implicit none
      double complex z,y,ris
      double precision coeff(60)
      integer n, Nmax

c     array starts with former bsh2m2_inside_coeff(1) (zeroth entry is 0)
c     coeff(n) = Convulute[BernoulliB,Car[Convolute[l,b]]][k]/(k+1)!
      data coeff /0.25d0,-0.3333333333333333d0,0.32465277777777777d0,-0
     &.30625d0,0.30379243827160493d0,-0.31872354497354497d0,0.349427536
     &06072058d0,-0.39610441205803308d0,0.460888390823150247d0,-0.54763
     &208903089855d0,0.662009111119538934d0,-0.81188925838826538d0,1.00
     &79564336222981625d0,-1.264599649091272633d0,1.6011498028257904542
     &d0,-2.043574582641123829d0,2.6267914678514710834d0,-3.39782126623
     &4530619d0,4.4200891674751389159d0,-5.779296054536628214d0,7.59144
     &20822867649412d0,-10.01380417551731684d0,13.259972434270600399d0,
     &-17.62046979451679106d0,23.491059571267035452d0,-31.4116491068060
     &3666d0,42.119811446777543608d0,-56.62449156061816696d0,76.3076073
     &22111233016d0,-103.0642326566170164d0,139.49618771598667273d0,-18
     &9.1796144758557257d0,257.03512190956252314d0,-349.840231144116085
     &0d0,476.93937865631952298d0,-651.2283789979552836d0,890.520432944
     &17179905d0,-1219.442885038297044d0,1672.0727345109416041d0,-2295.
     &601033726150628d0,3155.4310810264688627d0,-4342.275779476606369d0
     &,5982.0439827309674170d0,-8249.619738145092382d0,11388.0780556475
     &20056d0,-15735.49668840664380d0,21762.386319547507958d0,-30123.97
     &109708559602d0,41733.247558206646695d0,-57863.12925898186469d0,80
     &289.32329366449d0,-111490.2720614106d0,154927.0759149843d0,-21543
     &5.5591115531d0,299775.6357407621d0,-417401.3987681358d0,581541.04
     &13407259d0,-810711.850057715d0,1.130846347545417d6,-1.57827721432
     &5144d6/

      parameter (Nmax = 60)
      
      y = dcmplx(1d0,0d0)+z
      y = -log(y)
      ris = dcmplx(0d0,0d0)
      do n=1,Nmax
         ris = ris + y**(n+1)*coeff(n)
      enddo
      bsh2m2_inside = ris
      return 
      end

C---- expansion of H2m2(z) = -Li22(-1,z) in y = log(z) (and Re(z) >= 0)
C---- requires the routine bsh2m2_outside_coeff in li22coeff.F for the coefficients
C---- Nmax is the highest order of the taylor expansion

      double complex function bsh2m2_outside(z)
      implicit none
      double complex z,y,ris,cli2,cli4
      double precision coeff(61),pi,zeta3,ll2
      integer n, Nmax

c     coeff(k+1) = Convulute[Car[del],Bar[BernoulliB]][k]/(k+2)!
      data coeff /-0.34657359027997265d0,-0.09942893171332877d0,-0.0187
     &02410976110731d0,-0.00208333333333333333d0,-0.0000489283712703769
     &565d0,0.00002066798941798941798d0,1.2884145995370693d-6,-4.592886
     &537330981d-7,-3.841355103101167d-8,1.5187840708674042d-8,1.311812
     &0316900653d-9,-6.557442900035492d-10,-5.153945830168646d-11,3.365
     &5258621402486d-11,2.3077743717648296d-12,-1.931471133735372d-12,-
     &1.155139350926238d-13,1.1964852873441280d-13,6.3245372250613078d-
     &15,-7.839892377577508d-15,-3.716509741676115d-16,5.36648237346418
     &43d-16,2.3098164428978878d-17,-3.805635847779211d-17,-1.501916986
     &436523d-18,2.7792170059732893d-18,1.0136168766044270d-19,-2.08071
     &4243170466d-19,-7.057701198743988d-21,1.59134829672210430d-20,5.0
     &4683282290655840d-22,-1.2398157254237306d-21,-3.6929281941419875d
     &-23,9.81732678916664963d-23,2.75713292850179439d-24,-7.8859343489
     &266344d-24,-2.0953089297831929d-25,6.41581981843306157d-25,1.6176
     &5642120919963d-26,-5.2797430490927965d-26,-1.2666416605710754d-27
     &,4.38978751038320894d-27,1.00447851365803680d-28,-3.6840236598483
     &991d-28,-8.0579510486468311d-30,3.11806182115730762d-29,6.5321288
     &5795057049d-31,-2.6595787179672540d-30,-5.3461220898892793d-32,2.
     &28469700745414128d-31,4.41401503099778189d-33,-1.975545204274158d
     &-32,-3.673985220125961d-34,1.718584764348089d-33,3.08093293007346
     &2d-35,-1.503444991124448d-34,-2.601539822964359d-36,1.32209980521
     &4306d-35,2.210902587052212d-37,-1.168278254554942d-36,-1.89020535
     &3734545d-38/

      parameter (Nmax = 61)
      parameter (pi=3.1415926535897932385D0)
      parameter (zeta3=1.20205690315959428539973816151d0)
      
      ll2 = dlog(2d0)
      y = log(z)
      ris = dcmplx(0d0,0d0)
      do n=1,Nmax
         ris = ris + y**(n+1)*coeff(n)
      enddo
      ! additional pieces (not part of the sum)
      ris = ris + 71d0/1440d0*pi**4 + 1d0/6d0*pi**2*ll2**2 
     &     - ll2**4/6d0 - 4d0*cli4(dcmplx(0.5d0,0d0)) 
     &     - 7d0/2d0*ll2*zeta3 
     &     - 5d0/8d0*zeta3*log(z) + pi**2/12d0*cli2(z) - pi**4/72d0
      
      bsh2m2_outside = ris
      return 
      end


C---- expansion of H_2,1,-1(z) in y = log(1-z)
C---- requires the routine bsh21m1_inside_coeff in li22coeff.F for the coefficients
C---- Nmax is the highest order of the taylor expansion

      double complex function bsh21m1_inside(z)
      implicit none
      double complex z,y,ris
      double precision coeff(61)
      integer n, Nmax

c     coeff(k+1) = Convolute[Car[L},BernoulliB][k]/(k+1)!
      data coeff /0d0,0d0,0.055555555555555555d0,-0.041666666666666667d0
     &,0.021111111111111111d0,-0.011342592592592592d0,0.0073932350718065
     &003d0,-0.0055762235449735449d0,0.0046271188516558886d0,-0.00411593
     &36419753086d0,0.0038629405982625679d0,-0.0037835075429114780d0,0.0
     &038369243287143673d0,-0.0040054452150998234d0,0.004285098530729088
     &1d0,-0.0046816317089567285d0,0.0052089445088451841d0,-0.0058889221
     &403273633d0,0.0067522208567755897d0,-0.0078398552418544069d0,0.009
     &2055983909726064d0,-0.010919316080588861d0,0.013071451431842447d0,
     &-0.015778977915698980d0,0.019193260007986854d0,-0.0235104153561822
     &82d0,0.028984974459923778d0,-0.035947901615873312d0,0.044830397737
     &630808d0,-0.056195382909281664d0,0.070779196153853163d0,-0.0895469
     &08495008895d0,0.113765799431433695d0,-0.14510309989558134d0,0.1857
     &5619678426856d0,-0.23862631525569750d0,0.30755050306441750d0,-0.39
     &761188650052154d0,0.51555512609759617d0,-0.67034341981460265d0,0.8
     &7390616355980603d0,-1.1421436848979772d0,1.4962789532328810d0,-1.9
     &646780730203453d0,2.5853047361479493d0,-3.4090328131501249d0,4.504
     &1215980597537d0,-5.9622676828419232d0,7.9067966940867543d0,-10.503
     &761791350370d0,13.9769939331314087d0,-18.62852893150884d0,24.86635
     &593168189d0,-33.24214291737328d0,44.50256824694052d0,-59.659220779
     &7288d0,80.0838592330268d0,-107.6383289103828d0,144.8518753970011d0
     &,-195.1633208123604d0,263.2520618632051d0/
      
      parameter (Nmax= 61)
      
      y = dcmplx(1d0,0d0)-z
      y = -log(y)
      ris = dcmplx(0d0,0d0)
      do n=1,Nmax
         ris = ris + y**(n)*coeff(n)
      enddo
      bsh21m1_inside = ris
      return
      end

C---- expansion of H_2,1,-1(z) in y = log(z) for Re(z) >= 0
C---- requires the routine bsh21m1_outside_1_coeff in li22coeff.F for the coefficients
C---- Nmax is the highest order of the taylor expansion

      double complex function bsh21m1_outside_1(z)
      implicit none
      double complex z,y,ris,cli2,cli4,cli3
      double precision coeff(61),pi,zeta3,ll2
      integer n, Nmax

c     coeff(k+1) = Convolute[Car[Convolute[Car[eta],Bar[BernoulliB]]],Bar[BernoulliB]][k]/(k+2)!
      data coeff /0.25d0,0.072916666666666667d0,0.014178240740740740d0,0
     &.0017144097222222222d0,0.00007301311728395061d0,-0.000012504133597
     &88359d0,-1.4500984752366555d-6,2.23621289275741160d-7,3.6408877787
     &2260015d-8,-6.75558808154294265d-9,-1.10997077282832335d-9,2.92072
     &909317950879d-10,4.07040664729507986d-11,-1.54232004675897065d-11,
     &-1.75868046912178768d-12,9.08205825364570960d-13,8.669694078292288
     &52d-14,-5.72693997823042031d-14,-4.72265764364769275d-15,3.7961006
     &7543463165d-15,2.77257823938750085d-16,-2.61818556127472709d-16,-1
     &.72389951229310029d-17,1.86623221122622925d-17,1.12174938342744883
     &d-18,-1.36784448627205447d-18,-7.57565156220199884d-20,1.026798362
     &51631059d-19,5.27780444854412377d-21,-7.86896678494559070d-21,-3.7
     &7574503246258978d-22,6.14038472358867460d-22,2.76380954800098593d-
     &23,-4.86830965140488328d-23,-2.06403852393105725d-24,3.91454487572
     &201673d-24,1.56894524286430607d-25,-3.18746217044925224d-25,-1.211
     &51253199583864d-26,2.62488166036095086d-26,9.48775732188417713d-28
     &,-2.18371877690768050d-27,-7.52502306611488092d-29,1.8335540084562
     &6932d-28,6.03726988240999908d-30,-1.55254146773965411d-29,-4.89454
     &744652263780d-31,1.32474791719770610d-30,4.00620187106289264d-32,-
     &1.13838724502036771d-31,-3.30795228544275656d-33,9.8462750253901d-
     &33,2.753533737449097d-34,-8.56771228709784d-34,-2.309189265362991d
     &-35,7.496835279005074d-35,1.949976046480221d-36,-6.59387311821575d
     &-36,-1.657248185332491d-37,5.827730590842361d-37,1.4169158195933d-
     &38/

      parameter (pi=3.1415926535897932385D0)
      parameter (zeta3=1.20205690315959428539973816151d0)
      parameter (Nmax = 61)

      ll2 = dlog(2d0)
      y = log(z)
      ris = dcmplx(0d0,0d0)
      do n=1,Nmax
         ris = ris + y**(n+1)*coeff(n)
      enddo
!     additional pieces (not part of the sum)
      ris = ris - pi**4/80d0 + pi**2/12d0*ll2**2+ll2**4/24d0 
     &     + cli4(dcmplx(0.5d0,0d0)) + 7d0/8d0*ll2*zeta3 
     &     + y*(7d0*zeta3/8d0 + ll2**3/6d0 
     &     - pi**2/12d0*ll2) - (0.5d0*ll2**2 - pi**2/12d0)
     &     *(pi**2/6d0 - cli2(z)) + ll2*(-cli3(1d0 - z) 
     &     + cli2(1d0 - z)*log(1d0 - z) + 0.5d0*y*log(1d0-z)**2)
   
      bsh21m1_outside_1 = ris
      return
      end

C---- expansion of H_2,1,-1(z) in y = log(-z) for Re(z) <= 0
C---- requires the routine bsh21m1_outside_2_coeff in li22coeff.F for the coefficients
C---- Nmax is the highest order of the taylor expansion

      double complex function bsh21m1_outside_2(z)
      implicit none
      double complex z,y,ris,cli2,cli4
      double precision coeff1(61),coeff2(61),coeff3(61),pi,zeta3,ll2
      integer n, Nmax

c     coeff1(k+1) = 1/4*Convolute[Car[gam],Bar[g]][k]/(k+2)!
      data coeff1 /0.d0,0.041666666666666664d0,0.015625d0,0.0015625d0,-0
     &.00043402777777777775d0,-0.00009300595238095238d0,0.00002170138888
     &888889d0,6.329571759259259d-6,-1.3175843253968255d-6,-4.6385544432
     &41943d-7,8.898717666078777d-8,3.573011792412834d-8,-6.440106731525
     &382d-9,-2.8542441194350915d-9,4.893886071031408d-10,2.344413033238
     &938d-10,-3.856644066310772d-11,-1.968447326291694d-11,3.1260778689
     &1506d-12,1.6824596996814433d-12,-2.591491983625591d-13,-1.45933070
     &17835898d-13,2.1881086259469504d-14,1.2815326053750633d-14,-1.8759
     &379346457062d-15,-1.1373075263231321d-15,1.6291907312337193d-16,1.
     &0185168333506942d-16,-1.4306199549864602d-17,-9.19374092648735d-18
     &,1.2683309378387418d-18,8.3566882828284155d-19,-1.133901104753092d
     &-19,-7.642747026760038d-20,1.0212284628172237d-20,7.02828640619486
     &2d-21,-9.258027908274031d-22,-6.495160220257945d-22,8.442309112741
     &091d-23,6.029235640418441d-23,-7.739195375033816d-24,-5.6193649124
     &97749d-24,7.128585780343022d-25,5.2566533649945075d-25,-6.59470075
     &137659d-26,-4.93392967748237d-26,6.1250097873844946d-27,4.64535674
     &1659573d-27,-5.709457821604045d-28,-4.386141875296571d-28,5.339898
     &687147003d-29,4.152318311693282d-29,-5.0096745880731456d-30,-3.940
     &5785329265186d-30,4.713300047174754d-31,3.748145794905431d-31,-4.4
     &46221644935614d-32,-3.572674119099468d-32,4.204633438143212d-33,3.
     &412169733875996d-33,-3.985333763003178d-34/

c     coeff2(k+1) = 1/4*Convolute[Car[Car[gam]],Bar[g]][k]/(k+2)!
      data coeff2 /0.d0,0.041666666666666664d0,0.013020833333333334d0,0.
     &00078125d0,-0.0003689236111111111d0,-0.000038752480158730157d0,0.0
     &00019117890211640213d0,2.3861480838477368d-6,-1.1894858493165786d-
     &6,-1.6288763440456148d-7,8.170640766126877d-8,1.1871052215290496d-
     &8,-5.985996641481926d-9,-9.062814050215188d-10,4.5909312190151783d
     &-10,7.164860053216312d-11,-3.6439614891245165d-11,-5.8209617103446
     &45d-12,2.970688033325715d-12,4.833759556181296d-13,-2.474257822461
     &5757d-13,-4.08667204085739d-14,2.0972977936448042d-14,3.5072146563
     &23402d-15,-1.8040269804842876d-15,-3.0483852733626177d-16,1.571171
     &1182980456d-16,2.6786142826729357d-17,-1.3830500796359254d-17,-2.3
     &760539988433582d-18,1.2287808333254907d-18,2.1251906928701036d-19,
     &-1.1006142352006812d-19,-1.914758950919985d-20,9.929086819155697d-
     &21,1.736415980225799d-21,-9.014761409182746d-22,-1.583875009223077
     &d-22,8.231536213165235d-23,1.4523246009621054d-23,-7.5551535337982
     &56d-24,-1.3380308874662655d-24,6.966751883890051d-25,1.23805967731
     &85416d-25,-6.451482502735583d-26,-1.150078830970023d-26,5.99752345
     &9709979d-27,1.0722154394963392d-27,-5.59536576479478d-28,-1.002951
     &0870647065d-28,5.237288869237119d-29,9.410440074501714d-30,-4.9169
     &70159919398d-30,-8.854705143231399d-31,4.629190652393522d-31,8.353
     &804322882449d-32,-4.369610683008465d-32,-7.900628258973509d-33,4.1
     &34597171055496d-33,7.489193910055058d-34,-3.921089311632362d-34/

c     coeff3(k+1) = 1/4*Convolute[Car[Convolute[Car[zeta0],Bar[g]]],Bar[g]][k]/(k+2)!
      data coeff3 /0.d0,0.d0,-0.005208333333333333d0,-0.0027777777777777
     &78d0,-0.00044849537037037037d0,0.00007171792328042328d0,0.00002896
     &102017195767d0,-3.346929371028247d-6,-2.0962631290942304d-6,1.8625
     &329226457698d-7,1.60459779406642d-7,-1.1441181000133152d-8,-1.2764
     &225130425063d-8,7.505939634206319d-10,1.0452440048137503d-9,-5.160
     &657975459795d-11,-8.75648538072804d-11,3.673361812319829d-12,7.471
     &673863414197d-12,-2.683984455617326d-13,-6.47253582686209d-13,2.00
     &03913648526046d-14,5.678497945410287d-14,-1.513304268111153d-15,-5
     &.035798991109289d-15,1.1572839941957358d-16,4.507394710505403d-16,
     &-8.914410489716942d-18,-4.0670492949357045d-17,6.892965976882595d-
     &19,3.695736121536782d-18,-5.331738527522966d-20,-3.379368570142838
     &5d-19,4.109583231716237d-21,3.107325998255573d-20,-3.1415554138242
     &515d-22,-2.8714707772078157d-21,2.366884399780646d-23,2.6654773878
     &106823d-22,-1.7413132417682244d-24,-2.4843649215594016d-23,1.23209
     &56980933767d-25,2.3241652071346157d-24,-8.146550067702229d-27,-2.1
     &81683231858151d-25,4.703074070673616d-28,2.0543215403071973d-26,-1
     &.8468288618937404d-29,-1.939949730058493d-27,-5.3388442090463d-31,
     &1.8368062259364445d-28,2.526733885848279d-31,-1.7434234837116958d-
     &29,-4.201139772094111d-32,1.658570420649873d-30,5.61252287055187d-
     &33,-1.5812076638267927d-31,-6.805461760217498d-34,1.51045244111749
     &38d-32,7.816064202735668d-35,-1.4455508007435913d-33/

      parameter (pi=3.1415926535897932385D0)
      parameter (zeta3=1.20205690315959428539973816151d0)
      parameter (Nmax = 61)
      
      ll2 = dlog(2d0)
      y = log(-z)
      ris = dcmplx(0d0,0d0)
      
      do n=1,Nmax
         ris = ris + y**(n+1)*(-coeff2(n)-coeff3(n)
     &        +coeff1(n)*(log(-y) - 1d0/n - 1d0/(n+1)) )
      enddo
      
!     additional pieces (not part of the sum)
      ris = ris + pi**4/80d0 - pi**2*ll2**2/24d0 
     &     - ll2**4/12d0 - (pi**2/12d0 
     &     - ll2**2/2d0)*(-pi**2/12d0 -ll2*y - cli2(z)) 
     &     - 2d0*cli4(dcmplx(0.5d0,0d0)) 
     &     - y*(-ll2**3/6d0 + zeta3/8d0)
      
      bsh21m1_outside_2 = ris 
      return
      end
c-source-file basis.f
C=============================================================================
C---  basis mappings
C=============================================================================

c ---------------------------------------------------------
      double complex function basis2_1(x)
      implicit none
      double complex x,cli2            
      basis2_1=cli2(x)
      return
      end
c ---------------------------------------------------------
      double complex function basis2_2(x)
      implicit none
      double complex x,cli2
      basis2_2=cli2(-x)
      return
      end
c ---------------------------------------------------------
      double complex function basis2_3(x)
      implicit none
      double complex x,cli2
      basis2_3=cli2((1d0-x)/2d0)
      return
      end
c ---------------------------------------------------------
c ---------------------------------------------------------
c     basis3_1(z) = cli3(z) 
      double complex function basis3_1(x)
      implicit none
      double complex x,cli3           
      basis3_1=cli3(x)
      return
      end
c ---------------------------------------------------------
c     basis3_2(z) = cli3(-z)
      double complex function basis3_2(x)
      implicit none
      double complex x,cli3
      basis3_2=cli3(-x)
      return
      end
c ---------------------------------------------------------
c     basis3_3(z) = cli3(1-z)
      double complex function basis3_3(x)
      implicit none
      double complex x,cli3
      basis3_3 = cli3(1d0-x)
      return
      end
c ---------------------------------------------------------
c     basis3_4(z) = cli3(1/(1+z)) 
      double complex function basis3_4(x)
      implicit none
      double complex x,cli3
      basis3_4 = cli3(1d0/(1d0+x))
      return
      end
c ---------------------------------------------------------
c     basis3_5(z) = cli3((1+z)/2) 
      double complex function basis3_5(x)
      implicit none
      double complex x,cli3
      basis3_5 = cli3((1d0+x)/2d0)
      return
      end
c ---------------------------------------------------------
c     basis3_6(z) = cli3((1-z)/2) 
      double complex function basis3_6(x)
      implicit none
      double complex x,cli3
      basis3_6 = cli3((1d0-x)/2d0)
      return
      end
c ---------------------------------------------------------
c     basis3_7(z) = cli3((1-z)/(1+z)) 
      double complex function basis3_7(x)
      implicit none
      double complex x,cli3
      basis3_7 = cli3((1d0-x)/(1d0+x))
      return
      end
c ---------------------------------------------------------
c     basis3_8(z) = cli3(2z/(z-1))
      double complex function basis3_8(x)
      implicit none
      double complex x,cli3
      basis3_8 = cli3(2d0*x/(x-1d0))
      return
      end
c ---------------------------------------------------------
c ---------------------------------------------------------
c     basis1(x) = cli4(x) 
      double complex function basis1(x)
      implicit none
      double complex x,cli4
      basis1=cli4(x)
      return
      end
c ---------------------------------------------------------
c     basis2(x) = cli4(-x)
      double complex function basis2(x)
      implicit none
      double complex x,cli4
      basis2=cli4(-x)
      return
      end
c ---------------------------------------------------------
c     basis3(x) = cli4(1-x)
      double complex function basis3(x)
      implicit none
      double complex x,cli4
      basis3 = cli4(1d0-x)
      return
      end
c ---------------------------------------------------------
c     basis4(x) = cli4(1/(1+x)) 
      double complex function basis4(x)
      implicit none
      double complex x,cli4
      basis4 = cli4(1d0/(1d0+x))
      return
      end
c ---------------------------------------------------------
c     basis5(x) = cli4(x/(x-1))
      double complex function basis5(x)
      implicit none
      double complex x,cli4
      basis5 = cli4(x/(x-1d0))
      return
      end
c ---------------------------------------------------------
c     basis6(x) = cli4(x/(x+1))
      double complex function basis6(x)
      implicit none
      double complex x,cli4
      basis6 = cli4(x/(1d0+x))
      return
      end
c ---------------------------------------------------------
c     basis7(x) = cli4((1+x)/2) 
      double complex function basis7(x)
      implicit none
      double complex x,cli4
      basis7 = cli4((1d0+x)/2d0)
      return
      end
c ---------------------------------------------------------
c     basis8(x) = cli4((1-x)/2)
      double complex function basis8(x)
      implicit none
      double complex x,cli4
      basis8 = cli4((1d0-x)/2d0)
      return
      end
c ---------------------------------------------------------
c     basis9(x) = cli4((1-x)/(1+x))
      double complex function basis9(x)
      implicit none
      double complex x,cli4
      basis9 = cli4((1d0-x)/(1d0+x))
      return
      end
c ---------------------------------------------------------
c     basis10(x) = cli4((x-1)/(x+1))
      double complex function basis10(x)
      implicit none
      double complex x,cli4
      basis10 = cli4((x-1d0)/(1d0+x))
      return
      end
c ---------------------------------------------------------
c     basis11(x) = cli4(2x/(1+x))
      double complex function basis11(x)
      implicit none
      double complex x,cli4
      basis11 = cli4(2d0*x/(1d0+x))
      return
      end
c ---------------------------------------------------------
c     basis12(x) = cli4(2x/(x-1)) 
      double complex function basis12(x)
      implicit none
      double complex x,cli4
      basis12 = cli4(2d0*x/(x-1d0))
      return
      end
c ---------------------------------------------------------
c     basis13(x) = cli4(1-x^2) = cli4_sbc 
      double complex function basis13(x)
      implicit none
      double complex x,cli4_sbc !,cli4
      basis13=cli4_sbc(x)
      return
      end
c ---------------------------------------------------------
c     basis14(x) = cli4(x^2/(x^2-1)) 
      double complex function basis14(x)
      implicit none
      double complex x,cli4
      basis14 = cli4(x**2/(x**2-1d0))
      return
      end
c ---------------------------------------------------------
c     basis15(x) = cli4(4x/(1+x)^2) = cli4_sbc_2  
      double complex function basis15(x)
      implicit none      
      double complex x,cli4_sbc_2
      basis15=cli4_sbc_2(x)
      return
      end
c ---------------------------------------------------------
c     basis16(x) = ch2m2(x) 
      double complex function basis16(x)
      implicit none
      double complex x,ch2m2
      basis16=ch2m2(x)
      return
      end
c ---------------------------------------------------------
c     basis17(x) = ch21m1(x) 
      double complex function basis17(x)
      implicit none
      double complex x,ch21m1
      basis17=ch21m1(x)
      return
      end
c ---------------------------------------------------------
c     basis18(x) = ch21m1(-x)
      double complex function basis18(x)
      implicit none
      double complex x,ch21m1
      basis18=ch21m1(-x)
      return
      end
c ---------------------------------------------------------
c-source-file basisfct.f
C=============================================================================
C---  basis functions
C=============================================================================
c---  Li2

      double complex  function cli2(z)
      implicit none
      double complex ris, z, bsli2_inside,bsli2_outside, wcli2
      double complex zlocal
      double precision zabs, pi, zeta2, border, tiny, arg

      pi=3.1415926535897932385D0
      zeta2=pi**2/6d0

      border = 0.3d0 
      tiny = 1d-14
      zabs = abs(z)
      zlocal=z

      if (zabs.gt.1d0+tiny) then
         ris=-wcli2(1d0/z)-zeta2-0.5d0*log(-z)**2
      elseif (zabs.le.border) then 
         ris=bsli2_inside(z)
      else
         if (zabs.gt.1d0) then
            arg=atan2(dimag(zlocal),dreal(zlocal))
            zlocal=dcmplx(cos(arg),sin(arg))
         endif
         ris=bsli2_outside(zlocal)
      endif

      cli2=ris
      return
      end
      
c---  recursion
      
      double complex  function wcli2(z)
      implicit none
      double complex z, cli2
      wcli2 =  cli2(z)
      return
      end

c--- Li3

      double complex  function cli3(z)
      implicit none
      double complex ris, z, bsli3_inside,bsli3_outside, wcli3
      double complex zlocal
      double precision zabs,border, pi, zeta2, zeta3,tiny,arg
      
      pi=3.1415926535897932385D0
      zeta2=pi**2/6d0
      zeta3=1.20205690315959428539973816151d0
     
      border = 0.3d0
      zabs = abs(z)
      tiny = 1d-14
      zlocal=z

      if (zabs.gt.1d0+tiny) then
         ris=wcli3(1d0/z)-log(-z)**3/6d0-zeta2*log(-z)
      elseif (zabs.le.border) then 
         ris=bsli3_inside(z)
      else
         if (zabs.gt.1d0) then
            arg=atan2(dimag(zlocal),dreal(zlocal))
            zlocal=dcmplx(cos(arg),sin(arg))
         endif
         ris=bsli3_outside(zlocal)
      endif

      cli3=ris
      return
      end
      
c---  recursion

      double complex  function wcli3(z)
      implicit none
      double complex z, cli3
      wcli3 =  cli3(z)
      return
      end

c--- Li4

      double complex  function cli4(z)
      implicit none
      double complex ris, z, bsli4_outside, bsli4_inside, wcli4
      double complex zlocal
      double precision zabs, pi, zeta2, zeta3, zeta4, border,tiny,arg
     
      pi=3.1415926535897932385D0
      zeta2=pi**2/6d0
      zeta3=1.20205690315959428539973816151d0
      zeta4=pi**4/90d0
      
      border = 0.3d0
      zabs = abs(z)
      tiny = 1d-14
      zlocal=z

      if (zabs.gt.1d0+tiny) then
         ris=-wcli4(1d0/z) -log(-z)**4/24d0 - 7d0*zeta4/4d0 
     &           - zeta2*log(-z)**2/2d0
      elseif (zabs.le.border) then 
         ris=bsli4_inside(z)
      else
         if (zabs.gt.1d0) then
            arg=atan2(dimag(zlocal),dreal(zlocal))
            zlocal=dcmplx(cos(arg),sin(arg))
         endif
         ris=bsli4_outside(zlocal)
      endif

      
      cli4=ris
      return
      end
      
c     --- recursion for li4
      
      double complex  function wcli4(z)
      implicit none
      double complex z, cli4
      wcli4 = cli4(z)
      return
      end
c --- the case Li4(1-z^2) needs some special treatment because of its branch cut structure
c --- (that's what 'sbc' stands for: special branch cut)

      double complex function cli4_sbc(z)
      implicit none
      double complex ris, z, cli4, myi,basis14
      double complex ll1,ll2,ll3
      double precision pi,zabs,zreal
      integer s
      
      pi=3.1415926535897932385D0
      zabs = abs(z)
      zreal = dreal(z)
      myi = dcmplx(0d0,1d0)
           
      if (zabs.le.1d0) then !normal li4
         if (zreal.gt.0d0) then
            ris = cli4(1d0 - z**2)
         else if (zreal.eq.0d0 .and. s(z).eq.1) then !also normal li4
            ris = cli4(1d0 - z**2 - dcmplx(0d0,1d-60))
         else                   ! special branch cut configuration
            ris = cli4(1d0 - z**2- dcmplx(0d0,1d-60))
     &           - myi*pi*s(z)/3d0*(log(1d0 - z)+log(1d0+z))**3 
         endif
      else 
         ll1=log(1d0/z)
         ll2=log(1d0 - 1d0/z)
         ll3=log(1d0 + 1d0/z)
         ris = -2d0/3d0*ll1**4 + 4d0/3d0*ll2
     &        *ll1**3 + 4d0/3d0*ll3*ll1**3 
     &        - ll2**2*ll1**2 - ll3**2
     &        *ll1**2 - 2d0*ll2*ll3
     &        *ll1**2 - pi**2/3d0*ll1**2 + 1d0/3d0
     &        *ll2**3*ll1 + 1d0/3d0
     &        *ll3**3*ll1 + ll2
     &        *ll3**2*ll1 + pi**2/3d0*ll1
     &        *ll2 + ll3*ll2**2
     &        *ll1 + pi**2/3d0*ll1*ll3 
     &        - 1d0/24d0*ll2**4 - 1d0/24d0
     &        *ll3**4 - 1d0/6d0*ll2
     &        *ll3**3 - pi**2/12d0*ll2**2 
     &        - pi**2/12d0*ll3**2 
     &        - 1d0/4d0*ll2**2*ll3**2 
     &        - 1d0/6d0*ll2**3*ll3 
     &        - pi**2/6d0*ll2*ll3 
     &        - 7*pi**4/360d0 
     &        -basis14(1d0/z)
      endif
      
      cli4_sbc = ris
      return
      end

c --- the case Li4(4z/(1+z)^2) also needs some special treatment because of its branch cut structure

      double complex function cli4_sbc_2(z)
      implicit none
      double complex ris, z, cli4, myi, wcli4_sbc_2
      double complex arg,llx,ll1px,wcli4sbc2,zlocal
      double precision pi,zabs,zreal,ll2,tiny,arg2
      integer s
      
      pi=3.1415926535897932385D0
      ll2 = dlog(2d0)
      zabs = abs(z)
      zreal = dreal(z)
      myi = dcmplx(0d0,1d0)
      llx = 1
      ll1px = 1
      wcli4sbc2 = 1
      tiny = 1d-14
      zlocal=z
     
      ris = dcmplx(0d0,0d0)

      if (zabs.lt.1d0) then
         ris = cli4(4d0*z/(1d0+z)**2)
      elseif (zabs.lt.1d0+tiny) then 
         arg2=atan2(dimag(zlocal),dreal(zlocal))
         zlocal=dcmplx(cos(arg2),sin(arg2))
         arg = dcmplx(dreal(4d0*zlocal/(1d0+zlocal)**2),s(zlocal)*1d-60)
         ris = cli4(arg)
      else    
         wcli4sbc2 =  wcli4_sbc_2(1d0/z)
         llx = log(1d0/z)    
         ll1px = log(1d0+1d0/z)          
         ris = wcli4sbc2 + myi*pi*s(z)*
     &(4d0*ll2**2*llx - 8d0*ll2**2* ll1px 
     &+ 2d0*ll2*llx**2 - 8d0*ll2*llx
     &* ll1px + 8d0*ll2*ll1px**2 
     &+ 1d0/3d0*llx**3 - 2d0* ll1px*llx**2 
     &+ 4d0* ll1px**2*llx - 8d0/3d0* ll1px**3 
     &+ 8d0/3d0*ll2**3)


c             ris = wcli4_sbc_2(1d0/z) + myi*pi*s(z)*
c     &(4d0*ll2**2*log(1d0/z) - 8d0*ll2**2*log(1d0+1d0/z) 
c     &+ 2d0*ll2*log(1d0/z)**2 - 8d0*ll2*log(1d0/z)
c     &*log(1d0+1d0/z) + 8d0*ll2*log(1d0+1d0/z)**2 
c     &+ 1d0/3d0*log(1d0/z)**3 - 2d0*log(1d0+1d0/z)*log(1d0/z)**2 
c     &+ 4d0*log(1d0+1d0/z)**2*log(1d0/z) - 8d0/3d0*log(1d0+1d0/z)**3 
c     &+ 8d0/3d0*ll2**3)
      endif
      
      cli4_sbc_2 = ris
      return
      end

c     --- recursion for cli4_sbc_2
      
      double complex  function wcli4_sbc_2(z)
      implicit none
      double complex z, cli4_sbc_2
      wcli4_sbc_2 =  cli4_sbc_2(z)
      return
      end

C-----------------------------------------------------------------------
C     mapping of H_2-2(z) into convergent region
      
      double complex  function ch2m2(z)
      implicit none
      double complex ris,z,bsh2m2_inside,bsh2m2_outside,cli4,cli2
      double complex HPL4,wch2m2,myi,zlocal !,cli4_sbc,cli3
      double precision pi,zeta2,zeta3,zeta4,zabs,zreal,border,tiny,arg
      integer s

      pi=3.1415926535897932385D0
      zeta2=pi**2/6d0
      zeta3=1.20205690315959428539973816151d0
      zeta4=pi**4/90d0
      myi = dcmplx(0d0,1d0)
      tiny=1d-14
      
      border = 0.3d0
      zabs = abs(z)
      zreal = dreal(z)
      zlocal=z


      if (zabs.lt.border) then ! inside circle of |z| = 0.3, we employ the log(1+z) expansion
         ris = bsh2m2_inside(z)            
      elseif (zabs.lt.1d0+tiny) then
         if (zabs.gt.1d0) then
            arg=atan2(dimag(zlocal),dreal(zlocal))
            zlocal=dcmplx(cos(arg),sin(arg))
         endif
         if (zreal.ge.0d0) then ! on the half annulus 0.3 < |z| < 1 ; Re(z) >= 0, we have the log(x) exp.
            ris = bsh2m2_outside(zlocal)
         else                ! for Re(z) < 0, we map back to Re(z) > 0 by using the fact that HPL4(n1,n2,n3,n4,z) = (+-) HPL4(-n1,-n2,-n3,-n4,-z) (if n4 =/= 0):
            ris = HPL4(0,-1,0,1,-zlocal) 
         endif
      else                      ! For |z| > 1, we use the inversion formula to map into the unit circle. 
         ris = dcmplx(0d0,0d0) +
     &        wch2m2(1d0/z) 
     &        + 37d0*pi**4/720d0 
     &        - HPL4(0,1,0,0,1d0/z) 
     &        - log(1d0/z)**4/24d0 
     &        - pi**2/12d0*log(1d0/z)**2 
     &        - pi**2/6d0*cli2(1d0/z) 
     &        - cli4(-1d0/z) 
     &        + 3d0*zeta3*log(1d0/z)/2d0 
     &        - pi**3*myi*s(z)*log(1d0/z)/12d0
      endif
      
      ch2m2=ris
      return
      end
      
      
c     --- recursion for H_2-2(z)
      
      double complex  function wch2m2(z)
      implicit none
      double complex z, ch2m2
           
      wch2m2 = ch2m2(z)
      return
      end

C------------------------------------------------------------------------------
C     mapping of H21-1(z) into convergent region 
      
      double complex  function ch21m1(z)
      implicit none
      double complex ris,z,bsh21m1_inside,bsh21m1_outside_1,zlocal
      double complex bsh21m1_outside_2,cli4,cli2,HPL4,wch21m1,myi,ch2m2
      double precision pi,zeta2,zeta3,zeta4,border,zreal,zabs,ll2,tiny
      double precision arg
      integer s

      pi=3.1415926535897932385D0
      zeta2=pi**2/6d0
      zeta3=1.20205690315959428539973816151d0
      zeta4=pi**4/90d0
      ll2 = dlog(2d0)
      border = 0.3d0
      myi = dcmplx(0d0,1d0)
      tiny=1d-14

      zabs = abs(z)
      zreal = dreal(z)
      zlocal=z


      if (zabs.lt.border) then ! inside circle of |z| = 0.3, we employ the log(1+z) expansion
         ris = bsh21m1_inside(z)           
      elseif (zabs.lt.1d0+tiny) then
         if (zabs.gt.1d0) then
            arg=atan2(dimag(zlocal),dreal(zlocal))
            zlocal=dcmplx(cos(arg),sin(arg))
         endif
         if (zreal.ge.0d0) then ! on the half annulus 0.3 < |z| < 1 ; Re(z) >= 0, we have the log(x) exp.
            ris = bsh21m1_outside_1(zlocal)
         else                ! for Re(z) < 0, we map back to Re(z) > 0 by using the fact that HPL4(n1,n2,n3,n4,z) = (+-) HPL4(-n1,-n2,-n3,-n4,-z) (if n4 =/= 0):
            ris = bsh21m1_outside_2(zlocal)
         endif
      else                      ! For |z| > 1, we use the inversion formula to map into the unit circle. 
         ris = -wch21m1(1d0/z)-pi**4/144d0 -ch2m2(1d0/z) 
     &        + log(1d0/z)**4/24d0 
     &        + pi**2*ll2**2/3d0 - ll2**4/12d0 
     &        + 3d0*pi**2*ll2*log(1d0/z)/4d0 
     &        + pi**2*log(1d0/z)**2/8d0 
     &        + pi**2*cli2(1d0/z)/4d0 
     &        - 2*cli4(dcmplx(0.5d0,0d0)) 
     &        + cli4(-1d0/z) 
     &        - 7d0*zeta3*log(1d0/z)/8d0 
     &        + myi*s(z)*(pi**3*ll2/6d0 
     &        + pi**3*log(1d0/z)/12d0 
     &        - 0.5d0*pi*ll2**2*log(1d0/z) 
     &        - 0.5d0*pi*ll2*log(1d0/z)**2 
     &        - pi*ll2*cli2(1d0/z))
     &        - HPL4(0,0,1,-1,1d0/z) 
     &        + HPL4(0,0,1,0,1d0/z) 
     &        + HPL4(0,1,0,0,1d0/z) 
     &        + HPL4(0,1,1,0,1d0/z)
      endif

      ch21m1=ris
      return
      end
      

c     --- recursion for H21-1(z)
      
      double complex  function wch21m1(z)
      implicit none
      double complex z, ch21m1
            
      wch21m1 =  ch21m1(z)
      return
      end
c-source-file HPL1.f
C=============================================================================
C---  HPLs  of Rank 1  
C=============================================================================
c --- main forking function
      double  complex function HPL1(n1, x)
      implicit none 
      integer n1
      double complex x,ris
      double complex HPL1at0,HPL1at1,HPL1atm1
      double complex HPL1ar1,HPL1arm1,HPL1ar0
      double complex HPL1else
      double precision rad1,radm1,rad0
      
      rad0 = 0.025d0
      rad1 = 0.01d0
      radm1 = 0.025d0      

      if (abs(n1).gt.1) then
         print*, ""
         print*, "****************"
         print*, "Error in HPL1:"
         print*, "Index",n1," out of range !"
         print*, "Aborting..."
         print*,"****************"
         stop
      endif

      ris = dcmplx(0d0,0d0)
      
      if (x.eq.dcmplx(0d0,0d0)) then
         ris = HPL1at0(n1)
      elseif (x.eq.dcmplx(1d0,0d0)) then
         ris = HPL1at1(n1)
      elseif (x.eq.dcmplx(-1d0,0d0)) then
         ris = HPL1atm1(n1)
      elseif (abs(x-dcmplx(1d0,0d0)).lt.rad1) then
         ris = HPL1ar1(n1,x)
      elseif (abs(x+dcmplx(1d0,0d0)).lt.radm1) then
         ris = HPL1arm1(n1,x)
      elseif (abs(x-dcmplx(0d0,0d0)).lt.rad0) then
         ris = HPL1ar0(n1,x)
      else 
         ris = HPL1else(n1,x)
      endif

      HPL1=ris 
      return
      end
c ------------------------------------------------
      double complex function HPL1at0(n1)
      implicit none
      integer n1
      double complex ris

      ris = dcmplx(0d0,0d0)

      if (n1.eq.0) then             
         print*, ""
         print*, "****************"
         print*, "ERROR in HPL1: "
         print*, "HPL1(",n1
     &        ,",0) is divergent!"
         print*, "Aborting..."
         print*,"****************"
         stop
      endif
      HPL1at0=ris
      return
      end function
c ------------------------------------------------
      double complex function HPL1at1(n1)
      implicit none
      integer n1
      double complex ris
      double precision ll2

      ll2 = dlog(2d0)
      ris = dcmplx(0d0,0d0)

      if(n1.ne.1) then
         select case (n1)
         case(-1)
            ris = ll2
         case(0)
            ris = 0d0
         end select
      else
         print*, ""
         print*, "****************"
         print*, "ERROR in HPL1: "
         print*, "HPL1(",n1
     &        ,",1) is divergent!"
         print*, "Aborting..."
         print*,"****************"
         stop
      endif
      HPL1at1=ris
      return
      end function
c ------------------------------------------------
      double complex function HPL1atm1(n1)
      implicit none
      integer n1
      double complex ris,myi
      double precision ll2,pi

      pi=3.1415926535897932385D0
      myi = dcmplx(0d0,1d0)
      ll2 = dlog(2d0)

      ris = dcmplx(0d0,0d0)

      if(n1.ne.-1) then
         select case (n1)
         case(0)
            ris = myi*pi
         case(1)
            ris = ll2
         end select
      else
         print*, ""
         print*, "****************"
         print*, "ERROR in HPL1: "
         print*, "HPL1(",n1
     &        ,",-1) is divergent!"
         print*, "Aborting..."
         print*,"****************"
         stop
      endif
      HPL1atm1=ris
      return
      end function
c ------------------------------------------------
      double complex function HPL1ar1(n1,x)
      implicit none
      integer n1,bcflag
      double complex x,ris,zp,llzp
      double precision pi,ll2,xre

      pi=3.1415926535897932385D0
      ll2 = dlog(2d0)

      ris = dcmplx(0d0,0d0)
      bcflag = 0
      
c---  +i*epsilon to get branch cuts right ---
      if (dimag(x).eq.0d0) then
         x = x + dcmplx(0d0,1d-60)
         bcflag = 1
      endif
c---  
      select case(n1)
         case(-1)            !-1

            zp = 1d0-x

            ris = -((zp)/2d0) - (zp**2)/8d0 - (zp**3)/2
     &4d0 - (zp**4)/64d0 - (zp**5)/160d0 - (zp**6)/384d0 + ll
     &2

         case(0)            !0

            zp = 1d0-x

            ris = -zp - (zp**2)/2d0 - (zp**3)/3d0 - (zp
     &**4)/4d0 - (zp**5)/5d0 - (zp**6)/6d0

         case(1)            !1

            zp = 1d0-x
            llzp = log(zp)

            ris = -llzp
c End of expansions around x = +1
      end select
c --- set the imaginary part back to zero if it has been modified to
c --- get the branch cuts right (and should be zero).
      if (bcflag.eq.1) then
         x = x - dcmplx(0d0,1d-60)
         xre = dreal(x)
         if (n1.eq.0.and.xre.gt.0d0) then
            ris = dcmplx(dreal(ris),0d0)
c     
         else if (n1.eq.1.and.xre.lt.1d0) then
            ris = dcmplx(dreal(ris),0d0)
c            
         else if (n1.eq.-1.and.xre.gt.-1d0) then
            ris = dcmplx(dreal(ris),0d0)
         endif
      endif        
      
      HPL1ar1=ris
      return
      end function
c ------------------------------------------------
      double complex function HPL1arm1(n1,x)
      implicit none
      integer n1,s,szp,bcflag
      double complex x,ris,zp,llzp,myi
      double precision pi,ll2,xre
      
      pi=3.1415926535897932385D0
      ll2 = dlog(2d0)
      myi = dcmplx(0d0,1d0)
      
      ris = dcmplx(0d0,0d0)
      bcflag = 0
      
c---  +i*epsilon to get branch cuts right ---
      if (dimag(x).eq.0d0) then
         x = x + dcmplx(0d0,1d-60)
         bcflag = 1
      endif
c---  
      select case(n1)
      case(-1)                  !-1
         
         zp = x+1d0
         llzp = log(zp)
         
         ris = llzp
         
      case(0)                   !0
         
         zp = x+1d0
         szp = s(zp)
         
         ris = myi*pi*szp - zp - (zp**2)/2d0 - (zp**
     &        3)/3d0 - (zp**4)/4d0 - (zp**5)/5d0 - (zp**6)/6d0 - (zp*
     &        *7)/7d0 - (zp**8)/8d0 - (zp**9)/9d0
         
      case(1)                   !1
         
         zp = x+1d0
         
         ris = (zp)/2d0 + (zp**2)/8d0 + (zp**3)/24d0
     &        + (zp**4)/64d0 + (zp**5)/160d0 + (zp**6)/384d0 + (zp**
     &        7)/896d0 + (zp**8)/2048d0 + (zp**9)/4608d0 - ll2
c     End of expansions around x = -1
      end select
c     --- set the imaginary part back to zero if it has been modified to
c     --- get the branch cuts right (and should be zero).
      if (bcflag.eq.1) then
         x = x - dcmplx(0d0,1d-60)
         xre = dreal(x)
         if (n1.eq.0.and.xre.gt.0d0) then
            ris = dcmplx(dreal(ris),0d0)
c              
         else if (n1.eq.1.and.xre.lt.1d0) then
            ris = dcmplx(dreal(ris),0d0)
c     
         elseif (n1.eq.-1.and.xre.gt.-1d0) then
            ris = dcmplx(dreal(ris),0d0)
         endif
      endif        
      HPL1arm1=ris
      return
      end function
c ------------------------------------------------
      double complex function HPL1ar0(n1,x)
      implicit none
      integer n1,bcflag
      double complex x,ris,llx
      double precision pi,ll2,xre
      
      pi=3.1415926535897932385D0
      ll2 = dlog(2d0)

      ris = dcmplx(0d0,0d0)
      bcflag = 0

c---  +i*epsilon to get branch cuts right ---
      if (dimag(x).eq.0d0) then
         x = x + dcmplx(0d0,1d-60)
         bcflag = 1
      endif
c---  
      select case(n1)
      case(-1)                  !-1
         
         
         ris = x - (x**2)/2d0 + (x**3)/3d0 - (x**4)/
     &        4d0 + (x**5)/5d0 - (x**6)/6d0 + (x**7)/7d0 - (x**8)/8d0
     &        + (x**9)/9d0 - (x**10)/10d0
         
      case(0)                   !0
         
         llx = log(x)
         
         ris = llx
         
      case(1)                   !1
         
         
         ris = x + (x**2)/2d0 + (x**3)/3d0 + (x**4)/
     &        4d0 + (x**5)/5d0 + (x**6)/6d0 + (x**7)/7d0 + (x**8)/8d0
     &        + (x**9)/9d0 + (x**10)/10d0
c     End of expansions around x = 0
      end select
c     --- set the imaginary part back to zero if it has been modified to
c     --- get the branch cuts right (and should be zero).
      if (bcflag.eq.1) then
         x = x - dcmplx(0d0,1d-60)
         xre = dreal(x)
         if (n1.eq.0.and.xre.gt.0d0) then
            ris = dcmplx(dreal(ris),0d0)
c              
         else if (n1.eq.1.and.xre.lt.1d0) then
            ris = dcmplx(dreal(ris),0d0)
c     
         elseif (n1.eq.-1.and.xre.gt.-1d0) then
            ris = dcmplx(dreal(ris),0d0)
         endif
      endif        
      HPL1ar0=ris
      return
      end function
c ------------------------------------------------
      double  complex function HPL1else(n1, x)
      implicit none
      double complex x, ris
      integer n1,bcflag
      double precision xre
      
      bcflag = 0
      ris = dcmplx(0d0,0d0)
      
c---  +i*epsilon to get branch cuts right ---
      if (dimag(x).eq.0d0) then
         x = x + dcmplx(0d0,1d-60)
         bcflag = 1
      endif
c---  

      select case(n1)
      case(-1)
         ris =log(1.0d0 + x)
      case(0)
         ris=log(x)
      case(1)
         ris=-log(1.0d0 - x)
      end select

c     --- set the imaginary part back to zero if it has been modified to
c     --- get the branch cuts right (and should be zero).

      if (bcflag.eq.1) then
         x = x - dcmplx(0d0,1d-60)
         xre = dreal(x)
         if (n1.eq.0.and.xre.gt.0d0) then
            ris = dcmplx(dreal(ris),0d0)
c              
         else if (n1.eq.1.and.xre.lt.1d0) then
            ris = dcmplx(dreal(ris),0d0)
c     
         elseif (n1.eq.-1.and.xre.gt.-1d0) then
            ris = dcmplx(dreal(ris),0d0)
         endif
      endif        

      HPL1else=ris 
      return
      end
c ------------------------------------------------
c --- Real part of HPL1     
      double precision function HPL1real(n1,xr,xi)
      implicit none
      double precision xr,xi
      integer n1
      double complex x,HPL1
      x=dcmplx(xr,xi)
      HPL1real = dreal(HPL1(n1,x))
      return
      end

c --- Imaginary part of HPL1     
      double precision function HPL1im(n1,xr,xi)
      implicit none
      double precision xr,xi
      integer n1
      double complex x,HPL1
      x=dcmplx(xr,xi)
      HPL1im = dimag(HPL1(n1,x))
      return
      end
c-source-file HPL2ar0.f
c ------------------------------------------------
      double complex function HPL2ar0(n1,n2,x)
      implicit none
      integer n1,n2,j,bcflag
      double complex x,ris,myi,llx
      double precision pi, zeta2,ll2,xre

      pi=3.1415926535897932385D0
      zeta2=pi**2/6d0
      myi = dcmplx(0d0,1d0)
      ll2 = dlog(2d0)
      bcflag = 0

      j=1+(n2+1)+(n1+1)*3
      ris = dcmplx(0d0,0d0)

c---  +i*epsilon to get branch cuts right ---
      if (dimag(x).eq.0d0) then
         x = x + dcmplx(0d0,1d-60)
         bcflag = 1
      endif
c---  
      select case(j)
      case(1)                   !-1-1


         ris = (x**2)/2d0 - (x**3)/2d0 + (11d0*x**4)
     &        /24d0 - (5d0*x**5)/12d0 + (137d0*x**6)/360d0 - (7d0*x**
     &        7)/20d0 + (363d0*x**8)/1120d0 - (761d0*x**9)/2520d0 + (
     &        7129d0*x**10)/25200d0
         
      case(2)                   !-10
         
         llx = log(x)
         
         ris = -x + (x**2)/4d0 - (x**3)/9d0 + (x**4)
     &        /16d0 - (x**5)/25d0 + (x**6)/36d0 - (x**7)/49d0 + (x**8
     &        )/64d0 - (x**9)/81d0 + (x**10)/100d0 + x*llx - (x**2*ll
     &        x)/2d0 + (x**3*llx)/3d0 - (x**4*llx)/4d0 + (x**5*llx)/5
     &        d0 - (x**6*llx)/6d0 + (x**7*llx)/7d0 - (x**8*llx)/8d0 +
     &        (x**9*llx)/9d0 - (x**10*llx)/10d0
         
      case(3)                   !-11
         
         
         ris = (x**2)/2d0 - (x**3)/6d0 + (5d0*x**4)/
     &        24d0 - (7d0*x**5)/60d0 + (47d0*x**6)/360d0 - (37d0*x**7
     &        )/420d0 + (319d0*x**8)/3360d0 - (533d0*x**9)/7560d0 + (
     &        1879d0*x**10)/25200d0
         
      case(4)                   !0-1
         
         
         ris = x - (x**2)/4d0 + (x**3)/9d0 - (x**4)/
     &        16d0 + (x**5)/25d0 - (x**6)/36d0 + (x**7)/49d0 - (x**8)
     &        /64d0 + (x**9)/81d0 - (x**10)/100d0
         
      case(5)                   !00
         
         llx = log(x)
         
         ris = (llx**2)/2d0
         
      case(6)                   !01
         
         
         ris = x + (x**2)/4d0 + (x**3)/9d0 + (x**4)/
     &        16d0 + (x**5)/25d0 + (x**6)/36d0 + (x**7)/49d0 + (x**8)
     &        /64d0 + (x**9)/81d0 + (x**10)/100d0
         
      case(7)                   !1-1
         
         
         ris = (x**2)/2d0 + (x**3)/6d0 + (5d0*x**4)/
     &        24d0 + (7d0*x**5)/60d0 + (47d0*x**6)/360d0 + (37d0*x**7
     &        )/420d0 + (319d0*x**8)/3360d0 + (533d0*x**9)/7560d0 + (
     &        1879d0*x**10)/25200d0
         
      case(8)                   !10
         
         llx = log(x)
         
         ris = -x - (x**2)/4d0 - (x**3)/9d0 - (x**4)
     &        /16d0 - (x**5)/25d0 - (x**6)/36d0 - (x**7)/49d0 - (x**8
     &        )/64d0 - (x**9)/81d0 - (x**10)/100d0 + x*llx + (x**2*ll
     &        x)/2d0 + (x**3*llx)/3d0 + (x**4*llx)/4d0 + (x**5*llx)/5
     &        d0 + (x**6*llx)/6d0 + (x**7*llx)/7d0 + (x**8*llx)/8d0 +
     &        (x**9*llx)/9d0 + (x**10*llx)/10d0
         
      case(9)                   !11
         
         
         ris = (x**2)/2d0 + (x**3)/2d0 + (11d0*x**4)
     &        /24d0 + (5d0*x**5)/12d0 + (137d0*x**6)/360d0 + (7d0*x**
     &        7)/20d0 + (363d0*x**8)/1120d0 + (761d0*x**9)/2520d0 + (
     &        7129d0*x**10)/25200d0
c     End of expansions around x = 0
      end select
c     --- set the imaginary part back to zero if it has been modified to
c     --- get the branch cuts right (and should be zero).
      if (bcflag.eq.1) then
         xre = dreal(x)
         if (n2.eq.0.and.xre.gt.0d0) then
            if (xre.lt.1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
c
         else if (n2.eq.1.and.xre.lt.1d0) then
            if (n1.ne.-1) then
               ris = dcmplx(dreal(ris),0d0)
            else if (xre.gt.-1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
c            
         else if (n2.eq.-1.and.xre.gt.-1d0) then
            if (n1.ne.1) then
               ris = dcmplx(dreal(ris),0d0)
            else if (xre.lt.1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
         endif
      endif        
      
      HPL2ar0=ris
      return
      end function
c-source-file HPL2ar1.f
c ------------------------------------------------
      double complex function HPL2ar1(n1,n2,x)
      implicit none
      integer n1,n2,j,bcflag
      double complex x,ris,myi,zp,llzp
      double precision pi, zeta2,ll2,xre

      pi=3.1415926535897932385D0
      zeta2=pi**2/6d0
      myi = dcmplx(0d0,1d0)
      ll2 = dlog(2d0)
      bcflag = 0

      j=1+(n2+1)+(n1+1)*3
      ris = dcmplx(0d0,0d0)

c---  +i*epsilon to get branch cuts right ---
      if (dimag(x).eq.0d0) then
         x = x + dcmplx(0d0,1d-60)
         bcflag = 1
      endif
c---  
      select case(j)
      case(1)                   !-1-1
         
         zp = 1d0-x
         
         ris = -((zp*ll2)/2d0) + (ll2**2)/2d0 + zp**
     &        5*(5d0/384d0 - (ll2)/160d0) + zp**3*(1d0/16d0 - (ll2)/2
     &        4d0) + zp**6*(137d0/23040d0 - (ll2)/384d0) + zp**4*(11d
     &        0/384d0 - (ll2)/64d0) + zp**2*(1d0/8d0 - (ll2)/8d0)
         
      case(2)                   !-10
         
         zp = 1d0-x
         
         ris = -((pi**2)/12d0) + (zp**2)/4d0 + (zp**
     &        3)/6d0 + (5d0*zp**4)/48d0 + (zp**5)/15d0 + (2d0*zp**6)/
     &        45d0
         
      case(3)                   !-11
         
         zp = 1d0-x
         llzp = log(zp)
         
         ris = (pi**2)/12d0 - (ll2**2)/2d0 + zp**5*(
     &        -(1d0/800d0) + (llzp)/160d0) + zp**3*(-(1d0/72d0) + (ll
     &        zp)/24d0) + zp*(-(1d0/2d0) + (llzp)/2d0) + zp**6*(-(1d0
     &        /2304d0) + (llzp)/384d0) + zp**4*(-(1d0/256d0) + (llzp)
     &        /64d0) + zp**2*(-(1d0/16d0) + (llzp)/8d0)
         
      case(4)                   !0-1
         
         zp = 1d0-x
         
         ris = (pi**2)/12d0 - zp*ll2 + zp**2*(1d0/4d
     &        0 - (ll2)/2d0) + zp**3*(5d0/24d0 - (ll2)/3d0) + zp**4*(
     &        1d0/6d0 - (ll2)/4d0) + zp**5*(131d0/960d0 - (ll2)/5d0) 
     &        + zp**6*(661d0/5760d0 - (ll2)/6d0)
         
      case(5)                   !00
         
         zp = 1d0-x
         
         ris = (zp**2)/2d0 + (zp**3)/2d0 + (11d0*zp*
     &        *4)/24d0 + (5d0*zp**5)/12d0 + (137d0*zp**6)/360d0
         
      case(6)                   !01
         
         zp = 1d0-x
         llzp = log(zp)
         
         ris = (pi**2)/6d0 + zp*(-1 + llzp) + zp**2*
     &        (-(1d0/4d0) + (llzp)/2d0) + zp**3*(-(1d0/9d0) + (llzp)/
     &        3d0) + zp**4*(-(1d0/16d0) + (llzp)/4d0) + zp**5*(-(1d0/
     &        25d0) + (llzp)/5d0) + zp**6*(-(1d0/36d0) + (llzp)/6d0)
         
      case(7)                   !1-1
         
         zp = 1d0-x
         llzp = log(zp)
         
         ris = -((pi**2)/12d0) + (zp)/2d0 + (zp**2)/
     &        16d0 + (zp**3)/72d0 + (zp**4)/256d0 + (zp**5)/800d0 + (
     &        zp**6)/2304d0 + (ll2**2)/2d0 - ll2*llzp
         
      case(8)                   !10
         
         zp = 1d0-x
         
         ris = -((pi**2)/6d0) + zp + (zp**2)/4d0 + (
     &        zp**3)/9d0 + (zp**4)/16d0 + (zp**5)/25d0 + (zp**6)/36d0
         
      case(9)                   !11
         
         zp = 1d0-x
         llzp = log(zp)
         
         ris = (llzp**2)/2d0
c     End of expansions around x = +1
      end select
c --- set the imaginary part back to zero if it has been modified to
c --- get the branch cuts right (and should be zero).
      if (bcflag.eq.1) then
         xre = dreal(x)
         if (n2.eq.0.and.xre.gt.0d0) then
            if (xre.lt.1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
c
         else if (n2.eq.1.and.xre.lt.1d0) then
            if (n1.ne.-1) then
               ris = dcmplx(dreal(ris),0d0)
            else if (xre.gt.-1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
c            
         else if (n2.eq.-1.and.xre.gt.-1d0) then
            if (n1.ne.1) then
               ris = dcmplx(dreal(ris),0d0)
            else if (xre.lt.1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
         endif
      endif        
      
      HPL2ar1=ris
      return
      end function
      
c-source-file HPL2arm1.f
c ------------------------------------------------
      double complex function HPL2arm1(n1,n2,x)
      implicit none
      integer n1,n2,j,bcflag,s,szp
      double complex x,ris,myi,zp,llzp
      double precision pi, zeta2,ll2,xre

      pi=3.1415926535897932385D0
      zeta2=pi**2/6d0
      myi = dcmplx(0d0,1d0)
      ll2 = dlog(2d0)
      bcflag = 0

      j=1+(n2+1)+(n1+1)*3
      ris = dcmplx(0d0,0d0)

c---  +i*epsilon to get branch cuts right ---
      if (dimag(x).eq.0d0) then
         x = x + dcmplx(0d0,1d-60)
         bcflag = 1
      endif
c---  
      select case(j)
      case(1)                   !-1-1

         zp = x+1d0
         llzp = log(zp)
         
         ris = (llzp**2)/2d0
         
      case(2)                   !-10
         
         zp = x+1d0
         llzp = log(zp)
         szp = s(zp)
         
         ris = (pi**2)/6d0 - zp - (zp**2)/4d0 - (zp*
     &        *3)/9d0 - (zp**4)/16d0 - (zp**5)/25d0 - (zp**6)/36d0 - 
     &        (zp**7)/49d0 - (zp**8)/64d0 - (zp**9)/81d0 + myi*pi*szp
     &        *llzp
         
      case(3)                   !-11
         
         zp = x+1d0
         llzp = log(zp)
         
         ris = -((pi**2)/12d0) + (zp)/2d0 + (zp**2)/
     &        16d0 + (zp**3)/72d0 + (zp**4)/256d0 + (zp**5)/800d0 + (
     &        zp**6)/2304d0 + (zp**7)/6272d0 + (zp**8)/16384d0 + (zp*
     &        *9)/41472d0 + (ll2**2)/2d0 - ll2*llzp
         
      case(4)                   !0-1
         
         zp = x+1d0
         llzp = log(zp)
         
         ris = -((pi**2)/6d0) + zp*(1 - llzp) + zp**
     &        2*(1d0/4d0 - (llzp)/2d0) + zp**3*(1d0/9d0 - (llzp)/3d0)
     &        + zp**4*(1d0/16d0 - (llzp)/4d0) + zp**5*(1d0/25d0 - (l
     &        lzp)/5d0) + zp**6*(1d0/36d0 - (llzp)/6d0) + zp**7*(1d0/
     &        49d0 - (llzp)/7d0) + zp**8*(1d0/64d0 - (llzp)/8d0) + zp
     &        **9*(1d0/81d0 - (llzp)/9d0)
         
      case(5)                   !00
         
         zp = x+1d0
         szp = s(zp)
         
         ris = -((pi**2)/2d0) - myi*pi*szp*zp + (1d0
     &        /2d0 - (myi*pi*szp)/2d0)*zp**2 + (1d0/2d0 - (myi*pi*szp
     &        )/3d0)*zp**3 + (11d0/24d0 - (myi*pi*szp)/4d0)*zp**4 + (
     &        5d0/12d0 - (myi*pi*szp)/5d0)*zp**5 + (137d0/360d0 - (my
     &        i*pi*szp)/6d0)*zp**6 + (7d0/20d0 - (myi*pi*szp)/7d0)*zp
     &        **7 + (363d0/1120d0 - (myi*pi*szp)/8d0)*zp**8 + (761d0/
     &        2520d0 - (myi*pi*szp)/9d0)*zp**9
         
      case(6)                   !01
         
         zp = x+1d0
         
         ris = -((pi**2)/12d0) + zp*ll2 + zp**2*(-(1
     &        d0/4d0) + (ll2)/2d0) + zp**3*(-(5d0/24d0) + (ll2)/3d0) 
     &        + zp**4*(-(1d0/6d0) + (ll2)/4d0) + zp**5*(-(131d0/960d0
     &        ) + (ll2)/5d0) + zp**6*(-(661d0/5760d0) + (ll2)/6d0) + 
     &        zp**7*(-(1327d0/13440d0) + (ll2)/7d0) + zp**8*(-(1163d0
     &        /13440d0) + (ll2)/8d0) + zp**9*(-(148969d0/1935360d0) +
     &        (ll2)/9d0)
         
      case(7)                   !1-1
         
         zp = x+1d0
         llzp = log(zp)
         
         ris = (pi**2)/12d0 - (ll2**2)/2d0 + zp**5*(
     &        -(1d0/800d0) + (llzp)/160d0) + zp**8*(-(1d0/16384d0) + 
     &        (llzp)/2048d0) + zp**3*(-(1d0/72d0) + (llzp)/24d0) + zp
     &        *(-(1d0/2d0) + (llzp)/2d0) + zp**6*(-(1d0/2304d0) + (ll
     &        zp)/384d0) + zp**9*(-(1d0/41472d0) + (llzp)/4608d0) + z
     &        p**4*(-(1d0/256d0) + (llzp)/64d0) + zp**7*(-(1d0/6272d0
     &        ) + (llzp)/896d0) + zp**2*(-(1d0/16d0) + (llzp)/8d0)
         
      case(8)                   !10
         
         zp = x+1d0
         szp = s(zp)
         
         ris = (pi**2)/12d0 + (myi*pi*szp*zp)/2d0 + 
     &        (-(1d0/4d0) + (myi*pi*szp)/8d0)*zp**2 + (-(1d0/6d0) + (
     &        myi*pi*szp)/24d0)*zp**3 + (-(5d0/48d0) + (myi*pi*szp)/6
     &        4d0)*zp**4 + (-(1d0/15d0) + (myi*pi*szp)/160d0)*zp**5 +
     &        (-(2d0/45d0) + (myi*pi*szp)/384d0)*zp**6 + (-(13d0/420
     &        d0) + (myi*pi*szp)/896d0)*zp**7 + (-(151d0/6720d0) + (m
     &        yi*pi*szp)/2048d0)*zp**8 + (-(16d0/945d0) + (myi*pi*szp
     &        )/4608d0)*zp**9 - myi*pi*szp*ll2
         
      case(9)                   !11
         
         zp = x+1d0
         
         ris = -((zp*ll2)/2d0) + (ll2**2)/2d0 + zp**
     &        5*(5d0/384d0 - (ll2)/160d0) + zp**8*(363d0/286720d0 - (
     &        ll2)/2048d0) + zp**3*(1d0/16d0 - (ll2)/24d0) + zp**6*(1
     &        37d0/23040d0 - (ll2)/384d0) + zp**9*(761d0/1290240d0 - 
     &        (ll2)/4608d0) + zp**4*(11d0/384d0 - (ll2)/64d0) + zp**7
     &        *(7d0/2560d0 - (ll2)/896d0) + zp**2*(1d0/8d0 - (ll2)/8d
     &        0)
c     End of expansions around x = -1
      end select
c     --- set the imaginary part back to zero if it has been modified to
c     --- get the branch cuts right (and should be zero).
      if (bcflag.eq.1) then
         xre = dreal(x)
         if (n2.eq.0.and.xre.gt.0d0) then
            if (xre.lt.1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
c     
         else if (n2.eq.1.and.xre.lt.1d0) then
            if (n1.ne.-1) then
               ris = dcmplx(dreal(ris),0d0)
            else if (xre.gt.-1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
c     
         else if (n2.eq.-1.and.xre.gt.-1d0) then
            if (n1.ne.1) then
               ris = dcmplx(dreal(ris),0d0)
            else if (xre.lt.1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
         endif
      endif        
      
      HPL2arm1=ris
      return
      end function
c-source-file HPL2at1.f
      double complex function HPL2at1(n1,n2)
      implicit none
      integer n1,n2,j
      double complex ris,myi
      double precision pi,ll2

      pi=3.1415926535897932385D0
      myi = dcmplx(0d0,1d0)
      ll2 = dlog(2d0)

      j=1+(n2+1)+(n1+1)*3
      ris = dcmplx(0d0,0d0)

      if ((j.eq.7).or.(j.eq.9)) then
         print*, ""
         print*, "****************"
         print*, "ERROR in HPL2: "
         print*, "HPL2(",n1,",",n2
     &        ,",1) is divergent!"
         print*, "Aborting..."
         print*,"****************"
         stop
      else
         select case (j)
         case(1)
            ris = ll2**2/2d0
         case(2)
            ris = -pi**2/12d0
         case(3)
            ris = pi**2/12d0 - ll2**2/2d0
         case(4)
            ris = pi**2/12d0
         case(5)
            ris = 0d0
         case(6)
            ris = pi**2/6d0
         case(8)
            ris = -pi**2/6d0
         end select
      endif
      HPL2at1=ris
      return
      end function
c-source-file HPL2atm1.f
      double complex function HPL2atm1(n1,n2)
      implicit none
      integer n1,n2,j
      double complex ris,myi
      double precision pi,ll2

      pi=3.1415926535897932385D0
      myi = dcmplx(0d0,1d0)
      ll2 = dlog(2d0)

      j=1+(n2+1)+(n1+1)*3
      ris = dcmplx(0d0,0d0)

      if(j.gt.3) then
         select case (j)
         case(4)
            ris = -pi**2/6d0
         case(5)
            ris = -pi**2/2d0
         case(6)
            ris = -pi**2/12d0
         case(7)
            ris = pi**2/12d0 - ll2**2/2d0
         case(8)
            ris = pi**2/12d0 - myi*pi*ll2
         case(9)
            ris = ll2**2/2d0
         end select
      else
         print*, ""
         print*, "****************"
         print*, "ERROR in HPL2: "
         print*, "HPL2(",n1,",",n2
     &        ,",-1) is divergent!"
         print*, "Aborting..."
         print*,"****************"
         stop
      endif
      HPL2atm1=ris
      return
      end function
c-source-file HPL2else.f
C=============================================================================
C---  HPLs  of Rank 2  
C=============================================================================

      double  complex function HPL2else(n1,n2,x)
      implicit none 
      double precision pi,ll2,xre
      double complex x, ris,myi,ll1x,ll1mx,llx
      double complex cli2
      integer n1,n2,j,bcflag

      pi=3.1415926535897932385D0
      myi=dcmplx(0d0,1d0)

      ll2 = dlog(2d0)
      j = 3*(n1+1) + (n2+1) +1

      ris=dcmplx(0d0,0d0)
      bcflag = 0

c---  +i*epsilon to get branch cuts right ---
      if (dimag(x).eq.0d0) then
         x = x + dcmplx(0d0,1d-60)
         bcflag = 1
      endif
c---  
      ll1x = log(1d0+x)
      ll1mx = log(1d0-x)
      llx = log(x)
      
      select case(j)
      case(1)
         ris=ll1x**2/2d0
      case(2)
         ris=cli2(-x) + llx*ll1x
      case(3)
         ris= pi**2/12d0 - ll2**2/2d0 + ll2*ll1mx 
     &        - ll1mx*ll1x - cli2((1d0-x)/2d0)
      case(4)
         ris=-cli2(-x)
      case(5)
         ris=llx**2/2d0
      case(6)
         ris=cli2(x)
      case(7)
         ris=-pi**2/12d0 + cli2((1d0-x)/2d0) + ll2**2/2d0 
     &        -ll2*ll1mx
      case(8)
         ris=-cli2(x)-ll1mx*llx
      case(9)
         ris=ll1mx**2/2d0
      end select

c --- set the imaginary part back to zero if it has been modified to
c --- get the branch cuts right (and should be zero). Also, set imaginary
c --- part of result to zero if x is between 0 and 1.

      if (bcflag.eq.1) then
         x = x - dcmplx(0d0,1d-60)
         xre = dreal(x)
         if (xre.ge.0d0.and.xre.le.1d0) then
            ris = dcmplx(dreal(ris),0d0)
         endif
      endif

      HPL2else=ris 
      return
      end
c-source-file HPL2.f
C=============================================================================
C---  HPLs  of Rank 2  
C=============================================================================
c --- main forking function
      double  complex function HPL2(n1,n2, x)
      implicit none 
      integer n1,n2
      double complex x,ris
      double complex HPL2at0,HPL2at1,HPL2atm1
      double complex HPL2ar1,HPL2arm1,HPL2ar0
      double complex HPL2else
      double precision rad1,radm1,rad0
      
      rad0 = 0.025d0
      rad1 = 0.01d0
      radm1 = 0.025d0      

      if ((abs(n1).gt.1).or.(abs(n2).gt.1)) then
         print*, ""
         print*, "****************"
         print*, "Error in HPL2:"
         print*, "Indices",n1,n2," out of range !"
         print*, "Aborting..."
         print*,"****************"
         stop
      endif

      ris = dcmplx(0d0,0d0)
      
      if (x.eq.dcmplx(0d0,0d0)) then
         ris = HPL2at0(n1,n2)
      elseif (x.eq.dcmplx(1d0,0d0)) then
         ris = HPL2at1(n1,n2)
      elseif (x.eq.dcmplx(-1d0,0d0)) then
         ris = HPL2atm1(n1,n2)
      elseif (abs(x-dcmplx(1d0,0d0)).lt.rad1) then
         ris = HPL2ar1(n1,n2,x)
      elseif (abs(x+dcmplx(1d0,0d0)).lt.radm1) then
         ris = HPL2arm1(n1,n2,x)
      elseif (abs(x-dcmplx(0d0,0d0)).lt.rad0) then
         ris = HPL2ar0(n1,n2,x)
      else 
         ris = HPL2else(n1,n2,x)
      endif

      HPL2=ris 
      return
      end
c ------------------------------------------------
      double complex function HPL2at0(n1, n2)
      implicit none
      integer n1,n2,j
      double complex ris
    
      j=1+(n2+1)+(n1+1)*3
      ris = dcmplx(0d0,0d0)

      if (j.eq.5) then             
         print*, ""
         print*, "****************"
         print*, "ERROR in HPL2: "
         print*, "HPL2(",n1,",",n2
     &        ,",0) is divergent!"
         print*, "Aborting..."
         print*,"****************"
         stop
      endif
      HPL2at0=ris
      return
      end function

c --- Real part of HPL2     
      double precision function HPL2real(n1,n2,xr,xi)
      implicit none
      double precision xr,xi
      integer n1,n2
      double complex x,HPL2
      x=dcmplx(xr,xi)
      HPL2real = dreal(HPL2(n1,n2,x))
      return
      end

c --- Imaginary part of HPL2     
      double precision function HPL2im(n1,n2,xr,xi)
      implicit none
      double precision xr,xi
      integer n1,n2
      double complex x,HPL2
      x=dcmplx(xr,xi)
      HPL2im = dimag(HPL2(n1,n2,x))
      return
      end
c-source-file HPL3ar0.f
c ------------------------------------------------
      double complex function HPL3ar0(n1,n2,n3,x)
      implicit none
      integer n1,n2,n3,j,bcflag
      double complex x,ris,myi,llx
      double precision pi, zeta2, zeta3,ll2,xre

      pi=3.1415926535897932385D0
      zeta3=1.20205690315959428539973816151d0
      zeta2=pi**2/6d0
      myi = dcmplx(0d0,1d0)
      ll2 = dlog(2d0)
      bcflag = 0

      j=1+(n3+1)*1+(n2+1)*3+(n1+1)*9
      ris = dcmplx(0d0,0d0)

c---  +i*epsilon to get branch cuts right ---
      if (dimag(x).eq.0d0) then
         x = x + dcmplx(0d0,1d-60)
         bcflag = 1
      endif
c---  

      select case(j)
      case(1)                   !-1-1-1
         
         
         ris = (x**3)/6d0 - (x**4)/4d0 + (7d0*x**5)/
     &        24d0 - (5d0*x**6)/16d0 + (29d0*x**7)/90d0 - (469d0*x**8
     &        )/1440d0 + (29531d0*x**9)/90720d0 - (1303d0*x**10)/4032
     &        d0
         
      case(2)                   !-1-10
         
         llx = log(x)
         
         ris = -((3d0*x**2)/4d0) + (7d0*x**3)/12d0 -
     &        (131d0*x**4)/288d0 + (53d0*x**5)/144d0 - (2213d0*x**6)
     &        /7200d0 + (947d0*x**7)/3600d0 - (647707d0*x**8)/2822400
     &        d0 + (1290829d0*x**9)/6350400d0 - (11574649d0*x**10)/63
     &        504000d0 + (x**2*llx)/2d0 - (x**3*llx)/2d0 + (11d0*x**4
     &        *llx)/24d0 - (5d0*x**5*llx)/12d0 + (137d0*x**6*llx)/360
     &        d0 - (7d0*x**7*llx)/20d0 + (363d0*x**8*llx)/1120d0 - (7
     &        61d0*x**9*llx)/2520d0 + (7129d0*x**10*llx)/25200d0
         
      case(3)                   !-1-11
         
         
         ris = (x**3)/6d0 - (x**4)/6d0 + (7d0*x**5)/
     &        40d0 - (119d0*x**6)/720d0 + (101d0*x**7)/630d0 - (305d0
     &        *x**8)/2016d0 + (13157d0*x**9)/90720d0 - (41603d0*x**10
     &        )/302400d0
         
      case(4)                   !-10-1
         
         
         ris = (x**2)/2d0 - (5d0*x**3)/12d0 + (49d0*
     &        x**4)/144d0 - (41d0*x**5)/144d0 + (5269d0*x**6)/21600d0
     &        - (767d0*x**7)/3600d0 + (266681d0*x**8)/1411200d0 - (1
     &        077749d0*x**9)/6350400d0 + (9778141d0*x**10)/63504000d0
         
      case(5)                   !-100
         
         llx = log(x)
         
         ris = x - (x**2)/8d0 + (x**3)/27d0 - (x**4)
     &        /64d0 + (x**5)/125d0 - (x**6)/216d0 + (x**7)/343d0 - (x
     &        **8)/512d0 + (x**9)/729d0 - (x**10)/1000d0 - x*llx + (x
     &        **2*llx)/4d0 - (x**3*llx)/9d0 + (x**4*llx)/16d0 - (x**5
     &        *llx)/25d0 + (x**6*llx)/36d0 - (x**7*llx)/49d0 + (x**8*
     &        llx)/64d0 - (x**9*llx)/81d0 + (x**10*llx)/100d0 + (x*ll
     &        x**2)/2d0 - (x**2*llx**2)/4d0 + (x**3*llx**2)/6d0 - (x*
     &        *4*llx**2)/8d0 + (x**5*llx**2)/10d0 - (x**6*llx**2)/12d
     &        0 + (x**7*llx**2)/14d0 - (x**8*llx**2)/16d0 + (x**9*llx
     &        **2)/18d0 - (x**10*llx**2)/20d0
         
      case(6)                   !-101
         
         
         ris = (x**2)/2d0 - (x**3)/4d0 + (31d0*x**4)
     &        /144d0 - (23d0*x**5)/144d0 + (3019d0*x**6)/21600d0 - (1
     &        39d0*x**7)/1200d0 + (48877d0*x**8)/470400d0 - (191833d0
     &        *x**9)/2116800d0 + (5257891d0*x**10)/63504000d0
         
      case(7)                   !-11-1
         
         
         ris = (x**3)/6d0 - (x**4)/12d0 + (13d0*x**5
     &        )/120d0 - (17d0*x**6)/240d0 + (5d0*x**7)/63d0 - (589d0*
     &        x**8)/10080d0 + (5669d0*x**9)/90720d0 - (85d0*x**10)/17
     &        28d0
         
      case(8)                   !-110
         
         llx = log(x)
         
         ris = -((3d0*x**2)/4d0) + (11d0*x**3)/36d0 
     &        - (77d0*x**4)/288d0 + (659d0*x**5)/3600d0 - (1163d0*x**
     &        6)/7200d0 + (2517d0*x**7)/19600d0 - (108919d0*x**8)/940
     &        800d0 + (1875737d0*x**9)/19051200d0 - (5731399d0*x**10)
     &        /63504000d0 + (x**2*llx)/2d0 - (x**3*llx)/6d0 + (5d0*x*
     &        *4*llx)/24d0 - (7d0*x**5*llx)/60d0 + (47d0*x**6*llx)/36
     &        0d0 - (37d0*x**7*llx)/420d0 + (319d0*x**8*llx)/3360d0 -
     &        (533d0*x**9*llx)/7560d0 + (1879d0*x**10*llx)/25200d0
         
      case(9)                   !-111
         
         
         ris = (x**3)/6d0 + (11d0*x**5)/120d0 - (x**
     &        6)/144d0 + (19d0*x**7)/315d0 - (13d0*x**8)/1440d0 + (79
     &        9d0*x**9)/18144d0 - (317d0*x**10)/33600d0
         
      case(10)                  !0-1-1
         
         
         ris = (x**2)/4d0 - (x**3)/6d0 + (11d0*x**4)
     &        /96d0 - (x**5)/12d0 + (137d0*x**6)/2160d0 - (x**7)/20d0
     &        + (363d0*x**8)/8960d0 - (761d0*x**9)/22680d0 + (7129d0
     &        *x**10)/252000d0
         
      case(11)                  !0-10
         
         llx = log(x)
         
         ris = -2*x + (x**2)/4d0 - (2d0*x**3)/27d0 +
     & (x**4)/32d0 - (2d0*x**5)/125d0 + (x**6)/108d0 - (2d0*x
     &        **7)/343d0 + (x**8)/256d0 - (2d0*x**9)/729d0 + (x**10)/
     &        500d0 + x*llx - (x**2*llx)/4d0 + (x**3*llx)/9d0 - (x**4
     &        *llx)/16d0 + (x**5*llx)/25d0 - (x**6*llx)/36d0 + (x**7*
     &        llx)/49d0 - (x**8*llx)/64d0 + (x**9*llx)/81d0 - (x**10*
     &        llx)/100d0
         
      case(12)                  !0-11
         
         
         ris = (x**2)/4d0 - (x**3)/18d0 + (5d0*x**4)
     &        /96d0 - (7d0*x**5)/300d0 + (47d0*x**6)/2160d0 - (37d0*x
     &        **7)/2940d0 + (319d0*x**8)/26880d0 - (533d0*x**9)/68040
     &        d0 + (1879d0*x**10)/252000d0
         
      case(13)                  !00-1
         
         
         ris = x - (x**2)/8d0 + (x**3)/27d0 - (x**4)
     &        /64d0 + (x**5)/125d0 - (x**6)/216d0 + (x**7)/343d0 - (x
     &        **8)/512d0 + (x**9)/729d0 - (x**10)/1000d0
         
      case(14)                  !000
         
         llx = log(x)
         
         ris = (llx**3)/6d0
         
      case(15)                  !001
         
         
         ris = x + (x**2)/8d0 + (x**3)/27d0 + (x**4)
     &        /64d0 + (x**5)/125d0 + (x**6)/216d0 + (x**7)/343d0 + (x
     &        **8)/512d0 + (x**9)/729d0 + (x**10)/1000d0
         
      case(16)                  !01-1
         
         
         ris = (x**2)/4d0 + (x**3)/18d0 + (5d0*x**4)
     &        /96d0 + (7d0*x**5)/300d0 + (47d0*x**6)/2160d0 + (37d0*x
     &        **7)/2940d0 + (319d0*x**8)/26880d0 + (533d0*x**9)/68040
     &        d0 + (1879d0*x**10)/252000d0
         
      case(17)                  !010
         
         llx = log(x)
         
         ris = -2*x - (x**2)/4d0 - (2d0*x**3)/27d0 -
     &        (x**4)/32d0 - (2d0*x**5)/125d0 - (x**6)/108d0 - (2d0*x
     &        **7)/343d0 - (x**8)/256d0 - (2d0*x**9)/729d0 - (x**10)/
     &        500d0 + x*llx + (x**2*llx)/4d0 + (x**3*llx)/9d0 + (x**4
     &        *llx)/16d0 + (x**5*llx)/25d0 + (x**6*llx)/36d0 + (x**7*
     &        llx)/49d0 + (x**8*llx)/64d0 + (x**9*llx)/81d0 + (x**10*
     &        llx)/100d0
         
      case(18)                  !011
         
         
         ris = (x**2)/4d0 + (x**3)/6d0 + (11d0*x**4)
     &        /96d0 + (x**5)/12d0 + (137d0*x**6)/2160d0 + (x**7)/20d0
     &        + (363d0*x**8)/8960d0 + (761d0*x**9)/22680d0 + (7129d0
     &        *x**10)/252000d0
         
      case(19)                  !1-1-1
         
         
         ris = (x**3)/6d0 + (11d0*x**5)/120d0 + (x**
     &        6)/144d0 + (19d0*x**7)/315d0 + (13d0*x**8)/1440d0 + (79
     &        9d0*x**9)/18144d0 + (317d0*x**10)/33600d0
         
      case(20)                  !1-10

         llx = log(x)
         
         ris = -((3d0*x**2)/4d0) - (11d0*x**3)/36d0 
     &        - (77d0*x**4)/288d0 - (659d0*x**5)/3600d0 - (1163d0*x**
     &        6)/7200d0 - (2517d0*x**7)/19600d0 - (108919d0*x**8)/940
     &        800d0 - (1875737d0*x**9)/19051200d0 - (5731399d0*x**10)
     &        /63504000d0 + (x**2*llx)/2d0 + (x**3*llx)/6d0 + (5d0*x*
     &        *4*llx)/24d0 + (7d0*x**5*llx)/60d0 + (47d0*x**6*llx)/36
     &        0d0 + (37d0*x**7*llx)/420d0 + (319d0*x**8*llx)/3360d0 +
     &        (533d0*x**9*llx)/7560d0 + (1879d0*x**10*llx)/25200d0
         
      case(21)                  !1-11
         
         
         ris = (x**3)/6d0 + (x**4)/12d0 + (13d0*x**5
     &        )/120d0 + (17d0*x**6)/240d0 + (5d0*x**7)/63d0 + (589d0*
     &        x**8)/10080d0 + (5669d0*x**9)/90720d0 + (85d0*x**10)/17
     &        28d0
         
      case(22)                  !10-1
         
         
         ris = (x**2)/2d0 + (x**3)/4d0 + (31d0*x**4)
     &        /144d0 + (23d0*x**5)/144d0 + (3019d0*x**6)/21600d0 + (1
     &        39d0*x**7)/1200d0 + (48877d0*x**8)/470400d0 + (191833d0
     &        *x**9)/2116800d0 + (5257891d0*x**10)/63504000d0
         
      case(23)                  !100
         
         llx = log(x)
         
         ris = x + (x**2)/8d0 + (x**3)/27d0 + (x**4)
     &        /64d0 + (x**5)/125d0 + (x**6)/216d0 + (x**7)/343d0 + (x
     &        **8)/512d0 + (x**9)/729d0 + (x**10)/1000d0 - x*llx - (x
     &        **2*llx)/4d0 - (x**3*llx)/9d0 - (x**4*llx)/16d0 - (x**5
     &        *llx)/25d0 - (x**6*llx)/36d0 - (x**7*llx)/49d0 - (x**8*
     &        llx)/64d0 - (x**9*llx)/81d0 - (x**10*llx)/100d0 + (x*ll
     &        x**2)/2d0 + (x**2*llx**2)/4d0 + (x**3*llx**2)/6d0 + (x*
     &        *4*llx**2)/8d0 + (x**5*llx**2)/10d0 + (x**6*llx**2)/12d
     &        0 + (x**7*llx**2)/14d0 + (x**8*llx**2)/16d0 + (x**9*llx
     &        **2)/18d0 + (x**10*llx**2)/20d0
         
      case(24)                  !101
         
         
         ris = (x**2)/2d0 + (5d0*x**3)/12d0 + (49d0*
     &        x**4)/144d0 + (41d0*x**5)/144d0 + (5269d0*x**6)/21600d0
     &        + (767d0*x**7)/3600d0 + (266681d0*x**8)/1411200d0 + (1
     &        077749d0*x**9)/6350400d0 + (9778141d0*x**10)/63504000d0
         
      case(25)                  !11-1
         
         
         ris = (x**3)/6d0 + (x**4)/6d0 + (7d0*x**5)/
     &        40d0 + (119d0*x**6)/720d0 + (101d0*x**7)/630d0 + (305d0
     &        *x**8)/2016d0 + (13157d0*x**9)/90720d0 + (41603d0*x**10
     &        )/302400d0
         
      case(26)                  !110
         
         llx = log(x)
         
         ris = -((3d0*x**2)/4d0) - (7d0*x**3)/12d0 -
     &        (131d0*x**4)/288d0 - (53d0*x**5)/144d0 - (2213d0*x**6)
     &        /7200d0 - (947d0*x**7)/3600d0 - (647707d0*x**8)/2822400
     &        d0 - (1290829d0*x**9)/6350400d0 - (11574649d0*x**10)/63
     &        504000d0 + (x**2*llx)/2d0 + (x**3*llx)/2d0 + (11d0*x**4
     &        *llx)/24d0 + (5d0*x**5*llx)/12d0 + (137d0*x**6*llx)/360
     &        d0 + (7d0*x**7*llx)/20d0 + (363d0*x**8*llx)/1120d0 + (7
     &        61d0*x**9*llx)/2520d0 + (7129d0*x**10*llx)/25200d0
         
      case(27)                  !111
         
         
         ris = (x**3)/6d0 + (x**4)/4d0 + (7d0*x**5)/
     &        24d0 + (5d0*x**6)/16d0 + (29d0*x**7)/90d0 + (469d0*x**8
     &        )/1440d0 + (29531d0*x**9)/90720d0 + (1303d0*x**10)/4032
     &        d0
c     End of expansions around x = 0
      end select
c     --- set the imaginary part back to zero if it has been modified to
c     --- get the branch cuts right (and should be zero).
      if (bcflag.eq.1) then
         xre = dreal(x)
         if (n3.eq.0.and.xre.gt.0d0) then
            if (xre.lt.1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
c
         else if (n3.eq.1.and.xre.lt.1d0) then
            if (n1.ne.-1.and.n2.ne.-1) then
               ris = dcmplx(dreal(ris),0d0)
            else if (xre.gt.-1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
c            
         else if (n3.eq.-1.and.xre.gt.-1d0) then
            if (n1.ne.1.and.n2.ne.1) then
               ris = dcmplx(dreal(ris),0d0)
            else if (xre.lt.1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
         endif
      endif    

      HPL3ar0=ris
      return
      end function
c-source-file HPL3ar1.f
c ------------------------------------------------
      double complex function HPL3ar1(n1,n2,n3,x)
      implicit none
      integer n1,n2,n3,j,bcflag
      double complex x,ris,myi,zp,llzp
      double precision pi, zeta2, zeta3,ll2,xre

      pi=3.1415926535897932385D0
      zeta3=1.20205690315959428539973816151d0
      zeta2=pi**2/6d0
      myi = dcmplx(0d0,1d0)
      ll2 = dlog(2d0)
      bcflag = 0

      j=1+(n3+1)*1+(n2+1)*3+(n1+1)*9
      ris = dcmplx(0d0,0d0)

c---  +i*epsilon to get branch cuts right ---
      if (dimag(x).eq.0d0) then
         x = x + dcmplx(0d0,1d-60)
         bcflag = 1
      endif
c---  

c This was file contains the Taylor 
c expansions around x = +1
c The expansion parameter is zp = 1-x
      select case(j)
      
         case(1)            !-1-1-1

            zp = 1d0-x

            ris = -((zp*ll2**2)/4d0) + (ll2**3)/6d0 + z
     &p**4*(-(1d0/64d0) + (11d0*ll2)/384d0 - (ll2**2)/128d0) 
     &+ zp**2*((ll2)/8d0 - (ll2**2)/16d0) + zp**5*(-(7d0/768d
     &0) + (5d0*ll2)/384d0 - (ll2**2)/320d0) + zp**3*(-(1d0/4
     &8d0) + (ll2)/16d0 - (ll2**2)/48d0) + zp**6*(-(5d0/1024d
     &0) + (137d0*ll2)/23040d0 - (ll2**2)/768d0)

         case(2)            !-1-10

            zp = 1d0-x

            ris = (pi**2*zp)/24d0 + (pi**2*zp**2)/96d0 
     &+ (-(1d0/24d0) + (pi**2)/288d0)*zp**3 + ((pi**2)/768d0 
     &- 7d0/192d0)*zp**4 + ((pi**2)/1920d0 - 1d0/40d0)*zp**5 
     &+ (-(23d0/1440d0) + (pi**2)/4608d0)*zp**6 - (pi**2*ll2)
     &/12d0 + (zeta3)/8d0

         case(3)            !-1-11

            zp = 1d0-x
            llzp = log(zp)

            ris = -((ll2**3)/6d0) + zp*(-((pi**2)/24d0)
     & + (ll2**2)/4d0) + zp**3*(-((pi**2)/288d0) + 7d0/96d0 +
     & (ll2**2)/48d0 - (llzp)/16d0) + zp**6*(2213d0/460800d0 
     &- (pi**2)/4608d0 + (ll2**2)/768d0 - (137d0*llzp)/23040d
     &0) + zp**4*(131d0/4608d0 - (pi**2)/768d0 + (ll2**2)/128
     &d0 - (11d0*llzp)/384d0) + zp**5*(-((pi**2)/1920d0) + 53
     &d0/4608d0 + (ll2**2)/320d0 - (5d0*llzp)/384d0) + zp**2*
     &(3d0/16d0 - (pi**2)/96d0 + (ll2**2)/16d0 - (llzp)/8d0) 
     &+ (zeta3)/8d0

         case(4)            !-10-1

            zp = 1d0-x

            ris = -((pi**2*zp)/24d0) + (pi**2*ll2)/12d0
     & + zp**5*(-((pi**2)/1920d0) - 1d0/30d0 + (ll2)/15d0) + 
     &zp**6*(-((pi**2)/4608d0) - 97d0/3840d0 + (2d0*ll2)/45d0
     &) + zp**2*(-((pi**2)/96d0) + (ll2)/4d0) + zp**4*(-(1d0/
     &24d0) - (pi**2)/768d0 + (5d0*ll2)/48d0) + zp**3*(-(1d0/
     &24d0) - (pi**2)/288d0 + (ll2)/6d0) - (zeta3)/4d0

         case(5)            !-100

            zp = 1d0-x

            ris = -((zp**3)/12d0) - (3d0*zp**4)/32d0 - 
     &(zp**5)/12d0 - (5d0*zp**6)/72d0 + (3d0*zeta3)/4d0

         case(6)            !-101

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**2*zp)/12d0) + (pi**2*ll2)/6d0 
     &+ zp**5*(79d0/1800d0 - (pi**2)/960d0 - (llzp)/15d0) + z
     &p**6*(-((pi**2)/2304d0) + 169d0/7200d0 - (2d0*llzp)/45d
     &0) + zp**2*(-((pi**2)/48d0) + 3d0/8d0 - (llzp)/4d0) + z
     &p**4*(25d0/288d0 - (pi**2)/384d0 - (5d0*llzp)/48d0) + z
     &p**3*(-((pi**2)/144d0) + 13d0/72d0 - (llzp)/6d0) - (5d0
     &*zeta3)/8d0

         case(7)            !-11-1

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**2*ll2)/12d0 - (ll2**3)/6d0 + zp*
     &*5*((pi**2)/1920d0 - 41d0/4608d0 - (ll2)/800d0 - (ll2**
     &2)/320d0 + (ll2*llzp)/160d0) + zp**3*((pi**2)/288d0 - 5
     &d0/96d0 - (ll2)/72d0 - (ll2**2)/48d0 + (ll2*llzp)/24d0)
     & + zp*((pi**2)/24d0 - (ll2)/2d0 - (ll2**2)/4d0 + (ll2*l
     &lzp)/2d0) + zp**6*((pi**2)/4608d0 - 5269d0/1382400d0 - 
     &(ll2)/2304d0 - (ll2**2)/768d0 + (ll2*llzp)/384d0) + zp*
     &*4*(-(49d0/2304d0) + (pi**2)/768d0 - (ll2)/256d0 - (ll2
     &**2)/128d0 + (ll2*llzp)/64d0) + zp**2*(-(1d0/8d0) + (pi
     &**2)/96d0 - (ll2)/16d0 - (ll2**2)/16d0 + (ll2*llzp)/8d0
     &) - (zeta3)/4d0

         case(8)            !-110

            zp = 1d0-x

            ris = (pi**2*zp)/12d0 + ((pi**2)/48d0 - 1d0
     &/4d0)*zp**2 + ((pi**2)/144d0 - 1d0/8d0)*zp**3 + ((pi**2
     &)/384d0 - 35d0/576d0)*zp**4 + (-(11d0/360d0) + (pi**2)/
     &960d0)*zp**5 + ((pi**2)/2304d0 - 347d0/21600d0)*zp**6 +
     & (pi**2*ll2)/12d0 - zeta3

         case(9)            !-111

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**2*ll2)/12d0) + (ll2**3)/6d0 + 
     &zp**4*(-(1d0/1024d0) + (llzp)/256d0 - (llzp**2)/128d0) 
     &+ zp**2*(-(1d0/32d0) + (llzp)/16d0 - (llzp**2)/16d0) + 
     &zp**5*(-(1d0/4000d0) + (llzp)/800d0 - (llzp**2)/320d0) 
     &+ zp**3*(-(1d0/216d0) + (llzp)/72d0 - (llzp**2)/48d0) +
     & zp*(-(1d0/2d0) + (llzp)/2d0 - (llzp**2)/4d0) + zp**6*(
     &-(1d0/13824d0) + (llzp)/2304d0 - (llzp**2)/768d0) + (7d
     &0*zeta3)/8d0

         case(10)            !0-1-1

            zp = 1d0-x

            ris = -((zp*ll2**2)/2d0) + zp**5*(-(83d0/19
     &20d0) + (131d0*ll2)/960d0 - (ll2**2)/10d0) + zp**6*(-(1
     &1d0/288d0) + (661d0*ll2)/5760d0 - (ll2**2)/12d0) + zp**
     &2*((ll2)/4d0 - (ll2**2)/4d0) + zp**3*(-(1d0/24d0) + (5d
     &0*ll2)/24d0 - (ll2**2)/6d0) + zp**4*(-(3d0/64d0) + (ll2
     &)/6d0 - (ll2**2)/8d0) + (zeta3)/8d0

         case(11)            !0-10

            zp = 1d0-x

            ris = (pi**2*zp)/12d0 + (pi**2*zp**2)/24d0 
     &+ (-(1d0/12d0) + (pi**2)/36d0)*zp**3 + ((pi**2)/48d0 - 
     &5d0/48d0)*zp**4 + (-(5d0/48d0) + (pi**2)/60d0)*zp**5 + 
     &(-(47d0/480d0) + (pi**2)/72d0)*zp**6 - (3d0*zeta3)/2d0

         case(12)            !0-11

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**2*ll2)/4d0) + zp*(-((pi**2)/12
     &d0) + (ll2**2)/2d0) + zp**2*(-((pi**2)/24d0) + 3d0/8d0 
     &+ (ll2**2)/4d0 - (llzp)/4d0) + zp**3*(-((pi**2)/36d0) +
     & 37d0/144d0 + (ll2**2)/6d0 - (5d0*llzp)/24d0) + zp**6*(
     &13369d0/115200d0 - (pi**2)/72d0 + (ll2**2)/12d0 - (661d
     &0*llzp)/5760d0) + zp**4*(-((pi**2)/48d0) + 107d0/576d0 
     &+ (ll2**2)/8d0 - (llzp)/6d0) + zp**5*(-((pi**2)/60d0) +
     & 8257d0/57600d0 + (ll2**2)/10d0 - (131d0*llzp)/960d0) +
     & (13d0*zeta3)/8d0

         case(13)            !00-1

            zp = 1d0-x

            ris = -((pi**2*zp)/12d0) + zp**4*(-((pi**2)
     &/48d0) - 11d0/96d0 + (11d0*ll2)/24d0) + zp**2*(-((pi**2
     &)/24d0) + (ll2)/2d0) + zp**3*(-(1d0/12d0) - (pi**2)/36d
     &0 + (ll2)/2d0) + zp**6*(-((pi**2)/72d0) - 731d0/5760d0 
     &+ (137d0*ll2)/360d0) + zp**5*(-((pi**2)/60d0) - 1d0/8d0
     & + (5d0*ll2)/12d0) + (3d0*zeta3)/4d0

         case(14)            !000

            zp = 1d0-x

            ris = -((zp**3)/6d0) - (zp**4)/4d0 - (7d0*z
     &p**5)/24d0 - (5d0*zp**6)/16d0

         case(15)            !001

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**2*zp)/6d0) + zp**4*(-((pi**2)/
     &24d0) + 131d0/288d0 - (11d0*llzp)/24d0) + zp**2*(-((pi*
     &*2)/12d0) + 3d0/4d0 - (llzp)/2d0) + zp**3*(-((pi**2)/18
     &d0) + 7d0/12d0 - (llzp)/2d0) + zp**6*(-((pi**2)/36d0) +
     & 2213d0/7200d0 - (137d0*llzp)/360d0) + zp**5*(-((pi**2)
     &/30d0) + 53d0/144d0 - (5d0*llzp)/12d0) + zeta3

         case(16)            !01-1

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**2*ll2)/4d0 + zp*((pi**2)/12d0 - 
     &ll2 - (ll2**2)/2d0 + ll2*llzp) + zp**2*((pi**2)/24d0 - 
     &1d0/4d0 - (ll2)/4d0 - (ll2**2)/4d0 + (ll2*llzp)/2d0) + 
     &zp**3*((pi**2)/36d0 - 3d0/16d0 - (ll2)/9d0 - (ll2**2)/6
     &d0 + (ll2*llzp)/3d0) + zp**4*((pi**2)/48d0 - 83d0/576d0
     & - (ll2)/16d0 - (ll2**2)/8d0 + (ll2*llzp)/4d0) + zp**5*
     &(-(1337d0/11520d0) + (pi**2)/60d0 - (ll2)/25d0 - (ll2**
     &2)/10d0 + (ll2*llzp)/5d0) + zp**6*(-(33497d0/345600d0) 
     &+ (pi**2)/72d0 - (ll2)/36d0 - (ll2**2)/12d0 + (ll2*llzp
     &)/6d0) - zeta3

         case(17)            !010

            zp = 1d0-x

            ris = (pi**2*zp)/6d0 + ((pi**2)/12d0 - 1d0/
     &2d0)*zp**2 + ((pi**2)/18d0 - 5d0/12d0)*zp**3 + ((pi**2)
     &/24d0 - 49d0/144d0)*zp**4 + ((pi**2)/30d0 - 41d0/144d0)
     &*zp**5 + ((pi**2)/36d0 - 5269d0/21600d0)*zp**6 - 2*zeta
     &3

         case(18)            !011

            zp = 1d0-x
            llzp = log(zp)

            ris = zp**5*(-(1d0/125d0) + (llzp)/25d0 - (
     &llzp**2)/10d0) + zp**6*(-(1d0/216d0) + (llzp)/36d0 - (l
     &lzp**2)/12d0) + zp*(-1 + llzp - (llzp**2)/2d0) + zp**2*
     &(-(1d0/8d0) + (llzp)/4d0 - (llzp**2)/4d0) + zp**3*(-(1d
     &0/27d0) + (llzp)/9d0 - (llzp**2)/6d0) + zp**4*(-(1d0/64
     &d0) + (llzp)/16d0 - (llzp**2)/8d0) + zeta3

         case(19)            !1-1-1

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**2*ll2)/12d0) + (zp*ll2)/2d0 + 
     &(ll2**3)/3d0 + zp**2*(-(1d0/16d0) + (ll2)/16d0) + zp**6
     &*(-(137d0/138240d0) + (ll2)/2304d0) + zp**4*(-(11d0/153
     &6d0) + (ll2)/256d0) + zp**3*(-(1d0/48d0) + (ll2)/72d0) 
     &+ zp**5*(-(1d0/384d0) + (ll2)/800d0) - (ll2**2*llzp)/2d
     &0 + (zeta3)/8d0

         case(20)            !1-10

            zp = 1d0-x
            llzp = log(zp)

            ris = -((zp**2)/8d0) - (zp**3)/18d0 - (5d0*
     &zp**4)/192d0 - (zp**5)/75d0 - (zp**6)/135d0 - (pi**2*ll
     &2)/4d0 + (pi**2*llzp)/12d0 + (13d0*zeta3)/8d0

         case(21)            !1-11

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**2*ll2)/6d0 - (ll2**3)/3d0 - (pi*
     &*2*llzp)/12d0 + (ll2**2*llzp)/2d0 + zp**2*(1d0/16d0 - (
     &llzp)/16d0) + zp**6*(1d0/6912d0 - (llzp)/2304d0) + zp**
     &4*(1d0/512d0 - (llzp)/256d0) + zp*(1 - (llzp)/2d0) + zp
     &**3*(1d0/108d0 - (llzp)/72d0) + zp**5*(1d0/2000d0 - (ll
     &zp)/800d0) - (7d0*zeta3)/4d0

         case(22)            !10-1

            zp = 1d0-x
            llzp = log(zp)

            ris = zp*ll2 + zp**4*(-(1d0/24d0) + (ll2)/1
     &6d0) + zp**5*(-(131d0/4800d0) + (ll2)/25d0) + zp**6*(-(
     &661d0/34560d0) + (ll2)/36d0) + zp**2*(-(1d0/8d0) + (ll2
     &)/4d0) + zp**3*(-(5d0/72d0) + (ll2)/9d0) - (pi**2*llzp)
     &/12d0 - (5d0*zeta3)/8d0

         case(23)            !100

            zp = 1d0-x

            ris = -((zp**2)/4d0) - (zp**3)/6d0 - (11d0*
     &zp**4)/96d0 - (zp**5)/12d0 - (137d0*zp**6)/2160d0 + zet
     &a3

         case(24)            !101

            zp = 1d0-x
            llzp = log(zp)

            ris = zp*(2 - llzp) - (pi**2*llzp)/6d0 + zp
     &**4*(1d0/32d0 - (llzp)/16d0) + zp**5*(2d0/125d0 - (llzp
     &)/25d0) + zp**6*(1d0/108d0 - (llzp)/36d0) + zp**2*(1d0/
     &4d0 - (llzp)/4d0) + zp**3*(2d0/27d0 - (llzp)/9d0) - 2*z
     &eta3

         case(25)            !11-1

            zp = 1d0-x
            llzp = log(zp)

            ris = -((zp)/2d0) - (zp**2)/32d0 - (zp**3)/
     &216d0 - (zp**4)/1024d0 - (zp**5)/4000d0 - (zp**6)/13824
     &d0 - (pi**2*ll2)/12d0 + (ll2**3)/6d0 + (pi**2*llzp)/12d
     &0 - (ll2**2*llzp)/2d0 + (ll2*llzp**2)/2d0 + (7d0*zeta3)
     &/8d0

         case(26)            !110

            zp = 1d0-x
            llzp = log(zp)

            ris = -zp - (zp**2)/8d0 - (zp**3)/27d0 - (z
     &p**4)/64d0 - (zp**5)/125d0 - (zp**6)/216d0 + (pi**2*llz
     &p)/6d0 + zeta3

         case(27)            !111

            zp = 1d0-x
            llzp = log(zp)

            ris = -((llzp**3)/6d0)
c End of expansions around x = +1

      end select
c --- set the imaginary part back to zero if it has been modified to
c --- get the branch cuts right (and should be zero).
      if (bcflag.eq.1) then
         xre = dreal(x)
         if (n3.eq.0.and.xre.gt.0d0) then
            if (xre.lt.1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
c
         else if (n3.eq.1.and.xre.lt.1d0) then
            if (n1.ne.-1.and.n2.ne.-1) then
               ris = dcmplx(dreal(ris),0d0)
            else if (xre.gt.-1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
c            
         else if (n3.eq.-1.and.xre.gt.-1d0) then
            if (n1.ne.1.and.n2.ne.1) then
               ris = dcmplx(dreal(ris),0d0)
            else if (xre.lt.1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
         endif
      endif    

      HPL3ar1=ris
      return
      end function
c-source-file HPL3arm1.f
c ------------------------------------------------
      double complex function HPL3arm1(n1,n2,n3,x)
      implicit none
      integer n1,n2,n3,j,bcflag,s,szp
      double complex x,ris,myi,zp,llzp
      double precision pi, zeta2, zeta3,ll2,xre

      pi=3.1415926535897932385D0
      zeta3=1.20205690315959428539973816151d0
      zeta2=pi**2/6d0
      myi = dcmplx(0d0,1d0)
      ll2 = dlog(2d0)
      bcflag = 0

      j=1+(n3+1)*1+(n2+1)*3+(n1+1)*9
      ris = dcmplx(0d0,0d0)

c---  +i*epsilon to get branch cuts right ---
      if (dimag(x).eq.0d0) then
         x = x + dcmplx(0d0,1d-60)
         bcflag = 1
      endif
c---  

      select case (j)

c This was file contains the Taylor 
c expansions around x = -1
c The expansion parameter is zp = x+1

         case(1)            !-1-1-1

            zp = x+1d0
            llzp = log(zp)

            ris = (llzp**3)/6d0

         case(2)            !-1-10

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = -zp - (zp**2)/8d0 - (zp**3)/27d0 - (z
     &p**4)/64d0 - (zp**5)/125d0 - (zp**6)/216d0 - (zp**7)/34
     &3d0 - (zp**8)/512d0 - (zp**9)/729d0 + (pi**2*llzp)/6d0 
     &+ (myi*pi*szp*llzp**2)/2d0 + zeta3

         case(3)            !-1-11

            zp = x+1d0
            llzp = log(zp)

            ris = (zp)/2d0 + (zp**2)/32d0 + (zp**3)/216
     &d0 + (zp**4)/1024d0 + (zp**5)/4000d0 + (zp**6)/13824d0 
     &+ (zp**7)/43904d0 + (zp**8)/131072d0 + (zp**9)/373248d0
     & + (pi**2*ll2)/12d0 - (ll2**3)/6d0 - (pi**2*llzp)/12d0 
     &+ (ll2**2*llzp)/2d0 - (ll2*llzp**2)/2d0 - (7d0*zeta3)/8
     &d0

         case(4)            !-10-1

            zp = x+1d0
            llzp = log(zp)

            ris = zp*(2 - llzp) - (pi**2*llzp)/6d0 + zp
     &**4*(1d0/32d0 - (llzp)/16d0) + zp**5*(2d0/125d0 - (llzp
     &)/25d0) + zp**6*(1d0/108d0 - (llzp)/36d0) + zp**7*(2d0/
     &343d0 - (llzp)/49d0) + zp**2*(1d0/4d0 - (llzp)/4d0) + z
     &p**8*(1d0/256d0 - (llzp)/64d0) + zp**9*(2d0/729d0 - (ll
     &zp)/81d0) + zp**3*(2d0/27d0 - (llzp)/9d0) - 2*zeta3

         case(5)            !-100

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = (myi*pi**3*szp)/6d0 - myi*pi*szp*zp +
     & (1d0/4d0 - (myi*pi*szp)/4d0)*zp**2 + (1d0/6d0 - (myi*p
     &i*szp)/9d0)*zp**3 + (11d0/96d0 - (myi*pi*szp)/16d0)*zp*
     &*4 + (1d0/12d0 - (myi*pi*szp)/25d0)*zp**5 + (137d0/2160
     &d0 - (myi*pi*szp)/36d0)*zp**6 + (1d0/20d0 - (myi*pi*szp
     &)/49d0)*zp**7 + (363d0/8960d0 - (myi*pi*szp)/64d0)*zp**
     &8 + (761d0/22680d0 - (myi*pi*szp)/81d0)*zp**9 - (pi**2*
     &llzp)/2d0 - zeta3

         case(6)            !-101

            zp = x+1d0
            llzp = log(zp)

            ris = zp*ll2 + zp**4*(-(1d0/24d0) + (ll2)/1
     &6d0) + zp**5*(-(131d0/4800d0) + (ll2)/25d0) + zp**6*(-(
     &661d0/34560d0) + (ll2)/36d0) + zp**7*(-(1327d0/94080d0)
     & + (ll2)/49d0) + zp**2*(-(1d0/8d0) + (ll2)/4d0) + zp**8
     &*(-(1163d0/107520d0) + (ll2)/64d0) + zp**9*(-(148969d0/
     &17418240d0) + (ll2)/81d0) + zp**3*(-(5d0/72d0) + (ll2)/
     &9d0) - (pi**2*llzp)/12d0 - (5d0*zeta3)/8d0

         case(7)            !-11-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**2*ll2)/6d0) + (ll2**3)/3d0 + (
     &pi**2*llzp)/12d0 - (ll2**2*llzp)/2d0 + zp**8*(-(1d0/655
     &36d0) + (llzp)/16384d0) + zp**2*(-(1d0/16d0) + (llzp)/1
     &6d0) + zp**6*(-(1d0/6912d0) + (llzp)/2304d0) + zp**4*(-
     &(1d0/512d0) + (llzp)/256d0) + zp*(-1 + (llzp)/2d0) + zp
     &**9*(-(1d0/186624d0) + (llzp)/41472d0) + zp**7*(-(1d0/2
     &1952d0) + (llzp)/6272d0) + zp**3*(-(1d0/108d0) + (llzp)
     &/72d0) + zp**5*(-(1d0/2000d0) + (llzp)/800d0) + (7d0*ze
     &ta3)/4d0

         case(8)            !-110

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = -((myi*pi**3*szp)/12d0) + (myi*pi*szp
     &*zp)/2d0 + (-(1d0/8d0) + (myi*pi*szp)/16d0)*zp**2 + (-(
     &1d0/18d0) + (myi*pi*szp)/72d0)*zp**3 + (-(5d0/192d0) + 
     &(myi*pi*szp)/256d0)*zp**4 + (-(1d0/75d0) + (myi*pi*szp)
     &/800d0)*zp**5 + (-(1d0/135d0) + (myi*pi*szp)/2304d0)*zp
     &**6 + (-(13d0/2940d0) + (myi*pi*szp)/6272d0)*zp**7 + (-
     &(151d0/53760d0) + (myi*pi*szp)/16384d0)*zp**8 + (-(16d0
     &/8505d0) + (myi*pi*szp)/41472d0)*zp**9 - (pi**2*ll2)/4d
     &0 + (myi*pi*szp*ll2**2)/2d0 + (pi**2*llzp)/12d0 - myi*p
     &i*szp*ll2*llzp + (13d0*zeta3)/8d0

         case(9)            !-111

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**2*ll2)/12d0 - (zp*ll2)/2d0 - (ll
     &2**3)/3d0 + zp**8*(363d0/2293760d0 - (ll2)/16384d0) + z
     &p**2*(1d0/16d0 - (ll2)/16d0) + zp**6*(137d0/138240d0 - 
     &(ll2)/2304d0) + zp**4*(11d0/1536d0 - (ll2)/256d0) + zp*
     &*9*(761d0/11612160d0 - (ll2)/41472d0) + zp**7*(1d0/2560
     &d0 - (ll2)/6272d0) + zp**3*(1d0/48d0 - (ll2)/72d0) + zp
     &**5*(1d0/384d0 - (ll2)/800d0) + (ll2**2*llzp)/2d0 - (ze
     &ta3)/8d0

         case(10)            !0-1-1

            zp = x+1d0
            llzp = log(zp)

            ris = zp**5*(-(1d0/125d0) + (llzp)/25d0 - (
     &llzp**2)/10d0) + zp**6*(-(1d0/216d0) + (llzp)/36d0 - (l
     &lzp**2)/12d0) + zp**7*(-(1d0/343d0) + (llzp)/49d0 - (ll
     &zp**2)/14d0) + zp**8*(-(1d0/512d0) + (llzp)/64d0 - (llz
     &p**2)/16d0) + zp**9*(-(1d0/729d0) + (llzp)/81d0 - (llzp
     &**2)/18d0) + zp*(-1 + llzp - (llzp**2)/2d0) + zp**2*(-(
     &1d0/8d0) + (llzp)/4d0 - (llzp**2)/4d0) + zp**3*(-(1d0/2
     &7d0) + (llzp)/9d0 - (llzp**2)/6d0) + zp**4*(-(1d0/64d0)
     & + (llzp)/16d0 - (llzp**2)/8d0) + zeta3

         case(11)            !0-10

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = -((myi*pi**3*szp)/6d0) + zp*(-((pi**2
     &)/6d0) + myi*pi*szp - myi*pi*szp*llzp) + zp**2*(-((pi**
     &2)/12d0) + 1d0/2d0 + (myi*pi*szp)/4d0 - (myi*pi*szp*llz
     &p)/2d0) + zp**3*(-((pi**2)/18d0) + 5d0/12d0 + (myi*pi*s
     &zp)/9d0 - (myi*pi*szp*llzp)/3d0) + zp**4*(-((pi**2)/24d
     &0) + 49d0/144d0 + (myi*pi*szp)/16d0 - (myi*pi*szp*llzp)
     &/4d0) + zp**5*(-((pi**2)/30d0) + 41d0/144d0 + (myi*pi*s
     &zp)/25d0 - (myi*pi*szp*llzp)/5d0) + zp**6*(-((pi**2)/36
     &d0) + 5269d0/21600d0 + (myi*pi*szp)/36d0 - (myi*pi*szp*
     &llzp)/6d0) + zp**7*(-((pi**2)/42d0) + 767d0/3600d0 + (m
     &yi*pi*szp)/49d0 - (myi*pi*szp*llzp)/7d0) + zp**8*(26668
     &1d0/1411200d0 - (pi**2)/48d0 + (myi*pi*szp)/64d0 - (myi
     &*pi*szp*llzp)/8d0) + zp**9*(-((pi**2)/54d0) + 1077749d0
     &/6350400d0 + (myi*pi*szp)/81d0 - (myi*pi*szp*llzp)/9d0)
     & + 2*zeta3

         case(12)            !0-11

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**2*ll2)/4d0 + zp*((pi**2)/12d0 - 
     &ll2 - (ll2**2)/2d0 + ll2*llzp) + zp**2*((pi**2)/24d0 - 
     &1d0/4d0 - (ll2)/4d0 - (ll2**2)/4d0 + (ll2*llzp)/2d0) + 
     &zp**3*((pi**2)/36d0 - 3d0/16d0 - (ll2)/9d0 - (ll2**2)/6
     &d0 + (ll2*llzp)/3d0) + zp**4*((pi**2)/48d0 - 83d0/576d0
     & - (ll2)/16d0 - (ll2**2)/8d0 + (ll2*llzp)/4d0) + zp**5*
     &(-(1337d0/11520d0) + (pi**2)/60d0 - (ll2)/25d0 - (ll2**
     &2)/10d0 + (ll2*llzp)/5d0) + zp**6*(-(33497d0/345600d0) 
     &+ (pi**2)/72d0 - (ll2)/36d0 - (ll2**2)/12d0 + (ll2*llzp
     &)/6d0) + zp**7*(-(5587d0/67200d0) + (pi**2)/84d0 - (ll2
     &)/49d0 - (ll2**2)/14d0 + (ll2*llzp)/7d0) + zp**8*(-(136
     &919d0/1881600d0) + (pi**2)/96d0 - (ll2)/64d0 - (ll2**2)
     &/16d0 + (ll2*llzp)/8d0) + zp**9*((pi**2)/108d0 - 350549
     &39d0/541900800d0 - (ll2)/81d0 - (ll2**2)/18d0 + (ll2*ll
     &zp)/9d0) - zeta3

         case(13)            !00-1

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**2*zp)/6d0 + zp**4*((pi**2)/24d0 
     &- 131d0/288d0 + (11d0*llzp)/24d0) + zp**2*((pi**2)/12d0
     & - 3d0/4d0 + (llzp)/2d0) + zp**3*((pi**2)/18d0 - 7d0/12
     &d0 + (llzp)/2d0) + zp**6*((pi**2)/36d0 - 2213d0/7200d0 
     &+ (137d0*llzp)/360d0) + zp**8*((pi**2)/48d0 - 647707d0/
     &2822400d0 + (363d0*llzp)/1120d0) + zp**5*((pi**2)/30d0 
     &- 53d0/144d0 + (5d0*llzp)/12d0) + zp**9*((pi**2)/54d0 -
     & 1290829d0/6350400d0 + (761d0*llzp)/2520d0) + zp**7*((p
     &i**2)/42d0 - 947d0/3600d0 + (7d0*llzp)/20d0) - zeta3

         case(14)            !000

            zp = x+1d0
            szp = s(zp)

            ris = -((myi*pi**3*szp)/6d0) + (pi**2*zp)/2
     &d0 + ((pi**2)/4d0 + (myi*pi*szp)/2d0)*zp**2 + (-(1d0/6d
     &0) + (pi**2)/6d0 + (myi*pi*szp)/2d0)*zp**3 + (-(1d0/4d0
     &) + (pi**2)/8d0 + (myi*pi*11d0*szp)/24d0)*zp**4 + ((pi*
     &*2)/10d0 - 7d0/24d0 + (myi*pi*5d0*szp)/12d0)*zp**5 + ((
     &pi**2)/12d0 - 5d0/16d0 + (myi*pi*137d0*szp)/360d0)*zp**
     &6 + ((pi**2)/14d0 - 29d0/90d0 + (myi*pi*7d0*szp)/20d0)*
     &zp**7 + ((pi**2)/16d0 - 469d0/1440d0 + (myi*pi*363d0*sz
     &p)/1120d0)*zp**8 + ((pi**2)/18d0 - 29531d0/90720d0 + (m
     &yi*pi*761d0*szp)/2520d0)*zp**9

         case(15)            !001

            zp = x+1d0

            ris = (pi**2*zp)/12d0 + zp**4*((pi**2)/48d0
     & + 11d0/96d0 - (11d0*ll2)/24d0) + zp**2*((pi**2)/24d0 -
     & (ll2)/2d0) + zp**3*(1d0/12d0 + (pi**2)/36d0 - (ll2)/2d
     &0) + zp**6*((pi**2)/72d0 + 731d0/5760d0 - (137d0*ll2)/3
     &60d0) + zp**8*(3931d0/32256d0 + (pi**2)/96d0 - (363d0*l
     &l2)/1120d0) + zp**5*((pi**2)/60d0 + 1d0/8d0 - (5d0*ll2)
     &/12d0) + zp**9*((pi**2)/108d0 + 42799d0/362880d0 - (761
     &d0*ll2)/2520d0) + zp**7*(721d0/5760d0 + (pi**2)/84d0 - 
     &(7d0*ll2)/20d0) - (3d0*zeta3)/4d0

         case(16)            !01-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**2*ll2)/4d0) + zp*(-((pi**2)/12
     &d0) + (ll2**2)/2d0) + zp**8*(314543d0/3763200d0 - (pi**
     &2)/96d0 + (ll2**2)/16d0 - (1163d0*llzp)/13440d0) + zp**
     &7*(-((pi**2)/84d0) + 953d0/9800d0 + (ll2**2)/14d0 - (13
     &27d0*llzp)/13440d0) + zp**9*(-((pi**2)/108d0) + 3572057
     &71d0/4877107200d0 + (ll2**2)/18d0 - (148969d0*llzp)/193
     &5360d0) + zp**2*(-((pi**2)/24d0) + 3d0/8d0 + (ll2**2)/4
     &d0 - (llzp)/4d0) + zp**3*(-((pi**2)/36d0) + 37d0/144d0 
     &+ (ll2**2)/6d0 - (5d0*llzp)/24d0) + zp**6*(13369d0/1152
     &00d0 - (pi**2)/72d0 + (ll2**2)/12d0 - (661d0*llzp)/5760
     &d0) + zp**4*(-((pi**2)/48d0) + 107d0/576d0 + (ll2**2)/8
     &d0 - (llzp)/6d0) + zp**5*(-((pi**2)/60d0) + 8257d0/5760
     &0d0 + (ll2**2)/10d0 - (131d0*llzp)/960d0) + (13d0*zeta3
     &)/8d0

         case(17)            !010

            zp = x+1d0
            szp = s(zp)

            ris = -((myi*pi**3*szp)/12d0) + zp*(-((pi**
     &2)/12d0) + myi*pi*szp*ll2) + zp**2*(-((pi**2)/24d0) - (
     &myi*pi*szp)/4d0 + (myi*pi*szp*ll2)/2d0) + zp**3*(1d0/12
     &d0 - (pi**2)/36d0 - (myi*pi*5d0*szp)/24d0 + (myi*pi*szp
     &*ll2)/3d0) + zp**4*(-((pi**2)/48d0) + 5d0/48d0 - (myi*p
     &i*szp)/6d0 + (myi*pi*szp*ll2)/4d0) + zp**5*(5d0/48d0 - 
     &(pi**2)/60d0 - (myi*pi*131d0*szp)/960d0 + (myi*pi*szp*l
     &l2)/5d0) + zp**6*(47d0/480d0 - (pi**2)/72d0 - (myi*pi*6
     &61d0*szp)/5760d0 + (myi*pi*szp*ll2)/6d0) + zp**7*(13d0/
     &144d0 - (pi**2)/84d0 - (myi*pi*1327d0*szp)/13440d0 + (m
     &yi*pi*szp*ll2)/7d0) + zp**8*(3341d0/40320d0 - (pi**2)/9
     &6d0 - (myi*pi*1163d0*szp)/13440d0 + (myi*pi*szp*ll2)/8d
     &0) + zp**9*(13817d0/181440d0 - (pi**2)/108d0 - (myi*pi*
     &148969d0*szp)/1935360d0 + (myi*pi*szp*ll2)/9d0) + (3d0*
     &zeta3)/2d0

         case(18)            !011

            zp = x+1d0

            ris = -((zp*ll2**2)/2d0) + zp**5*(-(83d0/19
     &20d0) + (131d0*ll2)/960d0 - (ll2**2)/10d0) + zp**6*(-(1
     &1d0/288d0) + (661d0*ll2)/5760d0 - (ll2**2)/12d0) + zp**
     &7*(-(5417d0/161280d0) + (1327d0*ll2)/13440d0 - (ll2**2)
     &/14d0) + zp**8*(-(137d0/4608d0) + (1163d0*ll2)/13440d0 
     &- (ll2**2)/16d0) + zp**9*(-(617027d0/23224320d0) + (148
     &969d0*ll2)/1935360d0 - (ll2**2)/18d0) + zp**2*((ll2)/4d
     &0 - (ll2**2)/4d0) + zp**3*(-(1d0/24d0) + (5d0*ll2)/24d0
     & - (ll2**2)/6d0) + zp**4*(-(3d0/64d0) + (ll2)/6d0 - (ll
     &2**2)/8d0) + (zeta3)/8d0

         case(19)            !1-1-1

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**2*ll2)/12d0 - (ll2**3)/6d0 + zp*
     &*4*(1d0/1024d0 - (llzp)/256d0 + (llzp**2)/128d0) + zp**
     &2*(1d0/32d0 - (llzp)/16d0 + (llzp**2)/16d0) + zp**7*(1d
     &0/43904d0 - (llzp)/6272d0 + (llzp**2)/1792d0) + zp**5*(
     &1d0/4000d0 - (llzp)/800d0 + (llzp**2)/320d0) + zp**8*(1
     &d0/131072d0 - (llzp)/16384d0 + (llzp**2)/4096d0) + zp**
     &3*(1d0/216d0 - (llzp)/72d0 + (llzp**2)/48d0) + zp*(1d0/
     &2d0 - (llzp)/2d0 + (llzp**2)/4d0) + zp**6*(1d0/13824d0 
     &- (llzp)/2304d0 + (llzp**2)/768d0) + zp**9*(1d0/373248d
     &0 - (llzp)/41472d0 + (llzp**2)/9216d0) - (7d0*zeta3)/8d
     &0

         case(20)            !1-10

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = (myi*pi**3*szp)/12d0 + (pi**2*ll2)/12
     &d0 - (myi*pi*szp*ll2**2)/2d0 + zp**5*(-(11d0/360d0) + (
     &pi**2)/960d0 - (myi*pi*szp)/800d0 + (myi*pi*szp*llzp)/1
     &60d0) + zp**8*((pi**2)/12288d0 - 9701d0/1881600d0 - (my
     &i*pi*szp)/16384d0 + (myi*pi*szp*llzp)/2048d0) + zp**3*(
     &(pi**2)/144d0 - 1d0/8d0 - (myi*pi*szp)/72d0 + (myi*pi*s
     &zp*llzp)/24d0) + zp*((pi**2)/12d0 - (myi*pi*szp)/2d0 + 
     &(myi*pi*szp*llzp)/2d0) + zp**6*((pi**2)/2304d0 - 347d0/
     &21600d0 - (myi*pi*szp)/2304d0 + (myi*pi*szp*llzp)/384d0
     &) + zp**9*((pi**2)/27648d0 - 209d0/66150d0 - (myi*pi*sz
     &p)/41472d0 + (myi*pi*szp*llzp)/4608d0) + zp**4*((pi**2)
     &/384d0 - 35d0/576d0 - (myi*pi*szp)/256d0 + (myi*pi*szp*
     &llzp)/64d0) + zp**7*(-(149d0/16800d0) + (pi**2)/5376d0 
     &- (myi*pi*szp)/6272d0 + (myi*pi*szp*llzp)/896d0) + zp**
     &2*((pi**2)/48d0 - 1d0/4d0 - (myi*pi*szp)/16d0 + (myi*pi
     &*szp*llzp)/8d0) - zeta3

         case(21)            !1-11

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**2*ll2)/12d0) + (ll2**3)/6d0 + 
     &zp**5*(-((pi**2)/1920d0) + 41d0/4608d0 + (ll2)/800d0 + 
     &(ll2**2)/320d0 - (ll2*llzp)/160d0) + zp**8*(-((pi**2)/2
     &4576d0) + 266681d0/361267200d0 + (ll2)/16384d0 + (ll2**
     &2)/4096d0 - (ll2*llzp)/2048d0) + zp**3*(-((pi**2)/288d0
     &) + 5d0/96d0 + (ll2)/72d0 + (ll2**2)/48d0 - (ll2*llzp)/
     &24d0) + zp*(-((pi**2)/24d0) + (ll2)/2d0 + (ll2**2)/4d0 
     &- (ll2*llzp)/2d0) + zp**6*(-((pi**2)/4608d0) + 5269d0/1
     &382400d0 + (ll2)/2304d0 + (ll2**2)/768d0 - (ll2*llzp)/3
     &84d0) + zp**9*(1077749d0/3251404800d0 - (pi**2)/55296d0
     & + (ll2)/41472d0 + (ll2**2)/9216d0 - (ll2*llzp)/4608d0)
     & + zp**4*(49d0/2304d0 - (pi**2)/768d0 + (ll2)/256d0 + (
     &ll2**2)/128d0 - (ll2*llzp)/64d0) + zp**7*(-((pi**2)/107
     &52d0) + 767d0/460800d0 + (ll2)/6272d0 + (ll2**2)/1792d0
     & - (ll2*llzp)/896d0) + zp**2*(1d0/8d0 - (pi**2)/96d0 + 
     &(ll2)/16d0 + (ll2**2)/16d0 - (ll2*llzp)/8d0) + (zeta3)/
     &4d0

         case(22)            !10-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**2*zp)/12d0) + (pi**2*ll2)/6d0 
     &+ zp**5*(79d0/1800d0 - (pi**2)/960d0 - (llzp)/15d0) + z
     &p**7*(521d0/39200d0 - (pi**2)/5376d0 - (13d0*llzp)/420d
     &0) + zp**6*(-((pi**2)/2304d0) + 169d0/7200d0 - (2d0*llz
     &p)/45d0) + zp**2*(-((pi**2)/48d0) + 3d0/8d0 - (llzp)/4d
     &0) + zp**4*(25d0/288d0 - (pi**2)/384d0 - (5d0*llzp)/48d
     &0) + zp**8*(-((pi**2)/12288d0) + 7493d0/940800d0 - (151
     &d0*llzp)/6720d0) + zp**3*(-((pi**2)/144d0) + 13d0/72d0 
     &- (llzp)/6d0) + zp**9*(-((pi**2)/27648d0) + 3001d0/5953
     &50d0 - (16d0*llzp)/945d0) - (5d0*zeta3)/8d0

         case(23)            !100

            zp = x+1d0
            szp = s(zp)

            ris = (myi*pi**3*szp)/12d0 - (pi**2*zp)/4d0
     & + (-((pi**2)/16d0) - (myi*pi*szp)/4d0)*zp**2 + (1d0/12
     &d0 - (pi**2)/48d0 - (myi*pi*szp)/6d0)*zp**3 + (-((pi**2
     &)/128d0) + 3d0/32d0 - (myi*pi*5d0*szp)/48d0)*zp**4 + (1
     &d0/12d0 - (pi**2)/320d0 - (myi*pi*szp)/15d0)*zp**5 + (5
     &d0/72d0 - (pi**2)/768d0 - (myi*pi*2d0*szp)/45d0)*zp**6 
     &+ (-((pi**2)/1792d0) + 41d0/720d0 - (myi*pi*13d0*szp)/4
     &20d0)*zp**7 + (-((pi**2)/4096d0) + 539d0/11520d0 - (myi
     &*pi*151d0*szp)/6720d0)*zp**8 + (22d0/567d0 - (pi**2)/92
     &16d0 - (myi*pi*16d0*szp)/945d0)*zp**9 + (pi**2*ll2)/2d0
     & - (3d0*zeta3)/4d0

         case(24)            !101

            zp = x+1d0

            ris = -((pi**2*zp)/24d0) + (pi**2*ll2)/12d0
     & + zp**5*(-((pi**2)/1920d0) - 1d0/30d0 + (ll2)/15d0) + 
     &zp**7*(-((pi**2)/10752d0) - 767d0/40320d0 + (13d0*ll2)/
     &420d0) + zp**6*(-((pi**2)/4608d0) - 97d0/3840d0 + (2d0*
     &ll2)/45d0) + zp**2*(-((pi**2)/96d0) + (ll2)/4d0) + zp**
     &4*(-(1d0/24d0) - (pi**2)/768d0 + (5d0*ll2)/48d0) + zp**
     &8*(-((pi**2)/24576d0) - 935d0/64512d0 + (151d0*ll2)/672
     &0d0) + zp**3*(-(1d0/24d0) - (pi**2)/288d0 + (ll2)/6d0) 
     &+ zp**9*(-(2041d0/181440d0) - (pi**2)/55296d0 + (16d0*l
     &l2)/945d0) - (zeta3)/4d0

         case(25)            !11-1

            zp = x+1d0
            llzp = log(zp)

            ris = (ll2**3)/6d0 + zp*((pi**2)/24d0 - (ll
     &2**2)/4d0) + zp**3*((pi**2)/288d0 - 7d0/96d0 - (ll2**2)
     &/48d0 + (llzp)/16d0) + zp**6*(-(2213d0/460800d0) + (pi*
     &*2)/4608d0 - (ll2**2)/768d0 + (137d0*llzp)/23040d0) + z
     &p**8*((pi**2)/24576d0 - 647707d0/722534400d0 - (ll2**2)
     &/4096d0 + (363d0*llzp)/286720d0) + zp**4*(-(131d0/4608d
     &0) + (pi**2)/768d0 - (ll2**2)/128d0 + (11d0*llzp)/384d0
     &) + zp**5*((pi**2)/1920d0 - 53d0/4608d0 - (ll2**2)/320d
     &0 + (5d0*llzp)/384d0) + zp**9*(-(1290829d0/3251404800d0
     &) + (pi**2)/55296d0 - (ll2**2)/9216d0 + (761d0*llzp)/12
     &90240d0) + zp**7*((pi**2)/10752d0 - 947d0/460800d0 - (l
     &l2**2)/1792d0 + (7d0*llzp)/2560d0) + zp**2*(-(3d0/16d0)
     & + (pi**2)/96d0 - (ll2**2)/16d0 + (llzp)/8d0) - (zeta3)
     &/8d0

         case(26)            !110

            zp = x+1d0
            szp = s(zp)

            ris = -((pi**2*ll2)/12d0) + (myi*pi*szp*ll2
     &**2)/2d0 + zp**5*((pi**2)/1920d0 - 1d0/40d0 + (myi*pi*5
     &d0*szp)/384d0 - (myi*pi*szp*ll2)/160d0) + zp**8*(-(1019
     &d0/161280d0) + (pi**2)/24576d0 + (myi*pi*363d0*szp)/286
     &720d0 - (myi*pi*szp*ll2)/2048d0) + zp**3*(-(1d0/24d0) +
     & (pi**2)/288d0 + (myi*pi*szp)/16d0 - (myi*pi*szp*ll2)/2
     &4d0) + zp*((pi**2)/24d0 - (myi*pi*szp*ll2)/2d0) + zp**6
     &*(-(23d0/1440d0) + (pi**2)/4608d0 + (myi*pi*137d0*szp)/
     &23040d0 - (myi*pi*szp*ll2)/384d0) + zp**9*((pi**2)/5529
     &6d0 - 23d0/5670d0 + (myi*pi*761d0*szp)/1290240d0 - (myi
     &*pi*szp*ll2)/4608d0) + zp**4*((pi**2)/768d0 - 7d0/192d0
     & + (myi*pi*11d0*szp)/384d0 - (myi*pi*szp*ll2)/64d0) + z
     &p**7*(-(101d0/10080d0) + (pi**2)/10752d0 + (myi*pi*7d0*
     &szp)/2560d0 - (myi*pi*szp*ll2)/896d0) + zp**2*((pi**2)/
     &96d0 + (myi*pi*szp)/8d0 - (myi*pi*szp*ll2)/8d0) + (zeta
     &3)/8d0

         case(27)            !111

            zp = x+1d0

            ris = (zp*ll2**2)/4d0 - (ll2**3)/6d0 + zp**
     &4*(1d0/64d0 - (11d0*ll2)/384d0 + (ll2**2)/128d0) + zp**
     &2*(-((ll2)/8d0) + (ll2**2)/16d0) + zp**7*(29d0/11520d0 
     &- (7d0*ll2)/2560d0 + (ll2**2)/1792d0) + zp**5*(7d0/768d
     &0 - (5d0*ll2)/384d0 + (ll2**2)/320d0) + zp**8*(469d0/36
     &8640d0 - (363d0*ll2)/286720d0 + (ll2**2)/4096d0) + zp**
     &3*(1d0/48d0 - (ll2)/16d0 + (ll2**2)/48d0) + zp**6*(5d0/
     &1024d0 - (137d0*ll2)/23040d0 + (ll2**2)/768d0) + zp**9*
     &(29531d0/46448640d0 - (761d0*ll2)/1290240d0 + (ll2**2)/
     &9216d0)
c End of expansions around x = -1

      end select
c --- set the imaginary part back to zero if it has been modified to
c --- get the branch cuts right (and should be zero).
      if (bcflag.eq.1) then
         xre = dreal(x)
         if (n3.eq.0.and.xre.gt.0d0) then
            if (xre.lt.1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
c
         else if (n3.eq.1.and.xre.lt.1d0) then
            if (n1.ne.-1.and.n2.ne.-1) then
               ris = dcmplx(dreal(ris),0d0)
            else if (xre.gt.-1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
c            
         else if (n3.eq.-1.and.xre.gt.-1d0) then
            if (n1.ne.1.and.n2.ne.1) then
               ris = dcmplx(dreal(ris),0d0)
            else if (xre.lt.1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
         endif
      endif    

      HPL3arm1=ris
      return
      end function
c-source-file HPL3at1.f
      double complex function HPL3at1(n1,n2,n3)
      implicit none
      integer n1,n2,n3,j
      double complex ris,myi
      double precision pi, zeta2, zeta3,ll2

      pi=3.1415926535897932385D0
      zeta3=1.20205690315959428539973816151d0
      zeta2=pi**2/6d0
      myi = dcmplx(0d0,1d0)
      ll2 = dlog(2d0)

      j=1+(n3+1)*1+(n2+1)*3+(n1+1)*9
      ris = dcmplx(0d0,0d0)

      if(j.le.18.or.j.eq.23) then
         select case (j)
         case (1)
            ris =ll2**3d0/6d0
         case (2)
            ris =-(pi**2d0*ll2)/12d0 + zeta3/8d0
         case (3)
            ris =-ll2**3d0/6d0 + zeta3/8d0
         case (4)
            ris =(pi**2d0*ll2)/12d0 - zeta3/4d0
         case (5)
            ris =(3d0*zeta3)/4d0
         case (6)
            ris =(pi**2d0*ll2)/6d0 - (5d0*zeta3)/8d0
         case (7)
            ris =(pi**2d0*ll2)/12d0 
     &           - ll2**3d0/6d0-zeta3/4d0
         case (8)
            ris =(pi**2d0*ll2)/12d0 - zeta3
         case (9)
            ris =-(pi**2*ll2)/12d0+ll2**3/6d0
     &           +(7d0*zeta3)/8d0
         case (10)
            ris =zeta3/8d0
         case (11)
            ris =(-3d0*zeta3)/2d0
         case (12)
            ris =-(pi**2d0*ll2)/4d0 + (13d0*zeta3)/8d0
         case (13)
            ris =(3d0*zeta3)/4d0
         case (14)
            ris =0d0
         case (15)
            ris =zeta3
         case (16)
            ris =(pi**2d0*ll2)/4d0 - zeta3
         case (17)
            ris =-2d0*zeta3
         case (18)
            ris =zeta3
         case (23)
            ris =zeta3
         end select
      else
         print*, ""
         print*, "****************"
         print*, "ERROR in HPL3: "
         print*, "HPL3(",n1,",",n2,",",n3
     &        ,",1) is divergent!"
         print*, "Aborting..."
         print*,"****************"
         stop
      endif
      HPL3at1=ris
      return
      end function
c-source-file HPL3atm1.f
      double complex function HPL3atm1(n1,n2,n3)
      implicit none
      integer n1,n2,n3,j
      double complex ris,myi
      double precision pi, zeta2, zeta3,ll2

      pi=3.1415926535897932385D0
      zeta3=1.20205690315959428539973816151d0
      zeta2=pi**2/6d0
      myi = dcmplx(0d0,1d0)
      ll2 = dlog(2d0)

      j=1+(n3+1)*1+(n2+1)*3+(n1+1)*9
      ris = dcmplx(0d0,0d0)

      if(j.ge.10) then
         select case (j)
         case (10)
            ris = zeta3
         case (11)
            ris = 2*zeta3 - myi*pi*zeta2
         case(12)
            ris = pi**2*ll2/4d0 - zeta3
         case(13)
            ris = -zeta3
         case(14)
            ris = -myi*pi*zeta2
         case(15)
            ris = 3*zeta3/4d0
         case(16)
            ris = -pi**2*ll2/4d0 + 13d0*zeta3/8d0
         case(17)
            ris = 3*zeta3/2d0 - myi*pi**3/12d0
         case(18)
            ris = zeta3/8d0
         case(19)
            ris = 0.5d0*zeta2*ll2 - ll2**3/6d0 
     &           - 7d0*zeta3/8d0
         case(20)
            ris = 0.5d0*zeta2*ll2 + myi*pi*(0.5d0*zeta2 
     &           - 0.5d0*ll2**2) - zeta3
         case(21)
            ris = -0.5d0*zeta2*ll2 + ll2**3/6d0 
     &           + zeta3/4d0
         case(22)
            ris = zeta2*ll2 - 5d0*zeta3/8d0
         case(23)
            ris = 0.5d0*myi*pi*zeta2 + pi**2*ll2/2d0 
     &           - 3*zeta3/4d0
         case(24)
            ris = 0.5d0*zeta2*ll2 - zeta3/4d0
         case(25)
            ris = ll2**3/6d0 - zeta3/8d0
         case(26)
            ris = -0.5d0*zeta2*ll2 + zeta3/8d0 
     &           + 0.5d0*myi*pi*ll2**2
         case(27)
            ris = -1d0/6d0*ll2**3
         end select
      else
         print*, ""
         print*, "****************"
         print*, "ERROR in HPL3: "
         print*, "HPL3(",n1,",",n2,",",n3
     &        ,",-1) is divergent!"
         print*, "Aborting..."
         print*,"****************"
         stop
      endif
      HPL3atm1=ris
      return
      end function
c-source-file HPL3else.f
      double complex function HPL3else(n1, n2, n3, x)
      implicit none
      double precision pi, zeta2, zeta3,ll2,xre
      double complex x, ris
      double complex basis3_1,basis3_2,basis3_3,basis3_4
      double complex basis3_5,basis3_6,basis3_7,basis3_8
      double complex basis2_1,basis2_2,basis2_3
      double complex cli2,cli3
      double complex ll1px,ll1mx,llx
      integer n1,n2,n3,j,bcflag

      pi=3.1415926535897932385D0
      zeta3=1.20205690315959428539973816151d0
      zeta2=pi**2/6d0
      ll2 = dlog(2d0)

      j=1+(n3+1)*1+(n2+1)*3+(n1+1)*9

      ris=dcmplx(0d0,0d0)
      bcflag = 0

c---  +i*epsilon to get branch cuts right ---
      if (dimag(x).eq.0d0) then
         x = x + dcmplx(0d0,1d-60)
         bcflag = 1
      endif
c---  
c     #####################################################
      ll1px = log(1d0+x)
      ll1mx = log(1d0-x)
      llx = log(x)
      basis3_1 = cli3(x) 
      basis3_2 = cli3(-x)
      basis3_3 = cli3(1d0-x)
      basis3_4 = cli3(1d0/(1d0+x)) 
      basis3_5 = cli3((1d0+x)/2d0) 
      basis3_6 = cli3((1d0-x)/2d0) 
      basis3_7 = cli3((1d0-x)/(1d0+x)) 
      basis3_8 = cli3(2d0*x/(x-1d0))
      basis2_1 = cli2(x)
      basis2_2 = cli2(-x) 
      basis2_3 = cli2((1d0-x)/2d0)
c     #####################################################

      select case(j)
      case(1)
         ris = ll1px**3/6d0
      case(2)
         ris = -(pi**2*ll1px)/6d0 + ll1px**3/6d0 
     &- basis3_4 
     &        + zeta3
      case(3)
         ris = (pi**2*ll2)/12d0 - ll2**3/6d0 
     &- (pi**2*ll1px)/12d0 + (ll2**2*ll1px)/2d0 
     &- (ll2*ll1px**2)/2d0 + basis3_5 - (7*zeta3)/8d0
      case(4)
         ris = (pi**2*ll1px)/3d0 + llx*ll1px**2 
     &- ll1px**3/3d0 + ll1px*basis2_2+2*basis3_4
     &-2*zeta3
      case(5)
         ris = (llx**2*ll1px)/2d0+llx*basis2_2
     &-basis3_2
      case(6)
         ris = (pi**2*ll2)/6d0 - ll2**3/3d0 
     &- (pi**2*ll1mx)/12d0 + (ll2**2*ll1mx)/2d0 
     &- (pi**2*ll1px)/12d0 + (ll2**2*ll1px)/2d0 
     &- ll2*ll1mx*ll1px - ll1mx*llx*ll1px 
     &+(ll1mx*ll1px**2)/2d0-ll1mx*basis2_2
     &+basis3_6 
     &- basis3_3-basis3_4+ basis3_7 + basis3_5
     &-(3*zeta3)/4d0
      case(7)
         ris = -(pi**2*ll2)/6d0 + ll2**3/3d0 
     &+ (pi**2*ll1px)/4d0 - (3*ll2**2*ll1px)/2d0 
     &+ ll2*ll1mx*ll1px + ll2*ll1px**2 
     &- ll1mx*ll1px**2 - ll1px*basis2_3 
     &- 2*basis3_5 + (7*zeta3)/4d0
      case(8)
         ris = -(pi**2*ll2)/12d0 + ll2**3/6d0 
     &+ (pi**2*ll1mx)/6d0 - (ll2*ll1mx**2)/2d0 
     &+ ll1mx**3/6d0+(pi**2*llx)/12d0- (ll2**2*llx)/2d0 
     &+ ll2*ll1mx*llx - (ll1mx**2*llx)/2d0 
     &+ (pi**2*ll1px)/12d0 - (ll2**2*ll1px)/2d0 
     &+ ll2*ll1mx*ll1px - (ll1mx*ll1px**2)/2d0 
     &- llx*basis2_3 + basis3_2-basis3_1-basis3_8 
     &+ basis3_4 - basis3_7 - basis3_5 + (7*zeta3)/8d0
      case(9)
         ris = -(pi**2*ll2)/12d0 + ll2**3/6d0 
     &- (ll2*ll1mx**2)/2d0 + (ll1mx**2*ll1px)/2d0 
     &+ ll1mx*basis2_3 - basis3_6 + (7*zeta3)/8d0
      case(10)
         ris = -(pi**2*ll1px)/6d0 - (llx*ll1px**2)/2d0 
     &+ ll1px**3/6d0 - ll1px*basis2_2 - basis3_4 
     &+ zeta3
      case(11)
         ris = -(llx*basis2_2) + 2*basis3_2
      case(12)
         ris = -(pi**2*ll2)/12d0 + ll2**3/6d0 
     &- (pi**2*ll1mx)/12d0 - (ll2**2*ll1mx)/2d0 
     &+ (ll2*ll1mx**2)/2d0 - ll1mx**3/6d0 
     &+ (ll1mx**2*llx)/2d0 + ll1mx*basis2_2 
     &- basis3_6 + basis3_3 - basis3_2 +basis3_1
     &+basis3_8 
     &- zeta3/8d0
      case(13)
         ris = -basis3_2
      case(14)
         ris = llx**3/6d0
      case(15)
         ris = basis3_1
      case(16)
         ris = -(pi**2*ll2)/12d0 + ll2**3/6d0 
     &+ (pi**2*ll1mx)/6d0 - (ll2*ll1mx**2)/2d0 
     &+ ll1mx**3/6d0 - (ll1mx**2*llx)/2d0 
     &+ (pi**2*ll1px)/12d0 - (ll2**2*ll1px)/2d0 
     &+ ll2*ll1mx*ll1px + ll1mx*llx*ll1px 
     &-(ll1mx*ll1px**2)/2d0+ll1px*basis2_1
     &+basis3_2 
     &- basis3_1 - basis3_8 + basis3_4 -basis3_7
     &-basis3_5 
     &+ (7*zeta3)/8d0
      case(17)
         ris = llx*basis2_1 - 2*basis3_1
      case(18)
         ris = (pi**2*ll1mx)/6d0 - (ll1mx**2*llx)/2d0 
     &- ll1mx*basis2_1 - basis3_3 + zeta3
      case(19)
         ris = (pi**2*ll2)/12d0 - ll2**3/6d0 
     &- (pi**2*ll1px)/6d0 + ll2**2*ll1px 
     &- ll2*ll1mx*ll1px - (ll2*ll1px**2)/2d0 
     &+ (ll1mx*ll1px**2)/2d0 + ll1px*basis2_3 
     &+ basis3_5 - (7*zeta3)/8d0
      case(20)
         ris = -(pi**2*ll2)/12d0 + ll2**3/6d0 
     &- (pi**2*ll1mx)/12d0 - (ll2**2*ll1mx)/2d0 
     &+ (ll2*ll1mx**2)/2d0 - ll1mx**3/6d0 
     &- (pi**2*llx)/12d0 + (ll2**2*llx)/2d0 
     &- ll2*ll1mx*llx + (ll1mx**2*llx)/2d0 
     &+ llx*basis2_3 - basis3_6 + basis3_3 
     &- basis3_2 + basis3_1 + basis3_8 - zeta3/8d0
      case(21)
         ris = (pi**2*ll2)/6d0 - ll2**3/3d0 
     &- (pi**2*ll1mx)/12d0 + (ll2**2*ll1mx)/2d0 
     &- ll1mx*basis2_3 + 2*basis3_6 - (7*zeta3)/4d0
      case(22)
         ris = (pi**2*ll2)/6d0 - ll2**3/3d0 
     &- (pi**2*ll1mx)/12d0 + (ll2**2*ll1mx)/2d0 
     &- (pi**2*ll1px)/12d0 + (ll2**2*ll1px)/2d0 
     &- ll2*ll1mx*ll1px - ll1mx*llx*ll1px 
     &+(ll1mx*ll1px**2)/2d0-ll1px*basis2_1
     &+basis3_6 
     &- basis3_3 - basis3_4+basis3_7+basis3_5
     &- (3*zeta3)/4d0
      case(23)
         ris =-(ll1mx*llx**2)/2d0-llx*basis2_1
     &+basis3_1
      case(24)
         ris = -(pi**2*ll1mx)/3d0 + ll1mx**2*llx 
     &+ ll1mx*basis2_1 + 2*basis3_3 - 2*zeta3
      case(25)
         ris = -(pi**2*ll2)/12d0 + ll2**3/6d0 
     &+ (pi**2*ll1mx)/12d0 - (ll2**2*ll1mx)/2d0 
     &+ (ll2*ll1mx**2)/2d0 - basis3_6 + (7*zeta3)/8d0
      case(26)
         ris = (pi**2*ll1mx)/6d0 - basis3_3 + zeta3
      case(27)
         ris = -ll1mx**3/6d0
      end select
c --- set the imaginary part back to zero if it has been modified to
c --- get the branch cuts right (and should be zero). Also, set imaginary
c --- part of result to zero if x is between 0 and 1.

      if (bcflag.eq.1) then
         x = x - dcmplx(0d0,1d-60)
         xre = dreal(x)
         if (xre.ge.0d0.and.xre.le.1d0) then
            ris = dcmplx(dreal(ris),0d0)
         endif
      endif  
      
      HPL3else=ris
      return
      end
c-source-file HPL3.f
C=============================================================================
C---  HPLs  of Rank 3  
C=============================================================================
c--- main forking function
      double complex function HPL3(n1, n2, n3, x)
      implicit none
      integer n1,n2,n3
      double complex x,ris
      double complex HPL3at0,HPL3at1,HPL3atm1
      double complex HPL3ar1,HPL3arm1,HPL3ar0,HPL3else
      double precision rad1,radm1,rad0

      rad1 = 0.01d0
      radm1 = 0.025d0
      rad0 = 0.025d0
      
      if ((abs(n1).gt.1).or.(abs(n2).gt.1).or.(abs(n3).gt.1)) then
         print*, ""
         print*, "****************"
         print*, "Error in HPL3:"
         print*, "Indices",n1,n2,n3," out of range !"
         print*, "Aborting..."
         print*,"****************"
         stop
      endif

      ris = dcmplx(0d0,0d0)
      
      if (x.eq.dcmplx(0d0,0d0)) then
         ris = HPL3at0(n1,n2,n3)
      elseif (x.eq.dcmplx(1d0,0d0)) then
         ris = HPL3at1(n1,n2,n3)
      elseif (x.eq.dcmplx(-1d0,0d0)) then
         ris = HPL3atm1(n1,n2,n3)
      elseif (abs(x-dcmplx(1d0,0d0)).lt.rad1) then
         ris = HPL3ar1(n1,n2,n3,x)
      elseif (abs(x+dcmplx(1d0,0d0)).lt.radm1) then
         ris = HPL3arm1(n1,n2,n3,x)
      elseif (abs(x-dcmplx(0d0,0d0)).lt.rad0) then
         ris = HPL3ar0(n1,n2,n3,x)
      else 
         ris = HPL3else(n1,n2,n3,x)
      endif
      HPL3=ris
      return
      end function
c ------------------------------------------------
      double complex function HPL3at0(n1, n2, n3)
      implicit none
      integer n1,n2,n3,j
      double complex ris
    
      j=1+(n3+1)*1+(n2+1)*3+(n1+1)*9
      ris = dcmplx(0d0,0d0)

      if (j.eq.14) then             
         print*, ""
         print*, "****************"
         print*, "ERROR in HPL3: "
         print*, "HPL3(",n1,",",n2,",",n3
     &        ,",0) is divergent!"
         print*, "Aborting..."
         print*,"****************"
         stop
      endif
      HPL3at0=ris
      return
      end function
c ------------------------------------------------
c --- Real part of HPL3     
      double precision function HPL3real(n1,n2,n3,xr,xi)
      implicit none
      double precision xr,xi
      integer n1,n2,n3
      double complex x,HPL3
      x=dcmplx(xr,xi)
      HPL3real = dreal(HPL3(n1,n2,n3,x))
      return
      end

c --- Imaginary part of HPL3     
      double precision function HPL3im(n1,n2,n3,xr,xi)
      implicit none
      double precision xr,xi
      integer n1,n2,n3
      double complex x,HPL3
      x=dcmplx(xr,xi)
      HPL3im = dimag(HPL3(n1,n2,n3,x))
      return
      end
c-source-file HPL4ar0.f
      double complex function HPL4ar0(n1,n2,n3,n4,x)
      implicit none
      integer n1,n2,n3,n4,j,bcflag
      double complex x,ris,myi,cli4pt5,cli4,llx
      double precision pi, zeta2, zeta3,zeta4,ll2,xre

      pi=3.1415926535897932385D0
      zeta3=1.20205690315959428539973816151d0
      zeta2=pi**2/6d0
      zeta4=pi**4/90d0
      myi = dcmplx(0d0,1d0)
      bcflag = 0

      ll2 = dlog(2d0)
      cli4pt5 = cli4(dcmplx(0.5d0,0d0))

      j=1+(n4+1)*1+(n3+1)*3+(n2+1)*9+(n1+1)*27
      ris = dcmplx(0d0,0d0)

c---  +i*epsilon to get branch cuts right ---
      if (dimag(x).eq.0d0) then
         x = x + dcmplx(0d0,1d-60)
         bcflag = 1
      endif
c---  
      select case (j)
c This was file contains the Taylor 
c expansions around x = 0

         case(1)            !-1-1-1-1


            ris = (x**4)/24d0 - (x**5)/12d0 + (17d0*x**
     &6)/144d0 - (7d0*x**7)/48d0 + (967d0*x**8)/5760d0 - (89d
     &0*x**9)/480d0 + (4523d0*x**10)/22680d0

         case(2)            !-1-1-10

            llx = log(x)

            ris = -((11d0*x**3)/36d0) + (19d0*x**4)/48d
     &0 - (599d0*x**5)/1440d0 + (79d0*x**6)/192d0 - (3343d0*x
     &**7)/8400d0 + (21977d0*x**8)/57600d0 - (83359739d0*x**9
     &)/228614400d0 + (3538531d0*x**10)/10160640d0 + (x**3*ll
     &x)/6d0 - (x**4*llx)/4d0 + (7d0*x**5*llx)/24d0 - (5d0*x*
     &*6*llx)/16d0 + (29d0*x**7*llx)/90d0 - (469d0*x**8*llx)/
     &1440d0 + (29531d0*x**9*llx)/90720d0 - (1303d0*x**10*llx
     &)/4032d0

         case(3)            !-1-1-11


            ris = (x**4)/24d0 - (x**5)/15d0 + (61d0*x**
     &6)/720d0 - (97d0*x**7)/1008d0 + (467d0*x**8)/4480d0 - (
     &9931d0*x**9)/90720d0 + (1831d0*x**10)/16200d0

         case(4)            !-1-10-1


            ris = (x**3)/6d0 - (11d0*x**4)/48d0 + (181d
     &0*x**5)/720d0 - (37d0*x**6)/144d0 + (38569d0*x**7)/1512
     &00d0 - (43171d0*x**8)/172800d0 + (9261559d0*x**9)/38102
     &400d0 - (1197607d0*x**10)/5080320d0

         case(5)            !-1-100

            llx = log(x)

            ris = (7d0*x**2)/8d0 - (41d0*x**3)/72d0 + (
     &1397d0*x**4)/3456d0 - (2671d0*x**5)/8640d0 + (322493d0*
     &x**6)/1296000d0 - (104641d0*x**7)/504000d0 + (140539517
     &d0*x**8)/790272000d0 - (2486560891d0*x**9)/16003008000d
     &0 + (22064922487d0*x**10)/160030080000d0 - (3d0*x**2*ll
     &x)/4d0 + (7d0*x**3*llx)/12d0 - (131d0*x**4*llx)/288d0 +
     & (53d0*x**5*llx)/144d0 - (2213d0*x**6*llx)/7200d0 + (94
     &7d0*x**7*llx)/3600d0 - (647707d0*x**8*llx)/2822400d0 + 
     &(1290829d0*x**9*llx)/6350400d0 - (11574649d0*x**10*llx)
     &/63504000d0 + (x**2*llx**2)/4d0 - (x**3*llx**2)/4d0 + (
     &11d0*x**4*llx**2)/48d0 - (5d0*x**5*llx**2)/24d0 + (137d
     &0*x**6*llx**2)/720d0 - (7d0*x**7*llx**2)/40d0 + (363d0*
     &x**8*llx**2)/2240d0 - (761d0*x**9*llx**2)/5040d0 + (712
     &9d0*x**10*llx**2)/50400d0

         case(6)            !-1-101


            ris = (x**3)/6d0 - (3d0*x**4)/16d0 + (139d0
     &*x**5)/720d0 - (3d0*x**6)/16d0 + (27319d0*x**7)/151200d
     &0 - (29821d0*x**8)/172800d0 + (6284809d0*x**9)/38102400
     &d0 - (88913d0*x**10)/564480d0

         case(7)            !-1-11-1


            ris = (x**4)/24d0 - (x**5)/20d0 + (43d0*x**
     &6)/720d0 - (103d0*x**7)/1680d0 + (2563d0*x**8)/40320d0 
     &- (127d0*x**9)/2016d0 + (3569d0*x**10)/56700d0

         case(8)            !-1-110

            llx = log(x)

            ris = -((11d0*x**3)/36d0) + (11d0*x**4)/36d
     &0 - (719d0*x**5)/2400d0 + (1337d0*x**6)/4800d0 - (46061
     &d0*x**7)/176400d0 + (27479d0*x**8)/112896d0 - (52243973
     &d0*x**9)/228614400d0 + (163666943d0*x**10)/762048000d0 
     &+ (x**3*llx)/6d0 - (x**4*llx)/6d0 + (7d0*x**5*llx)/40d0
     & - (119d0*x**6*llx)/720d0 + (101d0*x**7*llx)/630d0 - (3
     &05d0*x**8*llx)/2016d0 + (13157d0*x**9*llx)/90720d0 - (4
     &1603d0*x**10*llx)/302400d0

         case(9)            !-1-111


            ris = (x**4)/24d0 - (x**5)/30d0 + (31d0*x**
     &6)/720d0 - (191d0*x**7)/5040d0 + (547d0*x**8)/13440d0 -
     & (3373d0*x**9)/90720d0 + (2147d0*x**10)/56700d0

         case(10)            !-10-1-1


            ris = (x**3)/12d0 - (5d0*x**4)/48d0 + (17d0
     &*x**5)/160d0 - (59d0*x**6)/576d0 + (2929d0*x**7)/30240d
     &0 - (629d0*x**8)/6912d0 + (185921d0*x**9)/2177280d0 - (
     &116423d0*x**10)/1451520d0

         case(11)            !-10-10

            llx = log(x)

            ris = -((5d0*x**2)/4d0) + (8d0*x**3)/9d0 - 
     &(1151d0*x**4)/1728d0 + (2281d0*x**5)/4320d0 - (17653d0*
     &x**6)/40500d0 + (93371d0*x**7)/252000d0 - (127203607d0*
     &x**8)/395136000d0 + (2276013631d0*x**9)/8001504000d0 - 
     &(4076031341d0*x**10)/16003008000d0 + (x**2*llx)/2d0 - (
     &5d0*x**3*llx)/12d0 + (49d0*x**4*llx)/144d0 - (41d0*x**5
     &*llx)/144d0 + (5269d0*x**6*llx)/21600d0 - (767d0*x**7*l
     &lx)/3600d0 + (266681d0*x**8*llx)/1411200d0 - (1077749d0
     &*x**9*llx)/6350400d0 + (9778141d0*x**10*llx)/63504000d0

         case(12)            !-10-11


            ris = (x**3)/12d0 - (11d0*x**4)/144d0 + (10
     &3d0*x**5)/1440d0 - (2743d0*x**6)/43200d0 + (8699d0*x**7
     &)/151200d0 - (439571d0*x**8)/8467200d0 + (3617053d0*x**
     &9)/76204800d0 - (33150437d0*x**10)/762048000d0

         case(13)            !-100-1


            ris = (x**2)/2d0 - (3d0*x**3)/8d0 + (251d0*
     &x**4)/864d0 - (407d0*x**5)/1728d0 + (256103d0*x**6)/129
     &6000d0 - (4081d0*x**7)/24000d0 + (9822481d0*x**8)/65856
     &000d0 - (78708473d0*x**9)/592704000d0 + (19148110939d0*
     &x**10)/160030080000d0

         case(14)            !-1000

            llx = log(x)

            ris = -x + (x**2)/16d0 - (x**3)/81d0 + (x**
     &4)/256d0 - (x**5)/625d0 + (x**6)/1296d0 - (x**7)/2401d0
     & + (x**8)/4096d0 - (x**9)/6561d0 + (x**10)/10000d0 + x*
     &llx - (x**2*llx)/8d0 + (x**3*llx)/27d0 - (x**4*llx)/64d
     &0 + (x**5*llx)/125d0 - (x**6*llx)/216d0 + (x**7*llx)/34
     &3d0 - (x**8*llx)/512d0 + (x**9*llx)/729d0 - (x**10*llx)
     &/1000d0 - (x*llx**2)/2d0 + (x**2*llx**2)/8d0 - (x**3*ll
     &x**2)/18d0 + (x**4*llx**2)/32d0 - (x**5*llx**2)/50d0 + 
     &(x**6*llx**2)/72d0 - (x**7*llx**2)/98d0 + (x**8*llx**2)
     &/128d0 - (x**9*llx**2)/162d0 + (x**10*llx**2)/200d0 + (
     &x*llx**3)/6d0 - (x**2*llx**3)/12d0 + (x**3*llx**3)/18d0
     & - (x**4*llx**3)/24d0 + (x**5*llx**3)/30d0 - (x**6*llx*
     &*3)/36d0 + (x**7*llx**3)/42d0 - (x**8*llx**3)/48d0 + (x
     &**9*llx**3)/54d0 - (x**10*llx**3)/60d0

         case(15)            !-1001


            ris = (x**2)/2d0 - (7d0*x**3)/24d0 + (197d0
     &*x**4)/864d0 - (1549d0*x**5)/8640d0 + (195353d0*x**6)/1
     &296000d0 - (194353d0*x**7)/1512000d0 + (66879079d0*x**8
     &)/592704000d0 - (533875007d0*x**9)/5334336000d0 + (1443
     &6577189d0*x**10)/160030080000d0

         case(16)            !-101-1


            ris = (x**3)/12d0 - (7d0*x**4)/144d0 + (71d
     &0*x**5)/1440d0 - (1607d0*x**6)/43200d0 + (5291d0*x**7)/
     &151200d0 - (245939d0*x**8)/8467200d0 + (2067997d0*x**9)
     &/76204800d0 - (18015013d0*x**10)/762048000d0

         case(17)            !-1010

            llx = log(x)

            ris = -((5d0*x**2)/4d0) + (2d0*x**3)/3d0 - 
     &(881d0*x**4)/1728d0 + (1687d0*x**5)/4320d0 - (13153d0*x
     &**6)/40500d0 + (206863d0*x**7)/756000d0 - (282912571d0*
     &x**8)/1185408000d0 + (560731627d0*x**9)/2667168000d0 - 
     &(3019814291d0*x**10)/16003008000d0 + (x**2*llx)/2d0 - (
     &x**3*llx)/4d0 + (31d0*x**4*llx)/144d0 - (23d0*x**5*llx)
     &/144d0 + (3019d0*x**6*llx)/21600d0 - (139d0*x**7*llx)/1
     &200d0 + (48877d0*x**8*llx)/470400d0 - (191833d0*x**9*ll
     &x)/2116800d0 + (5257891d0*x**10*llx)/63504000d0

         case(18)            !-1011


            ris = (x**3)/12d0 - (x**4)/48d0 + (19d0*x**
     &5)/480d0 - (11d0*x**6)/576d0 + (769d0*x**7)/30240d0 - (
     &553d0*x**8)/34560d0 + (40769d0*x**9)/2177280d0 - (19591
     &d0*x**10)/1451520d0

         case(19)            !-11-1-1


            ris = (x**4)/24d0 - (x**5)/30d0 + (31d0*x**
     &6)/720d0 - (181d0*x**7)/5040d0 + (1571d0*x**8)/40320d0 
     &- (113d0*x**9)/3360d0 + (15727d0*x**10)/453600d0

         case(20)            !-11-10

            llx = log(x)

            ris = -((11d0*x**3)/36d0) + (19d0*x**4)/144
     &d0 - (1181d0*x**5)/7200d0 + (1439d0*x**6)/14400d0 - (43
     &1d0*x**7)/3920d0 + (218839d0*x**8)/2822400d0 - (1880026
     &1d0*x**9)/228614400d0 + (54725d0*x**10)/870912d0 + (x**
     &3*llx)/6d0 - (x**4*llx)/12d0 + (13d0*x**5*llx)/120d0 - 
     &(17d0*x**6*llx)/240d0 + (5d0*x**7*llx)/63d0 - (589d0*x*
     &*8*llx)/10080d0 + (5669d0*x**9*llx)/90720d0 - (85d0*x**
     &10*llx)/1728d0

         case(21)            !-11-11


            ris = (x**4)/24d0 - (x**5)/60d0 + (23d0*x**
     &6)/720d0 - (29d0*x**7)/1680d0 + (1009d0*x**8)/40320d0 -
     & (1429d0*x**9)/90720d0 + (1853d0*x**10)/90720d0

         case(22)            !-110-1


            ris = (x**3)/6d0 - (x**4)/16d0 + (67d0*x**5
     &)/720d0 - (11d0*x**6)/216d0 + (9619d0*x**7)/151200d0 - 
     &(7117d0*x**8)/172800d0 + (73393d0*x**9)/1524096d0 - (14
     &51159d0*x**10)/42336000d0

         case(23)            !-1100

            llx = log(x)

            ris = (7d0*x**2)/8d0 - (85d0*x**3)/216d0 + 
     &(1019d0*x**4)/3456d0 - (46633d0*x**5)/216000d0 + (23024
     &3d0*x**6)/1296000d0 - (10882477d0*x**7)/74088000d0 + (3
     &01825801d0*x**8)/2370816000d0 - (5330081423d0*x**9)/480
     &09024000d0 + (15880889737d0*x**10)/160030080000d0 - (3d
     &0*x**2*llx)/4d0 + (11d0*x**3*llx)/36d0 - (77d0*x**4*llx
     &)/288d0 + (659d0*x**5*llx)/3600d0 - (1163d0*x**6*llx)/7
     &200d0 + (2517d0*x**7*llx)/19600d0 - (108919d0*x**8*llx)
     &/940800d0 + (1875737d0*x**9*llx)/19051200d0 - (5731399d
     &0*x**10*llx)/63504000d0 + (x**2*llx**2)/4d0 - (x**3*llx
     &**2)/12d0 + (5d0*x**4*llx**2)/48d0 - (7d0*x**5*llx**2)/
     &120d0 + (47d0*x**6*llx**2)/720d0 - (37d0*x**7*llx**2)/8
     &40d0 + (319d0*x**8*llx**2)/6720d0 - (533d0*x**9*llx**2)
     &/15120d0 + (1879d0*x**10*llx**2)/50400d0

         case(24)            !-1101


            ris = (x**3)/6d0 - (x**4)/48d0 + (61d0*x**5
     &)/720d0 - (5d0*x**6)/216d0 + (8269d0*x**7)/151200d0 - (
     &3667d0*x**8)/172800d0 + (60751d0*x**9)/1524096d0 - (240
     &0827d0*x**10)/127008000d0

         case(25)            !-111-1


            ris = (x**4)/24d0 + (7d0*x**6)/240d0 - (x**
     &7)/720d0 + (857d0*x**8)/40320d0 - (x**9)/480d0 + (7429d
     &0*x**10)/453600d0

         case(26)            !-1110

            llx = log(x)

            ris = -((11d0*x**3)/36d0) + (x**4)/24d0 - (
     &1027d0*x**5)/7200d0 + (25d0*x**6)/576d0 - (15653d0*x**7
     &)/176400d0 + (2209d0*x**8)/57600d0 - (2902399d0*x**9)/4
     &5722880d0 + (8469731d0*x**10)/254016000d0 + (x**3*llx)/
     &6d0 + (11d0*x**5*llx)/120d0 - (x**6*llx)/144d0 + (19d0*
     &x**7*llx)/315d0 - (13d0*x**8*llx)/1440d0 + (799d0*x**9*
     &llx)/18144d0 - (317d0*x**10*llx)/33600d0

         case(27)            !-1111


            ris = (x**4)/24d0 + (x**5)/60d0 + (5d0*x**6
     &)/144d0 + (5d0*x**7)/336d0 + (157d0*x**8)/5760d0 + (31d
     &0*x**9)/2592d0 + (9883d0*x**10)/453600d0

         case(28)            !0-1-1-1


            ris = (x**3)/18d0 - (x**4)/16d0 + (7d0*x**5
     &)/120d0 - (5d0*x**6)/96d0 + (29d0*x**7)/630d0 - (469d0*
     &x**8)/11520d0 + (29531d0*x**9)/816480d0 - (1303d0*x**10
     &)/40320d0

         case(29)            !0-1-10

            llx = log(x)

            ris = -((x**2)/2d0) + (x**3)/4d0 - (41d0*x*
     &*4)/288d0 + (13d0*x**5)/144d0 - (8009d0*x**6)/129600d0 
     &+ (161d0*x**7)/3600d0 - (190513d0*x**8)/5644800d0 + (16
     &7101d0*x**9)/6350400d0 - (13371157d0*x**10)/635040000d0
     & + (x**2*llx)/4d0 - (x**3*llx)/6d0 + (11d0*x**4*llx)/96
     &d0 - (x**5*llx)/12d0 + (137d0*x**6*llx)/2160d0 - (x**7*
     &llx)/20d0 + (363d0*x**8*llx)/8960d0 - (761d0*x**9*llx)/
     &22680d0 + (7129d0*x**10*llx)/252000d0

         case(30)            !0-1-11


            ris = (x**3)/18d0 - (x**4)/24d0 + (7d0*x**5
     &)/200d0 - (119d0*x**6)/4320d0 + (101d0*x**7)/4410d0 - (
     &305d0*x**8)/16128d0 + (13157d0*x**9)/816480d0 - (41603d
     &0*x**10)/3024000d0

         case(31)            !0-10-1


            ris = (x**2)/4d0 - (5d0*x**3)/36d0 + (49d0*
     &x**4)/576d0 - (41d0*x**5)/720d0 + (5269d0*x**6)/129600d
     &0 - (767d0*x**7)/25200d0 + (266681d0*x**8)/11289600d0 -
     & (1077749d0*x**9)/57153600d0 + (9778141d0*x**10)/635040
     &000d0

         case(32)            !0-100

            llx = log(x)

            ris = 3*x - (3d0*x**2)/16d0 + (x**3)/27d0 -
     & (3d0*x**4)/256d0 + (3d0*x**5)/625d0 - (x**6)/432d0 + (
     &3d0*x**7)/2401d0 - (3d0*x**8)/4096d0 + (x**9)/2187d0 - 
     &(3d0*x**10)/10000d0 - 2*x*llx + (x**2*llx)/4d0 - (2d0*x
     &**3*llx)/27d0 + (x**4*llx)/32d0 - (2d0*x**5*llx)/125d0 
     &+ (x**6*llx)/108d0 - (2d0*x**7*llx)/343d0 + (x**8*llx)/
     &256d0 - (2d0*x**9*llx)/729d0 + (x**10*llx)/500d0 + (x*l
     &lx**2)/2d0 - (x**2*llx**2)/8d0 + (x**3*llx**2)/18d0 - (
     &x**4*llx**2)/32d0 + (x**5*llx**2)/50d0 - (x**6*llx**2)/
     &72d0 + (x**7*llx**2)/98d0 - (x**8*llx**2)/128d0 + (x**9
     &*llx**2)/162d0 - (x**10*llx**2)/200d0

         case(33)            !0-101


            ris = (x**2)/4d0 - (x**3)/12d0 + (31d0*x**4
     &)/576d0 - (23d0*x**5)/720d0 + (3019d0*x**6)/129600d0 - 
     &(139d0*x**7)/8400d0 + (48877d0*x**8)/3763200d0 - (19183
     &3d0*x**9)/19051200d0 + (5257891d0*x**10)/635040000d0

         case(34)            !0-11-1


            ris = (x**3)/18d0 - (x**4)/48d0 + (13d0*x**
     &5)/600d0 - (17d0*x**6)/1440d0 + (5d0*x**7)/441d0 - (589
     &d0*x**8)/80640d0 + (5669d0*x**9)/816480d0 - (17d0*x**10
     &)/3456d0

         case(35)            !0-110

            llx = log(x)

            ris = -((x**2)/2d0) + (13d0*x**3)/108d0 - (
     &23d0*x**4)/288d0 + (743d0*x**5)/18000d0 - (3959d0*x**6)
     &/129600d0 + (8291d0*x**7)/411600d0 - (10007d0*x**8)/627
     &200d0 + (2024977d0*x**9)/171460800d0 - (6204907d0*x**10
     &)/635040000d0 + (x**2*llx)/4d0 - (x**3*llx)/18d0 + (5d0
     &*x**4*llx)/96d0 - (7d0*x**5*llx)/300d0 + (47d0*x**6*llx
     &)/2160d0 - (37d0*x**7*llx)/2940d0 + (319d0*x**8*llx)/26
     &880d0 - (533d0*x**9*llx)/68040d0 + (1879d0*x**10*llx)/2
     &52000d0

         case(36)            !0-111


            ris = (x**3)/18d0 + (11d0*x**5)/600d0 - (x*
     &*6)/864d0 + (19d0*x**7)/2205d0 - (13d0*x**8)/11520d0 + 
     &(799d0*x**9)/163296d0 - (317d0*x**10)/336000d0

         case(37)            !00-1-1


            ris = (x**2)/8d0 - (x**3)/18d0 + (11d0*x**4
     &)/384d0 - (x**5)/60d0 + (137d0*x**6)/12960d0 - (x**7)/1
     &40d0 + (363d0*x**8)/71680d0 - (761d0*x**9)/204120d0 + (
     &7129d0*x**10)/2520000d0

         case(38)            !00-10

            llx = log(x)

            ris = -3*x + (3d0*x**2)/16d0 - (x**3)/27d0 
     &+ (3d0*x**4)/256d0 - (3d0*x**5)/625d0 + (x**6)/432d0 - 
     &(3d0*x**7)/2401d0 + (3d0*x**8)/4096d0 - (x**9)/2187d0 +
     & (3d0*x**10)/10000d0 + x*llx - (x**2*llx)/8d0 + (x**3*l
     &lx)/27d0 - (x**4*llx)/64d0 + (x**5*llx)/125d0 - (x**6*l
     &lx)/216d0 + (x**7*llx)/343d0 - (x**8*llx)/512d0 + (x**9
     &*llx)/729d0 - (x**10*llx)/1000d0

         case(39)            !00-11


            ris = (x**2)/8d0 - (x**3)/54d0 + (5d0*x**4)
     &/384d0 - (7d0*x**5)/1500d0 + (47d0*x**6)/12960d0 - (37d
     &0*x**7)/20580d0 + (319d0*x**8)/215040d0 - (533d0*x**9)/
     &612360d0 + (1879d0*x**10)/2520000d0

         case(40)            !000-1


            ris = x - (x**2)/16d0 + (x**3)/81d0 - (x**4
     &)/256d0 + (x**5)/625d0 - (x**6)/1296d0 + (x**7)/2401d0 
     &- (x**8)/4096d0 + (x**9)/6561d0 - (x**10)/10000d0

         case(41)            !0000

            llx = log(x)

            ris = (llx**4)/24d0

         case(42)            !0001


            ris = x + (x**2)/16d0 + (x**3)/81d0 + (x**4
     &)/256d0 + (x**5)/625d0 + (x**6)/1296d0 + (x**7)/2401d0 
     &+ (x**8)/4096d0 + (x**9)/6561d0 + (x**10)/10000d0

         case(43)            !001-1


            ris = (x**2)/8d0 + (x**3)/54d0 + (5d0*x**4)
     &/384d0 + (7d0*x**5)/1500d0 + (47d0*x**6)/12960d0 + (37d
     &0*x**7)/20580d0 + (319d0*x**8)/215040d0 + (533d0*x**9)/
     &612360d0 + (1879d0*x**10)/2520000d0

         case(44)            !0010

            llx = log(x)

            ris = -3*x - (3d0*x**2)/16d0 - (x**3)/27d0 
     &- (3d0*x**4)/256d0 - (3d0*x**5)/625d0 - (x**6)/432d0 - 
     &(3d0*x**7)/2401d0 - (3d0*x**8)/4096d0 - (x**9)/2187d0 -
     & (3d0*x**10)/10000d0 + x*llx + (x**2*llx)/8d0 + (x**3*l
     &lx)/27d0 + (x**4*llx)/64d0 + (x**5*llx)/125d0 + (x**6*l
     &lx)/216d0 + (x**7*llx)/343d0 + (x**8*llx)/512d0 + (x**9
     &*llx)/729d0 + (x**10*llx)/1000d0

         case(45)            !0011


            ris = (x**2)/8d0 + (x**3)/18d0 + (11d0*x**4
     &)/384d0 + (x**5)/60d0 + (137d0*x**6)/12960d0 + (x**7)/1
     &40d0 + (363d0*x**8)/71680d0 + (761d0*x**9)/204120d0 + (
     &7129d0*x**10)/2520000d0

         case(46)            !01-1-1


            ris = (x**3)/18d0 + (11d0*x**5)/600d0 + (x*
     &*6)/864d0 + (19d0*x**7)/2205d0 + (13d0*x**8)/11520d0 + 
     &(799d0*x**9)/163296d0 + (317d0*x**10)/336000d0

         case(47)            !01-10

            llx = log(x)

            ris = -((x**2)/2d0) - (13d0*x**3)/108d0 - (
     &23d0*x**4)/288d0 - (743d0*x**5)/18000d0 - (3959d0*x**6)
     &/129600d0 - (8291d0*x**7)/411600d0 - (10007d0*x**8)/627
     &200d0 - (2024977d0*x**9)/171460800d0 - (6204907d0*x**10
     &)/635040000d0 + (x**2*llx)/4d0 + (x**3*llx)/18d0 + (5d0
     &*x**4*llx)/96d0 + (7d0*x**5*llx)/300d0 + (47d0*x**6*llx
     &)/2160d0 + (37d0*x**7*llx)/2940d0 + (319d0*x**8*llx)/26
     &880d0 + (533d0*x**9*llx)/68040d0 + (1879d0*x**10*llx)/2
     &52000d0

         case(48)            !01-11


            ris = (x**3)/18d0 + (x**4)/48d0 + (13d0*x**
     &5)/600d0 + (17d0*x**6)/1440d0 + (5d0*x**7)/441d0 + (589
     &d0*x**8)/80640d0 + (5669d0*x**9)/816480d0 + (17d0*x**10
     &)/3456d0

         case(49)            !010-1


            ris = (x**2)/4d0 + (x**3)/12d0 + (31d0*x**4
     &)/576d0 + (23d0*x**5)/720d0 + (3019d0*x**6)/129600d0 + 
     &(139d0*x**7)/8400d0 + (48877d0*x**8)/3763200d0 + (19183
     &3d0*x**9)/19051200d0 + (5257891d0*x**10)/635040000d0

         case(50)            !0100

            llx = log(x)

            ris = 3*x + (3d0*x**2)/16d0 + (x**3)/27d0 +
     & (3d0*x**4)/256d0 + (3d0*x**5)/625d0 + (x**6)/432d0 + (
     &3d0*x**7)/2401d0 + (3d0*x**8)/4096d0 + (x**9)/2187d0 + 
     &(3d0*x**10)/10000d0 - 2*x*llx - (x**2*llx)/4d0 - (2d0*x
     &**3*llx)/27d0 - (x**4*llx)/32d0 - (2d0*x**5*llx)/125d0 
     &- (x**6*llx)/108d0 - (2d0*x**7*llx)/343d0 - (x**8*llx)/
     &256d0 - (2d0*x**9*llx)/729d0 - (x**10*llx)/500d0 + (x*l
     &lx**2)/2d0 + (x**2*llx**2)/8d0 + (x**3*llx**2)/18d0 + (
     &x**4*llx**2)/32d0 + (x**5*llx**2)/50d0 + (x**6*llx**2)/
     &72d0 + (x**7*llx**2)/98d0 + (x**8*llx**2)/128d0 + (x**9
     &*llx**2)/162d0 + (x**10*llx**2)/200d0

         case(51)            !0101


            ris = (x**2)/4d0 + (5d0*x**3)/36d0 + (49d0*
     &x**4)/576d0 + (41d0*x**5)/720d0 + (5269d0*x**6)/129600d
     &0 + (767d0*x**7)/25200d0 + (266681d0*x**8)/11289600d0 +
     & (1077749d0*x**9)/57153600d0 + (9778141d0*x**10)/635040
     &000d0

         case(52)            !011-1


            ris = (x**3)/18d0 + (x**4)/24d0 + (7d0*x**5
     &)/200d0 + (119d0*x**6)/4320d0 + (101d0*x**7)/4410d0 + (
     &305d0*x**8)/16128d0 + (13157d0*x**9)/816480d0 + (41603d
     &0*x**10)/3024000d0

         case(53)            !0110

            llx = log(x)

            ris = -((x**2)/2d0) - (x**3)/4d0 - (41d0*x*
     &*4)/288d0 - (13d0*x**5)/144d0 - (8009d0*x**6)/129600d0 
     &- (161d0*x**7)/3600d0 - (190513d0*x**8)/5644800d0 - (16
     &7101d0*x**9)/6350400d0 - (13371157d0*x**10)/635040000d0
     & + (x**2*llx)/4d0 + (x**3*llx)/6d0 + (11d0*x**4*llx)/96
     &d0 + (x**5*llx)/12d0 + (137d0*x**6*llx)/2160d0 + (x**7*
     &llx)/20d0 + (363d0*x**8*llx)/8960d0 + (761d0*x**9*llx)/
     &22680d0 + (7129d0*x**10*llx)/252000d0

         case(54)            !0111


            ris = (x**3)/18d0 + (x**4)/16d0 + (7d0*x**5
     &)/120d0 + (5d0*x**6)/96d0 + (29d0*x**7)/630d0 + (469d0*
     &x**8)/11520d0 + (29531d0*x**9)/816480d0 + (1303d0*x**10
     &)/40320d0

         case(55)            !1-1-1-1


            ris = (x**4)/24d0 - (x**5)/60d0 + (5d0*x**6
     &)/144d0 - (5d0*x**7)/336d0 + (157d0*x**8)/5760d0 - (31d
     &0*x**9)/2592d0 + (9883d0*x**10)/453600d0

         case(56)            !1-1-10

            llx = log(x)

            ris = -((11d0*x**3)/36d0) - (x**4)/24d0 - (
     &1027d0*x**5)/7200d0 - (25d0*x**6)/576d0 - (15653d0*x**7
     &)/176400d0 - (2209d0*x**8)/57600d0 - (2902399d0*x**9)/4
     &5722880d0 - (8469731d0*x**10)/254016000d0 + (x**3*llx)/
     &6d0 + (11d0*x**5*llx)/120d0 + (x**6*llx)/144d0 + (19d0*
     &x**7*llx)/315d0 + (13d0*x**8*llx)/1440d0 + (799d0*x**9*
     &llx)/18144d0 + (317d0*x**10*llx)/33600d0

         case(57)            !1-1-11


            ris = (x**4)/24d0 + (7d0*x**6)/240d0 + (x**
     &7)/720d0 + (857d0*x**8)/40320d0 + (x**9)/480d0 + (7429d
     &0*x**10)/453600d0

         case(58)            !1-10-1


            ris = (x**3)/6d0 + (x**4)/48d0 + (61d0*x**5
     &)/720d0 + (5d0*x**6)/216d0 + (8269d0*x**7)/151200d0 + (
     &3667d0*x**8)/172800d0 + (60751d0*x**9)/1524096d0 + (240
     &0827d0*x**10)/127008000d0

         case(59)            !1-100

            llx = log(x)

            ris = (7d0*x**2)/8d0 + (85d0*x**3)/216d0 + 
     &(1019d0*x**4)/3456d0 + (46633d0*x**5)/216000d0 + (23024
     &3d0*x**6)/1296000d0 + (10882477d0*x**7)/74088000d0 + (3
     &01825801d0*x**8)/2370816000d0 + (5330081423d0*x**9)/480
     &09024000d0 + (15880889737d0*x**10)/160030080000d0 - (3d
     &0*x**2*llx)/4d0 - (11d0*x**3*llx)/36d0 - (77d0*x**4*llx
     &)/288d0 - (659d0*x**5*llx)/3600d0 - (1163d0*x**6*llx)/7
     &200d0 - (2517d0*x**7*llx)/19600d0 - (108919d0*x**8*llx)
     &/940800d0 - (1875737d0*x**9*llx)/19051200d0 - (5731399d
     &0*x**10*llx)/63504000d0 + (x**2*llx**2)/4d0 + (x**3*llx
     &**2)/12d0 + (5d0*x**4*llx**2)/48d0 + (7d0*x**5*llx**2)/
     &120d0 + (47d0*x**6*llx**2)/720d0 + (37d0*x**7*llx**2)/8
     &40d0 + (319d0*x**8*llx**2)/6720d0 + (533d0*x**9*llx**2)
     &/15120d0 + (1879d0*x**10*llx**2)/50400d0

         case(60)            !1-101


            ris = (x**3)/6d0 + (x**4)/16d0 + (67d0*x**5
     &)/720d0 + (11d0*x**6)/216d0 + (9619d0*x**7)/151200d0 + 
     &(7117d0*x**8)/172800d0 + (73393d0*x**9)/1524096d0 + (14
     &51159d0*x**10)/42336000d0

         case(61)            !1-11-1


            ris = (x**4)/24d0 + (x**5)/60d0 + (23d0*x**
     &6)/720d0 + (29d0*x**7)/1680d0 + (1009d0*x**8)/40320d0 +
     & (1429d0*x**9)/90720d0 + (1853d0*x**10)/90720d0

         case(62)            !1-110

            llx = log(x)

            ris = -((11d0*x**3)/36d0) - (19d0*x**4)/144
     &d0 - (1181d0*x**5)/7200d0 - (1439d0*x**6)/14400d0 - (43
     &1d0*x**7)/3920d0 - (218839d0*x**8)/2822400d0 - (1880026
     &1d0*x**9)/228614400d0 - (54725d0*x**10)/870912d0 + (x**
     &3*llx)/6d0 + (x**4*llx)/12d0 + (13d0*x**5*llx)/120d0 + 
     &(17d0*x**6*llx)/240d0 + (5d0*x**7*llx)/63d0 + (589d0*x*
     &*8*llx)/10080d0 + (5669d0*x**9*llx)/90720d0 + (85d0*x**
     &10*llx)/1728d0

         case(63)            !1-111


            ris = (x**4)/24d0 + (x**5)/30d0 + (31d0*x**
     &6)/720d0 + (181d0*x**7)/5040d0 + (1571d0*x**8)/40320d0 
     &+ (113d0*x**9)/3360d0 + (15727d0*x**10)/453600d0

         case(64)            !10-1-1


            ris = (x**3)/12d0 + (x**4)/48d0 + (19d0*x**
     &5)/480d0 + (11d0*x**6)/576d0 + (769d0*x**7)/30240d0 + (
     &553d0*x**8)/34560d0 + (40769d0*x**9)/2177280d0 + (19591
     &d0*x**10)/1451520d0

         case(65)            !10-10

            llx = log(x)

            ris = -((5d0*x**2)/4d0) - (2d0*x**3)/3d0 - 
     &(881d0*x**4)/1728d0 - (1687d0*x**5)/4320d0 - (13153d0*x
     &**6)/40500d0 - (206863d0*x**7)/756000d0 - (282912571d0*
     &x**8)/1185408000d0 - (560731627d0*x**9)/2667168000d0 - 
     &(3019814291d0*x**10)/16003008000d0 + (x**2*llx)/2d0 + (
     &x**3*llx)/4d0 + (31d0*x**4*llx)/144d0 + (23d0*x**5*llx)
     &/144d0 + (3019d0*x**6*llx)/21600d0 + (139d0*x**7*llx)/1
     &200d0 + (48877d0*x**8*llx)/470400d0 + (191833d0*x**9*ll
     &x)/2116800d0 + (5257891d0*x**10*llx)/63504000d0

         case(66)            !10-11


            ris = (x**3)/12d0 + (7d0*x**4)/144d0 + (71d
     &0*x**5)/1440d0 + (1607d0*x**6)/43200d0 + (5291d0*x**7)/
     &151200d0 + (245939d0*x**8)/8467200d0 + (2067997d0*x**9)
     &/76204800d0 + (18015013d0*x**10)/762048000d0

         case(67)            !100-1


            ris = (x**2)/2d0 + (7d0*x**3)/24d0 + (197d0
     &*x**4)/864d0 + (1549d0*x**5)/8640d0 + (195353d0*x**6)/1
     &296000d0 + (194353d0*x**7)/1512000d0 + (66879079d0*x**8
     &)/592704000d0 + (533875007d0*x**9)/5334336000d0 + (1443
     &6577189d0*x**10)/160030080000d0

         case(68)            !1000

            llx = log(x)

            ris = -x - (x**2)/16d0 - (x**3)/81d0 - (x**
     &4)/256d0 - (x**5)/625d0 - (x**6)/1296d0 - (x**7)/2401d0
     & - (x**8)/4096d0 - (x**9)/6561d0 - (x**10)/10000d0 + x*
     &llx + (x**2*llx)/8d0 + (x**3*llx)/27d0 + (x**4*llx)/64d
     &0 + (x**5*llx)/125d0 + (x**6*llx)/216d0 + (x**7*llx)/34
     &3d0 + (x**8*llx)/512d0 + (x**9*llx)/729d0 + (x**10*llx)
     &/1000d0 - (x*llx**2)/2d0 - (x**2*llx**2)/8d0 - (x**3*ll
     &x**2)/18d0 - (x**4*llx**2)/32d0 - (x**5*llx**2)/50d0 - 
     &(x**6*llx**2)/72d0 - (x**7*llx**2)/98d0 - (x**8*llx**2)
     &/128d0 - (x**9*llx**2)/162d0 - (x**10*llx**2)/200d0 + (
     &x*llx**3)/6d0 + (x**2*llx**3)/12d0 + (x**3*llx**3)/18d0
     & + (x**4*llx**3)/24d0 + (x**5*llx**3)/30d0 + (x**6*llx*
     &*3)/36d0 + (x**7*llx**3)/42d0 + (x**8*llx**3)/48d0 + (x
     &**9*llx**3)/54d0 + (x**10*llx**3)/60d0

         case(69)            !1001


            ris = (x**2)/2d0 + (3d0*x**3)/8d0 + (251d0*
     &x**4)/864d0 + (407d0*x**5)/1728d0 + (256103d0*x**6)/129
     &6000d0 + (4081d0*x**7)/24000d0 + (9822481d0*x**8)/65856
     &000d0 + (78708473d0*x**9)/592704000d0 + (19148110939d0*
     &x**10)/160030080000d0

         case(70)            !101-1


            ris = (x**3)/12d0 + (11d0*x**4)/144d0 + (10
     &3d0*x**5)/1440d0 + (2743d0*x**6)/43200d0 + (8699d0*x**7
     &)/151200d0 + (439571d0*x**8)/8467200d0 + (3617053d0*x**
     &9)/76204800d0 + (33150437d0*x**10)/762048000d0

         case(71)            !1010

            llx = log(x)

            ris = -((5d0*x**2)/4d0) - (8d0*x**3)/9d0 - 
     &(1151d0*x**4)/1728d0 - (2281d0*x**5)/4320d0 - (17653d0*
     &x**6)/40500d0 - (93371d0*x**7)/252000d0 - (127203607d0*
     &x**8)/395136000d0 - (2276013631d0*x**9)/8001504000d0 - 
     &(4076031341d0*x**10)/16003008000d0 + (x**2*llx)/2d0 + (
     &5d0*x**3*llx)/12d0 + (49d0*x**4*llx)/144d0 + (41d0*x**5
     &*llx)/144d0 + (5269d0*x**6*llx)/21600d0 + (767d0*x**7*l
     &lx)/3600d0 + (266681d0*x**8*llx)/1411200d0 + (1077749d0
     &*x**9*llx)/6350400d0 + (9778141d0*x**10*llx)/63504000d0

         case(72)            !1011


            ris = (x**3)/12d0 + (5d0*x**4)/48d0 + (17d0
     &*x**5)/160d0 + (59d0*x**6)/576d0 + (2929d0*x**7)/30240d
     &0 + (629d0*x**8)/6912d0 + (185921d0*x**9)/2177280d0 + (
     &116423d0*x**10)/1451520d0

         case(73)            !11-1-1


            ris = (x**4)/24d0 + (x**5)/30d0 + (31d0*x**
     &6)/720d0 + (191d0*x**7)/5040d0 + (547d0*x**8)/13440d0 +
     & (3373d0*x**9)/90720d0 + (2147d0*x**10)/56700d0

         case(74)            !11-10

            llx = log(x)

            ris = -((11d0*x**3)/36d0) - (11d0*x**4)/36d
     &0 - (719d0*x**5)/2400d0 - (1337d0*x**6)/4800d0 - (46061
     &d0*x**7)/176400d0 - (27479d0*x**8)/112896d0 - (52243973
     &d0*x**9)/228614400d0 - (163666943d0*x**10)/762048000d0 
     &+ (x**3*llx)/6d0 + (x**4*llx)/6d0 + (7d0*x**5*llx)/40d0
     & + (119d0*x**6*llx)/720d0 + (101d0*x**7*llx)/630d0 + (3
     &05d0*x**8*llx)/2016d0 + (13157d0*x**9*llx)/90720d0 + (4
     &1603d0*x**10*llx)/302400d0

         case(75)            !11-11


            ris = (x**4)/24d0 + (x**5)/20d0 + (43d0*x**
     &6)/720d0 + (103d0*x**7)/1680d0 + (2563d0*x**8)/40320d0 
     &+ (127d0*x**9)/2016d0 + (3569d0*x**10)/56700d0

         case(76)            !110-1


            ris = (x**3)/6d0 + (3d0*x**4)/16d0 + (139d0
     &*x**5)/720d0 + (3d0*x**6)/16d0 + (27319d0*x**7)/151200d
     &0 + (29821d0*x**8)/172800d0 + (6284809d0*x**9)/38102400
     &d0 + (88913d0*x**10)/564480d0

         case(77)            !1100

            llx = log(x)

            ris = (7d0*x**2)/8d0 + (41d0*x**3)/72d0 + (
     &1397d0*x**4)/3456d0 + (2671d0*x**5)/8640d0 + (322493d0*
     &x**6)/1296000d0 + (104641d0*x**7)/504000d0 + (140539517
     &d0*x**8)/790272000d0 + (2486560891d0*x**9)/16003008000d
     &0 + (22064922487d0*x**10)/160030080000d0 - (3d0*x**2*ll
     &x)/4d0 - (7d0*x**3*llx)/12d0 - (131d0*x**4*llx)/288d0 -
     & (53d0*x**5*llx)/144d0 - (2213d0*x**6*llx)/7200d0 - (94
     &7d0*x**7*llx)/3600d0 - (647707d0*x**8*llx)/2822400d0 - 
     &(1290829d0*x**9*llx)/6350400d0 - (11574649d0*x**10*llx)
     &/63504000d0 + (x**2*llx**2)/4d0 + (x**3*llx**2)/4d0 + (
     &11d0*x**4*llx**2)/48d0 + (5d0*x**5*llx**2)/24d0 + (137d
     &0*x**6*llx**2)/720d0 + (7d0*x**7*llx**2)/40d0 + (363d0*
     &x**8*llx**2)/2240d0 + (761d0*x**9*llx**2)/5040d0 + (712
     &9d0*x**10*llx**2)/50400d0

         case(78)            !1101


            ris = (x**3)/6d0 + (11d0*x**4)/48d0 + (181d
     &0*x**5)/720d0 + (37d0*x**6)/144d0 + (38569d0*x**7)/1512
     &00d0 + (43171d0*x**8)/172800d0 + (9261559d0*x**9)/38102
     &400d0 + (1197607d0*x**10)/5080320d0

         case(79)            !111-1


            ris = (x**4)/24d0 + (x**5)/15d0 + (61d0*x**
     &6)/720d0 + (97d0*x**7)/1008d0 + (467d0*x**8)/4480d0 + (
     &9931d0*x**9)/90720d0 + (1831d0*x**10)/16200d0

         case(80)            !1110

            llx = log(x)

            ris = -((11d0*x**3)/36d0) - (19d0*x**4)/48d
     &0 - (599d0*x**5)/1440d0 - (79d0*x**6)/192d0 - (3343d0*x
     &**7)/8400d0 - (21977d0*x**8)/57600d0 - (83359739d0*x**9
     &)/228614400d0 - (3538531d0*x**10)/10160640d0 + (x**3*ll
     &x)/6d0 + (x**4*llx)/4d0 + (7d0*x**5*llx)/24d0 + (5d0*x*
     &*6*llx)/16d0 + (29d0*x**7*llx)/90d0 + (469d0*x**8*llx)/
     &1440d0 + (29531d0*x**9*llx)/90720d0 + (1303d0*x**10*llx
     &)/4032d0

         case(81)            !1111


            ris = (x**4)/24d0 + (x**5)/12d0 + (17d0*x**
     &6)/144d0 + (7d0*x**7)/48d0 + (967d0*x**8)/5760d0 + (89d
     &0*x**9)/480d0 + (4523d0*x**10)/22680d0
c End of expansions around x = 0

      end select
c --- set the imaginary part back to zero if it has been modified to
c --- get the branch cuts right (and should be zero).
      if (bcflag.eq.1) then
         xre = dreal(x)

         if (n4.eq.0.and.xre.gt.0d0) then
            if (xre.lt.1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
c
         else if (n4.eq.1.and.xre.lt.1d0) then
            if (n1.ne.-1.and.n2.ne.-1.and.n3.ne.-1) then
               ris = dcmplx(dreal(ris),0d0)
            else if (xre.gt.-1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
c            
         else if (n4.eq.-1.and.xre.gt.-1d0) then
            if (n1.ne.1.and.n2.ne.1.and.n3.ne.1) then
               ris = dcmplx(dreal(ris),0d0)
            else if (xre.lt.1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
            
         endif
      endif

      HPL4ar0=ris
      return
      end function
c-source-file HPL4ar1.f
      double complex function HPL4ar1(n1,n2,n3,n4,x)
      implicit none
      integer n1,n2,n3,n4,j,bcflag
      double complex x,ris,myi,cli4pt5,cli4,zp,llzp
      double precision pi, zeta2, zeta3,zeta4,ll2,xre

      pi=3.1415926535897932385D0
      zeta3=1.20205690315959428539973816151d0
      zeta2=pi**2/6d0
      zeta4=pi**4/90d0
      myi = dcmplx(0d0,1d0)
      bcflag = 0

      ll2 = dlog(2d0)
      cli4pt5 = cli4(dcmplx(0.5d0,0d0))

      j=1+(n4+1)*1+(n3+1)*3+(n2+1)*9+(n1+1)*27
      ris = dcmplx(0d0,0d0)

c---  +i*epsilon to get branch cuts right ---
      if (dimag(x).eq.0d0) then
         x = x + dcmplx(0d0,1d-60)
         bcflag = 1
      endif
c---  
      select case (j)

c This was file contains the Taylor 
c expansions around x = +1
c The expansion parameter is zp = 1-x

         case(1)            !-1-1-1-1

            zp = 1d0-x

            ris = -((zp*ll2**3)/12d0) + (ll2**4)/24d0 +
     & zp**3*(-((ll2)/48d0) + (ll2**2)/32d0 - (ll2**3)/144d0)
     & + zp**6*(17d0/9216d0 - (5d0*ll2)/1024d0 + (137d0*ll2**
     &2)/46080d0 - (ll2**3)/2304d0) + zp**4*(1d0/384d0 - (ll2
     &)/64d0 + (11d0*ll2**2)/768d0 - (ll2**3)/384d0) + zp**2*
     &((ll2**2)/16d0 - (ll2**3)/48d0) + zp**5*(1d0/384d0 - (7
     &d0*ll2)/768d0 + (5d0*ll2**2)/768d0 - (ll2**3)/960d0)

         case(2)            !-1-1-10

            zp = 1d0-x

            ris = -((pi**4)/90d0) - (pi**2*ll2**2)/12d0
     & + (ll2**4)/24d0 + cli4pt5 + ll2*zeta3 + zp**5*(11d0/19
     &20d0 - (pi**2*5d0)/4608d0 + (pi**2*ll2)/1920d0 - (zeta3
     &)/1280d0) + zp*((pi**2*ll2)/24d0 - (zeta3)/16d0) + zp**
     &3*(-((pi**2)/192d0) + (pi**2*ll2)/288d0 - (zeta3)/192d0
     &) + zp**6*(103d0/23040d0 - (pi**2*137d0)/276480d0 + (pi
     &**2*ll2)/4608d0 - (zeta3)/3072d0) + zp**4*(1d0/192d0 - 
     &(pi**2*11d0)/4608d0 + (pi**2*ll2)/768d0 - (zeta3)/512d0
     &) + zp**2*(-((pi**2)/96d0) + (pi**2*ll2)/96d0 - (zeta3)
     &/64d0)

         case(3)            !-1-1-11

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**4)/90d0 + (pi**2*ll2**2)/24d0 - 
     &(ll2**4)/12d0 - cli4pt5 - (7d0*ll2*zeta3)/8d0 + zp**5*(
     &-(599d0/46080d0) + (pi**2*5d0)/4608d0 - (5d0*ll2**2)/76
     &8d0 + (ll2**3)/960d0 + (7d0*llzp)/768d0 - (zeta3)/1280d
     &0) + zp*((ll2**3)/12d0 - (zeta3)/16d0) + zp**3*((pi**2)
     &/192d0 - 11d0/288d0 - (ll2**2)/32d0 + (ll2**3)/144d0 + 
     &(llzp)/48d0 - (zeta3)/192d0) + zp**6*((pi**2*137d0)/276
     &480d0 - 79d0/12288d0 - (137d0*ll2**2)/46080d0 + (ll2**3
     &)/2304d0 + (5d0*llzp)/1024d0 - (zeta3)/3072d0) + zp**4*
     &((pi**2*11d0)/4608d0 - 19d0/768d0 - (11d0*ll2**2)/768d0
     & + (ll2**3)/384d0 + (llzp)/64d0 - (zeta3)/512d0) + zp**
     &2*((pi**2)/96d0 - (ll2**2)/16d0 + (ll2**3)/48d0 - (zeta
     &3)/64d0)

         case(4)            !-1-10-1

            zp = 1d0-x

            ris = (pi**4)/30d0 + (pi**2*ll2**2)/6d0 - (
     &ll2**4)/8d0 - 3*cli4pt5 - (23d0*ll2*zeta3)/8d0 + zp**6*
     &((pi**2*137d0)/276480d0 + 31d0/5760d0 - (23d0*ll2)/1440
     &d0 - (pi**2*ll2)/4608d0 + (zeta3)/1536d0) + zp**4*(1d0/
     &192d0 + (pi**2*11d0)/4608d0 - (pi**2*ll2)/768d0 - (7d0*
     &ll2)/192d0 + (zeta3)/256d0) + zp**2*((pi**2)/96d0 - (pi
     &**2*ll2)/96d0 + (zeta3)/32d0) + zp**5*(1d0/160d0 + (pi*
     &*2*5d0)/4608d0 - (pi**2*ll2)/1920d0 - (ll2)/40d0 + (zet
     &a3)/640d0) + zp*(-((pi**2*ll2)/24d0) + (zeta3)/8d0) + z
     &p**3*((pi**2)/192d0 - (ll2)/24d0 - (pi**2*ll2)/288d0 + 
     &(zeta3)/96d0)

         case(5)            !-1-100

            zp = 1d0-x

            ris = (pi**4)/48d0 + (pi**2*ll2**2)/12d0 - 
     &(ll2**4)/12d0 - 2*cli4pt5 - (3d0*zp*zeta3)/8d0 - (3d0*z
     &p**2*zeta3)/32d0 - (zp**3*zeta3)/32d0 - ll2*zeta3 + zp*
     &*4*(1d0/96d0 - (3d0*zeta3)/256d0) + zp**6*(29d0/2304d0 
     &- (zeta3)/512d0) + zp**5*(13d0/960d0 - (3d0*zeta3)/640d
     &0)

         case(6)            !-1-101

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**4)/480d0 + (pi**2*ll2**2)/12d0 -
     & (5d0*ll2*zeta3)/8d0 + zp**5*((pi**2*5d0)/2304d0 - 77d0
     &/2400d0 - (pi**2*ll2)/960d0 + (llzp)/40d0 + (zeta3)/256
     &d0) + zp*(-((pi**2*ll2)/12d0) + (5d0*zeta3)/16d0) + zp*
     &*3*(-(11d0/144d0) + (pi**2)/96d0 - (pi**2*ll2)/144d0 + 
     &(llzp)/24d0 + (5d0*zeta3)/192d0) + zp**6*((pi**2*137d0)
     &/138240d0 - 169d0/9600d0 - (pi**2*ll2)/2304d0 + (23d0*l
     &lzp)/1440d0 + (5d0*zeta3)/3072d0) + zp**4*((pi**2*11d0)
     &/2304d0 - 127d0/2304d0 - (pi**2*ll2)/384d0 + (7d0*llzp)
     &/192d0 + (5d0*zeta3)/512d0) + zp**2*((pi**2)/48d0 - (pi
     &**2*ll2)/48d0 + (5d0*zeta3)/64d0)

         case(7)            !-1-11-1

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**4)/30d0) - (pi**2*ll2**2)/8d0 
     &+ (ll2**4)/12d0 + 3*cli4pt5 + (11d0*ll2*zeta3)/4d0 + zp
     &**6*(-((pi**2*137d0)/276480d0) + 37d0/9216d0 + (2213d0*
     &ll2)/460800d0 - (pi**2*ll2)/4608d0 + (137d0*ll2**2)/460
     &80d0 + (ll2**3)/2304d0 - (137d0*ll2*llzp)/23040d0 + (ze
     &ta3)/1536d0) + zp**4*(-((pi**2*11d0)/4608d0) + 11d0/768
     &d0 + (131d0*ll2)/4608d0 - (pi**2*ll2)/768d0 + (11d0*ll2
     &**2)/768d0 + (ll2**3)/384d0 - (11d0*ll2*llzp)/384d0 + (
     &zeta3)/256d0) + zp**2*(-((pi**2)/96d0) + (3d0*ll2)/16d0
     & - (pi**2*ll2)/96d0 + (ll2**2)/16d0 + (ll2**3)/48d0 - (
     &ll2*llzp)/8d0 + (zeta3)/32d0) + zp**5*(181d0/23040d0 - 
     &(pi**2*5d0)/4608d0 - (pi**2*ll2)/1920d0 + (53d0*ll2)/46
     &08d0 + (5d0*ll2**2)/768d0 + (ll2**3)/960d0 - (5d0*ll2*l
     &lzp)/384d0 + (zeta3)/640d0) + zp*(-((pi**2*ll2)/24d0) +
     & (ll2**3)/12d0 + (zeta3)/8d0) + zp**3*(-((pi**2)/192d0)
     & + 1d0/48d0 - (pi**2*ll2)/288d0 + (7d0*ll2)/96d0 + (ll2
     &**2)/32d0 + (ll2**3)/144d0 - (ll2*llzp)/16d0 + (zeta3)/
     &96d0)

         case(8)            !-1-110

            zp = 1d0-x

            ris = -((pi**4*7d0)/288d0) - (pi**2*ll2**2)
     &/24d0 + (ll2**4)/12d0 + 2*cli4pt5 + (13d0*ll2*zeta3)/8d
     &0 + zp**5*(107d0/5760d0 - (pi**2*5d0)/2304d0 - (pi**2*l
     &l2)/1920d0 + (zeta3)/160d0) + zp**3*(1d0/24d0 - (pi**2)
     &/96d0 - (pi**2*ll2)/288d0 + (zeta3)/24d0) + zp*(-((pi**
     &2*ll2)/24d0) + (zeta3)/2d0) + zp**6*(-((pi**2*137d0)/13
     &8240d0) + 79d0/7680d0 - (pi**2*ll2)/4608d0 + (zeta3)/38
     &4d0) + zp**4*(-((pi**2*11d0)/2304d0) + 1d0/32d0 - (pi**
     &2*ll2)/768d0 + (zeta3)/64d0) + zp**2*(-((pi**2)/48d0) -
     & (pi**2*ll2)/96d0 + (zeta3)/8d0)

         case(9)            !-1-111

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**4)/720d0 + (ll2**4)/24d0 - (ll2*
     &zeta3)/8d0 + zp**5*(2671d0/276480d0 + (pi**2*ll2)/1920d
     &0 - (ll2**3)/960d0 - (53d0*llzp)/4608d0 + (5d0*llzp**2)
     &/768d0 - (7d0*zeta3)/1280d0) + zp*((pi**2*ll2)/24d0 - (
     &ll2**3)/12d0 - (7d0*zeta3)/16d0) + zp**3*(41d0/576d0 + 
     &(pi**2*ll2)/288d0 - (ll2**3)/144d0 - (7d0*llzp)/96d0 + 
     &(llzp**2)/32d0 - (7d0*zeta3)/192d0) + zp**6*(322493d0/8
     &2944000d0 + (pi**2*ll2)/4608d0 - (ll2**3)/2304d0 - (221
     &3d0*llzp)/460800d0 + (137d0*llzp**2)/46080d0 - (7d0*zet
     &a3)/3072d0) + zp**4*(1397d0/55296d0 + (pi**2*ll2)/768d0
     & - (ll2**3)/384d0 - (131d0*llzp)/4608d0 + (11d0*llzp**2
     &)/768d0 - (7d0*zeta3)/512d0) + zp**2*(7d0/32d0 + (pi**2
     &*ll2)/96d0 - (ll2**3)/48d0 - (3d0*llzp)/16d0 + (llzp**2
     &)/16d0 - (7d0*zeta3)/64d0)

         case(10)            !-10-1-1

            zp = 1d0-x

            ris = -((pi**4)/30d0) - (pi**2*ll2**2)/8d0 
     &+ (ll2**4)/8d0 + 3*cli4pt5 - (zp*zeta3)/16d0 + (11d0*ll
     &2*zeta3)/4d0 + zp**5*(13d0/1920d0 - (ll2)/30d0 + (ll2**
     &2)/30d0 - (zeta3)/1280d0) + zp**3*(-((ll2)/24d0) + (ll2
     &**2)/12d0 - (zeta3)/192d0) + zp**6*(37d0/5760d0 - (97d0
     &*ll2)/3840d0 + (ll2**2)/45d0 - (zeta3)/3072d0) + zp**4*
     &(1d0/192d0 - (ll2)/24d0 + (5d0*ll2**2)/96d0 - (zeta3)/5
     &12d0) + zp**2*((ll2**2)/8d0 - (zeta3)/64d0)

         case(11)            !-10-10

            zp = 1d0-x

            ris = -((pi**4*11d0)/288d0) - (pi**2*ll2**2
     &)/6d0 + (ll2**4)/6d0 + 4*cli4pt5 + (3d0*zp*zeta3)/4d0 +
     & 2*ll2*zeta3 + zp**3*(-((pi**2)/72d0) + (zeta3)/16d0) +
     & zp**6*(17d0/1152d0 - (pi**2)/270d0 + (zeta3)/256d0) + 
     &zp**4*(-((pi**2*5d0)/576d0) + 1d0/96d0 + (3d0*zeta3)/12
     &8d0) + zp**2*(-((pi**2)/48d0) + (3d0*zeta3)/16d0) + zp*
     &*5*(-((pi**2)/180d0) + 7d0/480d0 + (3d0*zeta3)/320d0)

         case(12)            !-10-11

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**4*5d0)/144d0 - (pi**2*ll2**2)/12
     &d0 - (ll2**4)/6d0 - 4*cli4pt5 - ll2*zeta3 + zp**5*((pi*
     &*2)/180d0 - 1367d0/28800d0 + (pi**2*ll2)/640d0 - (ll2**
     &2)/30d0 + (llzp)/30d0 - (13d0*zeta3)/1280d0) + zp*((pi*
     &*2*ll2)/8d0 - (13d0*zeta3)/16d0) + zp**3*(-(11d0/144d0)
     & + (pi**2)/72d0 + (pi**2*ll2)/96d0 - (ll2**2)/12d0 + (l
     &lzp)/24d0 - (13d0*zeta3)/192d0) + zp**6*((pi**2)/270d0 
     &- 7639d0/230400d0 + (pi**2*ll2)/1536d0 - (ll2**2)/45d0 
     &+ (97d0*llzp)/3840d0 - (13d0*zeta3)/3072d0) + zp**4*(-(
     &19d0/288d0) + (pi**2*5d0)/576d0 + (pi**2*ll2)/256d0 - (
     &5d0*ll2**2)/96d0 + (llzp)/24d0 - (13d0*zeta3)/512d0) + 
     &zp**2*((pi**2)/48d0 + (pi**2*ll2)/32d0 - (ll2**2)/8d0 -
     & (13d0*zeta3)/64d0)

         case(13)            !-100-1

            zp = 1d0-x

            ris = -((pi**4)/288d0) - (3d0*zp*zeta3)/8d0
     & + (3d0*ll2*zeta3)/4d0 + zp**3*((pi**2)/72d0 - (ll2)/12
     &d0 - (zeta3)/32d0) + zp**4*((pi**2*5d0)/576d0 + 1d0/96d
     &0 - (3d0*ll2)/32d0 - (3d0*zeta3)/256d0) + zp**2*((pi**2
     &)/48d0 - (3d0*zeta3)/32d0) + zp**6*((pi**2)/270d0 + 13d
     &0/768d0 - (5d0*ll2)/72d0 - (zeta3)/512d0) + zp**5*((pi*
     &*2)/180d0 + 1d0/64d0 - (ll2)/12d0 - (3d0*zeta3)/640d0)

         case(14)            !-1000

            zp = 1d0-x

            ris = -((pi**4*7d0)/720d0) + (zp**4)/48d0 +
     & (zp**5)/30d0 + (11d0*zp**6)/288d0

         case(15)            !-1001

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**4)/60d0 + (pi**2*ll2**2)/12d0 - 
     &(ll2**4)/12d0 - 2*cli4pt5 - (zp*zeta3)/2d0 - (3d0*ll2*z
     &eta3)/4d0 + zp**5*(-(317d0/2880d0) + (pi**2)/90d0 + (ll
     &zp)/12d0 - (zeta3)/160d0) + zp**3*((pi**2)/36d0 - 11d0/
     &72d0 + (llzp)/12d0 - (zeta3)/24d0) + zp**6*((pi**2)/135
     &d0 - 187d0/2304d0 + (5d0*llzp)/72d0 - (zeta3)/384d0) + 
     &zp**4*(-(55d0/384d0) + (pi**2*5d0)/288d0 + (3d0*llzp)/3
     &2d0 - (zeta3)/64d0) + zp**2*((pi**2)/24d0 - (zeta3)/8d0
     &)

         case(16)            !-101-1

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**4*7d0)/180d0) + (pi**2*ll2**2)
     &/12d0 + (ll2**4)/6d0 + 4*cli4pt5 + (13d0*ll2*zeta3)/8d0
     & + zp**5*(-((pi**2)/180d0) + 173d0/5760d0 - (pi**2*ll2)
     &/640d0 + (79d0*ll2)/1800d0 + (ll2**2)/30d0 - (ll2*llzp)
     &/15d0 + (zeta3)/160d0) + zp**3*(1d0/24d0 - (pi**2)/72d0
     & + (13d0*ll2)/72d0 - (pi**2*ll2)/96d0 + (ll2**2)/12d0 -
     & (ll2*llzp)/6d0 + (zeta3)/24d0) + zp*(-((pi**2*ll2)/8d0
     &) + (zeta3)/2d0) + zp**6*(-((pi**2)/270d0) + 3067d0/138
     &240d0 - (pi**2*ll2)/1536d0 + (169d0*ll2)/7200d0 + (ll2*
     &*2)/45d0 - (2d0*ll2*llzp)/45d0 + (zeta3)/384d0) + zp**4
     &*(5d0/128d0 - (pi**2*5d0)/576d0 - (pi**2*ll2)/256d0 + (
     &25d0*ll2)/288d0 + (5d0*ll2**2)/96d0 - (5d0*ll2*llzp)/48
     &d0 + (zeta3)/64d0) + zp**2*(-((pi**2)/48d0) - (pi**2*ll
     &2)/32d0 + (3d0*ll2)/8d0 + (ll2**2)/8d0 - (ll2*llzp)/4d0
     & + (zeta3)/8d0)

         case(17)            !-1010

            zp = 1d0-x

            ris = -((pi**4*17d0)/480d0) - (pi**2*ll2**2
     &)/6d0 + (ll2**4)/6d0 + 4*cli4pt5 + zp*zeta3 + (3d0*ll2*
     &zeta3)/2d0 + zp**3*(1d0/12d0 - (pi**2)/36d0 + (zeta3)/1
     &2d0) + zp**6*(-((pi**2)/135d0) + 179d0/3456d0 + (zeta3)
     &/192d0) + zp**4*(1d0/12d0 - (pi**2*5d0)/288d0 + (zeta3)
     &/32d0) + zp**2*(-((pi**2)/24d0) + (zeta3)/4d0) + zp**5*
     &(-((pi**2)/90d0) + 97d0/1440d0 + (zeta3)/80d0)

         case(18)            !-1011

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**4)/288d0 + (pi**2*ll2**2)/24d0 -
     & (ll2**4)/24d0 - cli4pt5 - (zp*zeta3)/2d0 + (ll2*zeta3)
     &/8d0 + zp**5*(12017d0/432000d0 - (79d0*llzp)/1800d0 + (
     &llzp**2)/30d0 - (zeta3)/160d0) + zp**3*(71d0/432d0 - (1
     &3d0*llzp)/72d0 + (llzp**2)/12d0 - (zeta3)/24d0) + zp**6
     &*(64861d0/5184000d0 - (169d0*llzp)/7200d0 + (llzp**2)/4
     &5d0 - (zeta3)/384d0) + zp**4*(113d0/1728d0 - (25d0*llzp
     &)/288d0 + (5d0*llzp**2)/96d0 - (zeta3)/64d0) + zp**2*(7
     &d0/16d0 - (3d0*llzp)/8d0 + (llzp**2)/8d0 - (zeta3)/8d0)

         case(19)            !-11-1-1

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**4)/30d0 + (pi**2*ll2**2)/6d0 - (
     &ll2**4)/6d0 - 3*cli4pt5 - (23d0*ll2*zeta3)/8d0 + zp**5*
     &(17d0/5120d0 + (pi**2*ll2)/1920d0 - (41d0*ll2)/4608d0 -
     & (ll2**2)/1600d0 - (ll2**3)/480d0 + (ll2**2*llzp)/320d0
     & - (zeta3)/1280d0) + zp*((pi**2*ll2)/24d0 - (ll2**2)/4d
     &0 - (ll2**3)/6d0 + (ll2**2*llzp)/4d0 - (zeta3)/16d0) + 
     &zp**3*(1d0/96d0 + (pi**2*ll2)/288d0 - (5d0*ll2)/96d0 - 
     &(ll2**2)/144d0 - (ll2**3)/72d0 + (ll2**2*llzp)/48d0 - (
     &zeta3)/192d0) + zp**6*(59d0/36864d0 + (pi**2*ll2)/4608d
     &0 - (5269d0*ll2)/1382400d0 - (ll2**2)/4608d0 - (ll2**3)
     &/1152d0 + (ll2**2*llzp)/768d0 - (zeta3)/3072d0) + zp**4
     &*(5d0/768d0 - (49d0*ll2)/2304d0 + (pi**2*ll2)/768d0 - (
     &ll2**2)/512d0 - (ll2**3)/192d0 + (ll2**2*llzp)/128d0 - 
     &(zeta3)/512d0) + zp**2*(-((ll2)/8d0) + (pi**2*ll2)/96d0
     & - (ll2**2)/32d0 - (ll2**3)/24d0 + (ll2**2*llzp)/16d0 -
     & (zeta3)/64d0)

         case(20)            !-11-10

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**4)/360d0 + (pi**2*ll2**2)/24d0 -
     & ll2*zeta3 + zp**5*(49d0/5760d0 + (pi**2)/9600d0 + (pi*
     &*2*ll2)/640d0 - (pi**2*llzp)/1920d0 - (13d0*zeta3)/1280
     &d0) + zp*((pi**2)/24d0 + (pi**2*ll2)/8d0 - (pi**2*llzp)
     &/24d0 - (13d0*zeta3)/16d0) + zp**3*(1d0/48d0 + (pi**2)/
     &864d0 + (pi**2*ll2)/96d0 - (pi**2*llzp)/288d0 - (13d0*z
     &eta3)/192d0) + zp**6*((pi**2)/27648d0 + 1609d0/345600d0
     & + (pi**2*ll2)/1536d0 - (pi**2*llzp)/4608d0 - (13d0*zet
     &a3)/3072d0) + zp**4*(17d0/1152d0 + (pi**2)/3072d0 + (pi
     &**2*ll2)/256d0 - (pi**2*llzp)/768d0 - (13d0*zeta3)/512d
     &0) + zp**2*((pi**2)/192d0 + (pi**2*ll2)/32d0 - (pi**2*l
     &lzp)/96d0 - (13d0*zeta3)/64d0)

         case(21)            !-11-11

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**4)/1440d0 - (pi**2*ll2**2)/24d0 
     &+ (ll2**4)/24d0 + (ll2*zeta3)/4d0 + zp**6*(-(17653d0/25
     &92000d0) - (pi**2)/27648d0 - (pi**2*ll2)/2304d0 + (ll2*
     &*2)/4608d0 + (ll2**3)/1152d0 + (pi**2*llzp)/4608d0 + (5
     &269d0*llzp)/1382400d0 - (ll2**2*llzp)/768d0 + (7d0*zeta
     &3)/1536d0) + zp**4*(-(1151d0/27648d0) - (pi**2)/3072d0 
     &- (pi**2*ll2)/384d0 + (ll2**2)/512d0 + (ll2**3)/192d0 +
     & (49d0*llzp)/2304d0 + (pi**2*llzp)/768d0 - (ll2**2*llzp
     &)/128d0 + (7d0*zeta3)/256d0) + zp**2*(-((pi**2)/192d0) 
     &- 5d0/16d0 - (pi**2*ll2)/48d0 + (ll2**2)/32d0 + (ll2**3
     &)/24d0 + (llzp)/8d0 + (pi**2*llzp)/96d0 - (ll2**2*llzp)
     &/16d0 + (7d0*zeta3)/32d0) + zp**5*(-(2281d0/138240d0) -
     & (pi**2)/9600d0 - (pi**2*ll2)/960d0 + (ll2**2)/1600d0 +
     & (ll2**3)/480d0 + (pi**2*llzp)/1920d0 + (41d0*llzp)/460
     &8d0 - (ll2**2*llzp)/320d0 + (7d0*zeta3)/640d0) + zp*(-(
     &(pi**2)/24d0) - (pi**2*ll2)/12d0 + (ll2**2)/4d0 + (ll2*
     &*3)/6d0 + (pi**2*llzp)/24d0 - (ll2**2*llzp)/4d0 + (7d0*
     &zeta3)/8d0) + zp**3*(-((pi**2)/864d0) - 1d0/9d0 - (pi**
     &2*ll2)/144d0 + (ll2**2)/144d0 + (ll2**3)/72d0 + (pi**2*
     &llzp)/288d0 + (5d0*llzp)/96d0 - (ll2**2*llzp)/48d0 + (7
     &d0*zeta3)/96d0)

         case(22)            !-110-1

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**4*11d0)/240d0 + (pi**2*ll2**2)/8
     &d0 - (ll2**4)/6d0 - 4*cli4pt5 - (13d0*ll2*zeta3)/4d0 + 
     &zp**5*(31d0/2880d0 - (pi**2)/9600d0 - (11d0*ll2)/360d0 
     &+ (pi**2*llzp)/1920d0 + (zeta3)/256d0) + zp*(-((pi**2)/
     &24d0) + (pi**2*llzp)/24d0 + (5d0*zeta3)/16d0) + zp**3*(
     &1d0/48d0 - (pi**2)/864d0 - (ll2)/8d0 + (pi**2*llzp)/288
     &d0 + (5d0*zeta3)/192d0) + zp**6*(-((pi**2)/27648d0) + 7
     &3d0/10800d0 - (347d0*ll2)/21600d0 + (pi**2*llzp)/4608d0
     & + (5d0*zeta3)/3072d0) + zp**4*(19d0/1152d0 - (pi**2)/3
     &072d0 - (35d0*ll2)/576d0 + (pi**2*llzp)/768d0 + (5d0*ze
     &ta3)/512d0) + zp**2*(-((pi**2)/192d0) - (ll2)/4d0 + (pi
     &**2*llzp)/96d0 + (5d0*zeta3)/64d0)

         case(23)            !-1100

            zp = 1d0-x

            ris = (pi**4*19d0)/1440d0 - (zp*zeta3)/2d0 
     &- (zp**2*zeta3)/8d0 - (3d0*ll2*zeta3)/4d0 + zp**5*(5d0/
     &192d0 - (zeta3)/160d0) + zp**3*(1d0/24d0 - (zeta3)/24d0
     &) + zp**6*(41d0/2304d0 - (zeta3)/384d0) + zp**4*(7d0/19
     &2d0 - (zeta3)/64d0)

         case(24)            !-1101

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**4*19d0)/1440d0 - (pi**2*ll2**2)/
     &24d0 - (ll2**4)/24d0 - cli4pt5 - (ll2*zeta3)/4d0 + zp*(
     &-((pi**2)/12d0) + (pi**2*llzp)/12d0 + zeta3) + zp**3*(-
     &((pi**2)/432d0) - 1d0/4d0 + (pi**2*llzp)/144d0 + (llzp)
     &/8d0 + (zeta3)/12d0) + zp**6*(-((pi**2)/13824d0) - 5152
     &1d0/2592000d0 + (pi**2*llzp)/2304d0 + (347d0*llzp)/2160
     &0d0 + (zeta3)/192d0) + zp**4*(-((pi**2)/1536d0) - 709d0
     &/6912d0 + (pi**2*llzp)/384d0 + (35d0*llzp)/576d0 + (zet
     &a3)/32d0) + zp**2*(-(5d0/8d0) - (pi**2)/96d0 + (pi**2*l
     &lzp)/48d0 + (llzp)/4d0 + (zeta3)/4d0) + zp**5*(-(1909d0
     &/43200d0) - (pi**2)/4800d0 + (11d0*llzp)/360d0 + (pi**2
     &*llzp)/960d0 + (zeta3)/80d0)

         case(25)            !-111-1

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**4)/288d0) - (pi**2*ll2**2)/24d
     &0 + (ll2**4)/24d0 + (7d0*ll2*zeta3)/8d0 + zp**5*(407d0/
     &55296d0 + (pi**2)/9600d0 + (pi**2*ll2)/1920d0 - (ll2)/4
     &000d0 - (ll2**2)/1600d0 - (ll2**3)/960d0 - (pi**2*llzp)
     &/1920d0 + (ll2*llzp)/800d0 + (ll2**2*llzp)/320d0 - (ll2
     &*llzp**2)/320d0 - (7d0*zeta3)/1280d0) + zp*((pi**2)/24d
     &0 + (pi**2*ll2)/24d0 - (ll2)/2d0 - (ll2**2)/4d0 - (ll2*
     &*3)/12d0 - (pi**2*llzp)/24d0 + (ll2*llzp)/2d0 + (ll2**2
     &*llzp)/4d0 - (ll2*llzp**2)/4d0 - (7d0*zeta3)/16d0) + zp
     &**3*(3d0/64d0 + (pi**2)/864d0 - (ll2)/216d0 + (pi**2*ll
     &2)/288d0 - (ll2**2)/144d0 - (ll2**3)/144d0 - (pi**2*llz
     &p)/288d0 + (ll2*llzp)/72d0 + (ll2**2*llzp)/48d0 - (ll2*
     &llzp**2)/48d0 - (7d0*zeta3)/192d0) + zp**6*((pi**2)/276
     &48d0 + 256103d0/82944000d0 - (ll2)/13824d0 + (pi**2*ll2
     &)/4608d0 - (ll2**2)/4608d0 - (ll2**3)/2304d0 - (pi**2*l
     &lzp)/4608d0 + (ll2*llzp)/2304d0 + (ll2**2*llzp)/768d0 -
     & (ll2*llzp**2)/768d0 - (7d0*zeta3)/3072d0) + zp**4*(251
     &d0/13824d0 + (pi**2)/3072d0 - (ll2)/1024d0 + (pi**2*ll2
     &)/768d0 - (ll2**2)/512d0 - (ll2**3)/384d0 - (pi**2*llzp
     &)/768d0 + (ll2*llzp)/256d0 + (ll2**2*llzp)/128d0 - (ll2
     &*llzp**2)/128d0 - (7d0*zeta3)/512d0) + zp**2*((pi**2)/1
     &92d0 + 1d0/8d0 - (ll2)/32d0 + (pi**2*ll2)/96d0 - (ll2**
     &2)/32d0 - (ll2**3)/48d0 - (pi**2*llzp)/96d0 + (ll2*llzp
     &)/16d0 + (ll2**2*llzp)/16d0 - (ll2*llzp**2)/16d0 - (7d0
     &*zeta3)/64d0)

         case(26)            !-1110

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**4)/720d0) - (ll2**4)/24d0 - cl
     &i4pt5 + (ll2*zeta3)/8d0 + zp**5*(329d0/17280d0 + (pi**2
     &)/4800d0 - (pi**2*llzp)/960d0 - (zeta3)/160d0) + zp**3*
     &((pi**2)/432d0 + 5d0/48d0 - (pi**2*llzp)/144d0 - (zeta3
     &)/24d0) + zp*((pi**2)/12d0 - (pi**2*llzp)/12d0 - (zeta3
     &)/2d0) + zp**6*((pi**2)/13824d0 + 44581d0/5184000d0 - (
     &pi**2*llzp)/2304d0 - (zeta3)/384d0) + zp**4*((pi**2)/15
     &36d0 + 151d0/3456d0 - (pi**2*llzp)/384d0 - (zeta3)/64d0
     &) + zp**2*(1d0/4d0 + (pi**2)/96d0 - (pi**2*llzp)/48d0 -
     & (zeta3)/8d0)

         case(27)            !-1111

            zp = 1d0-x
            llzp = log(zp)

            ris = zp*(-(1d0/2d0) + (llzp)/2d0 - (llzp**
     &2)/4d0 + (llzp**3)/12d0) + zp**3*(-(1d0/648d0) + (llzp)
     &/216d0 - (llzp**2)/144d0 + (llzp**3)/144d0) + zp**6*(-(
     &1d0/82944d0) + (llzp)/13824d0 - (llzp**2)/4608d0 + (llz
     &p**3)/2304d0) + zp**4*(-(1d0/4096d0) + (llzp)/1024d0 - 
     &(llzp**2)/512d0 + (llzp**3)/384d0) + zp**2*(-(1d0/64d0)
     & + (llzp)/32d0 - (llzp**2)/32d0 + (llzp**3)/48d0) + zp*
     &*5*(-(1d0/20000d0) + (llzp)/4000d0 - (llzp**2)/1600d0 +
     & (llzp**3)/960d0) + cli4pt5

         case(28)            !0-1-1-1

            zp = 1d0-x

            ris = (pi**4)/90d0 + (pi**2*ll2**2)/24d0 - 
     &(zp*ll2**3)/6d0 - (ll2**4)/24d0 + zp**2*((ll2**2)/8d0 -
     & (ll2**3)/12d0) + zp**3*(-((ll2)/24d0) + (5d0*ll2**2)/4
     &8d0 - (ll2**3)/18d0) + zp**4*(1d0/192d0 - (3d0*ll2)/64d
     &0 + (ll2**2)/12d0 - (ll2**3)/24d0) + zp**5*(7d0/960d0 -
     & (83d0*ll2)/1920d0 + (131d0*ll2**2)/1920d0 - (ll2**3)/3
     &0d0) + zp**6*(35d0/4608d0 - (11d0*ll2)/288d0 + (661d0*l
     &l2**2)/11520d0 - (ll2**3)/36d0) - cli4pt5 - (7d0*ll2*ze
     &ta3)/8d0

         case(29)            !0-1-10

            zp = 1d0-x

            ris = -((pi**4)/288d0) + zp**2*(-((pi**2)/4
     &8d0) + (pi**2*ll2)/24d0 - (zeta3)/16d0) + zp**3*(-((pi*
     &*2*5d0)/288d0) + (pi**2*ll2)/36d0 - (zeta3)/24d0) + zp*
     &*4*(-((pi**2)/72d0) + 1d0/96d0 + (pi**2*ll2)/48d0 - (ze
     &ta3)/32d0) + zp**5*(-((pi**2*131d0)/11520d0) + 1d0/64d0
     & + (pi**2*ll2)/60d0 - (zeta3)/40d0) + zp**6*(11d0/640d0
     & - (pi**2*661d0)/69120d0 + (pi**2*ll2)/72d0 - (zeta3)/4
     &8d0) + zp*((pi**2*ll2)/12d0 - (zeta3)/8d0)

         case(30)            !0-1-11

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**4)/80d0) + (pi**2*ll2**2)/24d0
     & + (ll2**4)/12d0 + 2*cli4pt5 + zp**2*((pi**2)/48d0 - (l
     &l2**2)/8d0 + (ll2**3)/12d0 - (zeta3)/16d0) + zp**3*(-(1
     &1d0/144d0) + (pi**2*5d0)/288d0 - (5d0*ll2**2)/48d0 + (l
     &l2**3)/18d0 + (llzp)/24d0 - (zeta3)/24d0) + zp**4*((pi*
     &*2)/72d0 - 59d0/768d0 - (ll2**2)/12d0 + (ll2**3)/24d0 +
     & (3d0*llzp)/64d0 - (zeta3)/32d0) + zp**5*((pi**2*131d0)
     &/11520d0 - 7651d0/115200d0 - (131d0*ll2**2)/1920d0 + (l
     &l2**3)/30d0 + (83d0*llzp)/1920d0 - (zeta3)/40d0) + zp**
     &6*(-(65d0/1152d0) + (pi**2*661d0)/69120d0 - (661d0*ll2*
     &*2)/11520d0 + (ll2**3)/36d0 + (11d0*llzp)/288d0 - (zeta
     &3)/48d0) + zp*((ll2**3)/6d0 - (zeta3)/8d0)

         case(31)            !0-10-1

            zp = 1d0-x

            ris = (pi**4*13d0)/288d0 + (pi**2*ll2**2)/6
     &d0 - (ll2**4)/6d0 - 4*cli4pt5 - (7d0*ll2*zeta3)/2d0 + z
     &p**3*((pi**2*5d0)/288d0 - (ll2)/12d0 - (pi**2*ll2)/36d0
     & + (zeta3)/12d0) + zp**4*((pi**2)/72d0 + 1d0/96d0 - (pi
     &**2*ll2)/48d0 - (5d0*ll2)/48d0 + (zeta3)/16d0) + zp**5*
     &((pi**2*131d0)/11520d0 + 1d0/60d0 - (5d0*ll2)/48d0 - (p
     &i**2*ll2)/60d0 + (zeta3)/20d0) + zp**6*((pi**2*661d0)/6
     &9120d0 + 7d0/360d0 - (47d0*ll2)/480d0 - (pi**2*ll2)/72d
     &0 + (zeta3)/24d0) + zp*(-((pi**2*ll2)/12d0) + (zeta3)/4
     &d0) + zp**2*((pi**2)/48d0 - (pi**2*ll2)/24d0 + (zeta3)/
     &8d0)

         case(32)            !0-100

            zp = 1d0-x

            ris = (pi**4*7d0)/240d0 - (3d0*zp*zeta3)/4d
     &0 - (3d0*zp**2*zeta3)/8d0 - (zp**3*zeta3)/4d0 + zp**4*(
     &1d0/48d0 - (3d0*zeta3)/16d0) + zp**5*(17d0/480d0 - (3d0
     &*zeta3)/20d0) + zp**6*(25d0/576d0 - (zeta3)/8d0)

         case(33)            !0-101

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**4)/480d0 + zp**2*((pi**2)/24d0 -
     & (pi**2*ll2)/12d0 + (5d0*zeta3)/16d0) + zp**3*((pi**2*5
     &d0)/144d0 - 11d0/72d0 - (pi**2*ll2)/18d0 + (llzp)/12d0 
     &+ (5d0*zeta3)/24d0) + zp**4*((pi**2)/36d0 - 95d0/576d0 
     &- (pi**2*ll2)/24d0 + (5d0*llzp)/48d0 + (5d0*zeta3)/32d0
     &) + zp**6*((pi**2*661d0)/34560d0 - 941d0/7200d0 - (pi**
     &2*ll2)/36d0 + (47d0*llzp)/480d0 + (5d0*zeta3)/48d0) + z
     &p**5*(-(43d0/288d0) + (pi**2*131d0)/5760d0 - (pi**2*ll2
     &)/30d0 + (5d0*llzp)/48d0 + (zeta3)/8d0) + zp*(-((pi**2*
     &ll2)/6d0) + (5d0*zeta3)/8d0)

         case(34)            !0-11-1

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**4*7d0)/720d0) - (pi**2*ll2**2)
     &/4d0 + (21d0*ll2*zeta3)/8d0 + zp**3*(1d0/24d0 - (pi**2*
     &5d0)/288d0 - (pi**2*ll2)/36d0 + (37d0*ll2)/144d0 + (5d0
     &*ll2**2)/48d0 + (ll2**3)/18d0 - (5d0*ll2*llzp)/24d0 + (
     &zeta3)/12d0) + zp**4*(17d0/384d0 - (pi**2)/72d0 - (pi**
     &2*ll2)/48d0 + (107d0*ll2)/576d0 + (ll2**2)/12d0 + (ll2*
     &*3)/24d0 - (ll2*llzp)/6d0 + (zeta3)/16d0) + zp**5*(-((p
     &i**2*131d0)/11520d0) + 457d0/11520d0 - (pi**2*ll2)/60d0
     & + (8257d0*ll2)/57600d0 + (131d0*ll2**2)/1920d0 + (ll2*
     &*3)/30d0 - (131d0*ll2*llzp)/960d0 + (zeta3)/20d0) + zp*
     &*6*(-((pi**2*661d0)/69120d0) + 955d0/27648d0 + (13369d0
     &*ll2)/115200d0 - (pi**2*ll2)/72d0 + (661d0*ll2**2)/1152
     &0d0 + (ll2**3)/36d0 - (661d0*ll2*llzp)/5760d0 + (zeta3)
     &/24d0) + zp*(-((pi**2*ll2)/12d0) + (ll2**3)/6d0 + (zeta
     &3)/4d0) + zp**2*(-((pi**2)/48d0) - (pi**2*ll2)/24d0 + (
     &3d0*ll2)/8d0 + (ll2**2)/8d0 + (ll2**3)/12d0 - (ll2*llzp
     &)/4d0 + (zeta3)/8d0)

         case(35)            !0-110

            zp = 1d0-x

            ris = (pi**4*13d0)/1440d0 + (pi**2*ll2**2)/
     &6d0 - (ll2**4)/6d0 - 4*cli4pt5 + zp*(-((pi**2*ll2)/12d0
     &) + zeta3) + zp**2*(-((pi**2)/24d0) - (pi**2*ll2)/24d0 
     &+ (zeta3)/2d0) + zp**3*(1d0/12d0 - (pi**2*5d0)/144d0 - 
     &(pi**2*ll2)/36d0 + (zeta3)/3d0) + zp**4*(-((pi**2)/36d0
     &) + 3d0/32d0 - (pi**2*ll2)/48d0 + (zeta3)/4d0) + zp**5*
     &(251d0/2880d0 - (pi**2*131d0)/5760d0 - (pi**2*ll2)/60d0
     & + (zeta3)/5d0) + zp**6*(1343d0/17280d0 - (pi**2*661d0)
     &/34560d0 - (pi**2*ll2)/72d0 + (zeta3)/6d0)

         case(36)            !0-111

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**4*11d0)/720d0) + (ll2**4)/8d0 
     &+ 3*cli4pt5 + zp**2*(7d0/16d0 + (pi**2*ll2)/24d0 - (ll2
     &**3)/12d0 - (3d0*llzp)/8d0 + (llzp**2)/8d0 - (7d0*zeta3
     &)/16d0) + zp**3*(227d0/864d0 + (pi**2*ll2)/36d0 - (ll2*
     &*3)/18d0 - (37d0*llzp)/144d0 + (5d0*llzp**2)/48d0 - (7d
     &0*zeta3)/24d0) + zp**4*(1247d0/6912d0 + (pi**2*ll2)/48d
     &0 - (ll2**3)/24d0 - (107d0*llzp)/576d0 + (llzp**2)/12d0
     & - (7d0*zeta3)/32d0) + zp**5*(470159d0/3456000d0 + (pi*
     &*2*ll2)/60d0 - (ll2**3)/30d0 - (8257d0*llzp)/57600d0 + 
     &(131d0*llzp**2)/1920d0 - (7d0*zeta3)/40d0) + zp**6*(225
     &7309d0/20736000d0 + (pi**2*ll2)/72d0 - (ll2**3)/36d0 - 
     &(13369d0*llzp)/115200d0 + (661d0*llzp**2)/11520d0 - (7d
     &0*zeta3)/48d0) + zp*((pi**2*ll2)/12d0 - (ll2**3)/6d0 - 
     &(7d0*zeta3)/8d0)

         case(37)            !00-1-1

            zp = 1d0-x

            ris = -((pi**4)/48d0) - (pi**2*ll2**2)/12d0
     & + (ll2**4)/12d0 + 2*cli4pt5 - (zp*zeta3)/8d0 + (7d0*ll
     &2*zeta3)/4d0 + zp**2*((ll2**2)/4d0 - (zeta3)/16d0) + zp
     &**3*(-((ll2)/12d0) + (ll2**2)/4d0 - (zeta3)/24d0) + zp*
     &*4*(1d0/96d0 - (11d0*ll2)/96d0 + (11d0*ll2**2)/48d0 - (
     &zeta3)/32d0) + zp**5*(17d0/960d0 - (ll2)/8d0 + (5d0*ll2
     &**2)/24d0 - (zeta3)/40d0) + zp**6*(253d0/11520d0 - (731
     &d0*ll2)/5760d0 + (137d0*ll2**2)/720d0 - (zeta3)/48d0)

         case(38)            !00-10

            zp = 1d0-x

            ris = -((pi**4*7d0)/240d0) + (3d0*zp*zeta3)
     &/2d0 + zp**3*(-((pi**2)/24d0) + (zeta3)/2d0) + zp**5*(-
     &((pi**2*5d0)/144d0) + 3d0/80d0 + (3d0*zeta3)/10d0) + zp
     &**6*(-((pi**2*137d0)/4320d0) + 7d0/144d0 + (zeta3)/4d0)
     & + zp**2*(-((pi**2)/24d0) + (3d0*zeta3)/4d0) + zp**4*(-
     &((pi**2*11d0)/288d0) + 1d0/48d0 + (3d0*zeta3)/8d0)

         case(39)            !00-11

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**4)/180d0) - (pi**2*ll2**2)/12d
     &0 + (ll2**4)/12d0 + 2*cli4pt5 + zp**2*((pi**2)/24d0 + (
     &pi**2*ll2)/8d0 - (ll2**2)/4d0 - (13d0*zeta3)/16d0) + zp
     &**3*((pi**2)/24d0 - 11d0/72d0 + (pi**2*ll2)/12d0 - (ll2
     &**2)/4d0 + (llzp)/12d0 - (13d0*zeta3)/24d0) + zp**4*(-(
     &215d0/1152d0) + (pi**2*11d0)/288d0 + (pi**2*ll2)/16d0 -
     & (11d0*ll2**2)/48d0 + (11d0*llzp)/96d0 - (13d0*zeta3)/3
     &2d0) + zp**5*((pi**2*5d0)/144d0 - 181d0/960d0 + (pi**2*
     &ll2)/20d0 - (5d0*ll2**2)/24d0 + (llzp)/8d0 - (13d0*zeta
     &3)/40d0) + zp**6*(-(2321d0/12800d0) + (pi**2*137d0)/432
     &0d0 + (pi**2*ll2)/24d0 - (137d0*ll2**2)/720d0 + (731d0*
     &llzp)/5760d0 - (13d0*zeta3)/48d0) + zp*((pi**2*ll2)/4d0
     & - (13d0*zeta3)/8d0)

         case(40)            !000-1

            zp = 1d0-x

            ris = (pi**4*7d0)/720d0 - (3d0*zp*zeta3)/4d
     &0 + zp**4*((pi**2*11d0)/288d0 + 1d0/48d0 - (ll2)/4d0 - 
     &(3d0*zeta3)/16d0) + zp**5*(19d0/480d0 + (pi**2*5d0)/144
     &d0 - (7d0*ll2)/24d0 - (3d0*zeta3)/20d0) + zp**3*((pi**2
     &)/24d0 - (ll2)/6d0 - (zeta3)/4d0) + zp**6*((pi**2*137d0
     &)/4320d0 + 31d0/576d0 - (5d0*ll2)/16d0 - (zeta3)/8d0) +
     & zp**2*((pi**2)/24d0 - (3d0*zeta3)/8d0)

         case(41)            !0000

            zp = 1d0-x

            ris = (zp**4)/24d0 + (zp**5)/12d0 + (17d0*z
     &p**6)/144d0

         case(42)            !0001

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**4)/90d0 - zp*zeta3 + zp**2*((pi*
     &*2)/12d0 - (zeta3)/2d0) + zp**3*((pi**2)/12d0 - 11d0/36
     &d0 + (llzp)/6d0 - (zeta3)/3d0) + zp**4*((pi**2*11d0)/14
     &4d0 - 19d0/48d0 + (llzp)/4d0 - (zeta3)/4d0) + zp**5*(-(
     &599d0/1440d0) + (pi**2*5d0)/72d0 + (7d0*llzp)/24d0 - (z
     &eta3)/5d0) + zp**6*((pi**2*137d0)/2160d0 - 79d0/192d0 +
     & (5d0*llzp)/16d0 - (zeta3)/6d0)

         case(43)            !001-1

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**4*19d0)/1440d0) + (7d0*ll2*zet
     &a3)/4d0 + zp*(-((pi**2*ll2)/4d0) + zeta3) + zp**2*(-((p
     &i**2)/24d0) + (3d0*ll2)/4d0 - (pi**2*ll2)/8d0 + (ll2**2
     &)/4d0 - (ll2*llzp)/2d0 + (zeta3)/2d0) + zp**3*(1d0/12d0
     & - (pi**2)/24d0 - (pi**2*ll2)/12d0 + (7d0*ll2)/12d0 + (
     &ll2**2)/4d0 - (ll2*llzp)/2d0 + (zeta3)/3d0) + zp**4*(-(
     &(pi**2*11d0)/288d0) + 7d0/64d0 - (pi**2*ll2)/16d0 + (13
     &1d0*ll2)/288d0 + (11d0*ll2**2)/48d0 - (11d0*ll2*llzp)/2
     &4d0 + (zeta3)/4d0) + zp**5*(-((pi**2*5d0)/144d0) + 67d0
     &/576d0 - (pi**2*ll2)/20d0 + (53d0*ll2)/144d0 + (5d0*ll2
     &**2)/24d0 - (5d0*ll2*llzp)/12d0 + (zeta3)/5d0) + zp**6*
     &(-((pi**2*137d0)/4320d0) + 893d0/7680d0 - (pi**2*ll2)/2
     &4d0 + (2213d0*ll2)/7200d0 + (137d0*ll2**2)/720d0 - (137
     &d0*ll2*llzp)/360d0 + (zeta3)/6d0)

         case(44)            !0010

            zp = 1d0-x

            ris = -((pi**4)/30d0) + 2*zp*zeta3 + zp**2*
     &(-((pi**2)/12d0) + zeta3) + zp**4*(-((pi**2*11d0)/144d0
     &) + 11d0/48d0 + (zeta3)/2d0) + zp**6*(-((pi**2*137d0)/2
     &160d0) + 37d0/144d0 + (zeta3)/3d0) + zp**3*(-((pi**2)/1
     &2d0) + 1d0/6d0 + (2d0*zeta3)/3d0) + zp**5*(181d0/720d0 
     &- (pi**2*5d0)/72d0 + (2d0*zeta3)/5d0)

         case(45)            !0011

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**4)/360d0 - zp*zeta3 + zp**2*(7d0
     &/8d0 - (3d0*llzp)/4d0 + (llzp**2)/4d0 - (zeta3)/2d0) + 
     &zp**3*(41d0/72d0 - (7d0*llzp)/12d0 + (llzp**2)/4d0 - (z
     &eta3)/3d0) + zp**4*(1397d0/3456d0 - (131d0*llzp)/288d0 
     &+ (11d0*llzp**2)/48d0 - (zeta3)/4d0) + zp**5*(2671d0/86
     &40d0 - (53d0*llzp)/144d0 + (5d0*llzp**2)/24d0 - (zeta3)
     &/5d0) + zp**6*(322493d0/1296000d0 - (2213d0*llzp)/7200d
     &0 + (137d0*llzp**2)/720d0 - (zeta3)/6d0)

         case(46)            !01-1-1

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**4*7d0)/288d0 + (pi**2*5d0*ll2**2
     &)/24d0 - (ll2**4)/12d0 - 2*cli4pt5 - (21d0*ll2*zeta3)/8
     &d0 + zp**2*((pi**2*ll2)/24d0 - (ll2)/4d0 - (ll2**2)/8d0
     & - (ll2**3)/6d0 + (ll2**2*llzp)/4d0 - (zeta3)/16d0) + z
     &p**3*(1d0/48d0 + (pi**2*ll2)/36d0 - (3d0*ll2)/16d0 - (l
     &l2**2)/18d0 - (ll2**3)/9d0 + (ll2**2*llzp)/6d0 - (zeta3
     &)/24d0) + zp**4*(1d0/48d0 + (pi**2*ll2)/48d0 - (83d0*ll
     &2)/576d0 - (ll2**2)/32d0 - (ll2**3)/12d0 + (ll2**2*llzp
     &)/8d0 - (zeta3)/32d0) + zp**5*(139d0/7680d0 - (1337d0*l
     &l2)/11520d0 + (pi**2*ll2)/60d0 - (ll2**2)/50d0 - (ll2**
     &3)/15d0 + (ll2**2*llzp)/10d0 - (zeta3)/40d0) + zp**6*(1
     &43d0/9216d0 - (33497d0*ll2)/345600d0 + (pi**2*ll2)/72d0
     & - (ll2**2)/72d0 - (ll2**3)/18d0 + (ll2**2*llzp)/12d0 -
     & (zeta3)/48d0) + zp*((pi**2*ll2)/12d0 - (ll2**2)/2d0 - 
     &(ll2**3)/3d0 + (ll2**2*llzp)/2d0 - (zeta3)/8d0)

         case(47)            !01-10

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**4*11d0)/480d0) - (pi**2*ll2**2
     &)/6d0 + (ll2**4)/6d0 + 4*cli4pt5 + zp**2*((pi**2)/48d0 
     &+ (pi**2*ll2)/8d0 - (pi**2*llzp)/24d0 - (13d0*zeta3)/16
     &d0) + zp**3*((pi**2)/108d0 + 1d0/24d0 + (pi**2*ll2)/12d
     &0 - (pi**2*llzp)/36d0 - (13d0*zeta3)/24d0) + zp**4*((pi
     &**2)/192d0 + 13d0/288d0 + (pi**2*ll2)/16d0 - (pi**2*llz
     &p)/48d0 - (13d0*zeta3)/32d0) + zp**5*(119d0/2880d0 + (p
     &i**2)/300d0 + (pi**2*ll2)/20d0 - (pi**2*llzp)/60d0 - (1
     &3d0*zeta3)/40d0) + zp**6*((pi**2)/432d0 + 3167d0/86400d
     &0 + (pi**2*ll2)/24d0 - (pi**2*llzp)/72d0 - (13d0*zeta3)
     &/48d0) + zp*((pi**2)/12d0 + (pi**2*ll2)/4d0 - (pi**2*ll
     &zp)/12d0 - (13d0*zeta3)/8d0)

         case(48)            !01-11

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**4*7d0)/288d0 - (pi**2*ll2**2)/8d
     &0 - (ll2**4)/8d0 - 3*cli4pt5 + zp**3*(-((pi**2)/108d0) 
     &- 5d0/12d0 - (pi**2*ll2)/18d0 + (ll2**2)/18d0 + (ll2**3
     &)/9d0 + (pi**2*llzp)/36d0 + (3d0*llzp)/16d0 - (ll2**2*l
     &lzp)/6d0 + (7d0*zeta3)/12d0) + zp**4*(-((pi**2)/192d0) 
     &- 2101d0/6912d0 - (pi**2*ll2)/24d0 + (ll2**2)/32d0 + (l
     &l2**3)/12d0 + (pi**2*llzp)/48d0 + (83d0*llzp)/576d0 - (
     &ll2**2*llzp)/8d0 + (7d0*zeta3)/16d0) + zp**5*(-((pi**2)
     &/300d0) - 82237d0/345600d0 - (pi**2*ll2)/30d0 + (ll2**2
     &)/50d0 + (ll2**3)/15d0 + (1337d0*llzp)/11520d0 + (pi**2
     &*llzp)/60d0 - (ll2**2*llzp)/10d0 + (7d0*zeta3)/20d0) + 
     &zp**6*(-((pi**2)/432d0) - 505931d0/2592000d0 - (pi**2*l
     &l2)/36d0 + (ll2**2)/72d0 + (ll2**3)/18d0 + (33497d0*llz
     &p)/345600d0 + (pi**2*llzp)/72d0 - (ll2**2*llzp)/12d0 + 
     &(7d0*zeta3)/24d0) + zp*(-((pi**2)/12d0) - (pi**2*ll2)/6
     &d0 + (ll2**2)/2d0 + (ll2**3)/3d0 + (pi**2*llzp)/12d0 - 
     &(ll2**2*llzp)/2d0 + (7d0*zeta3)/4d0) + zp**2*(-((pi**2)
     &/48d0) - 5d0/8d0 - (pi**2*ll2)/12d0 + (ll2**2)/8d0 + (l
     &l2**3)/6d0 + (pi**2*llzp)/24d0 + (llzp)/4d0 - (ll2**2*l
     &lzp)/4d0 + (7d0*zeta3)/8d0)

         case(49)            !010-1

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**4*71d0)/1440d0 + (pi**2*ll2**2)/
     &6d0 - (ll2**4)/6d0 - 4*cli4pt5 - (7d0*ll2*zeta3)/2d0 + 
     &zp**2*(-((pi**2)/48d0) - (ll2)/2d0 + (pi**2*llzp)/24d0 
     &+ (5d0*zeta3)/16d0) + zp**3*(-((pi**2)/108d0) + 1d0/24d
     &0 - (5d0*ll2)/12d0 + (pi**2*llzp)/36d0 + (5d0*zeta3)/24
     &d0) + zp**4*(-((pi**2)/192d0) + 7d0/144d0 - (49d0*ll2)/
     &144d0 + (pi**2*llzp)/48d0 + (5d0*zeta3)/32d0) + zp**6*(
     &-((pi**2)/432d0) + 3793d0/86400d0 - (5269d0*ll2)/21600d
     &0 + (pi**2*llzp)/72d0 + (5d0*zeta3)/48d0) + zp**5*(-((p
     &i**2)/300d0) + 17d0/360d0 - (41d0*ll2)/144d0 + (pi**2*l
     &lzp)/60d0 + (zeta3)/8d0) + zp*(-((pi**2)/12d0) + (pi**2
     &*llzp)/12d0 + (5d0*zeta3)/8d0)

         case(50)            !0100

            zp = 1d0-x

            ris = (pi**4)/30d0 - zp*zeta3 - (zp**2*zeta
     &3)/2d0 + zp**3*(1d0/12d0 - (zeta3)/3d0) + zp**4*(5d0/48
     &d0 - (zeta3)/4d0) + zp**5*(17d0/160d0 - (zeta3)/5d0) + 
     &zp**6*(59d0/576d0 - (zeta3)/6d0)

         case(51)            !0101

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**4)/120d0 + zp**2*(-((pi**2)/24d0
     &) - 5d0/4d0 + (pi**2*llzp)/12d0 + (llzp)/2d0 + zeta3) +
     & zp*(-((pi**2)/6d0) + (pi**2*llzp)/6d0 + 2*zeta3) + zp*
     &*4*(-(1151d0/1728d0) - (pi**2)/96d0 + (pi**2*llzp)/24d0
     & + (49d0*llzp)/144d0 + (zeta3)/2d0) + zp**6*(-((pi**2)/
     &216d0) - 17653d0/40500d0 + (pi**2*llzp)/36d0 + (5269d0*
     &llzp)/21600d0 + (zeta3)/3d0) + zp**3*(-((pi**2)/54d0) -
     & 8d0/9d0 + (pi**2*llzp)/18d0 + (5d0*llzp)/12d0 + (2d0*z
     &eta3)/3d0) + zp**5*(-((pi**2)/150d0) - 2281d0/4320d0 + 
     &(pi**2*llzp)/30d0 + (41d0*llzp)/144d0 + (2d0*zeta3)/5d0
     &)

         case(52)            !011-1

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**4)/80d0) + (pi**2*ll2**2)/12d0
     & + (ll2**4)/24d0 + cli4pt5 + (7d0*ll2*zeta3)/8d0 + zp**
     &2*((pi**2)/48d0 + 1d0/4d0 + (pi**2*ll2)/24d0 - (ll2)/8d
     &0 - (ll2**2)/8d0 - (ll2**3)/12d0 - (pi**2*llzp)/24d0 + 
     &(ll2*llzp)/4d0 + (ll2**2*llzp)/4d0 - (ll2*llzp**2)/4d0 
     &- (7d0*zeta3)/16d0) + zp**3*((pi**2)/108d0 + 17d0/96d0 
     &- (ll2)/27d0 + (pi**2*ll2)/36d0 - (ll2**2)/18d0 - (ll2*
     &*3)/18d0 - (pi**2*llzp)/36d0 + (ll2*llzp)/9d0 + (ll2**2
     &*llzp)/6d0 - (ll2*llzp**2)/6d0 - (7d0*zeta3)/24d0) + zp
     &**4*((pi**2)/192d0 + 463d0/3456d0 + (pi**2*ll2)/48d0 - 
     &(ll2)/64d0 - (ll2**2)/32d0 - (ll2**3)/24d0 - (pi**2*llz
     &p)/48d0 + (ll2*llzp)/16d0 + (ll2**2*llzp)/8d0 - (ll2*ll
     &zp**2)/8d0 - (7d0*zeta3)/32d0) + zp**5*(14843d0/138240d
     &0 + (pi**2)/300d0 - (ll2)/125d0 + (pi**2*ll2)/60d0 - (l
     &l2**2)/50d0 - (ll2**3)/30d0 - (pi**2*llzp)/60d0 + (ll2*
     &llzp)/25d0 + (ll2**2*llzp)/10d0 - (ll2*llzp**2)/10d0 - 
     &(7d0*zeta3)/40d0) + zp**6*(1856239d0/20736000d0 + (pi**
     &2)/432d0 - (ll2)/216d0 + (pi**2*ll2)/72d0 - (ll2**2)/72
     &d0 - (ll2**3)/36d0 - (pi**2*llzp)/72d0 + (ll2*llzp)/36d
     &0 + (ll2**2*llzp)/12d0 - (ll2*llzp**2)/12d0 - (7d0*zeta
     &3)/48d0) + zp*((pi**2)/12d0 - ll2 + (pi**2*ll2)/12d0 - 
     &(ll2**2)/2d0 - (ll2**3)/6d0 - (pi**2*llzp)/12d0 + ll2*l
     &lzp + (ll2**2*llzp)/2d0 - (ll2*llzp**2)/2d0 - (7d0*zeta
     &3)/8d0)

         case(53)            !0110

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**4)/72d0) + zp*((pi**2)/6d0 - (
     &pi**2*llzp)/6d0 - zeta3) + zp**2*((pi**2)/24d0 + 1d0/2d
     &0 - (pi**2*llzp)/12d0 - (zeta3)/2d0) + zp**3*((pi**2)/5
     &4d0 + 3d0/8d0 - (pi**2*llzp)/18d0 - (zeta3)/3d0) + zp**
     &4*(251d0/864d0 + (pi**2)/96d0 - (pi**2*llzp)/24d0 - (ze
     &ta3)/4d0) + zp**5*((pi**2)/150d0 + 407d0/1728d0 - (pi**
     &2*llzp)/30d0 - (zeta3)/5d0) + zp**6*((pi**2)/216d0 + 25
     &6103d0/1296000d0 - (pi**2*llzp)/36d0 - (zeta3)/6d0)

         case(54)            !0111

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**4)/90d0 + zp**2*(-(1d0/16d0) + (
     &llzp)/8d0 - (llzp**2)/8d0 + (llzp**3)/12d0) + zp**3*(-(
     &1d0/81d0) + (llzp)/27d0 - (llzp**2)/18d0 + (llzp**3)/18
     &d0) + zp**4*(-(1d0/256d0) + (llzp)/64d0 - (llzp**2)/32d
     &0 + (llzp**3)/24d0) + zp**5*(-(1d0/625d0) + (llzp)/125d
     &0 - (llzp**2)/50d0 + (llzp**3)/30d0) + zp**6*(-(1d0/129
     &6d0) + (llzp)/216d0 - (llzp**2)/72d0 + (llzp**3)/36d0) 
     &+ zp*(-1 + llzp - (llzp**2)/2d0 + (llzp**3)/6d0)

         case(55)            !1-1-1-1

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**4)/90d0) - (pi**2*ll2**2)/12d0
     & + (zp*ll2**2)/4d0 + (ll2**4)/6d0 + zp**3*(1d0/144d0 - 
     &(ll2)/48d0 + (ll2**2)/144d0) + zp**5*(7d0/3840d0 - (ll2
     &)/384d0 + (ll2**2)/1600d0) + zp**2*(-((ll2)/16d0) + (ll
     &2**2)/32d0) + zp**6*(5d0/6144d0 - (137d0*ll2)/138240d0 
     &+ (ll2**2)/4608d0) + zp**4*(1d0/256d0 - (11d0*ll2)/1536
     &d0 + (ll2**2)/512d0) - (ll2**3*llzp)/6d0 + cli4pt5 + ll
     &2*zeta3

         case(56)            !1-1-10

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**4*7d0)/360d0 - (pi**2*zp)/24d0 -
     & (pi**2*zp**2)/192d0 + (1d0/72d0 - (pi**2)/864d0)*zp**3
     & + (-((pi**2)/3072d0) + 7d0/768d0)*zp**4 + (1d0/200d0 -
     & (pi**2)/9600d0)*zp**5 + (-((pi**2)/27648d0) + 23d0/864
     &0d0)*zp**6 - (pi**2*ll2**2)/12d0 - (ll2**4)/12d0 + (pi*
     &*2*ll2*llzp)/12d0 - 2*cli4pt5 - (llzp*zeta3)/8d0

         case(57)            !1-1-11

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**4)/288d0) + (pi**2*ll2**2)/24d
     &0 - (ll2**4)/8d0 + zp*((pi**2)/24d0 - (ll2**2)/4d0) + (
     &ll2**3*llzp)/6d0 + zp**6*((pi**2)/27648d0 - 8009d0/8294
     &400d0 - (ll2**2)/4608d0 + (137d0*llzp)/138240d0) + zp**
     &4*((pi**2)/3072d0 - 41d0/4608d0 - (ll2**2)/512d0 + (11d
     &0*llzp)/1536d0) + zp**2*((pi**2)/192d0 - 1d0/8d0 - (ll2
     &**2)/32d0 + (llzp)/16d0) + zp**5*(-(13d0/4608d0) + (pi*
     &*2)/9600d0 - (ll2**2)/1600d0 + (llzp)/384d0) + zp**3*(-
     &(1d0/32d0) + (pi**2)/864d0 - (ll2**2)/144d0 + (llzp)/48
     &d0) - (llzp*zeta3)/8d0

         case(58)            !1-10-1

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**4)/24d0) + (pi**2*zp)/24d0 - (
     &pi**2*ll2**2)/8d0 + (ll2**4)/6d0 + zp**6*((pi**2)/27648
     &d0 + 97d0/23040d0 - (ll2)/135d0) + zp**3*(1d0/72d0 + (p
     &i**2)/864d0 - (ll2)/18d0) + zp**4*((pi**2)/3072d0 + 1d0
     &/96d0 - (5d0*ll2)/192d0) + zp**5*(1d0/150d0 + (pi**2)/9
     &600d0 - (ll2)/75d0) + zp**2*((pi**2)/192d0 - (ll2)/8d0)
     & - (pi**2*ll2*llzp)/12d0 + 4*cli4pt5 + (21d0*ll2*zeta3)
     &/8d0 + (llzp*zeta3)/4d0

         case(59)            !1-100

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**4)/180d0 + (zp**3)/36d0 + (3d0*z
     &p**4)/128d0 + (zp**5)/60d0 + (5d0*zp**6)/432d0 + (pi**2
     &*ll2**2)/12d0 - (ll2**4)/12d0 - 2*cli4pt5 - (3d0*llzp*z
     &eta3)/4d0

         case(60)            !1-101

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**4*29d0)/1440d0) + (pi**2*zp)/1
     &2d0 - (pi**2*ll2**2)/24d0 + (ll2**4)/8d0 - (pi**2*ll2*l
     &lzp)/6d0 + zp**6*((pi**2)/13824d0 - 667d0/129600d0 + (l
     &lzp)/135d0) + zp**3*(-(17d0/216d0) + (pi**2)/432d0 + (l
     &lzp)/18d0) + zp**4*((pi**2)/1536d0 - 65d0/2304d0 + (5d0
     &*llzp)/192d0) + zp**5*((pi**2)/4800d0 - 103d0/9000d0 + 
     &(llzp)/75d0) + zp**2*(-(1d0/4d0) + (pi**2)/96d0 + (llzp
     &)/8d0) + 3*cli4pt5 + (5d0*llzp*zeta3)/8d0

         case(61)            !1-11-1

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**4)/160d0 + (pi**2*ll2**2)/8d0 - 
     &(ll2**4)/8d0 - (pi**2*ll2*llzp)/12d0 + (ll2**3*llzp)/6d
     &0 + zp**2*(1d0/16d0 - (pi**2)/192d0 + (ll2)/16d0 + (ll2
     &**2)/32d0 - (ll2*llzp)/16d0) + zp**6*(-((pi**2)/27648d0
     &) + 5269d0/8294400d0 + (ll2)/6912d0 + (ll2**2)/4608d0 -
     & (ll2*llzp)/2304d0) + zp**4*(-((pi**2)/3072d0) + 49d0/9
     &216d0 + (ll2)/512d0 + (ll2**2)/512d0 - (ll2*llzp)/256d0
     &) + zp*(-((pi**2)/24d0) + ll2 + (ll2**2)/4d0 - (ll2*llz
     &p)/2d0) + zp**3*(5d0/288d0 - (pi**2)/864d0 + (ll2)/108d
     &0 + (ll2**2)/144d0 - (ll2*llzp)/72d0) + zp**5*(41d0/230
     &40d0 - (pi**2)/9600d0 + (ll2)/2000d0 + (ll2**2)/1600d0 
     &- (ll2*llzp)/800d0) - 2*ll2*zeta3 + (llzp*zeta3)/4d0

         case(62)            !1-110

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**4)/96d0) - (pi**2*zp)/12d0 + (
     &1d0/8d0 - (pi**2)/96d0)*zp**2 + (1d0/24d0 - (pi**2)/432
     &d0)*zp**3 + (-((pi**2)/1536d0) + 35d0/2304d0)*zp**4 + (
     &11d0/1800d0 - (pi**2)/4800d0)*zp**5 + (-((pi**2)/13824d
     &0) + 347d0/129600d0)*zp**6 + (pi**2*ll2**2)/24d0 + (ll2
     &**4)/8d0 - (pi**2*ll2*llzp)/12d0 + 3*cli4pt5 + llzp*zet
     &a3

         case(63)            !1-111

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**2*ll2*llzp)/12d0 - (ll2**3*llzp)
     &/6d0 + zp**3*(1d0/216d0 - (llzp)/108d0 + (llzp**2)/144d
     &0) + zp**5*(3d0/20000d0 - (llzp)/2000d0 + (llzp**2)/160
     &0d0) + zp**2*(3d0/64d0 - (llzp)/16d0 + (llzp**2)/32d0) 
     &+ zp**6*(1d0/27648d0 - (llzp)/6912d0 + (llzp**2)/4608d0
     &) + zp*(3d0/2d0 - llzp + (llzp**2)/4d0) + zp**4*(3d0/40
     &96d0 - (llzp)/512d0 + (llzp**2)/512d0) - 3*cli4pt5 - (7
     &d0*llzp*zeta3)/8d0

         case(64)            !10-1-1

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**4)/480d0) + (zp*ll2**2)/2d0 + 
     &zp**3*(1d0/72d0 - (5d0*ll2)/72d0 + (ll2**2)/18d0) + zp*
     &*4*(3d0/256d0 - (ll2)/24d0 + (ll2**2)/32d0) + zp**5*(83
     &d0/9600d0 - (131d0*ll2)/4800d0 + (ll2**2)/50d0) + zp**6
     &*(11d0/1728d0 - (661d0*ll2)/34560d0 + (ll2**2)/72d0) + 
     &zp**2*(-((ll2)/8d0) + (ll2**2)/8d0) - (llzp*zeta3)/8d0

         case(65)            !10-10

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**4*17d0)/1440d0 - (pi**2*zp)/12d0
     & - (pi**2*zp**2)/48d0 + (-((pi**2)/108d0) + 1d0/36d0)*z
     &p**3 + (-((pi**2)/192d0) + 5d0/192d0)*zp**4 + (-((pi**2
     &)/300d0) + 1d0/48d0)*zp**5 + (-((pi**2)/432d0) + 47d0/2
     &880d0)*zp**6 + (3d0*llzp*zeta3)/2d0

         case(66)            !10-11

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**4)/160d0 + (pi**2*ll2**2)/8d0 - 
     &(ll2**4)/8d0 + zp*((pi**2)/12d0 - (ll2**2)/2d0) + (pi**
     &2*ll2*llzp)/4d0 + zp**4*((pi**2)/192d0 - 131d0/2304d0 -
     & (ll2**2)/32d0 + (llzp)/24d0) + zp**5*((pi**2)/300d0 - 
     &9829d0/288000d0 - (ll2**2)/50d0 + (131d0*llzp)/4800d0) 
     &+ zp**6*((pi**2)/432d0 - 46717d0/2073600d0 - (ll2**2)/7
     &2d0 + (661d0*llzp)/34560d0) + zp**3*((pi**2)/108d0 - 47
     &d0/432d0 - (ll2**2)/18d0 + (5d0*llzp)/72d0) + zp**2*((p
     &i**2)/48d0 - 1d0/4d0 - (ll2**2)/8d0 + (llzp)/8d0) - 3*c
     &li4pt5 - (13d0*llzp*zeta3)/8d0

         case(67)            !100-1

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**4*11d0)/360d0) + (pi**2*zp)/12
     &d0 - (pi**2*ll2**2)/12d0 + (ll2**4)/12d0 + zp**5*((pi**
     &2)/300d0 + 1d0/40d0 - (ll2)/12d0) + zp**6*((pi**2)/432d
     &0 + 731d0/34560d0 - (137d0*ll2)/2160d0) + zp**2*((pi**2
     &)/48d0 - (ll2)/4d0) + zp**3*((pi**2)/108d0 + 1d0/36d0 -
     & (ll2)/6d0) + zp**4*((pi**2)/192d0 + 11d0/384d0 - (11d0
     &*ll2)/96d0) + 2*cli4pt5 + (7d0*ll2*zeta3)/4d0 - (3d0*ll
     &zp*zeta3)/4d0

         case(68)            !1000

            zp = 1d0-x

            ris = -((pi**4)/90d0) + (zp**3)/18d0 + (zp*
     &*4)/16d0 + (7d0*zp**5)/120d0 + (5d0*zp**6)/96d0

         case(69)            !1001

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**4)/72d0) + (pi**2*zp)/6d0 + zp
     &**5*(-(13d0/144d0) + (pi**2)/150d0 + (llzp)/12d0) + zp*
     &*6*((pi**2)/216d0 - 8009d0/129600d0 + (137d0*llzp)/2160
     &d0) + zp**2*((pi**2)/24d0 - 1d0/2d0 + (llzp)/4d0) + zp*
     &*3*(-(1d0/4d0) + (pi**2)/54d0 + (llzp)/6d0) + zp**4*(-(
     &41d0/288d0) + (pi**2)/96d0 + (11d0*llzp)/96d0) - llzp*z
     &eta3

         case(70)            !101-1

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**4)/1440d0 - (pi**2*ll2**2)/24d0 
     &+ (ll2**4)/24d0 - (pi**2*ll2*llzp)/4d0 + zp*(-((pi**2)/
     &12d0) + 2*ll2 + (ll2**2)/2d0 - ll2*llzp) + zp**4*(-((pi
     &**2)/192d0) + 83d0/2304d0 + (ll2)/32d0 + (ll2**2)/32d0 
     &- (ll2*llzp)/16d0) + zp**5*(-((pi**2)/300d0) + 1337d0/5
     &7600d0 + (2d0*ll2)/125d0 + (ll2**2)/50d0 - (ll2*llzp)/2
     &5d0) + zp**6*(33497d0/2073600d0 - (pi**2)/432d0 + (ll2)
     &/108d0 + (ll2**2)/72d0 - (ll2*llzp)/36d0) + zp**2*(-((p
     &i**2)/48d0) + 1d0/8d0 + (ll2)/4d0 + (ll2**2)/8d0 - (ll2
     &*llzp)/4d0) + zp**3*(-((pi**2)/108d0) + 1d0/16d0 + (2d0
     &*ll2)/27d0 + (ll2**2)/18d0 - (ll2*llzp)/9d0) + cli4pt5 
     &- (7d0*ll2*zeta3)/4d0 + llzp*zeta3

         case(71)            !1010

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**4*7d0)/360d0 - (pi**2*zp)/6d0 + 
     &(-((pi**2)/24d0) + 1d0/4d0)*zp**2 + (-((pi**2)/54d0) + 
     &5d0/36d0)*zp**3 + (49d0/576d0 - (pi**2)/96d0)*zp**4 + (
     &-((pi**2)/150d0) + 41d0/720d0)*zp**5 + (-((pi**2)/216d0
     &) + 5269d0/129600d0)*zp**6 + 2*llzp*zeta3

         case(72)            !1011

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**4)/30d0) + zp**3*(1d0/27d0 - (
     &2d0*llzp)/27d0 + (llzp**2)/18d0) + zp*(3 - 2*llzp + (ll
     &zp**2)/2d0) + zp**4*(3d0/256d0 - (llzp)/32d0 + (llzp**2
     &)/32d0) + zp**5*(3d0/625d0 - (2d0*llzp)/125d0 + (llzp**
     &2)/50d0) + zp**6*(1d0/432d0 - (llzp)/108d0 + (llzp**2)/
     &72d0) + zp**2*(3d0/16d0 - (llzp)/4d0 + (llzp**2)/8d0) -
     & llzp*zeta3

         case(73)            !11-1-1

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**4)/720d0) - (zp*ll2)/2d0 - (pi
     &**2*ll2**2)/12d0 + (ll2**4)/8d0 + zp**4*(11d0/6144d0 - 
     &(ll2)/1024d0) + zp**6*(137d0/829440d0 - (ll2)/13824d0) 
     &+ zp**3*(1d0/144d0 - (ll2)/216d0) + zp**2*(1d0/32d0 - (
     &ll2)/32d0) + zp**5*(1d0/1920d0 - (ll2)/4000d0) + (pi**2
     &*ll2*llzp)/12d0 - (ll2**3*llzp)/3d0 + (ll2**2*llzp**2)/
     &4d0 + ll2*zeta3 - (llzp*zeta3)/8d0

         case(74)            !11-10

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**4*11d0)/720d0 + (zp**2)/16d0 + (
     &zp**3)/54d0 + (5d0*zp**4)/768d0 + (zp**5)/375d0 + (zp**
     &6)/810d0 - (ll2**4)/8d0 + (pi**2*ll2*llzp)/4d0 - (pi**2
     &*llzp**2)/24d0 - 3*cli4pt5 - (13d0*llzp*zeta3)/8d0

         case(75)            !11-11

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**2*ll2*llzp)/6d0) + (ll2**3*llz
     &p)/3d0 + (pi**2*llzp**2)/24d0 - (ll2**2*llzp**2)/4d0 + 
     &zp**4*(-(3d0/4096d0) + (llzp)/1024d0) + zp**6*(-(1d0/27
     &648d0) + (llzp)/13824d0) + zp**3*(-(1d0/216d0) + (llzp)
     &/216d0) + zp*(-(3d0/2d0) + (llzp)/2d0) + zp**2*(-(3d0/6
     &4d0) + (llzp)/32d0) + zp**5*(-(3d0/20000d0) + (llzp)/40
     &00d0) + 3*cli4pt5 + (7d0*llzp*zeta3)/4d0

         case(76)            !110-1

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**4)/288d0) - zp*ll2 - (pi**2*ll
     &2**2)/24d0 + (ll2**4)/24d0 + zp**5*(131d0/24000d0 - (ll
     &2)/125d0) + zp**6*(661d0/207360d0 - (ll2)/216d0) + zp**
     &3*(5d0/216d0 - (ll2)/27d0) + zp**4*(1d0/96d0 - (ll2)/64
     &d0) + zp**2*(1d0/16d0 - (ll2)/8d0) + (pi**2*llzp**2)/24
     &d0 + cli4pt5 + (7d0*ll2*zeta3)/8d0 + (5d0*llzp*zeta3)/8
     &d0

         case(77)            !1100

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**4)/360d0) + (zp**2)/8d0 + (zp*
     &*3)/18d0 + (11d0*zp**4)/384d0 + (zp**5)/60d0 + (137d0*z
     &p**6)/12960d0 - llzp*zeta3

         case(78)            !1101

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**4)/30d0 + zp*(-3 + llzp) + (pi**
     &2*llzp**2)/12d0 + zp**5*(-(3d0/625d0) + (llzp)/125d0) +
     & zp**6*(-(1d0/432d0) + (llzp)/216d0) + zp**3*(-(1d0/27d
     &0) + (llzp)/27d0) + zp**4*(-(3d0/256d0) + (llzp)/64d0) 
     &+ zp**2*(-(3d0/16d0) + (llzp)/8d0) + 2*llzp*zeta3

         case(79)            !111-1

            zp = 1d0-x
            llzp = log(zp)

            ris = (zp)/2d0 + (zp**2)/64d0 + (zp**3)/648
     &d0 + (zp**4)/4096d0 + (zp**5)/20000d0 + (zp**6)/82944d0
     & + (pi**2*ll2*llzp)/12d0 - (ll2**3*llzp)/6d0 - (pi**2*l
     &lzp**2)/24d0 + (ll2**2*llzp**2)/4d0 - (ll2*llzp**3)/6d0
     & - cli4pt5 - (7d0*llzp*zeta3)/8d0

         case(80)            !1110

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**4)/90d0) + zp + (zp**2)/16d0 +
     & (zp**3)/81d0 + (zp**4)/256d0 + (zp**5)/625d0 + (zp**6)
     &/1296d0 - (pi**2*llzp**2)/12d0 - llzp*zeta3

         case(81)            !1111

            zp = 1d0-x
            llzp = log(zp)

            ris = (llzp**4)/24d0
c End of expansions around x = +1

      end select
c --- set the imaginary part back to zero if it has been modified to
c --- get the branch cuts right (and should be zero).
      if (bcflag.eq.1) then
         xre = dreal(x)

         if (n4.eq.0.and.xre.gt.0d0) then
            if (xre.lt.1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
c
         else if (n4.eq.1.and.xre.lt.1d0) then
            if (n1.ne.-1.and.n2.ne.-1.and.n3.ne.-1) then
               ris = dcmplx(dreal(ris),0d0)
            else if (xre.gt.-1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
c            
         else if (n4.eq.-1.and.xre.gt.-1d0) then
            if (n1.ne.1.and.n2.ne.1.and.n3.ne.1) then
               ris = dcmplx(dreal(ris),0d0)
            else if (xre.lt.1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
            
         endif
      endif

      HPL4ar1=ris
      return
      end function
c-source-file HPL4arm1.f
      double complex function HPL4arm1(n1,n2,n3,n4,x)
      implicit none
      integer n1,n2,n3,n4,j,bcflag,s,szp
      double complex x,ris,myi,zp,llzp,cli4,cli4pt5
      double precision pi, zeta2, zeta3,ll2,xre

      pi=3.1415926535897932385D0
      zeta3=1.20205690315959428539973816151d0
      zeta2=pi**2/6d0
      myi = dcmplx(0d0,1d0)
      ll2 = dlog(2d0)
      cli4pt5 = cli4(dcmplx(0.5d0,0d0))
      bcflag = 0

      j=1+(n4+1)*1+(n3+1)*3+(n2+1)*9+(n1+1)*27
      ris = dcmplx(0d0,0d0)

c---  +i*epsilon to get branch cuts right ---
      if (dimag(x).eq.0d0) then
         x = x + dcmplx(0d0,1d-60)
         bcflag = 1
      endif
c---  

      select case (j)

c This was file contains the Taylor 
c expansions around x = -1
c The expansion parameter is zp = x+1

         case(1)            !-1-1-1-1

            zp = x+1d0
            llzp = log(zp)

            ris = (llzp**4)/24d0

         case(2)            !-1-1-10

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = (pi**4)/90d0 - zp - (zp**2)/16d0 - (z
     &p**3)/81d0 - (zp**4)/256d0 - (zp**5)/625d0 - (zp**6)/12
     &96d0 - (zp**7)/2401d0 - (zp**8)/4096d0 - (zp**9)/6561d0
     & + (pi**2*llzp**2)/12d0 + (myi*pi*szp*llzp**3)/6d0 + ll
     &zp*zeta3

         case(3)            !-1-1-11

            zp = x+1d0
            llzp = log(zp)

            ris = (zp)/2d0 + (zp**2)/64d0 + (zp**3)/648
     &d0 + (zp**4)/4096d0 + (zp**5)/20000d0 + (zp**6)/82944d0
     & + (zp**7)/307328d0 + (zp**8)/1048576d0 + (zp**9)/33592
     &32d0 + (pi**2*ll2*llzp)/12d0 - (ll2**3*llzp)/6d0 - (pi*
     &*2*llzp**2)/24d0 + (ll2**2*llzp**2)/4d0 - (ll2*llzp**3)
     &/6d0 - cli4pt5 - (7d0*llzp*zeta3)/8d0

         case(4)            !-1-10-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4)/30d0) + zp*(3 - llzp) - (pi
     &**2*llzp**2)/12d0 + zp**5*(3d0/625d0 - (llzp)/125d0) + 
     &zp**6*(1d0/432d0 - (llzp)/216d0) + zp**3*(1d0/27d0 - (l
     &lzp)/27d0) + zp**7*(3d0/2401d0 - (llzp)/343d0) + zp**8*
     &(3d0/4096d0 - (llzp)/512d0) + zp**4*(3d0/256d0 - (llzp)
     &/64d0) + zp**9*(1d0/2187d0 - (llzp)/729d0) + zp**2*(3d0
     &/16d0 - (llzp)/8d0) - 2*llzp*zeta3

         case(5)            !-1-100

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = -((pi**4)/360d0) - myi*pi*szp*zp + (1
     &d0/8d0 - (myi*pi*szp)/8d0)*zp**2 + (1d0/18d0 - (myi*pi*
     &szp)/27d0)*zp**3 + (11d0/384d0 - (myi*pi*szp)/64d0)*zp*
     &*4 + (1d0/60d0 - (myi*pi*szp)/125d0)*zp**5 + (137d0/129
     &60d0 - (myi*pi*szp)/216d0)*zp**6 + (1d0/140d0 - (myi*pi
     &*szp)/343d0)*zp**7 + (363d0/71680d0 - (myi*pi*szp)/512d
     &0)*zp**8 + (761d0/204120d0 - (myi*pi*szp)/729d0)*zp**9 
     &+ (myi*pi**3*szp*llzp)/6d0 - (pi**2*llzp**2)/4d0 + myi*
     &pi*szp*zeta3 - llzp*zeta3

         case(6)            !-1-101

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4)/288d0 + zp*ll2 + (pi**2*ll2**
     &2)/24d0 - (ll2**4)/24d0 + zp**5*(-(131d0/24000d0) + (ll
     &2)/125d0) + zp**6*(-(661d0/207360d0) + (ll2)/216d0) + z
     &p**3*(-(5d0/216d0) + (ll2)/27d0) + zp**7*(-(1327d0/6585
     &60d0) + (ll2)/343d0) + zp**8*(-(1163d0/860160d0) + (ll2
     &)/512d0) + zp**4*(-(1d0/96d0) + (ll2)/64d0) + zp**9*(-(
     &148969d0/156764160d0) + (ll2)/729d0) + zp**2*(-(1d0/16d
     &0) + (ll2)/8d0) - (pi**2*llzp**2)/24d0 - cli4pt5 - (7d0
     &*ll2*zeta3)/8d0 - (5d0*llzp*zeta3)/8d0

         case(7)            !-1-11-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**2*ll2*llzp)/6d0) + (ll2**3*llz
     &p)/3d0 + (pi**2*llzp**2)/24d0 - (ll2**2*llzp**2)/4d0 + 
     &zp**4*(-(3d0/4096d0) + (llzp)/1024d0) + zp**8*(-(3d0/10
     &48576d0) + (llzp)/131072d0) + zp**6*(-(1d0/27648d0) + (
     &llzp)/13824d0) + zp**3*(-(1d0/216d0) + (llzp)/216d0) + 
     &zp*(-(3d0/2d0) + (llzp)/2d0) + zp**2*(-(3d0/64d0) + (ll
     &zp)/32d0) + zp**9*(-(1d0/1119744d0) + (llzp)/373248d0) 
     &+ zp**5*(-(3d0/20000d0) + (llzp)/4000d0) + zp**7*(-(3d0
     &/307328d0) + (llzp)/43904d0) + 3*cli4pt5 + (7d0*llzp*ze
     &ta3)/4d0

         case(8)            !-1-110

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = -((pi**4*11d0)/720d0) + (myi*pi*szp*z
     &p)/2d0 + (-(1d0/16d0) + (myi*pi*szp)/32d0)*zp**2 + (-(1
     &d0/54d0) + (myi*pi*szp)/216d0)*zp**3 + (-(5d0/768d0) + 
     &(myi*pi*szp)/1024d0)*zp**4 + (-(1d0/375d0) + (myi*pi*sz
     &p)/4000d0)*zp**5 + (-(1d0/810d0) + (myi*pi*szp)/13824d0
     &)*zp**6 + (-(13d0/20580d0) + (myi*pi*szp)/43904d0)*zp**
     &7 + (-(151d0/430080d0) + (myi*pi*szp)/131072d0)*zp**8 +
     & (-(16d0/76545d0) + (myi*pi*szp)/373248d0)*zp**9 + (myi
     &*pi**3*szp*ll2)/12d0 - (myi*pi*szp*ll2**3)/6d0 + (ll2**
     &4)/8d0 - (myi*pi**3*szp*llzp)/12d0 - (pi**2*ll2*llzp)/4
     &d0 + (myi*pi*szp*ll2**2*llzp)/2d0 + (pi**2*llzp**2)/24d
     &0 - (myi*pi*szp*ll2*llzp**2)/2d0 + 3*cli4pt5 - (myi*pi*
     &7d0*szp*zeta3)/8d0 + (13d0*llzp*zeta3)/8d0

         case(9)            !-1-111

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4)/720d0) - (zp*ll2)/2d0 - (pi
     &**2*ll2**2)/12d0 + (ll2**4)/8d0 + zp**4*(11d0/6144d0 - 
     &(ll2)/1024d0) + zp**8*(363d0/18350080d0 - (ll2)/131072d
     &0) + zp**6*(137d0/829440d0 - (ll2)/13824d0) + zp**3*(1d
     &0/144d0 - (ll2)/216d0) + zp**2*(1d0/32d0 - (ll2)/32d0) 
     &+ zp**9*(761d0/104509440d0 - (ll2)/373248d0) + zp**5*(1
     &d0/1920d0 - (ll2)/4000d0) + zp**7*(1d0/17920d0 - (ll2)/
     &43904d0) + (pi**2*ll2*llzp)/12d0 - (ll2**3*llzp)/3d0 + 
     &(ll2**2*llzp**2)/4d0 + ll2*zeta3 - (llzp*zeta3)/8d0

         case(10)            !-10-1-1

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4)/30d0 + zp**8*(-(3d0/4096d0) +
     & (llzp)/256d0 - (llzp**2)/128d0) + zp**9*(-(1d0/2187d0)
     & + (2d0*llzp)/729d0 - (llzp**2)/162d0) + zp**3*(-(1d0/2
     &7d0) + (2d0*llzp)/27d0 - (llzp**2)/18d0) + zp*(-3 + 2*l
     &lzp - (llzp**2)/2d0) + zp**4*(-(3d0/256d0) + (llzp)/32d
     &0 - (llzp**2)/32d0) + zp**5*(-(3d0/625d0) + (2d0*llzp)/
     &125d0 - (llzp**2)/50d0) + zp**6*(-(1d0/432d0) + (llzp)/
     &108d0 - (llzp**2)/72d0) + zp**2*(-(3d0/16d0) + (llzp)/4
     &d0 - (llzp**2)/8d0) + zp**7*(-(3d0/2401d0) + (2d0*llzp)
     &/343d0 - (llzp**2)/98d0) + llzp*zeta3

         case(11)            !-10-10

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = (pi**4*7d0)/360d0 - (myi*pi**3*szp*ll
     &zp)/6d0 + zp*(-((pi**2)/6d0) + 2*myi*pi*szp - myi*pi*sz
     &p*llzp) + zp**4*(49d0/576d0 - (pi**2)/96d0 + (myi*pi*sz
     &p)/32d0 - (myi*pi*szp*llzp)/16d0) + zp**5*(-((pi**2)/15
     &0d0) + 41d0/720d0 + (myi*pi*2d0*szp)/125d0 - (myi*pi*sz
     &p*llzp)/25d0) + zp**6*(-((pi**2)/216d0) + 5269d0/129600
     &d0 + (myi*pi*szp)/108d0 - (myi*pi*szp*llzp)/36d0) + zp*
     &*7*(-((pi**2)/294d0) + 767d0/25200d0 + (myi*pi*2d0*szp)
     &/343d0 - (myi*pi*szp*llzp)/49d0) + zp**2*(-((pi**2)/24d
     &0) + 1d0/4d0 + (myi*pi*szp)/4d0 - (myi*pi*szp*llzp)/4d0
     &) + zp**8*(266681d0/11289600d0 - (pi**2)/384d0 + (myi*p
     &i*szp)/256d0 - (myi*pi*szp*llzp)/64d0) + zp**9*(-((pi**
     &2)/486d0) + 1077749d0/57153600d0 + (myi*pi*2d0*szp)/729
     &d0 - (myi*pi*szp*llzp)/81d0) + zp**3*(-((pi**2)/54d0) +
     & 5d0/36d0 + (myi*pi*2d0*szp)/27d0 - (myi*pi*szp*llzp)/9
     &d0) - 2*myi*pi*szp*zeta3 + 2*llzp*zeta3

         case(12)            !-10-11

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4)/1440d0) + (pi**2*ll2**2)/24
     &d0 - (ll2**4)/24d0 + (pi**2*ll2*llzp)/4d0 + zp*((pi**2)
     &/12d0 - 2*ll2 - (ll2**2)/2d0 + ll2*llzp) + zp**4*((pi**
     &2)/192d0 - 83d0/2304d0 - (ll2)/32d0 - (ll2**2)/32d0 + (
     &ll2*llzp)/16d0) + zp**5*((pi**2)/300d0 - 1337d0/57600d0
     & - (2d0*ll2)/125d0 - (ll2**2)/50d0 + (ll2*llzp)/25d0) +
     & zp**6*(-(33497d0/2073600d0) + (pi**2)/432d0 - (ll2)/10
     &8d0 - (ll2**2)/72d0 + (ll2*llzp)/36d0) + zp**7*(-(5587d
     &0/470400d0) + (pi**2)/588d0 - (2d0*ll2)/343d0 - (ll2**2
     &)/98d0 + (ll2*llzp)/49d0) + zp**2*((pi**2)/48d0 - 1d0/8
     &d0 - (ll2)/4d0 - (ll2**2)/8d0 + (ll2*llzp)/4d0) + zp**8
     &*(-(136919d0/15052800d0) + (pi**2)/768d0 - (ll2)/256d0 
     &- (ll2**2)/128d0 + (ll2*llzp)/64d0) + zp**9*(-(35054939
     &d0/4877107200d0) + (pi**2)/972d0 - (2d0*ll2)/729d0 - (l
     &l2**2)/162d0 + (ll2*llzp)/81d0) + zp**3*((pi**2)/108d0 
     &- 1d0/16d0 - (2d0*ll2)/27d0 - (ll2**2)/18d0 + (ll2*llzp
     &)/9d0) - cli4pt5 + (7d0*ll2*zeta3)/4d0 - llzp*zeta3

         case(13)            !-100-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4)/72d0) + (pi**2*zp)/6d0 + zp
     &**5*(-(13d0/144d0) + (pi**2)/150d0 + (llzp)/12d0) + zp*
     &*7*((pi**2)/294d0 - 161d0/3600d0 + (llzp)/20d0) + zp**6
     &*((pi**2)/216d0 - 8009d0/129600d0 + (137d0*llzp)/2160d0
     &) + zp**2*((pi**2)/24d0 - 1d0/2d0 + (llzp)/4d0) + zp**3
     &*(-(1d0/4d0) + (pi**2)/54d0 + (llzp)/6d0) + zp**9*((pi*
     &*2)/486d0 - 167101d0/6350400d0 + (761d0*llzp)/22680d0) 
     &+ zp**8*((pi**2)/384d0 - 190513d0/5644800d0 + (363d0*ll
     &zp)/8960d0) + zp**4*(-(41d0/288d0) + (pi**2)/96d0 + (11
     &d0*llzp)/96d0) - llzp*zeta3

         case(14)            !-1000

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = -((pi**4*13d0)/180d0) + (pi**2*zp)/2d
     &0 + ((pi**2)/8d0 + (myi*pi*szp)/4d0)*zp**2 + (-(1d0/18d
     &0) + (pi**2)/18d0 + (myi*pi*szp)/6d0)*zp**3 + (-(1d0/16
     &d0) + (pi**2)/32d0 + (myi*pi*11d0*szp)/96d0)*zp**4 + ((
     &pi**2)/50d0 - 7d0/120d0 + (myi*pi*szp)/12d0)*zp**5 + ((
     &pi**2)/72d0 - 5d0/96d0 + (myi*pi*137d0*szp)/2160d0)*zp*
     &*6 + (-(29d0/630d0) + (pi**2)/98d0 + (myi*pi*szp)/20d0)
     &*zp**7 + ((pi**2)/128d0 - 469d0/11520d0 + (myi*pi*363d0
     &*szp)/8960d0)*zp**8 + ((pi**2)/162d0 - 29531d0/816480d0
     & + (myi*pi*761d0*szp)/22680d0)*zp**9 - (myi*pi**3*szp*l
     &lzp)/6d0 - myi*pi*szp*zeta3

         case(15)            !-1001

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4*11d0)/360d0) + (pi**2*zp)/12
     &d0 - (pi**2*ll2**2)/12d0 + (ll2**4)/12d0 + zp**5*((pi**
     &2)/300d0 + 1d0/40d0 - (ll2)/12d0) + zp**7*(103d0/5760d0
     & + (pi**2)/588d0 - (ll2)/20d0) + zp**6*((pi**2)/432d0 +
     & 731d0/34560d0 - (137d0*ll2)/2160d0) + zp**2*((pi**2)/4
     &8d0 - (ll2)/4d0) + zp**3*((pi**2)/108d0 + 1d0/36d0 - (l
     &l2)/6d0) + zp**9*(42799d0/3265920d0 + (pi**2)/972d0 - (
     &761d0*ll2)/22680d0) + zp**8*(3931d0/258048d0 + (pi**2)/
     &768d0 - (363d0*ll2)/8960d0) + zp**4*((pi**2)/192d0 + 11
     &d0/384d0 - (11d0*ll2)/96d0) + 2*cli4pt5 + (7d0*ll2*zeta
     &3)/4d0 - (3d0*llzp*zeta3)/4d0

         case(16)            !-101-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4)/160d0) - (pi**2*ll2**2)/8d0
     & + (ll2**4)/8d0 + zp*(-((pi**2)/12d0) + (ll2**2)/2d0) -
     & (pi**2*ll2*llzp)/4d0 + zp**8*(7401d0/627200d0 - (pi**2
     &)/768d0 + (ll2**2)/128d0 - (1163d0*llzp)/107520d0) + zp
     &**9*(398917091d0/43893964800d0 - (pi**2)/972d0 + (ll2**
     &2)/162d0 - (148969d0*llzp)/17418240d0) + zp**4*(-((pi**
     &2)/192d0) + 131d0/2304d0 + (ll2**2)/32d0 - (llzp)/24d0)
     & + zp**5*(-((pi**2)/300d0) + 9829d0/288000d0 + (ll2**2)
     &/50d0 - (131d0*llzp)/4800d0) + zp**6*(-((pi**2)/432d0) 
     &+ 46717d0/2073600d0 + (ll2**2)/72d0 - (661d0*llzp)/3456
     &0d0) + zp**3*(-((pi**2)/108d0) + 47d0/432d0 + (ll2**2)/
     &18d0 - (5d0*llzp)/72d0) + zp**2*(-((pi**2)/48d0) + 1d0/
     &4d0 + (ll2**2)/8d0 - (llzp)/8d0) + zp**7*(52379d0/32928
     &00d0 - (pi**2)/588d0 + (ll2**2)/98d0 - (1327d0*llzp)/94
     &080d0) + 3*cli4pt5 + (13d0*llzp*zeta3)/8d0

         case(17)            !-1010

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = (pi**4*17d0)/1440d0 + zp*(-((pi**2)/1
     &2d0) + myi*pi*szp*ll2) + zp**4*(-((pi**2)/192d0) + 5d0/
     &192d0 - (myi*pi*szp)/24d0 + (myi*pi*szp*ll2)/16d0) + zp
     &**5*(-((pi**2)/300d0) + 1d0/48d0 - (myi*pi*131d0*szp)/4
     &800d0 + (myi*pi*szp*ll2)/25d0) + zp**6*(-((pi**2)/432d0
     &) + 47d0/2880d0 - (myi*pi*661d0*szp)/34560d0 + (myi*pi*
     &szp*ll2)/36d0) + zp**7*(13d0/1008d0 - (pi**2)/588d0 - (
     &myi*pi*1327d0*szp)/94080d0 + (myi*pi*szp*ll2)/49d0) + z
     &p**2*(-((pi**2)/48d0) - (myi*pi*szp)/8d0 + (myi*pi*szp*
     &ll2)/4d0) + zp**8*(3341d0/322560d0 - (pi**2)/768d0 - (m
     &yi*pi*1163d0*szp)/107520d0 + (myi*pi*szp*ll2)/64d0) + z
     &p**9*(13817d0/1632960d0 - (pi**2)/972d0 - (myi*pi*14896
     &9d0*szp)/17418240d0 + (myi*pi*szp*ll2)/81d0) + zp**3*(-
     &((pi**2)/108d0) + 1d0/36d0 - (myi*pi*5d0*szp)/72d0 + (m
     &yi*pi*szp*ll2)/9d0) - (myi*pi**3*szp*llzp)/12d0 - (myi*
     &pi*5d0*szp*zeta3)/8d0 + (3d0*llzp*zeta3)/2d0

         case(18)            !-1011

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4)/480d0 - (zp*ll2**2)/2d0 + zp*
     &*8*(-(137d0/36864d0) + (1163d0*ll2)/107520d0 - (ll2**2)
     &/128d0) + zp**9*(-(617027d0/209018880d0) + (148969d0*ll
     &2)/17418240d0 - (ll2**2)/162d0) + zp**3*(-(1d0/72d0) + 
     &(5d0*ll2)/72d0 - (ll2**2)/18d0) + zp**4*(-(3d0/256d0) +
     & (ll2)/24d0 - (ll2**2)/32d0) + zp**5*(-(83d0/9600d0) + 
     &(131d0*ll2)/4800d0 - (ll2**2)/50d0) + zp**6*(-(11d0/172
     &8d0) + (661d0*ll2)/34560d0 - (ll2**2)/72d0) + zp**2*((l
     &l2)/8d0 - (ll2**2)/8d0) + zp**7*(-(5417d0/1128960d0) + 
     &(1327d0*ll2)/94080d0 - (ll2**2)/98d0) + (llzp*zeta3)/8d
     &0

         case(19)            !-11-1-1

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**2*ll2*llzp)/12d0 - (ll2**3*llzp)
     &/6d0 + zp**7*(3d0/307328d0 - (llzp)/21952d0 + (llzp**2)
     &/12544d0) + zp**3*(1d0/216d0 - (llzp)/108d0 + (llzp**2)
     &/144d0) + zp**5*(3d0/20000d0 - (llzp)/2000d0 + (llzp**2
     &)/1600d0) + zp**8*(3d0/1048576d0 - (llzp)/65536d0 + (ll
     &zp**2)/32768d0) + zp**2*(3d0/64d0 - (llzp)/16d0 + (llzp
     &**2)/32d0) + zp**6*(1d0/27648d0 - (llzp)/6912d0 + (llzp
     &**2)/4608d0) + zp*(3d0/2d0 - llzp + (llzp**2)/4d0) + zp
     &**4*(3d0/4096d0 - (llzp)/512d0 + (llzp**2)/512d0) + zp*
     &*9*(1d0/1119744d0 - (llzp)/186624d0 + (llzp**2)/82944d0
     &) - 3*cli4pt5 - (7d0*llzp*zeta3)/8d0

         case(20)            !-11-10

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = (pi**4)/96d0 - (myi*pi**3*szp*ll2)/6d
     &0 - (pi**2*ll2**2)/24d0 + (myi*pi*szp*ll2**3)/3d0 - (ll
     &2**4)/8d0 + (myi*pi**3*szp*llzp)/12d0 + (pi**2*ll2*llzp
     &)/12d0 - (myi*pi*szp*ll2**2*llzp)/2d0 + zp**8*(-(9701d0
     &/15052800d0) + (pi**2)/98304d0 - (myi*pi*szp)/65536d0 +
     & (myi*pi*szp*llzp)/16384d0) + zp**2*(-(1d0/8d0) + (pi**
     &2)/96d0 - (myi*pi*szp)/16d0 + (myi*pi*szp*llzp)/16d0) +
     & zp**6*((pi**2)/13824d0 - 347d0/129600d0 - (myi*pi*szp)
     &/6912d0 + (myi*pi*szp*llzp)/2304d0) + zp**4*((pi**2)/15
     &36d0 - 35d0/2304d0 - (myi*pi*szp)/512d0 + (myi*pi*szp*l
     &lzp)/256d0) + zp*((pi**2)/12d0 - myi*pi*szp + (myi*pi*s
     &zp*llzp)/2d0) + zp**9*((pi**2)/248832d0 - 209d0/595350d
     &0 - (myi*pi*szp)/186624d0 + (myi*pi*szp*llzp)/41472d0) 
     &+ zp**7*(-(149d0/117600d0) + (pi**2)/37632d0 - (myi*pi*
     &szp)/21952d0 + (myi*pi*szp*llzp)/6272d0) + zp**3*(-(1d0
     &/24d0) + (pi**2)/432d0 - (myi*pi*szp)/108d0 + (myi*pi*s
     &zp*llzp)/72d0) + zp**5*(-(11d0/1800d0) + (pi**2)/4800d0
     & - (myi*pi*szp)/2000d0 + (myi*pi*szp*llzp)/800d0) - 3*c
     &li4pt5 + (myi*pi*7d0*szp*zeta3)/4d0 - llzp*zeta3

         case(21)            !-11-11

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4)/160d0 + (pi**2*ll2**2)/8d0 - 
     &(ll2**4)/8d0 - (pi**2*ll2*llzp)/12d0 + (ll2**3*llzp)/6d
     &0 + zp**8*(-((pi**2)/196608d0) + 266681d0/2890137600d0 
     &+ (ll2)/65536d0 + (ll2**2)/32768d0 - (ll2*llzp)/16384d0
     &) + zp**2*(1d0/16d0 - (pi**2)/192d0 + (ll2)/16d0 + (ll2
     &**2)/32d0 - (ll2*llzp)/16d0) + zp**6*(-((pi**2)/27648d0
     &) + 5269d0/8294400d0 + (ll2)/6912d0 + (ll2**2)/4608d0 -
     & (ll2*llzp)/2304d0) + zp**4*(-((pi**2)/3072d0) + 49d0/9
     &216d0 + (ll2)/512d0 + (ll2**2)/512d0 - (ll2*llzp)/256d0
     &) + zp*(-((pi**2)/24d0) + ll2 + (ll2**2)/4d0 - (ll2*llz
     &p)/2d0) + zp**9*(1077749d0/29262643200d0 - (pi**2)/4976
     &64d0 + (ll2)/186624d0 + (ll2**2)/82944d0 - (ll2*llzp)/4
     &1472d0) + zp**7*(-((pi**2)/75264d0) + 767d0/3225600d0 +
     & (ll2)/21952d0 + (ll2**2)/12544d0 - (ll2*llzp)/6272d0) 
     &+ zp**3*(5d0/288d0 - (pi**2)/864d0 + (ll2)/108d0 + (ll2
     &**2)/144d0 - (ll2*llzp)/72d0) + zp**5*(41d0/23040d0 - (
     &pi**2)/9600d0 + (ll2)/2000d0 + (ll2**2)/1600d0 - (ll2*l
     &lzp)/800d0) - 2*ll2*zeta3 + (llzp*zeta3)/4d0

         case(22)            !-110-1

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4*29d0)/1440d0 - (pi**2*zp)/12d0
     & + (pi**2*ll2**2)/24d0 - (ll2**4)/8d0 + (pi**2*ll2*llzp
     &)/6d0 + zp**6*(-((pi**2)/13824d0) + 667d0/129600d0 - (l
     &lzp)/135d0) + zp**3*(17d0/216d0 - (pi**2)/432d0 - (llzp
     &)/18d0) + zp**7*(-((pi**2)/37632d0) + 2083d0/823200d0 -
     & (13d0*llzp)/2940d0) + zp**8*(6757d0/5017600d0 - (pi**2
     &)/98304d0 - (151d0*llzp)/53760d0) + zp**4*(-((pi**2)/15
     &36d0) + 65d0/2304d0 - (5d0*llzp)/192d0) + zp**5*(-((pi*
     &*2)/4800d0) + 103d0/9000d0 - (llzp)/75d0) + zp**9*(-((p
     &i**2)/248832d0) + 4121d0/5358150d0 - (16d0*llzp)/8505d0
     &) + zp**2*(1d0/4d0 - (pi**2)/96d0 - (llzp)/8d0) - 3*cli
     &4pt5 - (5d0*llzp*zeta3)/8d0

         case(23)            !-1100

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = (pi**4*17d0)/360d0 - (pi**2*zp)/4d0 +
     & (-((pi**2)/32d0) - (myi*pi*szp)/8d0)*zp**2 + (-((pi**2
     &)/144d0) + 1d0/36d0 - (myi*pi*szp)/18d0)*zp**3 + (3d0/1
     &28d0 - (pi**2)/512d0 - (myi*pi*5d0*szp)/192d0)*zp**4 + 
     &(-((pi**2)/1600d0) + 1d0/60d0 - (myi*pi*szp)/75d0)*zp**
     &5 + (-((pi**2)/4608d0) + 5d0/432d0 - (myi*pi*szp)/135d0
     &)*zp**6 + (-((pi**2)/12544d0) + 41d0/5040d0 - (myi*pi*1
     &3d0*szp)/2940d0)*zp**7 + (-((pi**2)/32768d0) + 539d0/92
     &160d0 - (myi*pi*151d0*szp)/53760d0)*zp**8 + (22d0/5103d
     &0 - (pi**2)/82944d0 - (myi*pi*16d0*szp)/8505d0)*zp**9 -
     & (myi*pi**3*szp*ll2)/4d0 - (pi**2*ll2**2)/6d0 - (ll2**4
     &)/12d0 + (myi*pi**3*szp*llzp)/12d0 + (pi**2*ll2*llzp)/2
     &d0 - 2*cli4pt5 + (myi*pi*13d0*szp*zeta3)/8d0 - (3d0*llz
     &p*zeta3)/4d0

         case(24)            !-1101

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4)/24d0 - (pi**2*zp)/24d0 + (pi*
     &*2*ll2**2)/8d0 - (ll2**4)/6d0 + zp**6*(-((pi**2)/27648d
     &0) - 97d0/23040d0 + (ll2)/135d0) + zp**3*(-(1d0/72d0) -
     & (pi**2)/864d0 + (ll2)/18d0) + zp**7*(-((pi**2)/75264d0
     &) - 767d0/282240d0 + (13d0*ll2)/2940d0) + zp**8*(-((pi*
     &*2)/196608d0) - 935d0/516096d0 + (151d0*ll2)/53760d0) +
     & zp**4*(-((pi**2)/3072d0) - 1d0/96d0 + (5d0*ll2)/192d0)
     & + zp**5*(-(1d0/150d0) - (pi**2)/9600d0 + (ll2)/75d0) +
     & zp**9*(-(2041d0/1632960d0) - (pi**2)/497664d0 + (16d0*
     &ll2)/8505d0) + zp**2*(-((pi**2)/192d0) + (ll2)/8d0) + (
     &pi**2*ll2*llzp)/12d0 - 4*cli4pt5 - (21d0*ll2*zeta3)/8d0
     & - (llzp*zeta3)/4d0

         case(25)            !-111-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4)/288d0) + (pi**2*ll2**2)/24d
     &0 - (ll2**4)/8d0 + zp*((pi**2)/24d0 - (ll2**2)/4d0) + (
     &ll2**3*llzp)/6d0 + zp**6*((pi**2)/27648d0 - 8009d0/8294
     &400d0 - (ll2**2)/4608d0 + (137d0*llzp)/138240d0) + zp**
     &4*((pi**2)/3072d0 - 41d0/4608d0 - (ll2**2)/512d0 + (11d
     &0*llzp)/1536d0) + zp**2*((pi**2)/192d0 - 1d0/8d0 - (ll2
     &**2)/32d0 + (llzp)/16d0) + zp**7*(-(161d0/460800d0) + (
     &pi**2)/75264d0 - (ll2**2)/12544d0 + (llzp)/2560d0) + zp
     &**8*(-(190513d0/1445068800d0) + (pi**2)/196608d0 - (ll2
     &**2)/32768d0 + (363d0*llzp)/2293760d0) + zp**5*(-(13d0/
     &4608d0) + (pi**2)/9600d0 - (ll2**2)/1600d0 + (llzp)/384
     &d0) + zp**3*(-(1d0/32d0) + (pi**2)/864d0 - (ll2**2)/144
     &d0 + (llzp)/48d0) + zp**9*(-(167101d0/3251404800d0) + (
     &pi**2)/497664d0 - (ll2**2)/82944d0 + (761d0*llzp)/11612
     &160d0) - (llzp*zeta3)/8d0

         case(26)            !-1110

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = -((pi**4*7d0)/360d0) + (myi*pi**3*szp
     &*ll2)/12d0 + (pi**2*ll2**2)/12d0 - (myi*pi*szp*ll2**3)/
     &3d0 + (ll2**4)/12d0 + zp**8*(-(1019d0/1290240d0) + (pi*
     &*2)/196608d0 + (myi*pi*363d0*szp)/2293760d0 - (myi*pi*s
     &zp*ll2)/16384d0) + zp**2*((pi**2)/192d0 + (myi*pi*szp)/
     &16d0 - (myi*pi*szp*ll2)/16d0) + zp**6*((pi**2)/27648d0 
     &- 23d0/8640d0 + (myi*pi*137d0*szp)/138240d0 - (myi*pi*s
     &zp*ll2)/2304d0) + zp**4*((pi**2)/3072d0 - 7d0/768d0 + (
     &myi*pi*11d0*szp)/1536d0 - (myi*pi*szp*ll2)/256d0) + zp*
     &((pi**2)/24d0 - (myi*pi*szp*ll2)/2d0) + zp**9*((pi**2)/
     &497664d0 - 23d0/51030d0 + (myi*pi*761d0*szp)/11612160d0
     & - (myi*pi*szp*ll2)/41472d0) + zp**7*(-(101d0/70560d0) 
     &+ (pi**2)/75264d0 + (myi*pi*szp)/2560d0 - (myi*pi*szp*l
     &l2)/6272d0) + zp**3*(-(1d0/72d0) + (pi**2)/864d0 + (myi
     &*pi*szp)/48d0 - (myi*pi*szp*ll2)/72d0) + zp**5*(-(1d0/2
     &00d0) + (pi**2)/9600d0 + (myi*pi*szp)/384d0 - (myi*pi*s
     &zp*ll2)/800d0) - (pi**2*ll2*llzp)/12d0 + (myi*pi*szp*ll
     &2**2*llzp)/2d0 + 2*cli4pt5 - (myi*pi*szp*zeta3)/8d0 + (
     &llzp*zeta3)/8d0

         case(27)            !-1111

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4)/90d0) - (pi**2*ll2**2)/12d0
     & + (zp*ll2**2)/4d0 + (ll2**4)/6d0 + zp**7*(29d0/80640d0
     & - (ll2)/2560d0 + (ll2**2)/12544d0) + zp**3*(1d0/144d0 
     &- (ll2)/48d0 + (ll2**2)/144d0) + zp**5*(7d0/3840d0 - (l
     &l2)/384d0 + (ll2**2)/1600d0) + zp**8*(469d0/2949120d0 -
     & (363d0*ll2)/2293760d0 + (ll2**2)/32768d0) + zp**2*(-((
     &ll2)/16d0) + (ll2**2)/32d0) + zp**6*(5d0/6144d0 - (137d
     &0*ll2)/138240d0 + (ll2**2)/4608d0) + zp**4*(1d0/256d0 -
     & (11d0*ll2)/1536d0 + (ll2**2)/512d0) + zp**9*(29531d0/4
     &18037760d0 - (761d0*ll2)/11612160d0 + (ll2**2)/82944d0)
     & - (ll2**3*llzp)/6d0 + cli4pt5 + ll2*zeta3

         case(28)            !0-1-1-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4)/90d0) + zp**2*(1d0/16d0 - (
     &llzp)/8d0 + (llzp**2)/8d0 - (llzp**3)/12d0) + zp**3*(1d
     &0/81d0 - (llzp)/27d0 + (llzp**2)/18d0 - (llzp**3)/18d0)
     & + zp**4*(1d0/256d0 - (llzp)/64d0 + (llzp**2)/32d0 - (l
     &lzp**3)/24d0) + zp**5*(1d0/625d0 - (llzp)/125d0 + (llzp
     &**2)/50d0 - (llzp**3)/30d0) + zp**6*(1d0/1296d0 - (llzp
     &)/216d0 + (llzp**2)/72d0 - (llzp**3)/36d0) + zp**7*(1d0
     &/2401d0 - (llzp)/343d0 + (llzp**2)/98d0 - (llzp**3)/42d
     &0) + zp**8*(1d0/4096d0 - (llzp)/512d0 + (llzp**2)/128d0
     & - (llzp**3)/48d0) + zp**9*(1d0/6561d0 - (llzp)/729d0 +
     & (llzp**2)/162d0 - (llzp**3)/54d0) + zp*(1 - llzp + (ll
     &zp**2)/2d0 - (llzp**3)/6d0)

         case(29)            !0-1-10

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = -((pi**4)/72d0) + zp*((pi**2)/6d0 - m
     &yi*pi*szp - (pi**2*llzp)/6d0 + myi*pi*szp*llzp - (myi*p
     &i*szp*llzp**2)/2d0 - zeta3) + myi*pi*szp*zeta3 + zp**2*
     &((pi**2)/24d0 + 1d0/2d0 - (myi*pi*szp)/8d0 - (pi**2*llz
     &p)/12d0 + (myi*pi*szp*llzp)/4d0 - (myi*pi*szp*llzp**2)/
     &4d0 - (zeta3)/2d0) + zp**3*((pi**2)/54d0 + 3d0/8d0 - (m
     &yi*pi*szp)/27d0 - (pi**2*llzp)/18d0 + (myi*pi*szp*llzp)
     &/9d0 - (myi*pi*szp*llzp**2)/6d0 - (zeta3)/3d0) + zp**4*
     &(251d0/864d0 + (pi**2)/96d0 - (myi*pi*szp)/64d0 - (pi**
     &2*llzp)/24d0 + (myi*pi*szp*llzp)/16d0 - (myi*pi*szp*llz
     &p**2)/8d0 - (zeta3)/4d0) + zp**5*((pi**2)/150d0 + 407d0
     &/1728d0 - (myi*pi*szp)/125d0 - (pi**2*llzp)/30d0 + (myi
     &*pi*szp*llzp)/25d0 - (myi*pi*szp*llzp**2)/10d0 - (zeta3
     &)/5d0) + zp**6*((pi**2)/216d0 + 256103d0/1296000d0 - (m
     &yi*pi*szp)/216d0 - (pi**2*llzp)/36d0 + (myi*pi*szp*llzp
     &)/36d0 - (myi*pi*szp*llzp**2)/12d0 - (zeta3)/6d0) + zp*
     &*7*((pi**2)/294d0 + 4081d0/24000d0 - (myi*pi*szp)/343d0
     & - (pi**2*llzp)/42d0 + (myi*pi*szp*llzp)/49d0 - (myi*pi
     &*szp*llzp**2)/14d0 - (zeta3)/7d0) + zp**8*((pi**2)/384d
     &0 + 9822481d0/65856000d0 - (myi*pi*szp)/512d0 - (pi**2*
     &llzp)/48d0 + (myi*pi*szp*llzp)/64d0 - (myi*pi*szp*llzp*
     &*2)/16d0 - (zeta3)/8d0) + zp**9*((pi**2)/486d0 + 787084
     &73d0/592704000d0 - (myi*pi*szp)/729d0 - (pi**2*llzp)/54
     &d0 + (myi*pi*szp*llzp)/81d0 - (myi*pi*szp*llzp**2)/18d0
     & - (zeta3)/9d0)

         case(30)            !0-1-11

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4)/80d0 - (pi**2*ll2**2)/12d0 - 
     &(ll2**4)/24d0 - cli4pt5 - (7d0*ll2*zeta3)/8d0 + zp**2*(
     &-((pi**2)/48d0) - 1d0/4d0 - (pi**2*ll2)/24d0 + (ll2)/8d
     &0 + (ll2**2)/8d0 + (ll2**3)/12d0 + (pi**2*llzp)/24d0 - 
     &(ll2*llzp)/4d0 - (ll2**2*llzp)/4d0 + (ll2*llzp**2)/4d0 
     &+ (7d0*zeta3)/16d0) + zp**3*(-((pi**2)/108d0) - 17d0/96
     &d0 + (ll2)/27d0 - (pi**2*ll2)/36d0 + (ll2**2)/18d0 + (l
     &l2**3)/18d0 + (pi**2*llzp)/36d0 - (ll2*llzp)/9d0 - (ll2
     &**2*llzp)/6d0 + (ll2*llzp**2)/6d0 + (7d0*zeta3)/24d0) +
     & zp**4*(-((pi**2)/192d0) - 463d0/3456d0 - (pi**2*ll2)/4
     &8d0 + (ll2)/64d0 + (ll2**2)/32d0 + (ll2**3)/24d0 + (pi*
     &*2*llzp)/48d0 - (ll2*llzp)/16d0 - (ll2**2*llzp)/8d0 + (
     &ll2*llzp**2)/8d0 + (7d0*zeta3)/32d0) + zp**5*(-(14843d0
     &/138240d0) - (pi**2)/300d0 + (ll2)/125d0 - (pi**2*ll2)/
     &60d0 + (ll2**2)/50d0 + (ll2**3)/30d0 + (pi**2*llzp)/60d
     &0 - (ll2*llzp)/25d0 - (ll2**2*llzp)/10d0 + (ll2*llzp**2
     &)/10d0 + (7d0*zeta3)/40d0) + zp**6*(-(1856239d0/2073600
     &0d0) - (pi**2)/432d0 + (ll2)/216d0 - (pi**2*ll2)/72d0 +
     & (ll2**2)/72d0 + (ll2**3)/36d0 + (pi**2*llzp)/72d0 - (l
     &l2*llzp)/36d0 - (ll2**2*llzp)/12d0 + (ll2*llzp**2)/12d0
     & + (7d0*zeta3)/48d0) + zp**8*(-((pi**2)/768d0) - 636802
     &727d0/9483264000d0 + (ll2)/512d0 - (pi**2*ll2)/96d0 + (
     &ll2**2)/128d0 + (ll2**3)/48d0 + (pi**2*llzp)/96d0 - (ll
     &2*llzp)/64d0 - (ll2**2*llzp)/16d0 + (ll2*llzp**2)/16d0 
     &+ (7d0*zeta3)/64d0) + zp**9*(-(81511906681d0/1365590016
     &000d0) - (pi**2)/972d0 - (pi**2*ll2)/108d0 + (ll2)/729d
     &0 + (ll2**2)/162d0 + (ll2**3)/54d0 + (pi**2*llzp)/108d0
     & - (ll2*llzp)/81d0 - (ll2**2*llzp)/18d0 + (ll2*llzp**2)
     &/18d0 + (7d0*zeta3)/72d0) + zp**7*(-(1856489d0/24192000
     &d0) - (pi**2)/588d0 + (ll2)/343d0 - (pi**2*ll2)/84d0 + 
     &(ll2**2)/98d0 + (ll2**3)/42d0 + (pi**2*llzp)/84d0 - (ll
     &2*llzp)/49d0 - (ll2**2*llzp)/14d0 + (ll2*llzp**2)/14d0 
     &+ (zeta3)/8d0) + zp*(-((pi**2)/12d0) + ll2 - (pi**2*ll2
     &)/12d0 + (ll2**2)/2d0 + (ll2**3)/6d0 + (pi**2*llzp)/12d
     &0 - ll2*llzp - (ll2**2*llzp)/2d0 + (ll2*llzp**2)/2d0 + 
     &(7d0*zeta3)/8d0)

         case(31)            !0-10-1

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4)/120d0 + zp**2*(-((pi**2)/24d0
     &) - 5d0/4d0 + (pi**2*llzp)/12d0 + (llzp)/2d0 + zeta3) +
     & zp*(-((pi**2)/6d0) + (pi**2*llzp)/6d0 + 2*zeta3) + zp*
     &*4*(-(1151d0/1728d0) - (pi**2)/96d0 + (pi**2*llzp)/24d0
     & + (49d0*llzp)/144d0 + (zeta3)/2d0) + zp**6*(-((pi**2)/
     &216d0) - 17653d0/40500d0 + (pi**2*llzp)/36d0 + (5269d0*
     &llzp)/21600d0 + (zeta3)/3d0) + zp**3*(-((pi**2)/54d0) -
     & 8d0/9d0 + (pi**2*llzp)/18d0 + (5d0*llzp)/12d0 + (2d0*z
     &eta3)/3d0) + zp**8*(-((pi**2)/384d0) - 127203607d0/3951
     &36000d0 + (266681d0*llzp)/1411200d0 + (pi**2*llzp)/48d0
     & + (zeta3)/4d0) + zp**5*(-((pi**2)/150d0) - 2281d0/4320
     &d0 + (pi**2*llzp)/30d0 + (41d0*llzp)/144d0 + (2d0*zeta3
     &)/5d0) + zp**7*(-((pi**2)/294d0) - 93371d0/252000d0 + (
     &pi**2*llzp)/42d0 + (767d0*llzp)/3600d0 + (2d0*zeta3)/7d
     &0) + zp**9*(-((pi**2)/486d0) - 2276013631d0/8001504000d
     &0 + (pi**2*llzp)/54d0 + (1077749d0*llzp)/6350400d0 + (2
     &d0*zeta3)/9d0)

         case(32)            !0-100

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = (pi**4)/20d0 + 2*myi*pi*szp*zeta3 + z
     &p*(-((pi**2)/2d0) - (myi*pi**3*szp)/6d0 + (pi**2*llzp)/
     &2d0 + zeta3) + zp**2*(-((pi**2)/8d0) - (myi*pi**3*szp)/
     &12d0 + (myi*pi*szp)/2d0 + (pi**2*llzp)/4d0 + (zeta3)/2d
     &0) + zp**3*(-(1d0/12d0) - (pi**2)/18d0 - (myi*pi**3*szp
     &)/18d0 + (myi*pi*5d0*szp)/12d0 + (pi**2*llzp)/6d0 + (ze
     &ta3)/3d0) + zp**4*(-((pi**2)/32d0) - 5d0/48d0 - (myi*pi
     &**3*szp)/24d0 + (myi*pi*49d0*szp)/144d0 + (pi**2*llzp)/
     &8d0 + (zeta3)/4d0) + zp**5*(-(17d0/160d0) - (pi**2)/50d
     &0 - (myi*pi**3*szp)/30d0 + (myi*pi*41d0*szp)/144d0 + (p
     &i**2*llzp)/10d0 + (zeta3)/5d0) + zp**6*(-(59d0/576d0) -
     & (pi**2)/72d0 - (myi*pi**3*szp)/36d0 + (myi*pi*5269d0*s
     &zp)/21600d0 + (pi**2*llzp)/12d0 + (zeta3)/6d0) + zp**7*
     &(-(2929d0/30240d0) - (pi**2)/98d0 - (myi*pi**3*szp)/42d
     &0 + (myi*pi*767d0*szp)/3600d0 + (pi**2*llzp)/14d0 + (ze
     &ta3)/7d0) + zp**8*(-((pi**2)/128d0) - 629d0/6912d0 + (m
     &yi*pi*266681d0*szp)/1411200d0 - (myi*pi**3*szp)/48d0 + 
     &(pi**2*llzp)/16d0 + (zeta3)/8d0) + zp**9*(-((pi**2)/162
     &d0) - 185921d0/2177280d0 - (myi*pi**3*szp)/54d0 + (myi*
     &pi*1077749d0*szp)/6350400d0 + (pi**2*llzp)/18d0 + (zeta
     &3)/9d0)

         case(33)            !0-101

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4*71d0)/1440d0 + (pi**2*ll2**2)/
     &6d0 - (ll2**4)/6d0 - 4*cli4pt5 - (7d0*ll2*zeta3)/2d0 + 
     &zp**2*(-((pi**2)/48d0) - (ll2)/2d0 + (pi**2*llzp)/24d0 
     &+ (5d0*zeta3)/16d0) + zp**3*(-((pi**2)/108d0) + 1d0/24d
     &0 - (5d0*ll2)/12d0 + (pi**2*llzp)/36d0 + (5d0*zeta3)/24
     &d0) + zp**4*(-((pi**2)/192d0) + 7d0/144d0 - (49d0*ll2)/
     &144d0 + (pi**2*llzp)/48d0 + (5d0*zeta3)/32d0) + zp**6*(
     &-((pi**2)/432d0) + 3793d0/86400d0 - (5269d0*ll2)/21600d
     &0 + (pi**2*llzp)/72d0 + (5d0*zeta3)/48d0) + zp**7*(4882
     &1d0/1209600d0 - (pi**2)/588d0 - (767d0*ll2)/3600d0 + (p
     &i**2*llzp)/84d0 + (5d0*zeta3)/56d0) + zp**8*(2511659d0/
     &67737600d0 - (pi**2)/768d0 - (266681d0*ll2)/1411200d0 +
     & (pi**2*llzp)/96d0 + (5d0*zeta3)/64d0) + zp**9*(1041298
     &1d0/304819200d0 - (pi**2)/972d0 - (1077749d0*ll2)/63504
     &00d0 + (pi**2*llzp)/108d0 + (5d0*zeta3)/72d0) + zp**5*(
     &-((pi**2)/300d0) + 17d0/360d0 - (41d0*ll2)/144d0 + (pi*
     &*2*llzp)/60d0 + (zeta3)/8d0) + zp*(-((pi**2)/12d0) + (p
     &i**2*llzp)/12d0 + (5d0*zeta3)/8d0)

         case(34)            !0-11-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4*7d0)/288d0) + (pi**2*ll2**2)
     &/8d0 + (ll2**4)/8d0 + 3*cli4pt5 + zp**7*((pi**2)/588d0 
     &+ 14001083d0/84672000d0 + (pi**2*ll2)/42d0 - (ll2**2)/9
     &8d0 - (ll2**3)/21d0 - (5587d0*llzp)/67200d0 - (pi**2*ll
     &zp)/84d0 + (ll2**2*llzp)/14d0 - (zeta3)/4d0) + zp**3*((
     &pi**2)/108d0 + 5d0/12d0 + (pi**2*ll2)/18d0 - (ll2**2)/1
     &8d0 - (ll2**3)/9d0 - (pi**2*llzp)/36d0 - (3d0*llzp)/16d
     &0 + (ll2**2*llzp)/6d0 - (7d0*zeta3)/12d0) + zp**4*((pi*
     &*2)/192d0 + 2101d0/6912d0 + (pi**2*ll2)/24d0 - (ll2**2)
     &/32d0 - (ll2**3)/12d0 - (pi**2*llzp)/48d0 - (83d0*llzp)
     &/576d0 + (ll2**2*llzp)/8d0 - (7d0*zeta3)/16d0) + zp**5*
     &((pi**2)/300d0 + 82237d0/345600d0 + (pi**2*ll2)/30d0 - 
     &(ll2**2)/50d0 - (ll2**3)/15d0 - (1337d0*llzp)/11520d0 -
     & (pi**2*llzp)/60d0 + (ll2**2*llzp)/10d0 - (7d0*zeta3)/2
     &0d0) + zp**6*((pi**2)/432d0 + 505931d0/2592000d0 + (pi*
     &*2*ll2)/36d0 - (ll2**2)/72d0 - (ll2**3)/18d0 - (33497d0
     &*llzp)/345600d0 - (pi**2*llzp)/72d0 + (ll2**2*llzp)/12d
     &0 - (7d0*zeta3)/24d0) + zp**8*(169983053d0/1185408000d0
     & + (pi**2)/768d0 + (pi**2*ll2)/48d0 - (ll2**2)/128d0 - 
     &(ll2**3)/24d0 - (136919d0*llzp)/1881600d0 - (pi**2*llzp
     &)/96d0 + (ll2**2*llzp)/16d0 - (7d0*zeta3)/32d0) + zp**9
     &*(86419598141d0/682795008000d0 + (pi**2)/972d0 + (pi**2
     &*ll2)/54d0 - (ll2**2)/162d0 - (ll2**3)/27d0 - (pi**2*ll
     &zp)/108d0 - (35054939d0*llzp)/541900800d0 + (ll2**2*llz
     &p)/18d0 - (7d0*zeta3)/36d0) + zp*((pi**2)/12d0 + (pi**2
     &*ll2)/6d0 - (ll2**2)/2d0 - (ll2**3)/3d0 - (pi**2*llzp)/
     &12d0 + (ll2**2*llzp)/2d0 - (7d0*zeta3)/4d0) + zp**2*((p
     &i**2)/48d0 + 5d0/8d0 + (pi**2*ll2)/12d0 - (ll2**2)/8d0 
     &- (ll2**3)/6d0 - (pi**2*llzp)/24d0 - (llzp)/4d0 + (ll2*
     &*2*llzp)/4d0 - (7d0*zeta3)/8d0)

         case(35)            !0-110

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = -((pi**4*11d0)/480d0) + (myi*pi**3*sz
     &p*ll2)/4d0 - (pi**2*ll2**2)/6d0 + (ll2**4)/6d0 + 4*cli4
     &pt5 - myi*pi*szp*zeta3 + zp**2*((pi**2)/48d0 + (myi*pi*
     &*3*szp)/24d0 - (myi*pi*szp)/4d0 + (pi**2*ll2)/8d0 - (my
     &i*pi*szp*ll2)/4d0 - (myi*pi*szp*ll2**2)/4d0 - (pi**2*ll
     &zp)/24d0 + (myi*pi*szp*ll2*llzp)/2d0 - (13d0*zeta3)/16d
     &0) + zp**3*((pi**2)/108d0 + 1d0/24d0 + (myi*pi**3*szp)/
     &36d0 - (myi*pi*3d0*szp)/16d0 + (pi**2*ll2)/12d0 - (myi*
     &pi*szp*ll2)/9d0 - (myi*pi*szp*ll2**2)/6d0 - (pi**2*llzp
     &)/36d0 + (myi*pi*szp*ll2*llzp)/3d0 - (13d0*zeta3)/24d0)
     & + zp**4*((pi**2)/192d0 + 13d0/288d0 + (myi*pi**3*szp)/
     &48d0 - (myi*pi*83d0*szp)/576d0 + (pi**2*ll2)/16d0 - (my
     &i*pi*szp*ll2)/16d0 - (myi*pi*szp*ll2**2)/8d0 - (pi**2*l
     &lzp)/48d0 + (myi*pi*szp*ll2*llzp)/4d0 - (13d0*zeta3)/32
     &d0) + zp**5*(119d0/2880d0 + (pi**2)/300d0 - (myi*pi*133
     &7d0*szp)/11520d0 + (myi*pi**3*szp)/60d0 + (pi**2*ll2)/2
     &0d0 - (myi*pi*szp*ll2)/25d0 - (myi*pi*szp*ll2**2)/10d0 
     &- (pi**2*llzp)/60d0 + (myi*pi*szp*ll2*llzp)/5d0 - (13d0
     &*zeta3)/40d0) + zp**6*((pi**2)/432d0 + 3167d0/86400d0 -
     & (myi*pi*33497d0*szp)/345600d0 + (myi*pi**3*szp)/72d0 +
     & (pi**2*ll2)/24d0 - (myi*pi*szp*ll2)/36d0 - (myi*pi*szp
     &*ll2**2)/12d0 - (pi**2*llzp)/72d0 + (myi*pi*szp*ll2*llz
     &p)/6d0 - (13d0*zeta3)/48d0) + zp**7*(1403d0/43200d0 + (
     &pi**2)/588d0 - (myi*pi*5587d0*szp)/67200d0 + (myi*pi**3
     &*szp)/84d0 + (pi**2*ll2)/28d0 - (myi*pi*szp*ll2)/49d0 -
     & (myi*pi*szp*ll2**2)/14d0 - (pi**2*llzp)/84d0 + (myi*pi
     &*szp*ll2*llzp)/7d0 - (13d0*zeta3)/56d0) + zp**8*(490589
     &d0/16934400d0 + (pi**2)/768d0 - (myi*pi*136919d0*szp)/1
     &881600d0 + (myi*pi**3*szp)/96d0 + (pi**2*ll2)/32d0 - (m
     &yi*pi*szp*ll2)/64d0 - (myi*pi*szp*ll2**2)/16d0 - (pi**2
     &*llzp)/96d0 + (myi*pi*szp*ll2*llzp)/8d0 - (13d0*zeta3)/
     &64d0) + zp**9*(3972277d0/152409600d0 + (pi**2)/972d0 + 
     &(myi*pi**3*szp)/108d0 - (myi*pi*35054939d0*szp)/5419008
     &00d0 + (pi**2*ll2)/36d0 - (myi*pi*szp*ll2)/81d0 - (myi*
     &pi*szp*ll2**2)/18d0 - (pi**2*llzp)/108d0 + (myi*pi*szp*
     &ll2*llzp)/9d0 - (13d0*zeta3)/72d0) + zp*((pi**2)/12d0 +
     & (myi*pi**3*szp)/12d0 + (pi**2*ll2)/4d0 - myi*pi*szp*ll
     &2 - (myi*pi*szp*ll2**2)/2d0 - (pi**2*llzp)/12d0 + myi*p
     &i*szp*ll2*llzp - (13d0*zeta3)/8d0)

         case(36)            !0-111

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4*7d0)/288d0) - (pi**2*5d0*ll2
     &**2)/24d0 + (ll2**4)/12d0 + 2*cli4pt5 + (21d0*ll2*zeta3
     &)/8d0 + zp**2*(-((pi**2*ll2)/24d0) + (ll2)/4d0 + (ll2**
     &2)/8d0 + (ll2**3)/6d0 - (ll2**2*llzp)/4d0 + (zeta3)/16d
     &0) + zp**3*(-(1d0/48d0) - (pi**2*ll2)/36d0 + (3d0*ll2)/
     &16d0 + (ll2**2)/18d0 + (ll2**3)/9d0 - (ll2**2*llzp)/6d0
     & + (zeta3)/24d0) + zp**4*(-(1d0/48d0) - (pi**2*ll2)/48d
     &0 + (83d0*ll2)/576d0 + (ll2**2)/32d0 + (ll2**3)/12d0 - 
     &(ll2**2*llzp)/8d0 + (zeta3)/32d0) + zp**5*(-(139d0/7680
     &d0) + (1337d0*ll2)/11520d0 - (pi**2*ll2)/60d0 + (ll2**2
     &)/50d0 + (ll2**3)/15d0 - (ll2**2*llzp)/10d0 + (zeta3)/4
     &0d0) + zp**6*(-(143d0/9216d0) + (33497d0*ll2)/345600d0 
     &- (pi**2*ll2)/72d0 + (ll2**2)/72d0 + (ll2**3)/18d0 - (l
     &l2**2*llzp)/12d0 + (zeta3)/48d0) + zp**7*(-(13007d0/967
     &680d0) + (5587d0*ll2)/67200d0 - (pi**2*ll2)/84d0 + (ll2
     &**2)/98d0 + (ll2**3)/21d0 - (ll2**2*llzp)/14d0 + (zeta3
     &)/56d0) + zp**8*(-(13061d0/1105920d0) + (136919d0*ll2)/
     &1881600d0 - (pi**2*ll2)/96d0 + (ll2**2)/128d0 + (ll2**3
     &)/24d0 - (ll2**2*llzp)/16d0 + (zeta3)/64d0) + zp**9*(-(
     &5861129d0/557383680d0) - (pi**2*ll2)/108d0 + (35054939d
     &0*ll2)/541900800d0 + (ll2**2)/162d0 + (ll2**3)/27d0 - (
     &ll2**2*llzp)/18d0 + (zeta3)/72d0) + zp*(-((pi**2*ll2)/1
     &2d0) + (ll2**2)/2d0 + (ll2**3)/3d0 - (ll2**2*llzp)/2d0 
     &+ (zeta3)/8d0)

         case(37)            !00-1-1

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4)/360d0 - zp*zeta3 + zp**2*(7d0
     &/8d0 - (3d0*llzp)/4d0 + (llzp**2)/4d0 - (zeta3)/2d0) + 
     &zp**3*(41d0/72d0 - (7d0*llzp)/12d0 + (llzp**2)/4d0 - (z
     &eta3)/3d0) + zp**4*(1397d0/3456d0 - (131d0*llzp)/288d0 
     &+ (11d0*llzp**2)/48d0 - (zeta3)/4d0) + zp**5*(2671d0/86
     &40d0 - (53d0*llzp)/144d0 + (5d0*llzp**2)/24d0 - (zeta3)
     &/5d0) + zp**6*(322493d0/1296000d0 - (2213d0*llzp)/7200d
     &0 + (137d0*llzp**2)/720d0 - (zeta3)/6d0) + zp**7*(10464
     &1d0/504000d0 - (947d0*llzp)/3600d0 + (7d0*llzp**2)/40d0
     & - (zeta3)/7d0) + zp**8*(140539517d0/790272000d0 - (647
     &707d0*llzp)/2822400d0 + (363d0*llzp**2)/2240d0 - (zeta3
     &)/8d0) + zp**9*(2486560891d0/16003008000d0 - (1290829d0
     &*llzp)/6350400d0 + (761d0*llzp**2)/5040d0 - (zeta3)/9d0
     &)

         case(38)            !00-10

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = (pi**4)/30d0 + zp*((myi*pi**3*szp)/6d
     &0 - 2*zeta3) + zp**2*((pi**2)/12d0 + (myi*pi**3*szp)/12
     &d0 - (myi*pi*3d0*szp)/4d0 + (myi*pi*szp*llzp)/2d0 - zet
     &a3) - myi*pi*szp*zeta3 + zp**4*((pi**2*11d0)/144d0 - 11
     &d0/48d0 + (myi*pi**3*szp)/24d0 - (myi*pi*131d0*szp)/288
     &d0 + (myi*pi*11d0*szp*llzp)/24d0 - (zeta3)/2d0) + zp**6
     &*((pi**2*137d0)/2160d0 - 37d0/144d0 + (myi*pi**3*szp)/3
     &6d0 - (myi*pi*2213d0*szp)/7200d0 + (myi*pi*137d0*szp*ll
     &zp)/360d0 - (zeta3)/3d0) + zp**3*((pi**2)/12d0 - 1d0/6d
     &0 + (myi*pi**3*szp)/18d0 - (myi*pi*7d0*szp)/12d0 + (myi
     &*pi*szp*llzp)/2d0 - (2d0*zeta3)/3d0) + zp**8*((pi**2*12
     &1d0)/2240d0 - 43171d0/172800d0 + (myi*pi**3*szp)/48d0 -
     & (myi*pi*647707d0*szp)/2822400d0 + (myi*pi*363d0*szp*ll
     &zp)/1120d0 - (zeta3)/4d0) + zp**5*(-(181d0/720d0) + (pi
     &**2*5d0)/72d0 + (myi*pi**3*szp)/30d0 - (myi*pi*53d0*szp
     &)/144d0 + (myi*pi*5d0*szp*llzp)/12d0 - (2d0*zeta3)/5d0)
     & + zp**7*(-(38569d0/151200d0) + (pi**2*7d0)/120d0 + (my
     &i*pi**3*szp)/42d0 - (myi*pi*947d0*szp)/3600d0 + (myi*pi
     &*7d0*szp*llzp)/20d0 - (2d0*zeta3)/7d0) + zp**9*((pi**2*
     &761d0)/15120d0 - 9261559d0/38102400d0 + (myi*pi**3*szp)
     &/54d0 - (myi*pi*1290829d0*szp)/6350400d0 + (myi*pi*761d
     &0*szp*llzp)/2520d0 - (2d0*zeta3)/9d0)

         case(39)            !00-11

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4*19d0)/1440d0) + (7d0*ll2*zet
     &a3)/4d0 + zp*(-((pi**2*ll2)/4d0) + zeta3) + zp**2*(-((p
     &i**2)/24d0) + (3d0*ll2)/4d0 - (pi**2*ll2)/8d0 + (ll2**2
     &)/4d0 - (ll2*llzp)/2d0 + (zeta3)/2d0) + zp**3*(1d0/12d0
     & - (pi**2)/24d0 - (pi**2*ll2)/12d0 + (7d0*ll2)/12d0 + (
     &ll2**2)/4d0 - (ll2*llzp)/2d0 + (zeta3)/3d0) + zp**4*(-(
     &(pi**2*11d0)/288d0) + 7d0/64d0 - (pi**2*ll2)/16d0 + (13
     &1d0*ll2)/288d0 + (11d0*ll2**2)/48d0 - (11d0*ll2*llzp)/2
     &4d0 + (zeta3)/4d0) + zp**5*(-((pi**2*5d0)/144d0) + 67d0
     &/576d0 - (pi**2*ll2)/20d0 + (53d0*ll2)/144d0 + (5d0*ll2
     &**2)/24d0 - (5d0*ll2*llzp)/12d0 + (zeta3)/5d0) + zp**6*
     &(-((pi**2*137d0)/4320d0) + 893d0/7680d0 - (pi**2*ll2)/2
     &4d0 + (2213d0*ll2)/7200d0 + (137d0*ll2**2)/720d0 - (137
     &d0*ll2*llzp)/360d0 + (zeta3)/6d0) + zp**7*(274607d0/241
     &9200d0 - (pi**2*7d0)/240d0 - (pi**2*ll2)/28d0 + (947d0*
     &ll2)/3600d0 + (7d0*ll2**2)/40d0 - (7d0*ll2*llzp)/20d0 +
     & (zeta3)/7d0) + zp**8*(2123381d0/19353600d0 - (pi**2*12
     &1d0)/4480d0 - (pi**2*ll2)/32d0 + (647707d0*ll2)/2822400
     &d0 + (363d0*ll2**2)/2240d0 - (363d0*ll2*llzp)/1120d0 + 
     &(zeta3)/8d0) + zp**9*(-((pi**2*761d0)/30240d0) + 804796
     &9d0/76204800d0 - (pi**2*ll2)/36d0 + (1290829d0*ll2)/635
     &0400d0 + (761d0*ll2**2)/5040d0 - (761d0*ll2*llzp)/2520d
     &0 + (zeta3)/9d0)

         case(40)            !000-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4)/90d0) + zp*zeta3 + zp**2*(-
     &((pi**2)/12d0) + (zeta3)/2d0) + zp**3*(-((pi**2)/12d0) 
     &+ 11d0/36d0 - (llzp)/6d0 + (zeta3)/3d0) + zp**4*(-((pi*
     &*2*11d0)/144d0) + 19d0/48d0 - (llzp)/4d0 + (zeta3)/4d0)
     & + zp**5*(599d0/1440d0 - (pi**2*5d0)/72d0 - (7d0*llzp)/
     &24d0 + (zeta3)/5d0) + zp**6*(-((pi**2*137d0)/2160d0) + 
     &79d0/192d0 - (5d0*llzp)/16d0 + (zeta3)/6d0) + zp**7*(-(
     &(pi**2*7d0)/120d0) + 3343d0/8400d0 - (29d0*llzp)/90d0 +
     & (zeta3)/7d0) + zp**8*(-((pi**2*121d0)/2240d0) + 21977d
     &0/57600d0 - (469d0*llzp)/1440d0 + (zeta3)/8d0) + zp**9*
     &(-((pi**2*761d0)/15120d0) + 83359739d0/228614400d0 - (2
     &9531d0*llzp)/90720d0 + (zeta3)/9d0)

         case(41)            !0000

            zp = x+1d0
            szp = s(zp)

            ris = (pi**4)/24d0 + (myi*pi**3*szp*zp)/6d0
     & + (-((pi**2)/4d0) + (myi*pi**3*szp)/12d0)*zp**2 + (-((
     &pi**2)/4d0) + (myi*pi**3*szp)/18d0 - (myi*pi*szp)/6d0)*
     &zp**3 + (1d0/24d0 - (pi**2*11d0)/48d0 + (myi*pi**3*szp)
     &/24d0 - (myi*pi*szp)/4d0)*zp**4 + (1d0/12d0 - (pi**2*5d
     &0)/24d0 + (myi*pi**3*szp)/30d0 - (myi*pi*7d0*szp)/24d0)
     &*zp**5 + (17d0/144d0 - (pi**2*137d0)/720d0 + (myi*pi**3
     &*szp)/36d0 - (myi*pi*5d0*szp)/16d0)*zp**6 + (-((pi**2*7
     &d0)/40d0) + 7d0/48d0 + (myi*pi**3*szp)/42d0 - (myi*pi*2
     &9d0*szp)/90d0)*zp**7 + (-((pi**2*363d0)/2240d0) + 967d0
     &/5760d0 - (myi*pi*469d0*szp)/1440d0 + (myi*pi**3*szp)/4
     &8d0)*zp**8 + (-((pi**2*761d0)/5040d0) + 89d0/480d0 + (m
     &yi*pi**3*szp)/54d0 - (myi*pi*29531d0*szp)/90720d0)*zp**
     &9

         case(42)            !0001

            zp = x+1d0

            ris = -((pi**4*7d0)/720d0) + (3d0*zp*zeta3)
     &/4d0 + zp**9*(-(112391d0/1451520d0) - (pi**2*761d0)/302
     &40d0 + (29531d0*ll2)/90720d0 + (zeta3)/12d0) + zp**4*(-
     &((pi**2*11d0)/288d0) - 1d0/48d0 + (ll2)/4d0 + (3d0*zeta
     &3)/16d0) + zp**5*(-(19d0/480d0) - (pi**2*5d0)/144d0 + (
     &7d0*ll2)/24d0 + (3d0*zeta3)/20d0) + zp**7*(-(2591d0/403
     &20d0) - (pi**2*7d0)/240d0 + (29d0*ll2)/90d0 + (3d0*zeta
     &3)/28d0) + zp**8*(-(23d0/320d0) - (pi**2*121d0)/4480d0 
     &+ (469d0*ll2)/1440d0 + (3d0*zeta3)/32d0) + zp**3*(-((pi
     &**2)/24d0) + (ll2)/6d0 + (zeta3)/4d0) + zp**6*(-((pi**2
     &*137d0)/4320d0) - 31d0/576d0 + (5d0*ll2)/16d0 + (zeta3)
     &/8d0) + zp**2*(-((pi**2)/24d0) + (3d0*zeta3)/8d0)

         case(43)            !001-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4)/180d0) - (pi**2*ll2**2)/12d
     &0 + (ll2**4)/12d0 + 2*cli4pt5 + zp**2*((pi**2)/24d0 + (
     &pi**2*ll2)/8d0 - (ll2**2)/4d0 - (13d0*zeta3)/16d0) + zp
     &**3*((pi**2)/24d0 - 11d0/72d0 + (pi**2*ll2)/12d0 - (ll2
     &**2)/4d0 + (llzp)/12d0 - (13d0*zeta3)/24d0) + zp**4*(-(
     &215d0/1152d0) + (pi**2*11d0)/288d0 + (pi**2*ll2)/16d0 -
     & (11d0*ll2**2)/48d0 + (11d0*llzp)/96d0 - (13d0*zeta3)/3
     &2d0) + zp**5*((pi**2*5d0)/144d0 - 181d0/960d0 + (pi**2*
     &ll2)/20d0 - (5d0*ll2**2)/24d0 + (llzp)/8d0 - (13d0*zeta
     &3)/40d0) + zp**6*(-(2321d0/12800d0) + (pi**2*137d0)/432
     &0d0 + (pi**2*ll2)/24d0 - (137d0*ll2**2)/720d0 + (731d0*
     &llzp)/5760d0 - (13d0*zeta3)/48d0) + zp**7*((pi**2*7d0)/
     &240d0 - 138503d0/806400d0 + (pi**2*ll2)/28d0 - (7d0*ll2
     &**2)/40d0 + (721d0*llzp)/5760d0 - (13d0*zeta3)/56d0) + 
     &zp**8*(-(182923d0/1128960d0) + (pi**2*121d0)/4480d0 + (
     &pi**2*ll2)/32d0 - (363d0*ll2**2)/2240d0 + (3931d0*llzp)
     &/32256d0 - (13d0*zeta3)/64d0) + zp**9*((pi**2*761d0)/30
     &240d0 - 139798291d0/914457600d0 + (pi**2*ll2)/36d0 - (7
     &61d0*ll2**2)/5040d0 + (42799d0*llzp)/362880d0 - (13d0*z
     &eta3)/72d0) + zp*((pi**2*ll2)/4d0 - (13d0*zeta3)/8d0)

         case(44)            !0010

            zp = x+1d0
            szp = s(zp)

            ris = (pi**4*7d0)/240d0 - (myi*pi*3d0*szp*z
     &eta3)/4d0 + zp**3*((pi**2)/24d0 + (myi*pi*szp)/12d0 + (
     &myi*pi**3*szp)/36d0 - (myi*pi*szp*ll2)/2d0 - (zeta3)/2d
     &0) + zp**5*((pi**2*5d0)/144d0 - 3d0/80d0 + (myi*pi**3*s
     &zp)/60d0 + (myi*pi*szp)/8d0 - (myi*pi*5d0*szp*ll2)/12d0
     & - (3d0*zeta3)/10d0) + zp**7*(-(187d0/3360d0) + (pi**2*
     &7d0)/240d0 + (myi*pi*721d0*szp)/5760d0 + (myi*pi**3*szp
     &)/84d0 - (myi*pi*7d0*szp*ll2)/20d0 - (3d0*zeta3)/14d0) 
     &+ zp**8*((pi**2*121d0)/4480d0 - 691d0/11520d0 + (myi*pi
     &*3931d0*szp)/32256d0 + (myi*pi**3*szp)/96d0 - (myi*pi*3
     &63d0*szp*ll2)/1120d0 - (3d0*zeta3)/16d0) + zp*((myi*pi*
     &*3*szp)/12d0 - (3d0*zeta3)/2d0) + zp**6*((pi**2*137d0)/
     &4320d0 - 7d0/144d0 + (myi*pi**3*szp)/72d0 + (myi*pi*731
     &d0*szp)/5760d0 - (myi*pi*137d0*szp*ll2)/360d0 - (zeta3)
     &/4d0) + zp**2*((pi**2)/24d0 + (myi*pi**3*szp)/24d0 - (m
     &yi*pi*szp*ll2)/2d0 - (3d0*zeta3)/4d0) + zp**9*(-(2521d0
     &/40320d0) + (pi**2*761d0)/30240d0 + (myi*pi**3*szp)/108
     &d0 + (myi*pi*42799d0*szp)/362880d0 - (myi*pi*761d0*szp*
     &ll2)/2520d0 - (zeta3)/6d0) + zp**4*((pi**2*11d0)/288d0 
     &- 1d0/48d0 + (myi*pi**3*szp)/48d0 + (myi*pi*11d0*szp)/9
     &6d0 - (myi*pi*11d0*szp*ll2)/24d0 - (3d0*zeta3)/8d0)

         case(45)            !0011

            zp = x+1d0

            ris = -((pi**4)/48d0) - (pi**2*ll2**2)/12d0
     & + (ll2**4)/12d0 + 2*cli4pt5 - (zp*zeta3)/8d0 + (7d0*ll
     &2*zeta3)/4d0 + zp**2*((ll2**2)/4d0 - (zeta3)/16d0) + zp
     &**3*(-((ll2)/12d0) + (ll2**2)/4d0 - (zeta3)/24d0) + zp*
     &*4*(1d0/96d0 - (11d0*ll2)/96d0 + (11d0*ll2**2)/48d0 - (
     &zeta3)/32d0) + zp**5*(17d0/960d0 - (ll2)/8d0 + (5d0*ll2
     &**2)/24d0 - (zeta3)/40d0) + zp**6*(253d0/11520d0 - (731
     &d0*ll2)/5760d0 + (137d0*ll2**2)/720d0 - (zeta3)/48d0) +
     & zp**7*(979d0/40320d0 - (721d0*ll2)/5760d0 + (7d0*ll2**
     &2)/40d0 - (zeta3)/56d0) + zp**8*(10943d0/430080d0 - (39
     &31d0*ll2)/32256d0 + (363d0*ll2**2)/2240d0 - (zeta3)/64d
     &0) + zp**9*(4703d0/181440d0 - (42799d0*ll2)/362880d0 + 
     &(761d0*ll2**2)/5040d0 - (zeta3)/72d0)

         case(46)            !01-1-1

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4*11d0)/720d0 - (ll2**4)/8d0 - 3
     &*cli4pt5 + zp**2*(-(7d0/16d0) - (pi**2*ll2)/24d0 + (ll2
     &**3)/12d0 + (3d0*llzp)/8d0 - (llzp**2)/8d0 + (7d0*zeta3
     &)/16d0) + zp**3*(-(227d0/864d0) - (pi**2*ll2)/36d0 + (l
     &l2**3)/18d0 + (37d0*llzp)/144d0 - (5d0*llzp**2)/48d0 + 
     &(7d0*zeta3)/24d0) + zp**4*(-(1247d0/6912d0) - (pi**2*ll
     &2)/48d0 + (ll2**3)/24d0 + (107d0*llzp)/576d0 - (llzp**2
     &)/12d0 + (7d0*zeta3)/32d0) + zp**5*(-(470159d0/3456000d
     &0) - (pi**2*ll2)/60d0 + (ll2**3)/30d0 + (8257d0*llzp)/5
     &7600d0 - (131d0*llzp**2)/1920d0 + (7d0*zeta3)/40d0) + z
     &p**6*(-(2257309d0/20736000d0) - (pi**2*ll2)/72d0 + (ll2
     &**3)/36d0 + (13369d0*llzp)/115200d0 - (661d0*llzp**2)/1
     &1520d0 + (7d0*zeta3)/48d0) + zp**8*(-(183970943d0/23708
     &16000d0) - (pi**2*ll2)/96d0 + (ll2**3)/48d0 + (314543d0
     &*llzp)/3763200d0 - (1163d0*llzp**2)/26880d0 + (7d0*zeta
     &3)/64d0) + zp**9*(-(833624776009d0/12290310144000d0) - 
     &(pi**2*ll2)/108d0 + (ll2**3)/54d0 + (357205771d0*llzp)/
     &4877107200d0 - (148969d0*llzp**2)/3870720d0 + (7d0*zeta
     &3)/72d0) + zp**7*(-(107435801d0/1185408000d0) - (pi**2*
     &ll2)/84d0 + (ll2**3)/42d0 + (953d0*llzp)/9800d0 - (1327
     &d0*llzp**2)/26880d0 + (zeta3)/8d0) + zp*(-((pi**2*ll2)/
     &12d0) + (ll2**3)/6d0 + (7d0*zeta3)/8d0)

         case(47)            !01-10

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = (pi**4*13d0)/1440d0 - (myi*pi**3*szp*
     &ll2)/4d0 + (pi**2*ll2**2)/6d0 - (ll2**4)/6d0 - 4*cli4pt
     &5 + (myi*pi*13d0*szp*zeta3)/8d0 + zp*(-((myi*pi**3*szp)
     &/12d0) - (pi**2*ll2)/12d0 + (myi*pi*szp*ll2**2)/2d0 + z
     &eta3) + zp**2*(-((pi**2)/24d0) - (myi*pi**3*szp)/24d0 +
     & (myi*pi*3d0*szp)/8d0 - (pi**2*ll2)/24d0 + (myi*pi*szp*
     &ll2**2)/4d0 - (myi*pi*szp*llzp)/4d0 + (zeta3)/2d0) + zp
     &**3*(1d0/12d0 - (pi**2*5d0)/144d0 - (myi*pi**3*szp)/36d
     &0 + (myi*pi*37d0*szp)/144d0 - (pi**2*ll2)/36d0 + (myi*p
     &i*szp*ll2**2)/6d0 - (myi*pi*5d0*szp*llzp)/24d0 + (zeta3
     &)/3d0) + zp**4*(-((pi**2)/36d0) + 3d0/32d0 - (myi*pi**3
     &*szp)/48d0 + (myi*pi*107d0*szp)/576d0 - (pi**2*ll2)/48d
     &0 + (myi*pi*szp*ll2**2)/8d0 - (myi*pi*szp*llzp)/6d0 + (
     &zeta3)/4d0) + zp**5*(251d0/2880d0 - (pi**2*131d0)/5760d
     &0 - (myi*pi**3*szp)/60d0 + (myi*pi*8257d0*szp)/57600d0 
     &- (pi**2*ll2)/60d0 + (myi*pi*szp*ll2**2)/10d0 - (myi*pi
     &*131d0*szp*llzp)/960d0 + (zeta3)/5d0) + zp**6*(1343d0/1
     &7280d0 - (pi**2*661d0)/34560d0 + (myi*pi*13369d0*szp)/1
     &15200d0 - (myi*pi**3*szp)/72d0 - (pi**2*ll2)/72d0 + (my
     &i*pi*szp*ll2**2)/12d0 - (myi*pi*661d0*szp*llzp)/5760d0 
     &+ (zeta3)/6d0) + zp**7*(2977d0/43200d0 - (pi**2*1327d0)
     &/80640d0 - (myi*pi**3*szp)/84d0 + (myi*pi*953d0*szp)/98
     &00d0 - (pi**2*ll2)/84d0 + (myi*pi*szp*ll2**2)/14d0 - (m
     &yi*pi*1327d0*szp*llzp)/13440d0 + (zeta3)/7d0) + zp**8*(
     &29711d0/483840d0 - (pi**2*1163d0)/80640d0 + (myi*pi*314
     &543d0*szp)/3763200d0 - (myi*pi**3*szp)/96d0 - (pi**2*ll
     &2)/96d0 + (myi*pi*szp*ll2**2)/16d0 - (myi*pi*1163d0*szp
     &*llzp)/13440d0 + (zeta3)/8d0) + zp**9*(-((pi**2*148969d
     &0)/11612160d0) + 8406389d0/152409600d0 - (myi*pi**3*szp
     &)/108d0 + (myi*pi*357205771d0*szp)/4877107200d0 - (pi**
     &2*ll2)/108d0 + (myi*pi*szp*ll2**2)/18d0 - (myi*pi*14896
     &9d0*szp*llzp)/1935360d0 + (zeta3)/9d0)

         case(48)            !01-11

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4*7d0)/720d0 + (pi**2*ll2**2)/4d
     &0 - (21d0*ll2*zeta3)/8d0 + zp**3*(-(1d0/24d0) + (pi**2*
     &5d0)/288d0 + (pi**2*ll2)/36d0 - (37d0*ll2)/144d0 - (5d0
     &*ll2**2)/48d0 - (ll2**3)/18d0 + (5d0*ll2*llzp)/24d0 - (
     &zeta3)/12d0) + zp**4*(-(17d0/384d0) + (pi**2)/72d0 + (p
     &i**2*ll2)/48d0 - (107d0*ll2)/576d0 - (ll2**2)/12d0 - (l
     &l2**3)/24d0 + (ll2*llzp)/6d0 - (zeta3)/16d0) + zp**5*((
     &pi**2*131d0)/11520d0 - 457d0/11520d0 + (pi**2*ll2)/60d0
     & - (8257d0*ll2)/57600d0 - (131d0*ll2**2)/1920d0 - (ll2*
     &*3)/30d0 + (131d0*ll2*llzp)/960d0 - (zeta3)/20d0) + zp*
     &*6*((pi**2*661d0)/69120d0 - 955d0/27648d0 - (13369d0*ll
     &2)/115200d0 + (pi**2*ll2)/72d0 - (661d0*ll2**2)/11520d0
     & - (ll2**3)/36d0 + (661d0*ll2*llzp)/5760d0 - (zeta3)/24
     &d0) + zp**7*((pi**2*1327d0)/161280d0 - 291769d0/9676800
     &d0 + (pi**2*ll2)/84d0 - (953d0*ll2)/9800d0 - (1327d0*ll
     &2**2)/26880d0 - (ll2**3)/42d0 + (1327d0*ll2*llzp)/13440
     &d0 - (zeta3)/28d0) + zp**8*((pi**2*1163d0)/161280d0 - 2
     &9407d0/1105920d0 - (314543d0*ll2)/3763200d0 + (pi**2*ll
     &2)/96d0 - (1163d0*ll2**2)/26880d0 - (ll2**3)/48d0 + (11
     &63d0*ll2*llzp)/13440d0 - (zeta3)/32d0) + zp**9*((pi**2*
     &148969d0)/23224320d0 - 231350923d0/9754214400d0 + (pi**
     &2*ll2)/108d0 - (357205771d0*ll2)/4877107200d0 - (148969
     &d0*ll2**2)/3870720d0 - (ll2**3)/54d0 + (148969d0*ll2*ll
     &zp)/1935360d0 - (zeta3)/36d0) + zp*((pi**2*ll2)/12d0 - 
     &(ll2**3)/6d0 - (zeta3)/4d0) + zp**2*((pi**2)/48d0 + (pi
     &**2*ll2)/24d0 - (3d0*ll2)/8d0 - (ll2**2)/8d0 - (ll2**3)
     &/12d0 + (ll2*llzp)/4d0 - (zeta3)/8d0)

         case(49)            !010-1

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4)/480d0 + zp**2*((pi**2)/24d0 -
     & (pi**2*ll2)/12d0 + (5d0*zeta3)/16d0) + zp**3*((pi**2*5
     &d0)/144d0 - 11d0/72d0 - (pi**2*ll2)/18d0 + (llzp)/12d0 
     &+ (5d0*zeta3)/24d0) + zp**4*((pi**2)/36d0 - 95d0/576d0 
     &- (pi**2*ll2)/24d0 + (5d0*llzp)/48d0 + (5d0*zeta3)/32d0
     &) + zp**6*((pi**2*661d0)/34560d0 - 941d0/7200d0 - (pi**
     &2*ll2)/36d0 + (47d0*llzp)/480d0 + (5d0*zeta3)/48d0) + z
     &p**7*(-(4d0/35d0) + (pi**2*1327d0)/80640d0 - (pi**2*ll2
     &)/42d0 + (13d0*llzp)/144d0 + (5d0*zeta3)/56d0) + zp**8*
     &(-(1137251d0/11289600d0) + (pi**2*1163d0)/80640d0 - (pi
     &**2*ll2)/48d0 + (3341d0*llzp)/40320d0 + (5d0*zeta3)/64d
     &0) + zp**9*((pi**2*148969d0)/11612160d0 - 20502379d0/22
     &8614400d0 - (pi**2*ll2)/54d0 + (13817d0*llzp)/181440d0 
     &+ (5d0*zeta3)/72d0) + zp**5*(-(43d0/288d0) + (pi**2*131
     &d0)/5760d0 - (pi**2*ll2)/30d0 + (5d0*llzp)/48d0 + (zeta
     &3)/8d0) + zp*(-((pi**2*ll2)/6d0) + (5d0*zeta3)/8d0)

         case(50)            !0100

            zp = x+1d0
            szp = s(zp)

            ris = (pi**4)/80d0 + (myi*pi*3d0*szp*zeta3)
     &/2d0 + zp**9*((pi**2*148969d0)/3870720d0 - 37d0/768d0 +
     & (myi*pi*13817d0*szp)/181440d0 - (myi*pi**3*szp)/108d0 
     &- (pi**2*ll2)/18d0 + (zeta3)/12d0) + zp**4*((pi**2)/12d
     &0 - 1d0/48d0 - (myi*pi**3*szp)/48d0 + (myi*pi*5d0*szp)/
     &48d0 - (pi**2*ll2)/8d0 + (3d0*zeta3)/16d0) + zp**5*((pi
     &**2*131d0)/1920d0 - 17d0/480d0 + (myi*pi*5d0*szp)/48d0 
     &- (myi*pi**3*szp)/60d0 - (pi**2*ll2)/10d0 + (3d0*zeta3)
     &/20d0) + zp**7*((pi**2*1327d0)/26880d0 - 95d0/2016d0 + 
     &(myi*pi*13d0*szp)/144d0 - (myi*pi**3*szp)/84d0 - (pi**2
     &*ll2)/14d0 + (3d0*zeta3)/28d0) + zp**8*((pi**2*1163d0)/
     &26880d0 - 557d0/11520d0 + (myi*pi*3341d0*szp)/40320d0 -
     & (myi*pi**3*szp)/96d0 - (pi**2*ll2)/16d0 + (3d0*zeta3)/
     &32d0) + zp**3*((pi**2*5d0)/48d0 + (myi*pi*szp)/12d0 - (
     &myi*pi**3*szp)/36d0 - (pi**2*ll2)/6d0 + (zeta3)/4d0) + 
     &zp*(-((myi*pi**3*szp)/12d0) - (pi**2*ll2)/2d0 + (3d0*ze
     &ta3)/4d0) + zp**6*(-(25d0/576d0) + (pi**2*661d0)/11520d
     &0 + (myi*pi*47d0*szp)/480d0 - (myi*pi**3*szp)/72d0 - (p
     &i**2*ll2)/12d0 + (zeta3)/8d0) + zp**2*((pi**2)/8d0 - (m
     &yi*pi**3*szp)/24d0 - (pi**2*ll2)/4d0 + (3d0*zeta3)/8d0)

         case(51)            !0101

            zp = x+1d0

            ris = (pi**4*13d0)/288d0 + (pi**2*ll2**2)/6
     &d0 - (ll2**4)/6d0 - 4*cli4pt5 - (7d0*ll2*zeta3)/2d0 + z
     &p**3*((pi**2*5d0)/288d0 - (ll2)/12d0 - (pi**2*ll2)/36d0
     & + (zeta3)/12d0) + zp**4*((pi**2)/72d0 + 1d0/96d0 - (pi
     &**2*ll2)/48d0 - (5d0*ll2)/48d0 + (zeta3)/16d0) + zp**5*
     &((pi**2*131d0)/11520d0 + 1d0/60d0 - (5d0*ll2)/48d0 - (p
     &i**2*ll2)/60d0 + (zeta3)/20d0) + zp**6*((pi**2*661d0)/6
     &9120d0 + 7d0/360d0 - (47d0*ll2)/480d0 - (pi**2*ll2)/72d
     &0 + (zeta3)/24d0) + zp**7*((pi**2*1327d0)/161280d0 + 10
     &9d0/5376d0 - (13d0*ll2)/144d0 - (pi**2*ll2)/84d0 + (zet
     &a3)/28d0) + zp**8*((pi**2*1163d0)/161280d0 + 12979d0/64
     &5120d0 - (3341d0*ll2)/40320d0 - (pi**2*ll2)/96d0 + (zet
     &a3)/32d0) + zp**9*((pi**2*148969d0)/23224320d0 + 56591d
     &0/2903040d0 - (13817d0*ll2)/181440d0 - (pi**2*ll2)/108d
     &0 + (zeta3)/36d0) + zp*(-((pi**2*ll2)/12d0) + (zeta3)/4
     &d0) + zp**2*((pi**2)/48d0 - (pi**2*ll2)/24d0 + (zeta3)/
     &8d0)

         case(52)            !011-1

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4)/80d0 - (pi**2*ll2**2)/24d0 - 
     &(ll2**4)/12d0 - 2*cli4pt5 + zp**2*(-((pi**2)/48d0) + (l
     &l2**2)/8d0 - (ll2**3)/12d0 + (zeta3)/16d0) + zp**3*(11d
     &0/144d0 - (pi**2*5d0)/288d0 + (5d0*ll2**2)/48d0 - (ll2*
     &*3)/18d0 - (llzp)/24d0 + (zeta3)/24d0) + zp**4*(-((pi**
     &2)/72d0) + 59d0/768d0 + (ll2**2)/12d0 - (ll2**3)/24d0 -
     & (3d0*llzp)/64d0 + (zeta3)/32d0) + zp**5*(-((pi**2*131d
     &0)/11520d0) + 7651d0/115200d0 + (131d0*ll2**2)/1920d0 -
     & (ll2**3)/30d0 - (83d0*llzp)/1920d0 + (zeta3)/40d0) + z
     &p**6*(65d0/1152d0 - (pi**2*661d0)/69120d0 + (661d0*ll2*
     &*2)/11520d0 - (ll2**3)/36d0 - (11d0*llzp)/288d0 + (zeta
     &3)/48d0) + zp**7*(-((pi**2*1327d0)/161280d0) + 1092631d
     &0/22579200d0 + (1327d0*ll2**2)/26880d0 - (ll2**3)/42d0 
     &- (5417d0*llzp)/161280d0 + (zeta3)/56d0) + zp**8*(-((pi
     &**2*1163d0)/161280d0) + 7763d0/184320d0 + (1163d0*ll2**
     &2)/26880d0 - (ll2**3)/48d0 - (137d0*llzp)/4608d0 + (zet
     &a3)/64d0) + zp**9*(-((pi**2*148969d0)/23224320d0) + 217
     &6291643d0/58525286400d0 + (148969d0*ll2**2)/3870720d0 -
     & (ll2**3)/54d0 - (617027d0*llzp)/23224320d0 + (zeta3)/7
     &2d0) + zp*(-((ll2**3)/6d0) + (zeta3)/8d0)

         case(53)            !0110

            zp = x+1d0
            szp = s(zp)

            ris = -((pi**4)/288d0) + (myi*pi*szp*zeta3)
     &/8d0 + zp**2*(-((pi**2)/48d0) + (pi**2*ll2)/24d0 + (myi
     &*pi*szp*ll2)/4d0 - (myi*pi*szp*ll2**2)/4d0 - (zeta3)/16
     &d0) + zp**3*(-((pi**2*5d0)/288d0) - (myi*pi*szp)/24d0 +
     & (pi**2*ll2)/36d0 + (myi*pi*5d0*szp*ll2)/24d0 - (myi*pi
     &*szp*ll2**2)/6d0 - (zeta3)/24d0) + zp**4*(-((pi**2)/72d
     &0) + 1d0/96d0 - (myi*pi*3d0*szp)/64d0 + (pi**2*ll2)/48d
     &0 + (myi*pi*szp*ll2)/6d0 - (myi*pi*szp*ll2**2)/8d0 - (z
     &eta3)/32d0) + zp**5*(-((pi**2*131d0)/11520d0) + 1d0/64d
     &0 - (myi*pi*83d0*szp)/1920d0 + (pi**2*ll2)/60d0 + (myi*
     &pi*131d0*szp*ll2)/960d0 - (myi*pi*szp*ll2**2)/10d0 - (z
     &eta3)/40d0) + zp**6*(11d0/640d0 - (pi**2*661d0)/69120d0
     & - (myi*pi*11d0*szp)/288d0 + (pi**2*ll2)/72d0 + (myi*pi
     &*661d0*szp*ll2)/5760d0 - (myi*pi*szp*ll2**2)/12d0 - (ze
     &ta3)/48d0) + zp**7*(-((pi**2*1327d0)/161280d0) + 49d0/2
     &880d0 - (myi*pi*5417d0*szp)/161280d0 + (pi**2*ll2)/84d0
     & + (myi*pi*1327d0*szp*ll2)/13440d0 - (myi*pi*szp*ll2**2
     &)/14d0 - (zeta3)/56d0) + zp**8*(-((pi**2*1163d0)/161280
     &d0) + 2603d0/161280d0 - (myi*pi*137d0*szp)/4608d0 + (pi
     &**2*ll2)/96d0 + (myi*pi*1163d0*szp*ll2)/13440d0 - (myi*
     &pi*szp*ll2**2)/16d0 - (zeta3)/64d0) + zp**9*(-((pi**2*1
     &48969d0)/23224320d0) + 809d0/53760d0 - (myi*pi*617027d0
     &*szp)/23224320d0 + (pi**2*ll2)/108d0 + (myi*pi*148969d0
     &*szp*ll2)/1935360d0 - (myi*pi*szp*ll2**2)/18d0 - (zeta3
     &)/72d0) + zp*((pi**2*ll2)/12d0 - (myi*pi*szp*ll2**2)/2d
     &0 - (zeta3)/8d0)

         case(54)            !0111

            zp = x+1d0

            ris = -((pi**4)/90d0) - (pi**2*ll2**2)/24d0
     & + (zp*ll2**3)/6d0 + (ll2**4)/24d0 + zp**2*(-((ll2**2)/
     &8d0) + (ll2**3)/12d0) + zp**3*((ll2)/24d0 - (5d0*ll2**2
     &)/48d0 + (ll2**3)/18d0) + zp**4*(-(1d0/192d0) + (3d0*ll
     &2)/64d0 - (ll2**2)/12d0 + (ll2**3)/24d0) + zp**5*(-(7d0
     &/960d0) + (83d0*ll2)/1920d0 - (131d0*ll2**2)/1920d0 + (
     &ll2**3)/30d0) + zp**6*(-(35d0/4608d0) + (11d0*ll2)/288d
     &0 - (661d0*ll2**2)/11520d0 + (ll2**3)/36d0) + zp**7*(-(
     &155d0/21504d0) + (5417d0*ll2)/161280d0 - (1327d0*ll2**2
     &)/26880d0 + (ll2**3)/42d0) + zp**8*(-(2441d0/368640d0) 
     &+ (137d0*ll2)/4608d0 - (1163d0*ll2**2)/26880d0 + (ll2**
     &3)/48d0) + zp**9*(-(19997d0/3317760d0) + (617027d0*ll2)
     &/23224320d0 - (148969d0*ll2**2)/3870720d0 + (ll2**3)/54
     &d0) + cli4pt5 + (7d0*ll2*zeta3)/8d0

         case(55)            !1-1-1-1

            zp = x+1d0
            llzp = log(zp)

            ris = zp**8*(-(1d0/1048576d0) + (llzp)/1310
     &72d0 - (llzp**2)/32768d0 + (llzp**3)/12288d0) + zp*(-(1
     &d0/2d0) + (llzp)/2d0 - (llzp**2)/4d0 + (llzp**3)/12d0) 
     &+ zp**3*(-(1d0/648d0) + (llzp)/216d0 - (llzp**2)/144d0 
     &+ (llzp**3)/144d0) + zp**6*(-(1d0/82944d0) + (llzp)/138
     &24d0 - (llzp**2)/4608d0 + (llzp**3)/2304d0) + zp**9*(-(
     &1d0/3359232d0) + (llzp)/373248d0 - (llzp**2)/82944d0 + 
     &(llzp**3)/27648d0) + zp**4*(-(1d0/4096d0) + (llzp)/1024
     &d0 - (llzp**2)/512d0 + (llzp**3)/384d0) + zp**2*(-(1d0/
     &64d0) + (llzp)/32d0 - (llzp**2)/32d0 + (llzp**3)/48d0) 
     &+ zp**7*(-(1d0/307328d0) + (llzp)/43904d0 - (llzp**2)/1
     &2544d0 + (llzp**3)/5376d0) + zp**5*(-(1d0/20000d0) + (l
     &lzp)/4000d0 - (llzp**2)/1600d0 + (llzp**3)/960d0) + cli
     &4pt5

         case(56)            !1-1-10

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = (pi**4)/720d0 + (myi*pi**3*szp*ll2)/1
     &2d0 - (myi*pi*szp*ll2**3)/6d0 + (ll2**4)/24d0 + cli4pt5
     & - (myi*pi*7d0*szp*zeta3)/8d0 - (ll2*zeta3)/8d0 + zp**5
     &*(-(329d0/17280d0) - (pi**2)/4800d0 + (myi*pi*szp)/4000
     &d0 + (pi**2*llzp)/960d0 - (myi*pi*szp*llzp)/800d0 + (my
     &i*pi*szp*llzp**2)/320d0 + (zeta3)/160d0) + zp**8*(-(183
     &91283d0/9483264000d0) - (pi**2)/98304d0 + (myi*pi*szp)/
     &131072d0 + (pi**2*llzp)/12288d0 - (myi*pi*szp*llzp)/163
     &84d0 + (myi*pi*szp*llzp**2)/4096d0 + (zeta3)/2048d0) + 
     &zp**3*(-((pi**2)/432d0) - 5d0/48d0 + (myi*pi*szp)/216d0
     & + (pi**2*llzp)/144d0 - (myi*pi*szp*llzp)/72d0 + (myi*p
     &i*szp*llzp**2)/48d0 + (zeta3)/24d0) + zp*(-((pi**2)/12d
     &0) + (myi*pi*szp)/2d0 + (pi**2*llzp)/12d0 - (myi*pi*szp
     &*llzp)/2d0 + (myi*pi*szp*llzp**2)/4d0 + (zeta3)/2d0) + 
     &zp**6*(-((pi**2)/13824d0) - 44581d0/5184000d0 + (myi*pi
     &*szp)/13824d0 + (pi**2*llzp)/2304d0 - (myi*pi*szp*llzp)
     &/2304d0 + (myi*pi*szp*llzp**2)/768d0 + (zeta3)/384d0) +
     & zp**9*(-(20706533d0/21337344000d0) - (pi**2)/248832d0 
     &+ (myi*pi*szp)/373248d0 + (pi**2*llzp)/27648d0 - (myi*p
     &i*szp*llzp)/41472d0 + (myi*pi*szp*llzp**2)/9216d0 + (ze
     &ta3)/4608d0) + zp**4*(-((pi**2)/1536d0) - 151d0/3456d0 
     &+ (myi*pi*szp)/1024d0 + (pi**2*llzp)/384d0 - (myi*pi*sz
     &p*llzp)/256d0 + (myi*pi*szp*llzp**2)/128d0 + (zeta3)/64
     &d0) + zp**7*(-((pi**2)/37632d0) - 48581d0/12096000d0 + 
     &(myi*pi*szp)/43904d0 + (pi**2*llzp)/5376d0 - (myi*pi*sz
     &p*llzp)/6272d0 + (myi*pi*szp*llzp**2)/1792d0 + (zeta3)/
     &896d0) + zp**2*(-(1d0/4d0) - (pi**2)/96d0 + (myi*pi*szp
     &)/32d0 + (pi**2*llzp)/48d0 - (myi*pi*szp*llzp)/16d0 + (
     &myi*pi*szp*llzp**2)/16d0 + (zeta3)/8d0)

         case(57)            !1-1-11

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4)/288d0) - (pi**2*ll2**2)/24d
     &0 + (ll2**4)/24d0 + (7d0*ll2*zeta3)/8d0 + zp**7*(4081d0
     &/3072000d0 + (pi**2)/75264d0 + (pi**2*ll2)/10752d0 - (l
     &l2)/43904d0 - (ll2**2)/12544d0 - (ll2**3)/5376d0 - (pi*
     &*2*llzp)/10752d0 + (ll2*llzp)/6272d0 + (ll2**2*llzp)/17
     &92d0 - (ll2*llzp**2)/1792d0 - (zeta3)/1024d0) + zp**5*(
     &407d0/55296d0 + (pi**2)/9600d0 + (pi**2*ll2)/1920d0 - (
     &ll2)/4000d0 - (ll2**2)/1600d0 - (ll2**3)/960d0 - (pi**2
     &*llzp)/1920d0 + (ll2*llzp)/800d0 + (ll2**2*llzp)/320d0 
     &- (ll2*llzp**2)/320d0 - (7d0*zeta3)/1280d0) + zp**8*((p
     &i**2)/196608d0 + 9822481d0/16859136000d0 - (ll2)/131072
     &d0 + (pi**2*ll2)/24576d0 - (ll2**2)/32768d0 - (ll2**3)/
     &12288d0 - (pi**2*llzp)/24576d0 + (ll2*llzp)/16384d0 + (
     &ll2**2*llzp)/4096d0 - (ll2*llzp**2)/4096d0 - (7d0*zeta3
     &)/16384d0) + zp*((pi**2)/24d0 + (pi**2*ll2)/24d0 - (ll2
     &)/2d0 - (ll2**2)/4d0 - (ll2**3)/12d0 - (pi**2*llzp)/24d
     &0 + (ll2*llzp)/2d0 + (ll2**2*llzp)/4d0 - (ll2*llzp**2)/
     &4d0 - (7d0*zeta3)/16d0) + zp**3*(3d0/64d0 + (pi**2)/864
     &d0 - (ll2)/216d0 + (pi**2*ll2)/288d0 - (ll2**2)/144d0 -
     & (ll2**3)/144d0 - (pi**2*llzp)/288d0 + (ll2*llzp)/72d0 
     &+ (ll2**2*llzp)/48d0 - (ll2*llzp**2)/48d0 - (7d0*zeta3)
     &/192d0) + zp**6*((pi**2)/27648d0 + 256103d0/82944000d0 
     &- (ll2)/13824d0 + (pi**2*ll2)/4608d0 - (ll2**2)/4608d0 
     &- (ll2**3)/2304d0 - (pi**2*llzp)/4608d0 + (ll2*llzp)/23
     &04d0 + (ll2**2*llzp)/768d0 - (ll2*llzp**2)/768d0 - (7d0
     &*zeta3)/3072d0) + zp**9*((pi**2)/497664d0 + 78708473d0/
     &303464448000d0 - (ll2)/373248d0 + (pi**2*ll2)/55296d0 -
     & (ll2**2)/82944d0 - (ll2**3)/27648d0 - (pi**2*llzp)/552
     &96d0 + (ll2*llzp)/41472d0 + (ll2**2*llzp)/9216d0 - (ll2
     &*llzp**2)/9216d0 - (7d0*zeta3)/36864d0) + zp**4*(251d0/
     &13824d0 + (pi**2)/3072d0 - (ll2)/1024d0 + (pi**2*ll2)/7
     &68d0 - (ll2**2)/512d0 - (ll2**3)/384d0 - (pi**2*llzp)/7
     &68d0 + (ll2*llzp)/256d0 + (ll2**2*llzp)/128d0 - (ll2*ll
     &zp**2)/128d0 - (7d0*zeta3)/512d0) + zp**2*((pi**2)/192d
     &0 + 1d0/8d0 - (ll2)/32d0 + (pi**2*ll2)/96d0 - (ll2**2)/
     &32d0 - (ll2**3)/48d0 - (pi**2*llzp)/96d0 + (ll2*llzp)/1
     &6d0 + (ll2**2*llzp)/16d0 - (ll2*llzp**2)/16d0 - (7d0*ze
     &ta3)/64d0)

         case(58)            !1-10-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4*19d0)/1440d0) + (pi**2*ll2**
     &2)/24d0 + (ll2**4)/24d0 + cli4pt5 + zp*((pi**2)/12d0 - 
     &(pi**2*llzp)/12d0 - zeta3) + (ll2*zeta3)/4d0 + zp**8*(1
     &0723549d0/2370816000d0 + (pi**2)/98304d0 - (pi**2*llzp)
     &/12288d0 - (9701d0*llzp)/1881600d0 - (zeta3)/1024d0) + 
     &zp**3*((pi**2)/432d0 + 1d0/4d0 - (pi**2*llzp)/144d0 - (
     &llzp)/8d0 - (zeta3)/12d0) + zp**6*((pi**2)/13824d0 + 51
     &521d0/2592000d0 - (pi**2*llzp)/2304d0 - (347d0*llzp)/21
     &600d0 - (zeta3)/192d0) + zp**9*(24451813d0/10668672000d
     &0 + (pi**2)/248832d0 - (pi**2*llzp)/27648d0 - (209d0*ll
     &zp)/66150d0 - (zeta3)/2304d0) + zp**4*((pi**2)/1536d0 +
     & 709d0/6912d0 - (pi**2*llzp)/384d0 - (35d0*llzp)/576d0 
     &- (zeta3)/32d0) + zp**7*((pi**2)/37632d0 + 393707d0/423
     &36000d0 - (149d0*llzp)/16800d0 - (pi**2*llzp)/5376d0 - 
     &(zeta3)/448d0) + zp**2*(5d0/8d0 + (pi**2)/96d0 - (pi**2
     &*llzp)/48d0 - (llzp)/4d0 - (zeta3)/4d0) + zp**5*(1909d0
     &/43200d0 + (pi**2)/4800d0 - (11d0*llzp)/360d0 - (pi**2*
     &llzp)/960d0 - (zeta3)/80d0)

         case(59)            !1-100

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = -((pi**4*41d0)/1440d0) + (myi*pi**3*s
     &zp*ll2)/12d0 + (pi**2*ll2**2)/4d0 - myi*pi*szp*zeta3 - 
     &(3d0*ll2*zeta3)/4d0 + zp**5*((pi**2)/1600d0 + 5d0/192d0
     & - (myi*pi*11d0*szp)/360d0 + (myi*pi**3*szp)/960d0 - (p
     &i**2*llzp)/320d0 - (zeta3)/160d0) + zp**8*((pi**2)/3276
     &8d0 + 4669d0/552960d0 + (myi*pi**3*szp)/12288d0 - (myi*
     &pi*9701d0*szp)/1881600d0 - (pi**2*llzp)/4096d0 - (zeta3
     &)/2048d0) + zp**3*((pi**2)/144d0 + 1d0/24d0 + (myi*pi**
     &3*szp)/144d0 - (myi*pi*szp)/8d0 - (pi**2*llzp)/48d0 - (
     &zeta3)/24d0) + zp*((pi**2)/4d0 + (myi*pi**3*szp)/12d0 -
     & (pi**2*llzp)/4d0 - (zeta3)/2d0) + zp**6*(41d0/2304d0 +
     & (pi**2)/4608d0 + (myi*pi**3*szp)/2304d0 - (myi*pi*347d
     &0*szp)/21600d0 - (pi**2*llzp)/768d0 - (zeta3)/384d0) + 
     &zp**9*(10457d0/1741824d0 + (pi**2)/82944d0 + (myi*pi**3
     &*szp)/27648d0 - (myi*pi*209d0*szp)/66150d0 - (pi**2*llz
     &p)/9216d0 - (zeta3)/4608d0) + zp**4*((pi**2)/512d0 + 7d
     &0/192d0 + (myi*pi**3*szp)/384d0 - (myi*pi*35d0*szp)/576
     &d0 - (pi**2*llzp)/128d0 - (zeta3)/64d0) + zp**7*((pi**2
     &)/12544d0 + 2941d0/241920d0 - (myi*pi*149d0*szp)/16800d
     &0 + (myi*pi**3*szp)/5376d0 - (pi**2*llzp)/1792d0 - (zet
     &a3)/896d0) + zp**2*((pi**2)/32d0 + (myi*pi**3*szp)/48d0
     & - (myi*pi*szp)/4d0 - (pi**2*llzp)/16d0 - (zeta3)/8d0)

         case(60)            !1-101

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4*11d0)/240d0) - (pi**2*ll2**2
     &)/8d0 + (ll2**4)/6d0 + 4*cli4pt5 + (13d0*ll2*zeta3)/4d0
     & + zp**5*(-(31d0/2880d0) + (pi**2)/9600d0 + (11d0*ll2)/
     &360d0 - (pi**2*llzp)/1920d0 - (zeta3)/256d0) + zp**8*((
     &pi**2)/196608d0 - 744197d0/270950400d0 + (9701d0*ll2)/1
     &881600d0 - (pi**2*llzp)/24576d0 - (5d0*zeta3)/16384d0) 
     &+ zp*((pi**2)/24d0 - (pi**2*llzp)/24d0 - (5d0*zeta3)/16
     &d0) + zp**3*(-(1d0/48d0) + (pi**2)/864d0 + (ll2)/8d0 - 
     &(pi**2*llzp)/288d0 - (5d0*zeta3)/192d0) + zp**6*((pi**2
     &)/27648d0 - 73d0/10800d0 + (347d0*ll2)/21600d0 - (pi**2
     &*llzp)/4608d0 - (5d0*zeta3)/3072d0) + zp**9*((pi**2)/49
     &7664d0 - 555271d0/304819200d0 + (209d0*ll2)/66150d0 - (
     &pi**2*llzp)/55296d0 - (5d0*zeta3)/36864d0) + zp**4*(-(1
     &9d0/1152d0) + (pi**2)/3072d0 + (35d0*ll2)/576d0 - (pi**
     &2*llzp)/768d0 - (5d0*zeta3)/512d0) + zp**2*((pi**2)/192
     &d0 + (ll2)/4d0 - (pi**2*llzp)/96d0 - (5d0*zeta3)/64d0) 
     &+ zp**7*(-(10313d0/2419200d0) + (pi**2)/75264d0 + (149d
     &0*ll2)/16800d0 - (pi**2*llzp)/10752d0 - (5d0*zeta3)/716
     &8d0)

         case(61)            !1-11-1

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4)/1440d0 - (pi**2*ll2**2)/24d0 
     &+ (ll2**4)/24d0 + (ll2*zeta3)/4d0 + zp**7*(-((pi**2)/75
     &264d0) - 93371d0/32256000d0 - (pi**2*ll2)/5376d0 + (ll2
     &**2)/12544d0 + (ll2**3)/2688d0 + (pi**2*llzp)/10752d0 +
     & (767d0*llzp)/460800d0 - (ll2**2*llzp)/1792d0 + (zeta3)
     &/512d0) + zp**6*(-(17653d0/2592000d0) - (pi**2)/27648d0
     & - (pi**2*ll2)/2304d0 + (ll2**2)/4608d0 + (ll2**3)/1152
     &d0 + (pi**2*llzp)/4608d0 + (5269d0*llzp)/1382400d0 - (l
     &l2**2*llzp)/768d0 + (7d0*zeta3)/1536d0) + zp**9*(-(2276
     &013631d0/4096770048000d0) - (pi**2)/497664d0 - (pi**2*l
     &l2)/27648d0 + (ll2**2)/82944d0 + (ll2**3)/13824d0 + (10
     &77749d0*llzp)/3251404800d0 + (pi**2*llzp)/55296d0 - (ll
     &2**2*llzp)/9216d0 + (7d0*zeta3)/18432d0) + zp**4*(-(115
     &1d0/27648d0) - (pi**2)/3072d0 - (pi**2*ll2)/384d0 + (ll
     &2**2)/512d0 + (ll2**3)/192d0 + (49d0*llzp)/2304d0 + (pi
     &**2*llzp)/768d0 - (ll2**2*llzp)/128d0 + (7d0*zeta3)/256
     &d0) + zp**2*(-((pi**2)/192d0) - 5d0/16d0 - (pi**2*ll2)/
     &48d0 + (ll2**2)/32d0 + (ll2**3)/24d0 + (llzp)/8d0 + (pi
     &**2*llzp)/96d0 - (ll2**2*llzp)/16d0 + (7d0*zeta3)/32d0)
     & + zp**5*(-(2281d0/138240d0) - (pi**2)/9600d0 - (pi**2*
     &ll2)/960d0 + (ll2**2)/1600d0 + (ll2**3)/480d0 + (pi**2*
     &llzp)/1920d0 + (41d0*llzp)/4608d0 - (ll2**2*llzp)/320d0
     & + (7d0*zeta3)/640d0) + zp**8*(-(127203607d0/1011548160
     &00d0) - (pi**2)/196608d0 - (pi**2*ll2)/12288d0 + (ll2**
     &2)/32768d0 + (ll2**3)/6144d0 + (pi**2*llzp)/24576d0 + (
     &266681d0*llzp)/361267200d0 - (ll2**2*llzp)/4096d0 + (7d
     &0*zeta3)/8192d0) + zp*(-((pi**2)/24d0) - (pi**2*ll2)/12
     &d0 + (ll2**2)/4d0 + (ll2**3)/6d0 + (pi**2*llzp)/24d0 - 
     &(ll2**2*llzp)/4d0 + (7d0*zeta3)/8d0) + zp**3*(-((pi**2)
     &/864d0) - 1d0/9d0 - (pi**2*ll2)/144d0 + (ll2**2)/144d0 
     &+ (ll2**3)/72d0 + (pi**2*llzp)/288d0 + (5d0*llzp)/96d0 
     &- (ll2**2*llzp)/48d0 + (7d0*zeta3)/96d0)

         case(62)            !1-110

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = -((pi**4)/360d0) - (myi*pi**3*szp*ll2
     &)/12d0 - (pi**2*ll2**2)/24d0 + (myi*pi*szp*ll2**3)/6d0 
     &+ (myi*pi*szp*zeta3)/4d0 + ll2*zeta3 + zp**5*(-(49d0/57
     &60d0) - (pi**2)/9600d0 - (myi*pi**3*szp)/1920d0 + (myi*
     &pi*41d0*szp)/4608d0 - (pi**2*ll2)/640d0 + (myi*pi*szp*l
     &l2)/800d0 + (myi*pi*szp*ll2**2)/320d0 + (pi**2*llzp)/19
     &20d0 - (myi*pi*szp*ll2*llzp)/160d0 + (13d0*zeta3)/1280d
     &0) + zp**8*(-((pi**2)/196608d0) - 374123d0/270950400d0 
     &- (myi*pi**3*szp)/24576d0 + (myi*pi*266681d0*szp)/36126
     &7200d0 - (pi**2*ll2)/8192d0 + (myi*pi*szp*ll2)/16384d0 
     &+ (myi*pi*szp*ll2**2)/4096d0 + (pi**2*llzp)/24576d0 - (
     &myi*pi*szp*ll2*llzp)/2048d0 + (13d0*zeta3)/16384d0) + z
     &p*(-((pi**2)/24d0) - (myi*pi**3*szp)/24d0 - (pi**2*ll2)
     &/8d0 + (myi*pi*szp*ll2)/2d0 + (myi*pi*szp*ll2**2)/4d0 +
     & (pi**2*llzp)/24d0 - (myi*pi*szp*ll2*llzp)/2d0 + (13d0*
     &zeta3)/16d0) + zp**3*(-(1d0/48d0) - (pi**2)/864d0 - (my
     &i*pi**3*szp)/288d0 + (myi*pi*5d0*szp)/96d0 - (pi**2*ll2
     &)/96d0 + (myi*pi*szp*ll2)/72d0 + (myi*pi*szp*ll2**2)/48
     &d0 + (pi**2*llzp)/288d0 - (myi*pi*szp*ll2*llzp)/24d0 + 
     &(13d0*zeta3)/192d0) + zp**6*(-((pi**2)/27648d0) - 1609d
     &0/345600d0 - (myi*pi**3*szp)/4608d0 + (myi*pi*5269d0*sz
     &p)/1382400d0 - (pi**2*ll2)/1536d0 + (myi*pi*szp*ll2)/23
     &04d0 + (myi*pi*szp*ll2**2)/768d0 + (pi**2*llzp)/4608d0 
     &- (myi*pi*szp*ll2*llzp)/384d0 + (13d0*zeta3)/3072d0) + 
     &zp**9*(-((pi**2)/497664d0) - 469253d0/609638400d0 + (my
     &i*pi*1077749d0*szp)/3251404800d0 - (myi*pi**3*szp)/5529
     &6d0 - (pi**2*ll2)/18432d0 + (myi*pi*szp*ll2)/41472d0 + 
     &(myi*pi*szp*ll2**2)/9216d0 + (pi**2*llzp)/55296d0 - (my
     &i*pi*szp*ll2*llzp)/4608d0 + (13d0*zeta3)/36864d0) + zp*
     &*4*(-(17d0/1152d0) - (pi**2)/3072d0 + (myi*pi*49d0*szp)
     &/2304d0 - (myi*pi**3*szp)/768d0 - (pi**2*ll2)/256d0 + (
     &myi*pi*szp*ll2)/256d0 + (myi*pi*szp*ll2**2)/128d0 + (pi
     &**2*llzp)/768d0 - (myi*pi*szp*ll2*llzp)/64d0 + (13d0*ze
     &ta3)/512d0) + zp**2*(-((pi**2)/192d0) + (myi*pi*szp)/8d
     &0 - (myi*pi**3*szp)/96d0 - (pi**2*ll2)/32d0 + (myi*pi*s
     &zp*ll2)/16d0 + (myi*pi*szp*ll2**2)/16d0 + (pi**2*llzp)/
     &96d0 - (myi*pi*szp*ll2*llzp)/8d0 + (13d0*zeta3)/64d0) +
     & zp**7*(-(6107d0/2419200d0) - (pi**2)/75264d0 - (myi*pi
     &**3*szp)/10752d0 + (myi*pi*767d0*szp)/460800d0 - (pi**2
     &*ll2)/3584d0 + (myi*pi*szp*ll2)/6272d0 + (myi*pi*szp*ll
     &2**2)/1792d0 + (pi**2*llzp)/10752d0 - (myi*pi*szp*ll2*l
     &lzp)/896d0 + (13d0*zeta3)/7168d0)

         case(63)            !1-111

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4)/30d0 + (pi**2*ll2**2)/6d0 - (
     &ll2**4)/6d0 - 3*cli4pt5 - (23d0*ll2*zeta3)/8d0 + zp**5*
     &(17d0/5120d0 + (pi**2*ll2)/1920d0 - (41d0*ll2)/4608d0 -
     & (ll2**2)/1600d0 - (ll2**3)/480d0 + (ll2**2*llzp)/320d0
     & - (zeta3)/1280d0) + zp**8*(629d0/1769472d0 + (pi**2*ll
     &2)/24576d0 - (266681d0*ll2)/361267200d0 - (ll2**2)/3276
     &8d0 - (ll2**3)/6144d0 + (ll2**2*llzp)/4096d0 - (zeta3)/
     &16384d0) + zp*((pi**2*ll2)/24d0 - (ll2**2)/4d0 - (ll2**
     &3)/6d0 + (ll2**2*llzp)/4d0 - (zeta3)/16d0) + zp**3*(1d0
     &/96d0 + (pi**2*ll2)/288d0 - (5d0*ll2)/96d0 - (ll2**2)/1
     &44d0 - (ll2**3)/72d0 + (ll2**2*llzp)/48d0 - (zeta3)/192
     &d0) + zp**6*(59d0/36864d0 + (pi**2*ll2)/4608d0 - (5269d
     &0*ll2)/1382400d0 - (ll2**2)/4608d0 - (ll2**3)/1152d0 + 
     &(ll2**2*llzp)/768d0 - (zeta3)/3072d0) + zp**9*(185921d0
     &/1114767360d0 - (1077749d0*ll2)/3251404800d0 + (pi**2*l
     &l2)/55296d0 - (ll2**2)/82944d0 - (ll2**3)/13824d0 + (ll
     &2**2*llzp)/9216d0 - (zeta3)/36864d0) + zp**4*(5d0/768d0
     & - (49d0*ll2)/2304d0 + (pi**2*ll2)/768d0 - (ll2**2)/512
     &d0 - (ll2**3)/192d0 + (ll2**2*llzp)/128d0 - (zeta3)/512
     &d0) + zp**2*(-((ll2)/8d0) + (pi**2*ll2)/96d0 - (ll2**2)
     &/32d0 - (ll2**3)/24d0 + (ll2**2*llzp)/16d0 - (zeta3)/64
     &d0) + zp**7*(2929d0/3870720d0 + (pi**2*ll2)/10752d0 - (
     &767d0*ll2)/460800d0 - (ll2**2)/12544d0 - (ll2**3)/2688d
     &0 + (ll2**2*llzp)/1792d0 - (zeta3)/7168d0)

         case(64)            !10-1-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4)/288d0) - (pi**2*ll2**2)/24d
     &0 + (ll2**4)/24d0 + cli4pt5 + (zp*zeta3)/2d0 - (ll2*zet
     &a3)/8d0 + zp**5*(-(12017d0/432000d0) + (79d0*llzp)/1800
     &d0 - (llzp**2)/30d0 + (zeta3)/160d0) + zp**8*(-(2783246
     &3d0/9483264000d0) + (7493d0*llzp)/940800d0 - (151d0*llz
     &p**2)/13440d0 + (zeta3)/2048d0) + zp**3*(-(71d0/432d0) 
     &+ (13d0*llzp)/72d0 - (llzp**2)/12d0 + (zeta3)/24d0) + z
     &p**6*(-(64861d0/5184000d0) + (169d0*llzp)/7200d0 - (llz
     &p**2)/45d0 + (zeta3)/384d0) + zp**9*(-(293914637d0/1920
     &36096000d0) + (3001d0*llzp)/595350d0 - (8d0*llzp**2)/94
     &5d0 + (zeta3)/4608d0) + zp**4*(-(113d0/1728d0) + (25d0*
     &llzp)/288d0 - (5d0*llzp**2)/96d0 + (zeta3)/64d0) + zp**
     &7*(-(3505829d0/592704000d0) + (521d0*llzp)/39200d0 - (1
     &3d0*llzp**2)/840d0 + (zeta3)/896d0) + zp**2*(-(7d0/16d0
     &) + (3d0*llzp)/8d0 - (llzp**2)/8d0 + (zeta3)/8d0)

         case(65)            !10-10

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = -((pi**4*17d0)/480d0) + (myi*pi**3*sz
     &p*ll2)/6d0 - (pi**2*ll2**2)/6d0 + (ll2**4)/6d0 + 4*cli4
     &pt5 - (myi*pi*5d0*szp*zeta3)/8d0 + (3d0*ll2*zeta3)/2d0 
     &+ zp*(-((myi*pi**3*szp)/12d0) + zeta3) + zp**8*(-((pi**
     &2*151d0)/40320d0) + 42371d0/1382400d0 - (myi*pi**3*szp)
     &/12288d0 + (myi*pi*7493d0*szp)/940800d0 - (myi*pi*151d0
     &*szp*llzp)/6720d0 + (zeta3)/1024d0) + zp**3*(1d0/12d0 -
     & (pi**2)/36d0 - (myi*pi**3*szp)/144d0 + (myi*pi*13d0*sz
     &p)/72d0 - (myi*pi*szp*llzp)/6d0 + (zeta3)/12d0) + zp**6
     &*(-((pi**2)/135d0) + 179d0/3456d0 - (myi*pi**3*szp)/230
     &4d0 + (myi*pi*169d0*szp)/7200d0 - (myi*pi*2d0*szp*llzp)
     &/45d0 + (zeta3)/192d0) + zp**9*(735253d0/30481920d0 - (
     &pi**2*8d0)/2835d0 - (myi*pi**3*szp)/27648d0 + (myi*pi*3
     &001d0*szp)/595350d0 - (myi*pi*16d0*szp*llzp)/945d0 + (z
     &eta3)/2304d0) + zp**4*(1d0/12d0 - (pi**2*5d0)/288d0 + (
     &myi*pi*25d0*szp)/288d0 - (myi*pi**3*szp)/384d0 - (myi*p
     &i*5d0*szp*llzp)/48d0 + (zeta3)/32d0) + zp**7*(-((pi**2*
     &13d0)/2520d0) + 23963d0/604800d0 + (myi*pi*521d0*szp)/3
     &9200d0 - (myi*pi**3*szp)/5376d0 - (myi*pi*13d0*szp*llzp
     &)/420d0 + (zeta3)/448d0) + zp**2*(-((pi**2)/24d0) - (my
     &i*pi**3*szp)/48d0 + (myi*pi*3d0*szp)/8d0 - (myi*pi*szp*
     &llzp)/4d0 + (zeta3)/4d0) + zp**5*(-((pi**2)/90d0) + 97d
     &0/1440d0 + (myi*pi*79d0*szp)/1800d0 - (myi*pi**3*szp)/9
     &60d0 - (myi*pi*szp*llzp)/15d0 + (zeta3)/80d0)

         case(66)            !10-11

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4*7d0)/180d0 - (pi**2*ll2**2)/12
     &d0 - (ll2**4)/6d0 - 4*cli4pt5 - (13d0*ll2*zeta3)/8d0 + 
     &zp**5*((pi**2)/180d0 - 173d0/5760d0 + (pi**2*ll2)/640d0
     & - (79d0*ll2)/1800d0 - (ll2**2)/30d0 + (ll2*llzp)/15d0 
     &- (zeta3)/160d0) + zp**8*(-(479389d0/38707200d0) + (pi*
     &*2*151d0)/80640d0 + (pi**2*ll2)/8192d0 - (7493d0*ll2)/9
     &40800d0 - (151d0*ll2**2)/13440d0 + (151d0*ll2*llzp)/672
     &0d0 - (zeta3)/2048d0) + zp**3*(-(1d0/24d0) + (pi**2)/72
     &d0 - (13d0*ll2)/72d0 + (pi**2*ll2)/96d0 - (ll2**2)/12d0
     & + (ll2*llzp)/6d0 - (zeta3)/24d0) + zp*((pi**2*ll2)/8d0
     & - (zeta3)/2d0) + zp**6*((pi**2)/270d0 - 3067d0/138240d
     &0 + (pi**2*ll2)/1536d0 - (169d0*ll2)/7200d0 - (ll2**2)/
     &45d0 + (2d0*ll2*llzp)/45d0 - (zeta3)/384d0) + zp**9*(-(
     &1164053d0/121927680d0) + (pi**2*4d0)/2835d0 + (pi**2*ll
     &2)/18432d0 - (3001d0*ll2)/595350d0 - (8d0*ll2**2)/945d0
     & + (16d0*ll2*llzp)/945d0 - (zeta3)/4608d0) + zp**4*(-(5
     &d0/128d0) + (pi**2*5d0)/576d0 + (pi**2*ll2)/256d0 - (25
     &d0*ll2)/288d0 - (5d0*ll2**2)/96d0 + (5d0*ll2*llzp)/48d0
     & - (zeta3)/64d0) + zp**7*(-(39751d0/2419200d0) + (pi**2
     &*13d0)/5040d0 + (pi**2*ll2)/3584d0 - (521d0*ll2)/39200d
     &0 - (13d0*ll2**2)/840d0 + (13d0*ll2*llzp)/420d0 - (zeta
     &3)/896d0) + zp**2*((pi**2)/48d0 + (pi**2*ll2)/32d0 - (3
     &d0*ll2)/8d0 - (ll2**2)/8d0 + (ll2*llzp)/4d0 - (zeta3)/8
     &d0)

         case(67)            !100-1

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4)/60d0 + (pi**2*ll2**2)/12d0 - 
     &(ll2**4)/12d0 - 2*cli4pt5 - (zp*zeta3)/2d0 - (3d0*ll2*z
     &eta3)/4d0 + zp**5*(-(317d0/2880d0) + (pi**2)/90d0 + (ll
     &zp)/12d0 - (zeta3)/160d0) + zp**8*((pi**2*151d0)/40320d
     &0 - 41419d0/921600d0 + (539d0*llzp)/11520d0 - (zeta3)/2
     &048d0) + zp**3*((pi**2)/36d0 - 11d0/72d0 + (llzp)/12d0 
     &- (zeta3)/24d0) + zp**6*((pi**2)/135d0 - 187d0/2304d0 +
     & (5d0*llzp)/72d0 - (zeta3)/384d0) + zp**9*(-(6297983d0/
     &182891520d0) + (pi**2*8d0)/2835d0 + (22d0*llzp)/567d0 -
     & (zeta3)/4608d0) + zp**4*(-(55d0/384d0) + (pi**2*5d0)/2
     &88d0 + (3d0*llzp)/32d0 - (zeta3)/64d0) + zp**7*((pi**2*
     &13d0)/2520d0 - 3451d0/57600d0 + (41d0*llzp)/720d0 - (ze
     &ta3)/896d0) + zp**2*((pi**2)/24d0 - (zeta3)/8d0)

         case(68)            !1000

            zp = x+1d0
            szp = s(zp)

            ris = -((pi**4*23d0)/720d0) - (myi*pi**3*sz
     &p*zp)/12d0 + ((pi**2)/8d0 - (myi*pi**3*szp)/48d0)*zp**2
     & + ((pi**2)/12d0 + (myi*pi*szp)/12d0 - (myi*pi**3*szp)/
     &144d0)*zp**3 + (-(1d0/48d0) + (pi**2*5d0)/96d0 - (myi*p
     &i**3*szp)/384d0 + (myi*pi*3d0*szp)/32d0)*zp**4 + (-(1d0
     &/30d0) + (pi**2)/30d0 + (myi*pi*szp)/12d0 - (myi*pi**3*
     &szp)/960d0)*zp**5 + (-(11d0/288d0) + (pi**2)/45d0 - (my
     &i*pi**3*szp)/2304d0 + (myi*pi*5d0*szp)/72d0)*zp**6 + (-
     &(13d0/336d0) + (pi**2*13d0)/840d0 - (myi*pi**3*szp)/537
     &6d0 + (myi*pi*41d0*szp)/720d0)*zp**7 + ((pi**2*151d0)/1
     &3440d0 - 427d0/11520d0 - (myi*pi**3*szp)/12288d0 + (myi
     &*pi*539d0*szp)/11520d0)*zp**8 + (-(14d0/405d0) + (pi**2
     &*8d0)/945d0 - (myi*pi**3*szp)/27648d0 + (myi*pi*22d0*sz
     &p)/567d0)*zp**9 + (myi*pi**3*szp*ll2)/6d0 - (myi*pi*3d0
     &*szp*zeta3)/4d0

         case(69)            !1001

            zp = x+1d0

            ris = -((pi**4)/288d0) - (3d0*zp*zeta3)/8d0
     & + (3d0*ll2*zeta3)/4d0 + zp**3*((pi**2)/72d0 - (ll2)/12
     &d0 - (zeta3)/32d0) + zp**4*((pi**2*5d0)/576d0 + 1d0/96d
     &0 - (3d0*ll2)/32d0 - (3d0*zeta3)/256d0) + zp**2*((pi**2
     &)/48d0 - (3d0*zeta3)/32d0) + zp**7*(47d0/2880d0 + (pi**
     &2*13d0)/5040d0 - (41d0*ll2)/720d0 - (3d0*zeta3)/3584d0)
     & + zp**6*((pi**2)/270d0 + 13d0/768d0 - (5d0*ll2)/72d0 -
     & (zeta3)/512d0) + zp**9*(481d0/35840d0 + (pi**2*4d0)/28
     &35d0 - (22d0*ll2)/567d0 - (zeta3)/6144d0) + zp**5*((pi*
     &*2)/180d0 + 1d0/64d0 - (ll2)/12d0 - (3d0*zeta3)/640d0) 
     &+ zp**8*((pi**2*151d0)/80640d0 + 1379d0/92160d0 - (539d
     &0*ll2)/11520d0 - (3d0*zeta3)/8192d0)

         case(70)            !101-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4*5d0)/144d0) + (pi**2*ll2**2)
     &/12d0 + (ll2**4)/6d0 + 4*cli4pt5 + ll2*zeta3 + zp**5*(-
     &((pi**2)/180d0) + 1367d0/28800d0 - (pi**2*ll2)/640d0 + 
     &(ll2**2)/30d0 - (llzp)/30d0 + (13d0*zeta3)/1280d0) + zp
     &**8*(306053d0/18063360d0 - (pi**2*151d0)/80640d0 - (pi*
     &*2*ll2)/8192d0 + (151d0*ll2**2)/13440d0 - (935d0*llzp)/
     &64512d0 + (13d0*zeta3)/16384d0) + zp*(-((pi**2*ll2)/8d0
     &) + (13d0*zeta3)/16d0) + zp**3*(11d0/144d0 - (pi**2)/72
     &d0 - (pi**2*ll2)/96d0 + (ll2**2)/12d0 - (llzp)/24d0 + (
     &13d0*zeta3)/192d0) + zp**6*(-((pi**2)/270d0) + 7639d0/2
     &30400d0 - (pi**2*ll2)/1536d0 + (ll2**2)/45d0 - (97d0*ll
     &zp)/3840d0 + (13d0*zeta3)/3072d0) + zp**9*(23078341d0/1
     &828915200d0 - (pi**2*4d0)/2835d0 - (pi**2*ll2)/18432d0 
     &+ (8d0*ll2**2)/945d0 - (2041d0*llzp)/181440d0 + (13d0*z
     &eta3)/36864d0) + zp**4*(19d0/288d0 - (pi**2*5d0)/576d0 
     &- (pi**2*ll2)/256d0 + (5d0*ll2**2)/96d0 - (llzp)/24d0 +
     & (13d0*zeta3)/512d0) + zp**2*(-((pi**2)/48d0) - (pi**2*
     &ll2)/32d0 + (ll2**2)/8d0 + (13d0*zeta3)/64d0) + zp**7*(
     &3671d0/156800d0 - (pi**2*13d0)/5040d0 - (pi**2*ll2)/358
     &4d0 + (13d0*ll2**2)/840d0 - (767d0*llzp)/40320d0 + (13d
     &0*zeta3)/7168d0)

         case(71)            !1010

            zp = x+1d0
            szp = s(zp)

            ris = -((pi**4*11d0)/288d0) + (myi*pi**3*sz
     &p*ll2)/12d0 - (pi**2*ll2**2)/6d0 + (ll2**4)/6d0 + 4*cli
     &4pt5 - (myi*pi*szp*zeta3)/4d0 + 2*ll2*zeta3 + zp**3*(-(
     &(pi**2)/72d0) - (myi*pi*szp)/24d0 - (myi*pi**3*szp)/288
     &d0 + (myi*pi*szp*ll2)/6d0 + (zeta3)/16d0) + zp**6*(17d0
     &/1152d0 - (pi**2)/270d0 - (myi*pi**3*szp)/4608d0 - (myi
     &*pi*97d0*szp)/3840d0 + (myi*pi*2d0*szp*ll2)/45d0 + (zet
     &a3)/256d0) + zp**9*(14081d0/1451520d0 - (pi**2*4d0)/283
     &5d0 - (myi*pi*2041d0*szp)/181440d0 - (myi*pi**3*szp)/55
     &296d0 + (myi*pi*16d0*szp*ll2)/945d0 + (zeta3)/3072d0) +
     & zp**4*(-((pi**2*5d0)/576d0) + 1d0/96d0 - (myi*pi*szp)/
     &24d0 - (myi*pi**3*szp)/768d0 + (myi*pi*5d0*szp*ll2)/48d
     &0 + (3d0*zeta3)/128d0) + zp**2*(-((pi**2)/48d0) - (myi*
     &pi**3*szp)/96d0 + (myi*pi*szp*ll2)/4d0 + (3d0*zeta3)/16
     &d0) + zp**7*(179d0/13440d0 - (pi**2*13d0)/5040d0 - (myi
     &*pi**3*szp)/10752d0 - (myi*pi*767d0*szp)/40320d0 + (myi
     &*pi*13d0*szp*ll2)/420d0 + (3d0*zeta3)/1792d0) + zp**5*(
     &-((pi**2)/180d0) + 7d0/480d0 - (myi*pi**3*szp)/1920d0 -
     & (myi*pi*szp)/30d0 + (myi*pi*szp*ll2)/15d0 + (3d0*zeta3
     &)/320d0) + zp**8*(-((pi**2*151d0)/80640d0) + 1057d0/921
     &60d0 - (myi*pi**3*szp)/24576d0 - (myi*pi*935d0*szp)/645
     &12d0 + (myi*pi*151d0*szp*ll2)/6720d0 + (3d0*zeta3)/4096
     &d0) + zp*(-((myi*pi**3*szp)/24d0) + (3d0*zeta3)/4d0)

         case(72)            !1011

            zp = x+1d0

            ris = (pi**4)/30d0 + (pi**2*ll2**2)/8d0 - (
     &ll2**4)/8d0 - 3*cli4pt5 + (zp*zeta3)/16d0 - (11d0*ll2*z
     &eta3)/4d0 + zp**5*(-(13d0/1920d0) + (ll2)/30d0 - (ll2**
     &2)/30d0 + (zeta3)/1280d0) + zp**8*(-(2321d0/516096d0) +
     & (935d0*ll2)/64512d0 - (151d0*ll2**2)/13440d0 + (zeta3)
     &/16384d0) + zp**3*((ll2)/24d0 - (ll2**2)/12d0 + (zeta3)
     &/192d0) + zp**6*(-(37d0/5760d0) + (97d0*ll2)/3840d0 - (
     &ll2**2)/45d0 + (zeta3)/3072d0) + zp**9*(-(157d0/43008d0
     &) + (2041d0*ll2)/181440d0 - (8d0*ll2**2)/945d0 + (zeta3
     &)/36864d0) + zp**4*(-(1d0/192d0) + (ll2)/24d0 - (5d0*ll
     &2**2)/96d0 + (zeta3)/512d0) + zp**2*(-((ll2**2)/8d0) + 
     &(zeta3)/64d0) + zp**7*(-(221d0/40320d0) + (767d0*ll2)/4
     &0320d0 - (13d0*ll2**2)/840d0 + (zeta3)/7168d0)

         case(73)            !11-1-1

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4)/720d0 + (ll2**4)/24d0 - (ll2*
     &zeta3)/8d0 + zp**7*(104641d0/64512000d0 + (pi**2*ll2)/1
     &0752d0 - (ll2**3)/5376d0 - (947d0*llzp)/460800d0 + (7d0
     &*llzp**2)/5120d0 - (zeta3)/1024d0) + zp**5*(2671d0/2764
     &80d0 + (pi**2*ll2)/1920d0 - (ll2**3)/960d0 - (53d0*llzp
     &)/4608d0 + (5d0*llzp**2)/768d0 - (7d0*zeta3)/1280d0) + 
     &zp**8*(140539517d0/202309632000d0 + (pi**2*ll2)/24576d0
     & - (ll2**3)/12288d0 - (647707d0*llzp)/722534400d0 + (36
     &3d0*llzp**2)/573440d0 - (7d0*zeta3)/16384d0) + zp*((pi*
     &*2*ll2)/24d0 - (ll2**3)/12d0 - (7d0*zeta3)/16d0) + zp**
     &3*(41d0/576d0 + (pi**2*ll2)/288d0 - (ll2**3)/144d0 - (7
     &d0*llzp)/96d0 + (llzp**2)/32d0 - (7d0*zeta3)/192d0) + z
     &p**6*(322493d0/82944000d0 + (pi**2*ll2)/4608d0 - (ll2**
     &3)/2304d0 - (2213d0*llzp)/460800d0 + (137d0*llzp**2)/46
     &080d0 - (7d0*zeta3)/3072d0) + zp**9*(2486560891d0/81935
     &40096000d0 + (pi**2*ll2)/55296d0 - (ll2**3)/27648d0 - (
     &1290829d0*llzp)/3251404800d0 + (761d0*llzp**2)/2580480d
     &0 - (7d0*zeta3)/36864d0) + zp**4*(1397d0/55296d0 + (pi*
     &*2*ll2)/768d0 - (ll2**3)/384d0 - (131d0*llzp)/4608d0 + 
     &(11d0*llzp**2)/768d0 - (7d0*zeta3)/512d0) + zp**2*(7d0/
     &32d0 + (pi**2*ll2)/96d0 - (ll2**3)/48d0 - (3d0*llzp)/16
     &d0 + (llzp**2)/16d0 - (7d0*zeta3)/64d0)

         case(74)            !11-10

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = (pi**4*7d0)/288d0 + (pi**2*ll2**2)/24
     &d0 + (myi*pi*szp*ll2**3)/6d0 - (ll2**4)/12d0 - 2*cli4pt
     &5 - (myi*pi*szp*zeta3)/8d0 - (13d0*ll2*zeta3)/8d0 + zp*
     &*5*(-(107d0/5760d0) + (pi**2*5d0)/2304d0 + (myi*pi**3*s
     &zp)/1920d0 - (myi*pi*53d0*szp)/4608d0 + (pi**2*ll2)/192
     &0d0 - (myi*pi*szp*ll2**2)/320d0 + (myi*pi*5d0*szp*llzp)
     &/384d0 - (zeta3)/160d0) + zp**8*(-(115543d0/38707200d0)
     & + (pi**2*121d0)/573440d0 + (myi*pi**3*szp)/24576d0 - (
     &myi*pi*647707d0*szp)/722534400d0 + (pi**2*ll2)/24576d0 
     &- (myi*pi*szp*ll2**2)/4096d0 + (myi*pi*363d0*szp*llzp)/
     &286720d0 - (zeta3)/2048d0) + zp**3*(-(1d0/24d0) + (pi**
     &2)/96d0 + (myi*pi**3*szp)/288d0 - (myi*pi*7d0*szp)/96d0
     & + (pi**2*ll2)/288d0 - (myi*pi*szp*ll2**2)/48d0 + (myi*
     &pi*szp*llzp)/16d0 - (zeta3)/24d0) + zp*((myi*pi**3*szp)
     &/24d0 + (pi**2*ll2)/24d0 - (myi*pi*szp*ll2**2)/4d0 - (z
     &eta3)/2d0) + zp**6*((pi**2*137d0)/138240d0 - 79d0/7680d
     &0 - (myi*pi*2213d0*szp)/460800d0 + (myi*pi**3*szp)/4608
     &d0 + (pi**2*ll2)/4608d0 - (myi*pi*szp*ll2**2)/768d0 + (
     &myi*pi*137d0*szp*llzp)/23040d0 - (zeta3)/384d0) + zp**9
     &*((pi**2*761d0)/7741440d0 - 983419d0/609638400d0 - (myi
     &*pi*1290829d0*szp)/3251404800d0 + (myi*pi**3*szp)/55296
     &d0 + (pi**2*ll2)/55296d0 - (myi*pi*szp*ll2**2)/9216d0 +
     & (myi*pi*761d0*szp*llzp)/1290240d0 - (zeta3)/4608d0) + 
     &zp**4*((pi**2*11d0)/2304d0 - 1d0/32d0 - (myi*pi*131d0*s
     &zp)/4608d0 + (myi*pi**3*szp)/768d0 + (pi**2*ll2)/768d0 
     &- (myi*pi*szp*ll2**2)/128d0 + (myi*pi*11d0*szp*llzp)/38
     &4d0 - (zeta3)/64d0) + zp**7*(-(13441d0/2419200d0) + (pi
     &**2*7d0)/15360d0 + (myi*pi**3*szp)/10752d0 - (myi*pi*94
     &7d0*szp)/460800d0 + (pi**2*ll2)/10752d0 - (myi*pi*szp*l
     &l2**2)/1792d0 + (myi*pi*7d0*szp*llzp)/2560d0 - (zeta3)/
     &896d0) + zp**2*((pi**2)/48d0 - (myi*pi*3d0*szp)/16d0 + 
     &(myi*pi**3*szp)/96d0 + (pi**2*ll2)/96d0 - (myi*pi*szp*l
     &l2**2)/16d0 + (myi*pi*szp*llzp)/8d0 - (zeta3)/8d0)

         case(75)            !11-11

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4)/30d0) - (pi**2*ll2**2)/8d0 
     &+ (ll2**4)/12d0 + 3*cli4pt5 + (11d0*ll2*zeta3)/4d0 + zp
     &**6*(-((pi**2*137d0)/276480d0) + 37d0/9216d0 + (2213d0*
     &ll2)/460800d0 - (pi**2*ll2)/4608d0 + (137d0*ll2**2)/460
     &80d0 + (ll2**3)/2304d0 - (137d0*ll2*llzp)/23040d0 + (ze
     &ta3)/1536d0) + zp**9*(-((pi**2*761d0)/15482880d0) + 926
     &1559d0/19508428800d0 + (1290829d0*ll2)/3251404800d0 - (
     &pi**2*ll2)/55296d0 + (761d0*ll2**2)/2580480d0 + (ll2**3
     &)/27648d0 - (761d0*ll2*llzp)/1290240d0 + (zeta3)/18432d
     &0) + zp**4*(-((pi**2*11d0)/4608d0) + 11d0/768d0 + (131d
     &0*ll2)/4608d0 - (pi**2*ll2)/768d0 + (11d0*ll2**2)/768d0
     & + (ll2**3)/384d0 - (11d0*ll2*llzp)/384d0 + (zeta3)/256
     &d0) + zp**2*(-((pi**2)/96d0) + (3d0*ll2)/16d0 - (pi**2*
     &ll2)/96d0 + (ll2**2)/16d0 + (ll2**3)/48d0 - (ll2*llzp)/
     &8d0 + (zeta3)/32d0) + zp**7*(38569d0/19353600d0 - (pi**
     &2*7d0)/30720d0 - (pi**2*ll2)/10752d0 + (947d0*ll2)/4608
     &00d0 + (7d0*ll2**2)/5120d0 + (ll2**3)/5376d0 - (7d0*ll2
     &*llzp)/2560d0 + (zeta3)/3584d0) + zp**5*(181d0/23040d0 
     &- (pi**2*5d0)/4608d0 - (pi**2*ll2)/1920d0 + (53d0*ll2)/
     &4608d0 + (5d0*ll2**2)/768d0 + (ll2**3)/960d0 - (5d0*ll2
     &*llzp)/384d0 + (zeta3)/640d0) + zp**8*(-((pi**2*121d0)/
     &1146880d0) + 43171d0/44236800d0 - (pi**2*ll2)/24576d0 +
     & (647707d0*ll2)/722534400d0 + (363d0*ll2**2)/573440d0 +
     & (ll2**3)/12288d0 - (363d0*ll2*llzp)/286720d0 + (zeta3)
     &/8192d0) + zp*(-((pi**2*ll2)/24d0) + (ll2**3)/12d0 + (z
     &eta3)/8d0) + zp**3*(-((pi**2)/192d0) + 1d0/48d0 - (pi**
     &2*ll2)/288d0 + (7d0*ll2)/96d0 + (ll2**2)/32d0 + (ll2**3
     &)/144d0 - (ll2*llzp)/16d0 + (zeta3)/96d0)

         case(76)            !110-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4)/480d0) - (pi**2*ll2**2)/12d
     &0 + (5d0*ll2*zeta3)/8d0 + zp**5*(-((pi**2*5d0)/2304d0) 
     &+ 77d0/2400d0 + (pi**2*ll2)/960d0 - (llzp)/40d0 - (zeta
     &3)/256d0) + zp**8*(232819d0/45158400d0 - (pi**2*121d0)/
     &573440d0 + (pi**2*ll2)/12288d0 - (1019d0*llzp)/161280d0
     & - (5d0*zeta3)/16384d0) + zp*((pi**2*ll2)/12d0 - (5d0*z
     &eta3)/16d0) + zp**3*(11d0/144d0 - (pi**2)/96d0 + (pi**2
     &*ll2)/144d0 - (llzp)/24d0 - (5d0*zeta3)/192d0) + zp**6*
     &(-((pi**2*137d0)/138240d0) + 169d0/9600d0 + (pi**2*ll2)
     &/2304d0 - (23d0*llzp)/1440d0 - (5d0*zeta3)/3072d0) + zp
     &**9*(40487d0/14288400d0 - (pi**2*761d0)/7741440d0 + (pi
     &**2*ll2)/27648d0 - (23d0*llzp)/5670d0 - (5d0*zeta3)/368
     &64d0) + zp**4*(-((pi**2*11d0)/2304d0) + 127d0/2304d0 + 
     &(pi**2*ll2)/384d0 - (7d0*llzp)/192d0 - (5d0*zeta3)/512d
     &0) + zp**2*(-((pi**2)/48d0) + (pi**2*ll2)/48d0 - (5d0*z
     &eta3)/64d0) + zp**7*(13423d0/1411200d0 - (pi**2*7d0)/15
     &360d0 + (pi**2*ll2)/5376d0 - (101d0*llzp)/10080d0 - (5d
     &0*zeta3)/7168d0)

         case(77)            !1100

            zp = x+1d0
            szp = s(zp)

            ris = (pi**4)/48d0 - (myi*pi**3*szp*ll2)/12
     &d0 - (pi**2*ll2**2)/6d0 - (ll2**4)/12d0 - 2*cli4pt5 + (
     &myi*pi*szp*zeta3)/8d0 - ll2*zeta3 + zp**3*(-((pi**2)/32
     &d0) - (myi*pi*szp)/24d0 + (myi*pi**3*szp)/288d0 + (pi**
     &2*ll2)/48d0 - (zeta3)/32d0) + zp**4*(-((pi**2*11d0)/768
     &d0) + 1d0/96d0 + (myi*pi**3*szp)/768d0 - (myi*pi*7d0*sz
     &p)/192d0 + (pi**2*ll2)/128d0 - (3d0*zeta3)/256d0) + zp*
     &*2*(-((pi**2)/16d0) + (myi*pi**3*szp)/96d0 + (pi**2*ll2
     &)/16d0 - (3d0*zeta3)/32d0) + zp**7*(167d0/16128d0 - (pi
     &**2*7d0)/5120d0 - (myi*pi*101d0*szp)/10080d0 + (myi*pi*
     &*3*szp)/10752d0 + (pi**2*ll2)/1792d0 - (3d0*zeta3)/3584
     &d0) + zp**6*(29d0/2304d0 - (pi**2*137d0)/46080d0 - (myi
     &*pi*23d0*szp)/1440d0 + (myi*pi**3*szp)/4608d0 + (pi**2*
     &ll2)/768d0 - (zeta3)/512d0) + zp**9*(2569d0/414720d0 - 
     &(pi**2*761d0)/2580480d0 + (myi*pi**3*szp)/55296d0 - (my
     &i*pi*23d0*szp)/5670d0 + (pi**2*ll2)/9216d0 - (zeta3)/61
     &44d0) + zp**5*(-((pi**2*5d0)/768d0) + 13d0/960d0 + (myi
     &*pi**3*szp)/1920d0 - (myi*pi*szp)/40d0 + (pi**2*ll2)/32
     &0d0 - (3d0*zeta3)/640d0) + zp**8*(-((pi**2*363d0)/57344
     &0d0) + 497d0/61440d0 - (myi*pi*1019d0*szp)/161280d0 + (
     &myi*pi**3*szp)/24576d0 + (pi**2*ll2)/4096d0 - (3d0*zeta
     &3)/8192d0) + zp*((myi*pi**3*szp)/24d0 + (pi**2*ll2)/4d0
     & - (3d0*zeta3)/8d0)

         case(78)            !1101

            zp = x+1d0

            ris = -((pi**4)/30d0) - (pi**2*ll2**2)/6d0 
     &+ (ll2**4)/8d0 + 3*cli4pt5 + (23d0*ll2*zeta3)/8d0 + zp*
     &*6*(-((pi**2*137d0)/276480d0) - 31d0/5760d0 + (23d0*ll2
     &)/1440d0 + (pi**2*ll2)/4608d0 - (zeta3)/1536d0) + zp**9
     &*(-(43d0/20160d0) - (pi**2*761d0)/15482880d0 + (pi**2*l
     &l2)/55296d0 + (23d0*ll2)/5670d0 - (zeta3)/18432d0) + zp
     &**4*(-(1d0/192d0) - (pi**2*11d0)/4608d0 + (pi**2*ll2)/7
     &68d0 + (7d0*ll2)/192d0 - (zeta3)/256d0) + zp**2*(-((pi*
     &*2)/96d0) + (pi**2*ll2)/96d0 - (zeta3)/32d0) + zp**7*(-
     &(221d0/53760d0) - (pi**2*7d0)/30720d0 + (101d0*ll2)/100
     &80d0 + (pi**2*ll2)/10752d0 - (zeta3)/3584d0) + zp**5*(-
     &(1d0/160d0) - (pi**2*5d0)/4608d0 + (pi**2*ll2)/1920d0 +
     & (ll2)/40d0 - (zeta3)/640d0) + zp**8*(-((pi**2*121d0)/1
     &146880d0) - 7709d0/2580480d0 + (1019d0*ll2)/161280d0 + 
     &(pi**2*ll2)/24576d0 - (zeta3)/8192d0) + zp*((pi**2*ll2)
     &/24d0 - (zeta3)/8d0) + zp**3*(-((pi**2)/192d0) + (ll2)/
     &24d0 + (pi**2*ll2)/288d0 - (zeta3)/96d0)

         case(79)            !111-1

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4)/90d0 + (pi**2*ll2**2)/24d0 - 
     &(ll2**4)/12d0 - cli4pt5 - (7d0*ll2*zeta3)/8d0 + zp**5*(
     &-(599d0/46080d0) + (pi**2*5d0)/4608d0 - (5d0*ll2**2)/76
     &8d0 + (ll2**3)/960d0 + (7d0*llzp)/768d0 - (zeta3)/1280d
     &0) + zp**8*((pi**2*121d0)/1146880d0 - 21977d0/14745600d
     &0 - (363d0*ll2**2)/573440d0 + (ll2**3)/12288d0 + (469d0
     &*llzp)/368640d0 - (zeta3)/16384d0) + zp*((ll2**3)/12d0 
     &- (zeta3)/16d0) + zp**3*((pi**2)/192d0 - 11d0/288d0 - (
     &ll2**2)/32d0 + (ll2**3)/144d0 + (llzp)/48d0 - (zeta3)/1
     &92d0) + zp**6*((pi**2*137d0)/276480d0 - 79d0/12288d0 - 
     &(137d0*ll2**2)/46080d0 + (ll2**3)/2304d0 + (5d0*llzp)/1
     &024d0 - (zeta3)/3072d0) + zp**9*((pi**2*761d0)/15482880
     &d0 - 83359739d0/117050572800d0 - (761d0*ll2**2)/2580480
     &d0 + (ll2**3)/27648d0 + (29531d0*llzp)/46448640d0 - (ze
     &ta3)/36864d0) + zp**4*((pi**2*11d0)/4608d0 - 19d0/768d0
     & - (11d0*ll2**2)/768d0 + (ll2**3)/384d0 + (llzp)/64d0 -
     & (zeta3)/512d0) + zp**2*((pi**2)/96d0 - (ll2**2)/16d0 +
     & (ll2**3)/48d0 - (zeta3)/64d0) + zp**7*(-(3343d0/107520
     &0d0) + (pi**2*7d0)/30720d0 - (7d0*ll2**2)/5120d0 + (ll2
     &**3)/5376d0 + (29d0*llzp)/11520d0 - (zeta3)/7168d0)

         case(80)            !1110

            zp = x+1d0
            szp = s(zp)

            ris = (pi**4)/90d0 + (pi**2*ll2**2)/12d0 - 
     &(myi*pi*szp*ll2**3)/6d0 - (ll2**4)/24d0 - cli4pt5 - ll2
     &*zeta3 + zp**5*(-(11d0/1920d0) + (pi**2*5d0)/4608d0 + (
     &myi*pi*7d0*szp)/768d0 - (pi**2*ll2)/1920d0 - (myi*pi*5d
     &0*szp*ll2)/384d0 + (myi*pi*szp*ll2**2)/320d0 + (zeta3)/
     &1280d0) + zp**8*((pi**2*121d0)/1146880d0 - 563d0/286720
     &d0 + (myi*pi*469d0*szp)/368640d0 - (pi**2*ll2)/24576d0 
     &- (myi*pi*363d0*szp*ll2)/286720d0 + (myi*pi*szp*ll2**2)
     &/4096d0 + (zeta3)/16384d0) + zp*(-((pi**2*ll2)/24d0) + 
     &(myi*pi*szp*ll2**2)/4d0 + (zeta3)/16d0) + zp**3*((pi**2
     &)/192d0 + (myi*pi*szp)/48d0 - (pi**2*ll2)/288d0 - (myi*
     &pi*szp*ll2)/16d0 + (myi*pi*szp*ll2**2)/48d0 + (zeta3)/1
     &92d0) + zp**6*(-(103d0/23040d0) + (pi**2*137d0)/276480d
     &0 + (myi*pi*5d0*szp)/1024d0 - (pi**2*ll2)/4608d0 - (myi
     &*pi*137d0*szp*ll2)/23040d0 + (myi*pi*szp*ll2**2)/768d0 
     &+ (zeta3)/3072d0) + zp**9*(-(203d0/165888d0) + (pi**2*7
     &61d0)/15482880d0 + (myi*pi*29531d0*szp)/46448640d0 - (p
     &i**2*ll2)/55296d0 - (myi*pi*761d0*szp*ll2)/1290240d0 + 
     &(myi*pi*szp*ll2**2)/9216d0 + (zeta3)/36864d0) + zp**4*(
     &-(1d0/192d0) + (pi**2*11d0)/4608d0 + (myi*pi*szp)/64d0 
     &- (pi**2*ll2)/768d0 - (myi*pi*11d0*szp*ll2)/384d0 + (my
     &i*pi*szp*ll2**2)/128d0 + (zeta3)/512d0) + zp**2*((pi**2
     &)/96d0 - (pi**2*ll2)/96d0 - (myi*pi*szp*ll2)/8d0 + (myi
     &*pi*szp*ll2**2)/16d0 + (zeta3)/64d0) + zp**7*(-(493d0/1
     &61280d0) + (pi**2*7d0)/30720d0 + (myi*pi*29d0*szp)/1152
     &0d0 - (pi**2*ll2)/10752d0 - (myi*pi*7d0*szp*ll2)/2560d0
     & + (myi*pi*szp*ll2**2)/1792d0 + (zeta3)/7168d0)

         case(81)            !1111

            zp = x+1d0

            ris = -((zp*ll2**3)/12d0) + (ll2**4)/24d0 +
     & zp**8*(967d0/1474560d0 - (469d0*ll2)/368640d0 + (363d0
     &*ll2**2)/573440d0 - (ll2**3)/12288d0) + zp**3*(-((ll2)/
     &48d0) + (ll2**2)/32d0 - (ll2**3)/144d0) + zp**6*(17d0/9
     &216d0 - (5d0*ll2)/1024d0 + (137d0*ll2**2)/46080d0 - (ll
     &2**3)/2304d0) + zp**9*(89d0/245760d0 - (29531d0*ll2)/46
     &448640d0 + (761d0*ll2**2)/2580480d0 - (ll2**3)/27648d0)
     & + zp**4*(1d0/384d0 - (ll2)/64d0 + (11d0*ll2**2)/768d0 
     &- (ll2**3)/384d0) + zp**2*((ll2**2)/16d0 - (ll2**3)/48d
     &0) + zp**7*(7d0/6144d0 - (29d0*ll2)/11520d0 + (7d0*ll2*
     &*2)/5120d0 - (ll2**3)/5376d0) + zp**5*(1d0/384d0 - (7d0
     &*ll2)/768d0 + (5d0*ll2**2)/768d0 - (ll2**3)/960d0)
c End of expansions around x = -1

      end select

c --- set the imaginary part back to zero if it has been modified to
c --- get the branch cuts right (and should be zero).
      if (bcflag.eq.1) then
         xre = dreal(x)

         if (n4.eq.0.and.xre.gt.0d0) then
            if (xre.lt.1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
c
         else if (n4.eq.1.and.xre.lt.1d0) then
            if (n1.ne.-1.and.n2.ne.-1.and.n3.ne.-1) then
               ris = dcmplx(dreal(ris),0d0)
            else if (xre.gt.-1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
c            
         else if (n4.eq.-1.and.xre.gt.-1d0) then
            if (n1.ne.1.and.n2.ne.1.and.n3.ne.1) then
               ris = dcmplx(dreal(ris),0d0)
            else if (xre.lt.1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
            
         endif
      endif

      HPL4arm1=ris
      return
      end function
c-source-file HPL4at1.f
      double complex function HPL4at1(n1,n2,n3,n4)
      implicit none
      integer n1,n2,n3,n4,j
      double complex ris,myi,cli4pt5,cli4
      double precision pi, zeta2, zeta3,zeta4,ll2

      pi=3.1415926535897932385D0
      zeta3=1.20205690315959428539973816151d0
      zeta2=pi**2/6d0
      zeta4=pi**4/90d0
      myi = dcmplx(0d0,1d0)

      ll2 = dlog(2d0)
      cli4pt5 = cli4(dcmplx(0.5d0,0d0))

      j=1+(n4+1)*1+(n3+1)*3+(n2+1)*9+(n1+1)*27
      ris = dcmplx(0d0,0d0)

      if ((j.le.54).or.(j.eq.68)) then
         select case (j)
         case(1)                !-1-1-1-1
            ris = ll2**4/24d0
         case(2)                !-1-1-10
            ris = -pi**4/90d0 - (pi**2*ll2**2)/12d0 
     &           + ll2**4/24d0 + cli4pt5 + ll2*zeta3
         case(3)                !-1-1-11
            ris = pi**4/90d0 + (pi**2*ll2**2)/24d0 
     &           - ll2**4/12d0 - cli4pt5 
     &           - (7*ll2*zeta3)/8d0
         case(4)                !-1-10-1
            ris = pi**4/30d0 + (pi**2*ll2**2)/6d0 
     &           - ll2**4/8d0 - 3*cli4pt5 
     &           - (23*ll2*zeta3)/8d0
         case(5)                !-1-100
            ris = pi**4/48d0 + (pi**2*ll2**2)/12d0 
     &           - ll2**4/12d0 - 2*cli4pt5 
     &           - ll2*zeta3
         case(6)                !-1-101
            ris = pi**4/480d0 + (pi**2*ll2**2)/12d0 
     &           - (5*ll2*zeta3)/8d0
         case(7)                !-1-11-1
            ris = -pi**4/30d0 - (pi**2*ll2**2)/8d0 
     &           + ll2**4/12d0 + 3*cli4pt5 
     &           + (11*ll2*zeta3)/4d0
         case(8)                !-1-110
            ris = (-7*pi**4)/288d0 - (pi**2*ll2**2)/24d0 
     &           + ll2**4/12d0 + 2*cli4pt5 
     &           + (13*ll2*zeta3)/8d0
         case(9)                !-1-111
            ris = pi**4/720d0 + ll2**4/24d0 
     &           -(ll2*zeta3)/8d0
         case(10)               !-10-1-1
            ris = -pi**4/30d0 - (pi**2*ll2**2)/8d0 
     &           + ll2**4/8d0 + 3*cli4pt5 
     &           + (11*ll2*zeta3)/4d0
         case(11)               !-10-10
            ris = (-11*pi**4)/288d0 - (pi**2*ll2**2)/6d0 
     &           + ll2**4/6d0 + 4*cli4pt5 
     &           + 2*ll2*zeta3
         case(12)               !-10-11
            ris = (5*pi**4)/144d0 - (pi**2*ll2**2)/12d0 
     &           - ll2**4/6d0 - 4*cli4pt5 - ll2*zeta3
         case(13)               !-100-1
            ris = -pi**4/288d0 + (3*ll2*zeta3)/4d0
         case(14)               !-1000
            ris = (-7*pi**4)/720d0
         case(15)               !-1001
            ris = pi**4/60d0 + (pi**2*ll2**2)/12d0 
     &           - ll2**4/12d0 - 2*cli4pt5 
     &           - (3*ll2*zeta3)/4d0
         case(16)               !-101-1
            ris = (-7*pi**4)/180d0 + (pi**2*ll2**2)/12d0 
     &           + ll2**4/6d0 + 4*cli4pt5 
     &           + (13*ll2*zeta3)/8d0
         case(17)               !-1010
            ris = (-17*pi**4)/480d0 - (pi**2*ll2**2)/6d0 
     &           + ll2**4/6d0 + 4*cli4pt5 
     &           + (3*ll2*zeta3)/2d0
         case(18)               !-1011
            ris = pi**4/288d0 + (pi**2*ll2**2)/24d0 
     &           - ll2**4/24d0 - cli4pt5 
     &           + (ll2*zeta3)/8d0
         case(19)               !-11-1-1
            ris = pi**4/30d0 + (pi**2*ll2**2)/6d0 
     &           - ll2**4/6d0 - 3*cli4pt5 
     &           - (23*ll2*zeta3)/8d0
         case(20)               !-11-10
            ris = pi**4/360d0 + (pi**2*ll2**2)/24d0 
     &           - ll2*zeta3
         case(21)               !-11-11
            ris = pi**4/1440d0 - (pi**2*ll2**2)/24d0 
     &           + ll2**4/24d0 + (ll2*zeta3)/4d0
         case(22)               !-110-1
            ris = (11*pi**4)/240d0 + (pi**2*ll2**2)/8d0 
     &           - ll2**4/6d0 - 4*cli4pt5 
     &           - (13*ll2*zeta3)/4d0
         case(23)               !-1100
            ris = (19*pi**4)/1440d0 - (3*ll2*zeta3)/4d0
         case(24)               !-1101
            ris = (19*pi**4)/1440d0 - (pi**2*ll2**2)/24d0 
     &           - ll2**4/24d0 - cli4pt5 
     &           - (ll2*zeta3)/4d0
         case(25)               !-111-1
            ris = -pi**4/288d0 - (pi**2*ll2**2)/24d0 
     &           + ll2**4/24d0 + (7*ll2*zeta3)/8d0
         case(26)               !-1110
            ris = -pi**4/720d0 - ll2**4/24d0 - cli4pt5 
     &           + (ll2*zeta3)/8d0
         case(27)               !-1111
            ris = cli4pt5
         case(28)               !0-1-1-1
            ris = pi**4/90d0 + (pi**2*ll2**2)/24d0 
     &           - ll2**4/24d0 - cli4pt5 
     &           - (7*ll2*zeta3)/8d0
         case(29)               !0-1-10
            ris = -pi**4/288d0
         case(30)               !0-1-11
            ris = -pi**4/80d0 + (pi**2*ll2**2)/24d0 
     &           + ll2**4/12d0 + 2*cli4pt5
         case(31)               !0-10-1
            ris = (13*pi**4)/288d0 + (pi**2*ll2**2)/6d0 
     &           - ll2**4/6d0 - 4*cli4pt5 
     &           - (7*ll2*zeta3)/2d0
         case(32)               !0-100
            ris = (7*pi**4)/240d0
         case(33)               !0-101
            ris = pi**4/480d0
         case(34)               !0-11-1
            ris = (-7*pi**4)/720d0 - (pi**2*ll2**2)/4d0 
     &           + (21*ll2*zeta3)/8d0
         case(35)               !0-110
            ris = (13*pi**4)/1440d0 + (pi**2*ll2**2)/6d0 
     &           - ll2**4/6d0 - 4*cli4pt5
         case(36)               !0-111
            ris = (-11*pi**4)/720d0 + ll2**4/8d0 
     &           + 3*cli4pt5
         case(37)               !00-1-1
            ris = -pi**4/48d0 - (pi**2*ll2**2)/12d0 
     &           + ll2**4/12d0 + 2*cli4pt5 
     &           + (7*ll2*zeta3)/4d0
         case(38)               !00-10
            ris = (-7*pi**4)/240d0
         case(39)               !00-11
            ris = -pi**4/180d0 - (pi**2*ll2**2)/12d0 
     &           + ll2**4/12d0 + 2*cli4pt5
         case(40)               !000-1
            ris = (7*pi**4)/720d0
         case(41)               !0000
            ris = 0d0
         case(42)               !0001
            ris = pi**4/90d0
         case(43)               !001-1
            ris = (-19*pi**4)/1440d0 + (7*ll2*zeta3)/4d0
         case(44)               !0010
            ris = -pi**4/30d0
         case(45)               !0011
            ris = pi**4/360d0
         case(46)               !01-1-1
            ris = (7*pi**4)/288d0 + (5*pi**2*ll2**2)/24d0 
     &           - ll2**4/12d0 - 2*cli4pt5 
     &           - (21*ll2*zeta3)/8d0
         case(47)               !01-10
            ris = (-11*pi**4)/480d0 - (pi**2*ll2**2)/6d0 
     &           + ll2**4/6d0 + 4*cli4pt5
         case(48)               !01-11
            ris = (7*pi**4)/288d0 - (pi**2*ll2**2)/8d0 
     &           - ll2**4/8d0 - 3*cli4pt5
         case(49)               !010-1
            ris = (71*pi**4)/1440d0 + (pi**2*ll2**2)/6d0 
     &           - ll2**4/6d0 - 4*cli4pt5 
     &           - (7*ll2*zeta3)/2d0
         case(50)               !0100
            ris = pi**4/30d0
         case(51)               !0101
            ris = pi**4/120d0
         case(52)               !011-1
            ris = -pi**4/80d0 + (pi**2*ll2**2)/12d0 
     &           + ll2**4/24d0 + cli4pt5 
     &           + (7*ll2*zeta3)/8d0
         case(53)               !0110
            ris = -pi**4/72d0
         case(54)               !0111
            ris = pi**4/90d0
         case(68)               !1000
            ris = -pi**4/90d0
         end select
      else
         print*, ""
         print*, "****************"
         print*, "ERROR in HPL4: "
         print*, "HPL4(",n1,",",n2,",",n3,",",n4
     &        ,",1) is divergent!"
         print*, "Aborting..."
         print*,"****************"
         stop
      endif
      HPL4at1=ris
      return
      end function
c-source-file HPL4atm1.f
      double complex function HPL4atm1(n1,n2,n3,n4)
      implicit none
      integer n1,n2,n3,n4,j
      double complex ris,myi,cli4pt5,cli4
      double precision pi, zeta2, zeta3,zeta4,ll2

      pi=3.1415926535897932385D0
      zeta3=1.20205690315959428539973816151d0
      zeta2=pi**2/6d0
      zeta4=pi**4/90d0
      myi = dcmplx(0d0,1d0)

      ll2 = dlog(2d0)
      cli4pt5 = cli4(dcmplx(0.5d0,0d0))

      j=1+(n4+1)*1+(n3+1)*3+(n2+1)*9+(n1+1)*27
      ris = dcmplx(0d0,0d0)

      if (j.ge.28) then
         select case (j)
         case(28)               !0-1-1-1
            ris = -pi**4/90d0
         case(29)               !0-1-10
            ris = -pi**4/72d0 + myi*pi*zeta3
         case(30)               !0-1-11
            ris = pi**4/80d0 - (pi**2*ll2**2)/12d0 
     &           - ll2**4/24d0 - cli4pt5 
     &           - (7*ll2*zeta3)/8d0
         case(31)               !0-10-1
            ris = pi**4/120d0
         case(32)               !0-100
            ris = pi**4/20d0 + 2d0*myi*pi*zeta3
         case(33)               !0-101
            ris = (71*pi**4)/1440d0 + (pi**2*ll2**2)/6d0 
     &           - ll2**4/6d0 - 4*cli4pt5 
     &           - (7*ll2*zeta3)/2d0
         case(34)               !0-11-1
            ris = (-7*pi**4)/288d0 + (pi**2*ll2**2)/8d0 
     &           + ll2**4/8d0 + 3*cli4pt5
         case(35)               !0-110
            ris = (-71*pi**4)/1440d0 - (pi**2*ll2**2)/6d0 
     &           + ll2**4/6d0 + 4*cli4pt5 
     &           + myi*pi*((pi**2*ll2)/4d0 
     &           - zeta3) + (7*ll2*zeta3)/2d0 
     &           - 2*((-19*pi**4)/1440d0 
     &           + (7*ll2*zeta3)/4d0)
         case(36)               !0-111
            ris = (-7*pi**4)/288d0 - (5*pi**2*ll2**2)/24d0 
     &           + ll2**4/12d0 + 2*cli4pt5 
     &           + (21*ll2*zeta3)/8d0
         case(37)               !00-1-1
            ris = pi**4/360d0
         case(38)               !00-10
            ris = pi**4/30d0 - myi*pi*zeta3
         case(39)               !00-11
            ris = (-19*pi**4)/1440d0 + (7*ll2*zeta3)/4d0
         case(40)               !000-1
            ris = -pi**4/90d0
         case(41)               !0000
            ris = pi**4/24d0
         case(42)               !0001
            ris = (-7*pi**4)/720d0
         case(43)               !001-1
            ris = -pi**4/180d0 - (pi**2*ll2**2)/12d0 
     &           + ll2**4/12d0 + 2*cli4pt5
         case(44)               !0010
            ris = (7*pi**4)/240d0 - myi*3d0/4d0*pi*zeta3
         case(45)               !0011
            ris = -pi**4/48d0 - (pi**2*ll2**2)/12d0 
     &           + ll2**4/12d0 + 2*cli4pt5 
     &           + (7*ll2*zeta3)/4d0
         case(46)               !01-1-1
            ris = (11*pi**4)/720d0 - ll2**4/8d0 - 3*cli4pt5
         case(47)               !01-10
            ris = -pi**4/480d0 - 2*(-pi**4/180d0 
     &           - (pi**2*ll2**2)/12d0 + ll2**4/12d0 
     &           + 2*cli4pt5) 
     &           + myi*pi*(-(pi**2*ll2)/4d0 + (13*zeta3)/8d0)
         case(48)               !01-11
            ris = (7*pi**4)/720d0 + (pi**2*ll2**2)/4d0 
     &           - (21*ll2*zeta3)/8d0
         case(49)               !010-1
            ris = pi**4/480d0
         case(50)               !0100
            ris = pi**4/80d0 + myi*3d0/2d0*pi*zeta3
         case(51)               !0101
            ris = (13*pi**4)/288d0 + (pi**2*ll2**2)/6d0 
     &           - ll2**4/6d0 - 4*cli4pt5 
     &           - (7*ll2*zeta3)/2d0
         case(52)               !011-1
            ris = pi**4/80d0 - (pi**2*ll2**2)/24d0 
     &           - ll2**4/12d0 - 2*cli4pt5
         case(53)               !0110
            ris = (-13*pi**4)/288d0 - (pi**2*ll2**2)/6d0 
     &           + ll2**4/6d0 + 4*cli4pt5 
     &           + myi/8d0*pi*zeta3 
     &           + (7*ll2*zeta3)/2d0 - 2*(-pi**4/48d0 
     &           - (pi**2*ll2**2)
     &           /12d0 + ll2**4/12d0 + 2*cli4pt5 
     &           + (7*ll2*zeta3)/4d0)
         case(54)               !0111
            ris = -pi**4/90d0 - (pi**2*ll2**2)/24d0 
     &           + ll2**4/24d0 + cli4pt5 
     &           + (7*ll2*zeta3)/8d0
         case(55)               !1-1-1-1
            ris = cli4pt5
         case(56)               !1-1-10
            ris = pi**4/720d0 + ll2**4/24d0 + cli4pt5 
     &           - myi*pi*(-(pi**2*ll2)/12d0 + ll2**3/6d0 
     &           + (7*zeta3)/8d0) - (ll2*zeta3)/8d0
         case(57)               !1-1-11
            ris = -pi**4/288d0 - (pi**2*ll2**2)/24d0 
     &           + ll2**4/24d0 + (7*ll2*zeta3)/8d0
         case(58)               !1-10-1
            ris = (-19*pi**4)/1440d0 + (pi**2*ll2**2)/24d0 
     &           + ll2**4/24d0 + cli4pt5 
     &           + (ll2*zeta3)/4d0
         case(59)               !1-100
            ris = (19*pi**4)/1440d0 - (pi**2*(pi**2/12d0 
     &           - ll2**2/2d0))/2d0 
     &           - myi*pi*((pi**2*ll2)/6d0 
     &           - (5*zeta3)/8d0) - (3*ll2*zeta3)/4d0 
     &           - myi*pi*(-(pi**2*ll2)/4d0 + (13*zeta3)/8d0)
         case(60)               !1-101
            ris = (-11*pi**4)/240d0 - (pi**2*ll2**2)/8d0 
     &           + ll2**4/6d0 + 4*cli4pt5 
     &           + (13*ll2*zeta3)/4d0
         case(61)               !1-11-1
            ris = pi**4/1440d0 - (pi**2*ll2**2)/24d0 
     &           + ll2**4/24d0 + (ll2*zeta3)/4d0
         case(62)               !1-110
            ris = -pi**4/360d0 - (pi**2*ll2**2)/24d0 
     &           - myi*pi*((pi**2*ll2)/12d0 - ll2**3/6d0 
     &           - zeta3/4d0) 
     &           + ll2*zeta3
         case(63)               !1-111
            ris = pi**4/30d0 + (pi**2*ll2**2)/6d0 
     &           - ll2**4/6d0 - 3*cli4pt5 
     &           - (23*ll2*zeta3)/8d0
         case(64)               !10-1-1
            ris = -pi**4/288d0 - (pi**2*ll2**2)/24d0 
     &           + ll2**4/24d0 + cli4pt5 
     &           - (ll2*zeta3)/8d0
         case(65)               !10-10
            ris = -pi**4/480d0 + myi*pi*((pi**2*ll2)/6d0 
     &           - (5*zeta3)/8d0) - 2*(pi**4/60d0 
     &           + (pi**2*ll2**2)/12d0 
     &           - ll2**4/12d0 - 2*cli4pt5 
     &           - (3*ll2*zeta3)/4d0)
         case(66)               !10-11
            ris = (7*pi**4)/180d0 - (pi**2*ll2**2)/12d0 
     &           - ll2**4/6d0 - 4*cli4pt5 
     &           - (13*ll2*zeta3)/8d0
         case(67)               !100-1
            ris = pi**4/60d0 + (pi**2*ll2**2)/12d0 
     &           - ll2**4/12d0 - 2*cli4pt5 
     &           - (3*ll2*zeta3)/4d0
         case(68)               !1000
            ris = (-23*pi**4)/720d0 + myi/6d0*pi**3*ll2 
     &           - myi*3d0/4d0*pi*zeta3
         case(69)               !1001
            ris = -pi**4/288d0 + (3*ll2*zeta3)/4d0
         case(70)               !101-1
            ris = (-5*pi**4)/144d0 + (pi**2*ll2**2)/12d0 
     &           + ll2**4/6d0 + 4*cli4pt5 + ll2*zeta3
         case(71)               !1010
            ris = (-13*pi**4)/288d0 - (pi**2*ll2**2)/6d0 
     &           + ll2**4/6d0 + 4*cli4pt5 
     &           + myi*pi*((pi**2*ll2)
     &           /12d0 - zeta3/4d0) + (7*ll2*zeta3)/2d0 
     &           - 2*(-pi**4/288d0 
     &           + (3*ll2*zeta3)/4d0)
         case(72)               !1011
            ris = pi**4/30d0 + (pi**2*ll2**2)/8d0 
     &           - ll2**4/8d0 - 3*cli4pt5 
     &           - (11*ll2*zeta3)/4d0
         case(73)               !11-1-1
            ris = pi**4/720d0 + ll2**4/24d0 
     &           - (ll2*zeta3)/8d0
         case(74)               !11-10
            ris = (7*pi**4)/288d0 + (pi**2*ll2**2)/24d0 
     &           - ll2**4/12d0 - 2*cli4pt5 
     &           - myi*pi*(-ll2**3/6d0 
     &           + zeta3/8d0) - (13*ll2*zeta3)/8d0
         case(75)               !11-11
            ris = -pi**4/30d0 - (pi**2*ll2**2)/8d0 
     &           + ll2**4/12d0 + 3*cli4pt5 
     &           + (11*ll2*zeta3)/4d0
         case(76)               !110-1
            ris = -pi**4/480d0 - (pi**2*ll2**2)/12d0 
     &           + (5*ll2*zeta3)/8d0
         case(77)               !1100
            ris = pi**4/48d0 - (pi**2*ll2**2)/6d0 
     &           - ll2**4/12d0 - 2*cli4pt5 
     &           - myi*pi*((pi**2*ll2)
     &           /12d0 - zeta3/4d0) - myi/8d0*pi*zeta3 
     &           - ll2*zeta3
         case(78)               !1101
            ris = -pi**4/30d0 - (pi**2*ll2**2)/6d0 
     &           + ll2**4/8d0 + 3*cli4pt5 
     &           + (23*ll2*zeta3)/8d0
         case(79)               !111-1
            ris = pi**4/90d0 + (pi**2*ll2**2)/24d0 
     &           - ll2**4/12d0 - cli4pt5 
     &           - (7*ll2*zeta3)/8d0
         case(80)               !1110
            ris = pi**4/90d0 + (pi**2*ll2**2)/12d0 
     &           - myi/6d0*pi*ll2**3 - ll2**4/24d0 
     &           - cli4pt5 
     &           - ll2*zeta3
         case(81)               !1111
            ris = ll2**4/24d0
         end select
      else
         print*, ""
         print*, "****************"
         print*, "ERROR in HPL4: "
         print*, "HPL4(",n1,",",n2,",",n3,",",n4
     &        ,",-1) is divergent!"
         print*, "Aborting..."
         print*,"****************"
         stop
      endif
      HPL4atm1=ris
      return
      end function
c-source-file HPL4else.f
      double complex function HPL4else(n1, n2, n3, n4, x)
      implicit  none 
      double precision pi, zeta2, zeta3, zeta4, xre
      double complex x, ris, myi
c      double complex cli2,cli3
      double complex cli4
      double complex basis1,basis2,basis3,basis4,basis5,basis6,basis7
      double complex basis8,basis9,basis10,basis11,basis12,basis13
      double complex basis14,basis15,basis16,basis17,basis18
      double complex basis3_1,basis3_2,basis3_3,basis3_4,basis3_5
      double complex basis3_6,basis3_7,basis3_8
      double complex basis2_1,basis2_2,basis2_3
      integer n1,n2,n3,n4,j,bcflag,bcflag_save
      double complex ll2,cli4pt5,ll1px,ll1mx,llx,b16,b17,b18

c     #####################################################
c     basis1(x) = cli4(x) 
c     basis2(x) = cli4(-x)
c     basis3(x) = cli4(1-x)
c     basis4(x) = cli4(1/(1+x)) 
c     basis5(x) = cli4(x/(x-1))
c     basis6(x) = cli4(x/(x+1)) 
c     basis7(x) = cli4((1+x)/2) 
c     basis8(x) = cli4((1-x)/2)
c     basis9(x) = cli4((1-x)/(1+x))
c     basis10(x) = cli4((x-1)/(x+1))
c     basis11(x) = cli4(2x/(1+x))
c     basis12(x) = cli4(2x/(x-1)) 
c     basis13(x) = cli4(1-x^2) = cli4_sbc 
c     basis14(x) = cli4(x^2/(x^2-1)) 
c     basis15(x) = cli4(4x/(1+x)^2) = cli4_sbc_2  
c     basis16(x) = ch2m2(x) 
c     basis17(x) = ch21m1(x) 
c     basis18(x) = ch21m1(-x) 
c     #####################################################
c     basis3_1(z) = cli3(z) 
c     basis3_2(z) = cli3(-z)
c     basis3_3(z) = cli3(1-z)
c     basis3_4(z) = cli3(1/(1+z)) 
c     basis3_5(z) = cli3((1+z)/2) 
c     basis3_6(z) = cli3((1-z)/2) 
c     basis3_7(z) = cli3((1-z)/(1+z)) 
c     basis3_8(z) = cli3(2z/(z-1))
c     basis2_1(z) = cli2(z)
c     basis2_2(z) = cli2(-z)) 
c     basis2_3(z) = cli2((1-z)/2)
c     #####################################################
c     #####################################################
    
      pi=3.1415926535897932385d0
      zeta3=1.20205690315959428539973816151d0
      zeta2=pi**2/6d0
      zeta4=pi**4/90d0
      myi = dcmplx(0d0,1d0)

      ll2 = dlog(2d0)
      cli4pt5 = cli4(dcmplx(0.5d0,0d0))

      j=1+(n4+1)*1+(n3+1)*3+(n2+1)*9+(n1+1)*27
      
      ris=dcmplx(0d0,0d0)
      bcflag = 0

c---  +i*epsilon to get branch cuts right ---
      if (dimag(x).eq.0d0) then
         x = x + dcmplx(0d0,1d-60)
         bcflag = 1
      endif
c---  
      ll1px = log(1d0+x)
      ll1mx = log(1d0-x)
      llx = log(x)

      select case (j)
            
      case(1)                !-1-1-1-1
            
         ris = ll1px**4/24d0
            
      case(2)                   !-1-1-10
         
         ris = -pi**4/90d0 + basis4(x) - (pi**2*ll1px**2)
     &/12d0 + ll1px**4/24d0 + ll1px*zeta3
         
      case(3)                   !-1-1-11
         
         ris = -cli4pt5 + basis7(x) 
     &+ (pi**2*ll2
     &*ll1px)/12d0 - (ll2**3*ll1px)/6d0 - (pi**2
     &*ll1px**2)/24d0 + (ll2**2*ll1px**2)/4d0 
     &- (ll2*ll1px**3)/6d0 - (7*ll1px*zeta3)/8d0

      case(4)                   !-1-10-1

         ris = pi**4/30d0 - 3*basis4(x) - basis3_4(x)
     &*ll1px + (pi**2*ll1px**2)/12d0 + ll1px**4/24d0 
     &- 2*ll1px*zeta3

      case(5)                   !-1-100

         ris = pi**4/90d0 - basis2(x) - basis4(x) 
     &- basis6(x) - basis3_4(x)*llx - (pi**2*llx
     &*ll1px)/6d0 + (pi**2*ll1px**2)/12d0 - (llx**2
     &*ll1px**2)/4d0 + (llx*ll1px**3)/3d0 
     &- ll1px**4/12d0 + llx*zeta3 - ll1px*zeta3

      case(6)                   !-1-101

         ris = pi**4/480d0 + basis3(x) - basis9(x)/2d0 
     &+ basis10(x)/2d0 - basis13(x)/4d0 
     &+ basis3_4(x)*ll1mx + (pi**2*ll1mx*ll1px)
     &/12d0 + (pi**2*ll1px**2)/12d0 - (ll1mx*ll1px**3)
     &/6d0 - (7d0*ll1mx*zeta3)/8d0 - (5d0*ll1px*zeta3)/8d0

      case(7)                   !-1-11-1

         ris = 3d0*cli4pt5 
     &- 3d0*basis7(x) 
     &+ basis3_5(x)*ll1px 
     &- (pi**2*ll2*ll1px)/6d0 
     &+ (ll2**3*ll1px)/3d0 + (pi**2*ll1px**2)/24d0 
     &- (ll2**2*ll1px**2)/4d0 + (7*ll1px*zeta3)/4d0

      case(8)                   !-1-110

         ris = pi**4/90d0 + 3*cli4pt5 + basis2(x)/2d0 
     &- basis1(x)/2d0 - basis15(x)/4d0 - basis4(x) 
     &+ basis11(x) - 3*basis7(x) 
     &+ basis3_5(x)*llx + (pi**2*ll2*llx)/12d0 
     &- (ll2**3*llx)/6d0 - (pi**2*ll2*ll1px)/4d0 
     &+ (ll2**3*ll1px)/2d0 - (pi**2*llx*ll1px)/12d0 
     &+ (ll2**2*llx*ll1px)/2d0 + (5*pi**2*ll1px**2)
     &/24d0 - (3*ll2**2*ll1px**2)/4d0 - (ll2*llx
     &*ll1px**2)/2d0 + (ll2*ll1px**3)/2d0 + (llx
     &*ll1px**3)/6d0 - ll1px**4/6d0 - (7*llx*zeta3)/8d0 
     &+ (13*ll1px*zeta3)/8d0

      case(9)                   !-1-111

         ris = (-7*pi**4)/720d0 - basis8(x) 
     &- basis10(x) + basis7(x) 
     &- basis3_5(x)*ll1mx - (pi**2*ll2*ll1mx)/6d0 
     &+ (ll2**3*ll1mx)/3d0 - (ll2**2*ll1mx**2)/4d0 
     &+ (pi**2*ll2*ll1px)/12d0 -(ll2**3*ll1px)/6d0 
     &+ (pi**2*ll1mx*ll1px)/6d0 - (ll2**2*ll1mx
     &*ll1px)/2d0 + (ll2*ll1mx**2*ll1px)/2d0 
     &- (pi**2*ll1px**2)/12d0 + (ll2**2*ll1px**2)/4d0 
     &- (ll1mx**2*ll1px**2)/4d0 + (ll1mx
     &*ll1px**3)/6d0 - ll1px**4/24d0 + ll1mx*zeta3 
     &- (ll1px*zeta3)/8d0

      case(10)                  !-10-1-1

         ris = -pi**4/30d0 + 3*basis4(x) + 2*basis3_4(x)
     &*ll1px + (pi**2*ll1px**2)/12d0 + (basis2_2(x)
     &*ll1px**2)/2d0 + (llx*ll1px**3)/2d0 
     &- (5*ll1px**4)/24d0 + ll1px*zeta3

      case(11)                  !-10-10

         ris = -pi**4/45d0 + basis2_2(x)**2/2d0 + 2*basis2(x) 
     &+ 2*basis4(x) + 2*basis6(x) + 2*basis3_4(x)
     &*llx + (pi**2*llx*ll1px)/3d0 + basis2_2(x)*llx
     &*ll1px - (pi**2*ll1px**2)/6d0 + llx**2
     &*ll1px**2 - (2*llx*ll1px**3)/3d0 + ll1px**4
     &/6d0 - 2*llx*zeta3 + 2*ll1px*zeta3

      case(12)                  !-10-11

         bcflag_save=bcflag
         b18 = basis18(x)
         bcflag = bcflag_save

         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = (-11*pi**4)/720d0 + (pi**2*basis2_2(x))/12d0 
     &- basis2_3(x)*basis2_2(x) + basis2_2(x)**2/2d0 
     &+ b18 
     &- 2*basis3(x) - (3*basis2(x))/2d0 
     &- basis1(x)/2d0 - basis15(x)/4d0 + basis4(x) 
     &+ basis9(x) - basis10(x) + basis3(x) 
     &+ basis11(x) + basis13(x)/2d0 
     &+ basis2(x) + basis1(x) + basis15(x)/2d0 
     &- basis9(x)/2d0 + basis10(x)/2d0 
     &- 2*basis11(x) 
     &+ 3*basis7(x) - basis13(x)/4d0 
     &+ pi**4/480d0 - (pi**2*basis2_2(x))/12d0 
     &+ basis2_3(x)*basis2_2(x) - basis2_2(x)**2/2d0 
     &- 3*cli4pt5 
     &+ (basis2_2(x)*ll2**2)/2d0
     &+ basis3_4(x)*ll1mx - basis2_2(x)
     &*ll2*ll1mx - 2*basis3_1(x)*ll1px 
     &- 2*basis3_8(x)*ll1px - 2*basis3_7(x)
     &*ll1px + (pi**2*ll2*ll1px)/4d0 - (ll2**3
     &*ll1px)/2d0 + (5*pi**2*ll1mx*ll1px)/12d0 
     &+ basis2_2(x)*ll1mx*ll1px - ll2*ll1mx**2
     &*ll1px + (ll1mx**3*ll1px)/3d0 - ll1mx**2
     &*llx*ll1px - (5*pi**2*ll1px**2)/12d0 
     &+ (basis2_3(x)*ll1px**2)/2d0 - (basis2_2(x)
     &*ll1px**2)/2d0 + (basis2_1(x)*ll1px**2)/2d0 
     &+ ll2**2*ll1px**2 + (3*ll2*ll1mx
     &*ll1px**2)/2d0 + 2*ll1mx*llx*ll1px**2 
     &- (3*ll2*ll1px**3)/2d0 - (2*ll1mx
     &*ll1px**3)/3d0 - llx*ll1px**3 
     &+ (5*ll1px**4)/8d0 - (7*ll1mx*zeta3)/8d0 
     &- (5*ll1px*zeta3)/4d0
     &- (basis2_2(x)*ll2**2)/2d0 - 2*basis3_4(x)*ll1mx 
     &+ basis2_2(x)*ll2*ll1mx + 2*basis3_1(x)*ll1px 
     &+ 2*basis3_8(x)*ll1px + 2*basis3_7(x)
     &*ll1px - (pi**2*ll1mx*ll1px)/2d0 - basis2_2(x)
     &*ll1mx*ll1px + ll2*ll1mx**2*ll1px 
     &- (ll1mx**3*ll1px)/3d0 + ll1mx**2*llx
     &*ll1px + (pi**2*ll1px**2)/8d0 - (basis2_3(x)
     &*ll1px**2)/2d0 + (basis2_2(x)*ll1px**2)/2d0
     &-(basis2_1(x)
     &*ll1px**2)/2d0 - (ll2**2*ll1px**2)/4d0 
     &- (3*ll2*ll1mx*ll1px**2)/2d0 - 2*ll1mx
     &*llx*ll1px**2 + ll2*ll1px**3 + (5*ll1mx
     &*ll1px**3)/6d0 + (5*llx*ll1px**3)/6d0 
     &- (11*ll1px**4)/24d0 + (7*ll1mx*zeta3)/4d0 
     &+ (ll1px*zeta3)/4d0

      case(13)                  !-100-1

         ris = -basis2_2(x)**2/2d0 - basis3_2(x)*ll1px

      case(14)                  !-1000

         ris = basis2(x) - basis3_2(x)*llx + (basis2_2(x)
     &*llx**2)/2d0 + (llx**3*ll1px)/6d0

      case(15)                  !-1001

         bcflag_save=bcflag
         b16 = basis16(x)
         bcflag = bcflag_save

         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = -pi**4/360d0 + basis2_2(x)*basis2_1(x) + b16 
     &- basis3(x) - basis2(x) - basis1(x) + basis5(x) 
     &+ basis4(x) + basis6(x) + basis13(x)/4d0 
     &- basis14(x)/4d0 + basis3_2(x)*ll1mx 
     &+ (pi**2*ll1mx**2)/16d0 + ll1mx**4/32d0 
     &- (ll1mx**3*llx)/12d0 + 2*basis3_1(x)*ll1px 
     &- (pi**2*ll1mx*ll1px)/24d0 - (ll1mx**3
     &*ll1px)/24d0 + (ll1mx**2*llx*ll1px)/4d0 
     &- (5*pi**2*ll1px**2)/48d0 - (ll1mx**2*ll1px**2)
     &/16d0 + (ll1mx*llx*ll1px**2)/4d0 - (ll1mx
     &*ll1px**3)/24d0 - (llx*ll1px**3)/12d0 
     &+ (7*ll1px**4)/96d0 + (3*ll1mx*zeta3)/4d0 
     &+ (3*ll1px*zeta3)/4d0

      case(16)                  !-101-1

         bcflag_save=bcflag
         b18 = basis18(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = pi**4/90d0 - (pi**2*basis2_2(x))/12d0 + basis2_3(x)
     &*basis2_2(x) - basis2_2(x)**2/2d0 
     &- b18 
     &- (pi**4/480d0 - (pi**2*basis2_2(x))/12d0 
     &+ basis2_3(x)*basis2_2(x) - basis2_2(x)**2/2d0 
     &- 3*cli4pt5 + basis3(x) 
     &+ basis2(x) + basis1(x) + basis15(x)/2d0 
     &- basis9(x)/2d0 + basis10(x)/2d0 
     &+ 2*basis6(x) - 2*basis11(x) 
     &+ 3*basis7(x) - basis13(x)/4d0 + (basis2_2(x)
     &*ll2**2)/2d0 + basis3_4(x)*ll1mx - basis2_2(x)
     &*ll2*ll1mx - 2*basis3_1(x)*ll1px 
     &- 2*basis3_8(x)*ll1px - 2*basis3_7(x)
     &*ll1px + (pi**2*ll2*ll1px)/4d0 - (ll2**3
     &*ll1px)/2d0 + (5*pi**2*ll1mx*ll1px)/12d0 
     &+ basis2_2(x)*ll1mx*ll1px - ll2*ll1mx**2
     &*ll1px + (ll1mx**3*ll1px)/3d0 - ll1mx**2
     &*llx*ll1px - (5*pi**2*ll1px**2)/12d0 
     &+ (basis2_3(x)*ll1px**2)/2d0 - (basis2_2(x)
     &*ll1px**2)/2d0 + (basis2_1(x)*ll1px**2)/2d0 
     &+ ll2**2*ll1px**2 + (3*ll2*ll1mx
     &*ll1px**2)/2d0 + 2*ll1mx*llx*ll1px**2 
     &- (3*ll2*ll1px**3)/2d0 - (2*ll1mx
     &*ll1px**3)/3d0 - llx*ll1px**3 
     &+ (5*ll1px**4)/8d0 - (7*ll1mx*zeta3)/8d0 
     &- (5*ll1px*zeta3)/4d0)
     &+ (3*basis2(x))/2d0 + basis1(x)/2d0 
     &+ basis15(x)/4d0 - basis4(x) 
     &+ 2*basis6(x) - basis11(x) + (basis2_2(x)
     &*ll2**2)/2d0 - basis2_2(x)*ll2*ll1mx 
     &+ basis3_6(x)*ll1px - basis3_3(x)*ll1px 
     &- 2*basis3_1(x)*ll1px - 2*basis3_8(x)*ll1px 
     &- basis3_4(x)*ll1px - basis3_7(x)
     &*ll1px + basis3_5(x)*ll1px + (pi**2*ll2
     &*ll1px)/6d0 - (ll2**3*ll1px)/3d0 
     &+ (pi**2*ll1mx*ll1px)/4d0 + (ll2**2*ll1mx
     &*ll1px)/2d0 - ll2*ll1mx**2*ll1px 
     &+ (ll1mx**3*ll1px)/3d0 - ll1mx**2*llx
     &*ll1px - (3*pi**2*ll1px**2)/8d0 + (basis2_3(x)
     &*ll1px**2)/2d0 - (basis2_2(x)*ll1px**2)/2d0
     &+(basis2_1(x)
     &*ll1px**2)/2d0 + (3*ll2**2*ll1px**2)/4d0 
     &+ (ll2*ll1mx*ll1px**2)/2d0 + ll1mx*llx
     &*ll1px**2 - ll2*ll1px**3 - (5*llx
     &*ll1px**3)/6d0 + (11*ll1px**4)/24d0 
     &+ (ll1px*zeta3)/4d0

      case(17)                  !-1010

         bcflag_save=bcflag
         b16 = basis16(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = -(basis2_2(x)*basis2_1(x)) - b16 
     &+ basis3_6(x)
     &*llx - basis3_3(x)*llx - basis3_4(x)*llx 
     &+ basis3_7(x)*llx + basis3_5(x)*llx 
     &+ (pi**2*ll2*llx)/6d0 - (ll2**3*llx)/3d0 
     &- (pi**2*ll1mx*llx)/12d0 - basis2_2(x)*ll1mx
     &*llx + (ll2**2*ll1mx*llx)/2d0 - 2*basis3_1(x)
     &*ll1px - (pi**2*llx*ll1px)/12d0 
     &+ (ll2**2*llx*ll1px)/2d0 - ll2*ll1mx
     &*llx*ll1px - ll1mx*llx**2*ll1px 
     &+ (ll1mx*llx*ll1px**2)/2d0 - (3*llx*zeta3)/4d0

      case(18)                  !-1011

         ris = pi**4/288d0 - basis4(x) + basis9(x)
     &/2d0 - basis10(x)/2d0 - basis13(x)/4d0 
     &- basis3_6(x)*ll1mx + basis3_3(x)*ll1mx 
     &+ basis3_4(x)*ll1mx - basis3_7(x)
     &*ll1mx - basis3_5(x)*ll1mx - (pi**2*ll2
     &*ll1mx)/6d0 + (ll2**3*ll1mx)/3d0 + (pi**2
     &*ll1mx**2)/24d0 + (basis2_2(x)*ll1mx**2)/2d0 
     &- (ll2**2*ll1mx**2)/2d0 + (pi**2*ll1mx
     &*ll1px)/4d0 - (ll2**2*ll1mx*ll1px)/2d0 
     &+ ll2*ll1mx**2*ll1px + (ll1mx**2*llx
     &*ll1px)/2d0 + (pi**2*ll1px**2)/24d0 
     &- (ll1mx**2*ll1px**2)/2d0 - ll1px**4/24d0 
     &+ (ll1mx*zeta3)/8d0 + (ll1px*zeta3)/8d0

      case(19)                  !-11-1-1

         ris = -3*cli4pt5 + 3*basis7(x) 
     &- 2*basis3_5(x)*ll1px + (pi**2*ll2*ll1px)
     &/12d0 - (ll2**3*ll1px)/6d0 + (pi**2*ll1px**2)/12d0 
     &- (basis2_3(x)*ll1px**2)/2d0 - (ll2**2
     &*ll1px**2)/2d0 + (ll2*ll1mx*ll1px**2)/2d0 
     &+ (ll2*ll1px**3)/2d0 - (ll1mx*ll1px**3)/2d0 
     &- (7*ll1px*zeta3)/8d0

      case(20)                  !-11-10

         bcflag_save=bcflag
         b18 = basis18(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = -pi**4/90d0 - basis2_2(x)**2/2d0 
     &- b18 
     &- (pi**4/480d0 - (pi**2*basis2_2(x))/12d0 
     &+ basis2_3(x)*basis2_2(x) - basis2_2(x)**2/2d0 
     &- 3*cli4pt5 + basis3(x) 
     &+ basis2(x) + basis1(x) + basis15(x)/2d0 
     &- basis9(x)/2d0 + basis10(x)/2d0 
     &+ 2*basis6(x) - 2*basis11(x) 
     &+ 3*basis7(x) - basis13(x)/4d0 + (basis2_2(x)
     &*ll2**2)/2d0 + basis3_4(x)*ll1mx - basis2_2(x)
     &*ll2*ll1mx - 2*basis3_1(x)*ll1px 
     &- 2*basis3_8(x)*ll1px - 2*basis3_7(x)
     &*ll1px + (pi**2*ll2*ll1px)/4d0 - (ll2**3
     &*ll1px)/2d0 + (5*pi**2*ll1mx*ll1px)/12d0 
     &+ basis2_2(x)*ll1mx*ll1px - ll2*ll1mx**2
     &*ll1px + (ll1mx**3*ll1px)/3d0 - ll1mx**2
     &*llx*ll1px - (5*pi**2*ll1px**2)/12d0 
     &+ (basis2_3(x)*ll1px**2)/2d0 - (basis2_2(x)
     &*ll1px**2)/2d0 + (basis2_1(x)*ll1px**2)/2d0 
     &+ ll2**2*ll1px**2 + (3*ll2*ll1mx
     &*ll1px**2)/2d0 + 2*ll1mx*llx*ll1px**2 
     &- (3*ll2*ll1px**3)/2d0 - (2*ll1mx
     &*ll1px**3)/3d0 - llx*ll1px**3 
     &+ (5*ll1px**4)/8d0 - (7*ll1mx*zeta3)/8d0 
     &- (5*ll1px*zeta3)/4d0)
     &- 6*cli4pt5 + basis2(x)/2d0 
     &+ (3*basis1(x))/2d0 + (3*basis15(x))/4d0 
     &+ basis4(x) + 2*basis6(x) - 3*basis11(x) 
     &+ 6*basis7(x) - 2*basis3_5(x)*llx 
     &- (pi**2*ll2*llx)/6d0 + (ll2**3*llx)/3d0 
     &- 2*basis3_1(x)*ll1px - 2*basis3_8(x)*ll1px 
     &- 2*basis3_7(x)*ll1px + (pi**2*ll2
     &*ll1px)/2d0 - ll2**3*ll1px 
     &+ (pi**2*ll1mx*ll1px)/3d0 - ll2*ll1mx**2
     &*ll1px + (ll1mx**3*ll1px)/3d0 + (pi**2*llx
     &*ll1px)/4d0 - basis2_3(x)*llx*ll1px 
     &- (3*ll2**2*llx*ll1px)/2d0 + ll2*ll1mx
     &*llx*ll1px - ll1mx**2*llx*ll1px 
     &- (17*pi**2*ll1px**2)/24d0 + (basis2_3(x)
     &*ll1px**2)/2d0 - (basis2_2(x)*ll1px**2)/2d0 
     &+ (basis2_1(x)*ll1px**2)/2d0 
     &+ (7*ll2**2*ll1px**2)
     &/4d0 + (3*ll2*ll1mx*ll1px**2)/2d0 + ll2
     &*llx*ll1px**2 + ll1mx*llx*ll1px**2 
     &- 2*ll2*ll1px**3 - (ll1mx*ll1px**3)/2d0 
     &- (7*llx*ll1px**3)/6d0 + (19*ll1px**4)/24d0 
     &+ (7*llx*zeta3)/4d0 - (9*ll1px*zeta3)/4d0

      case(21)                  !-11-11

         ris = (11*pi**4)/480d0 - (pi**2*basis2_3(x))/12d0 
     &+ basis2_3(x)**2/2d0 + 2*basis8(x) + 2
     &*basis10(x) - 2*basis7(x) 
     &- (pi**2*ll2**2)/24d0 + (basis2_3(x)*ll2**2)/2d0 
     &+ ll2**4/8d0 + 2*basis3_5(x)*ll1mx + (5*pi**2
     &*ll2*ll1mx)/12d0 - basis2_3(x)*ll2
     &*ll1mx - (7*ll2**3*ll1mx)/6d0 + ll2**2
     &*ll1mx**2 - (pi**2*ll2*ll1px)/6d0 + (ll2**3
     &*ll1px)/3d0 - (5*pi**2*ll1mx*ll1px)/12d0 
     &+ basis2_3(x)*ll1mx*ll1px + (3*ll2**2
     &*ll1mx*ll1px)/2d0 - 2*ll2*ll1mx**2
     &*ll1px + (pi**2*ll1px**2)/6d0 - (ll2**2
     &*ll1px**2)/2d0 + ll1mx**2*ll1px**2 
     &- (ll1mx*ll1px**3)/3d0 + ll1px**4/12d0 
     &- 2*ll1mx*zeta3 + (ll1px*zeta3)/4d0

      case(22)                  !-110-1

         bcflag_save=bcflag
         b18 = basis18(x)
         bcflag=bcflag_save         

         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = -pi**4/90d0 + basis2_2(x)**2/2d0 
     &+ b18 
     &+ pi**4/480d0 - (pi**2*basis2_2(x))/12d0 
     &+ basis2_3(x)*basis2_2(x) - basis2_2(x)**2/2d0 
     &- 3*cli4pt5 + basis3(x) 
     &+ basis2(x) + basis1(x) + basis15(x)/2d0 
     &- basis9(x)/2d0 + basis10(x)/2d0 
     &+ 2*basis6(x) - 2*basis11(x) 
     &+ 3*basis7(x) - basis13(x)/4d0 + (basis2_2(x)
     &*ll2**2)/2d0 + basis3_4(x)*ll1mx - basis2_2(x)
     &*ll2*ll1mx - 2*basis3_1(x)*ll1px 
     &- 2*basis3_8(x)*ll1px - 2*basis3_7(x)
     &*ll1px + (pi**2*ll2*ll1px)/4d0 - (ll2**3
     &*ll1px)/2d0 + (5*pi**2*ll1mx*ll1px)/12d0 
     &+ basis2_2(x)*ll1mx*ll1px - ll2*ll1mx**2
     &*ll1px + (ll1mx**3*ll1px)/3d0 - ll1mx**2
     &*llx*ll1px - (5*pi**2*ll1px**2)/12d0 
     &+ (basis2_3(x)*ll1px**2)/2d0 - (basis2_2(x)
     &*ll1px**2)/2d0 + (basis2_1(x)*ll1px**2)/2d0 
     &+ ll2**2*ll1px**2 + (3*ll2*ll1mx
     &*ll1px**2)/2d0 + 2*ll1mx*llx*ll1px**2 
     &- (3*ll2*ll1px**3)/2d0 - (2*ll1mx
     &*ll1px**3)/3d0 - llx*ll1px**3 
     &+ (5*ll1px**4)/8d0 - (7*ll1mx*zeta3)/8d0 
     &- (5*ll1px*zeta3)/4d0
     &- (3*basis2(x))/2d0 - basis1(x)/2d0 
     &- basis15(x)/4d0 + basis4(x) 
     &- 2*basis6(x) + basis11(x) + basis3_2(x)
     &*ll1px + basis3_1(x)*ll1px + basis3_8(x)
     &*ll1px + basis3_4(x)*ll1px 
     &+ basis3_7(x)*ll1px - basis3_5(x)
     &*ll1px - (pi**2*ll2*ll1px)/12d0 + (ll2**3
     &*ll1px)/6d0 - (pi**2*ll1mx*ll1px)/6d0 
     &+ (ll2*ll1mx**2*ll1px)/2d0 - (ll1mx**3
     &*ll1px)/6d0 + (ll1mx**2*llx*ll1px)/2d0 
     &+ (3*pi**2*ll1px**2)/8d0 - (basis2_3(x)
     &*ll1px**2)/2d0 + (basis2_2(x)*ll1px**2)/2d0 
     &- (basis2_1(x)*ll1px**2)/2d0 
     &- (3*ll2**2*ll1px**2)
     &/4d0 - (ll2*ll1mx*ll1px**2)/2d0 - ll1mx
     &*llx*ll1px**2 + ll2*ll1px**3 + (5*llx
     &*ll1px**3)/6d0 - (11*ll1px**4)/24d0 
     &- (ll1px*zeta3)/8d0

      case(23)                  !-1100

         ris = -pi**4/90d0 - 2*cli4pt5 - basis2(x)/2d0 
     &+ (3*basis1(x))/2d0 + basis15(x)/4d0 
     &+ basis4(x) + basis6(x) - 2*basis11(x) 
     &+ 2*basis7(x) + basis3_2(x)*llx - basis3_1(x)*llx 
     &- basis3_8(x)*llx + basis3_4(x)*llx 
     &- basis3_7(x)*llx - basis3_5(x)*llx 
     &- (pi**2*ll2*llx)/12d0 + (ll2**3*llx)/6d0 
     &+ (pi**2*ll1mx*llx)/6d0 - (ll2*ll1mx**2
     &*llx)/2d0 + (ll1mx**3*llx)/6d0 + (pi**2*llx**2)
     &/24d0 - (basis2_3(x)*llx**2)/2d0 - (ll2**2
     &*llx**2)/4d0 + (ll2*ll1mx*llx**2)/2d0 
     &- (ll1mx**2*llx**2)/2d0 + (pi**2*ll2*ll1px)
     &/6d0 - (ll2**3*ll1px)/3d0 + (pi**2*llx*ll1px)
     &/12d0 - (ll2**2*llx*ll1px)/2d0 +ll2*ll1mx
     &*llx*ll1px + (ll1mx*llx**2*ll1px)/2d0 
     &- (pi**2*ll1px**2)/6d0 + (ll2**2*ll1px**2)/2d0 
     &- (ll1mx*llx*ll1px**2)/2d0 - (ll2
     &*ll1px**3)/3d0 - (llx*ll1px**3)/6d0 
     &+ ll1px**4/6d0 + (7*llx*zeta3)/8d0 
     &        - (3*ll1px*zeta3)/4d0

      case(24)                  !-1101

         bcflag_save=bcflag
         b17= basis17(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = pi**4/144d0 - basis2_2(x)*basis2_1(x)
     &+ basis2_1(x)**2/2d0 
     &+ b17-( -pi**4/288d0 
     &- (pi**2*basis2_1(x))/12d0 + basis2_3(x)
     &*basis2_1(x) - basis2_2(x)*basis2_1(x) 
     &+ basis2_1(x)**2/2d0 
     &- 5*cli4pt5 
     &+ basis8(x) - 2*basis3(x) + basis2(x) + basis1(x) 
     &- 2*basis12(x) + basis15(x)/2d0 
     &+ 3*basis4(x) - basis9(x)/2d0 
     &+ basis10(x)/2d0 + 2*basis6(x) 
     &- 4*basis11(x) + 4*basis7(x) 
     &+ basis13(x)/4d0 + (basis2_1(x)*ll2**2)/2d0 
     &- 2*basis3_8(x)*ll1mx + (pi**2*ll2
     &*ll1mx)/12d0 - basis2_1(x)*ll2*ll1mx 
     &- (ll2**3*ll1mx)/6d0 + (pi**2*ll1mx**2)/8d0 
     &+ (basis2_3(x)*ll1mx**2)/2d0 - (basis2_2(x)
     &*ll1mx**2)/2d0 + (basis2_1(x)*ll1mx**2)/2d0 
     &+ (ll2**2*ll1mx**2)/2d0-(2*ll2*ll1mx**3)/3d0 
     &+ (7*ll1mx**4)/24d0 - (ll1mx**3*llx)/3d0 
     &- basis3_3(x)*ll1px - 2*basis3_1(x)*ll1px 
     &+ (pi**2*ll2*ll1px)/3d0-(2*ll2**3*ll1px)/3d0 
     &- (3*pi**2*ll1px**2)/8d0 + ll2**2*ll1px**2 
     &- (2*ll2*ll1px**3)/3d0 - (llx*ll1px**3)/3d0 
     &+ (3*ll1px**4)/8d0 + (7*ll1mx*zeta3)/4d0 
     &- (5*ll1px*zeta3)/8d0)
     &- 4*cli4pt5 + 2*basis8(x) - basis3(x) 
     &+ (3*basis2(x))/2d0 + basis1(x)/2d0 - 2*basis5(x) 
     &- basis12(x) + basis15(x)/4d0 
     &+ 2*basis4(x) - basis9(x) 
     &+ basis10(x) - 2*basis11(x) 
     &+ 2*basis7(x) + basis13(x)/4d0 
     &+ basis14(x)/4d0 - basis3_2(x)
     &*ll1mx + basis3_1(x)*ll1mx - basis3_8(x)
     &*ll1mx - basis3_4(x)*ll1mx 
     &+ basis3_7(x)*ll1mx + basis3_5(x)
     &*ll1mx + (pi**2*ll2*ll1mx)/4d0 - (ll2**3
     &*ll1mx)/2d0 - (5*pi**2*ll1mx**2)/48d0 
     &+ (basis2_3(x)*ll1mx**2)/2d0 - (basis2_2(x)
     &*ll1mx**2)/2d0 + (basis2_1(x)*ll1mx**2)/2d0 
     &+ (3*ll2**2*ll1mx**2)/4d0-(ll2*ll1mx**3)/3d0 
     &+ (3*ll1mx**4)/32d0 + (ll1mx**3*llx)/4d0 
     &- 2*basis3_1(x)*ll1px + (pi**2*ll2*ll1px)/6d0 
     &- (ll2**3*ll1px)/3d0 - (3*pi**2*ll1mx*ll1px)
     &/8d0 + (ll2**2*ll1mx*ll1px)/2d0 - ll2
     &*ll1mx**2*ll1px + (ll1mx**3*ll1px)/24d0 
     &- (ll1mx**2*llx*ll1px)/4d0 
     &- (7*pi**2*ll1px**2)/48d0 + (ll2**2*ll1px**2)/2d0 
     &+ (9*ll1mx**2*ll1px**2)/16d0 - (ll1mx*llx
     &*ll1px**2)/4d0 - (ll2*ll1px**3)/3d0 
     &+ (ll1mx*ll1px**3)/24d0 - (llx*ll1px**3)
     &/12d0 + (17*ll1px**4)/96d0 - (ll1mx*zeta3)/8d0 
     &- (7*ll1px*zeta3)/4d0

      case(25)                  !-111-1

         ris = -pi**4/288d0 + (pi**2*basis2_3(x))/12d0 
     &- basis2_3(x)**2/2d0 + (pi**2*ll2**2)/24d0 
     &- (basis2_3(x)*ll2**2)/2d0 - ll2**4/8d0 
     &- (pi**2*ll2*ll1mx)/12d0 + basis2_3(x)*ll2
     &*ll1mx + (ll2**3*ll1mx)/2d0 - (ll2**2
     &*ll1mx**2)/2d0 - basis3_6(x)*ll1px - (pi**2
     &*ll2*ll1px)/12d0 + (ll2**3*ll1px)/6d0 
     &+ (pi**2*ll1mx*ll1px)/12d0 - (ll2**2*ll1mx
     &*ll1px)/2d0 + (ll2*ll1mx**2*ll1px)/2d0 
     &+ (7*ll1px*zeta3)/8d0

      case(26)                  !-1110

         bcflag_save=bcflag
         b17= basis17(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = pi**4/288d0 + basis2_2(x)*basis2_1(x)
     &- basis2_1(x)**2/2d0 
     &- (b17-( -pi**4/288d0 
     &- (pi**2*basis2_1(x))/12d0 + basis2_3(x)
     &*basis2_1(x) - basis2_2(x)*basis2_1(x) 
     &+ basis2_1(x)**2/2d0 
     &- 5*cli4pt5 
     &+ basis8(x) - 2*basis3(x) + basis2(x) + basis1(x) 
     &- 2*basis12(x) + basis15(x)/2d0 
     &+ 3*basis4(x) - basis9(x)/2d0 
     &+ basis10(x)/2d0 + 2*basis6(x) 
     &- 4*basis11(x) + 4*basis7(x) 
     &+ basis13(x)/4d0 + (basis2_1(x)*ll2**2)/2d0 
     &- 2*basis3_8(x)*ll1mx + (pi**2*ll2
     &*ll1mx)/12d0 - basis2_1(x)*ll2*ll1mx 
     &- (ll2**3*ll1mx)/6d0 + (pi**2*ll1mx**2)/8d0 
     &+ (basis2_3(x)*ll1mx**2)/2d0 - (basis2_2(x)
     &*ll1mx**2)/2d0 + (basis2_1(x)*ll1mx**2)/2d0 
     &+ (ll2**2*ll1mx**2)/2d0-(2*ll2*ll1mx**3)/3d0 
     &+ (7*ll1mx**4)/24d0 - (ll1mx**3*llx)/3d0 
     &- basis3_3(x)*ll1px - 2*basis3_1(x)*ll1px 
     &+ (pi**2*ll2*ll1px)/3d0-(2*ll2**3*ll1px)/3d0 
     &- (3*pi**2*ll1px**2)/8d0 + ll2**2*ll1px**2 
     &- (2*ll2*ll1px**3)/3d0 - (llx*ll1px**3)/3d0 
     &+ (3*ll1px**4)/8d0 + (7*ll1mx*zeta3)/4d0 
     &- (5*ll1px*zeta3)/8d0))
     &+ 5*cli4pt5 
     &- basis8(x) + 2*basis3(x) - basis2(x) - basis1(x) 
     &+ 2*basis12(x) - basis15(x)/2d0 
     &- 3*basis4(x) + basis9(x)/2d0 
     &- basis10(x)/2d0 - 2*basis6(x) 
     &+ 4*basis11(x) - 4*basis7(x) 
     &- basis13(x)/4d0 + 2*basis3_8(x)*ll1mx 
     &- (pi**2*ll2*ll1mx)/12d0 +(ll2**3*ll1mx)/6d0 
     &- (pi**2*ll1mx**2)/8d0 - (basis2_3(x)
     &*ll1mx**2)/2d0 + (basis2_2(x)*ll1mx**2)/2d0 
     &- (basis2_1(x)*ll1mx**2)/2d0
     &-(ll2**2*ll1mx**2)/2d0 
     &+ (2*ll2*ll1mx**3)/3d0 - (7*ll1mx**4)/24d0 
     &- basis3_6(x)*llx - (pi**2*ll2*llx)/12d0 
     &+ (ll2**3*llx)/6d0 + basis2_3(x)*ll1mx
     &*llx - (ll2*ll1mx**2*llx)/2d0 + (ll1mx**3
     &*llx)/3d0 + 2*basis3_1(x)*ll1px - (pi**2*ll2
     &*ll1px)/3d0 + (2*ll2**3*ll1px)/3d0 + (pi**2
     &*ll1mx*ll1px)/6d0 + (3*pi**2*ll1px**2)/8d0 
     &- ll2**2*ll1px**2 + (2*ll2*ll1px**3)/3d0 
     &+ (llx*ll1px**3)/3d0 - (3*ll1px**4)/8d0 
     &- (7*ll1mx*zeta3)/4d0 + (7*llx*zeta3)/8d0 
     &+ (13*ll1px*zeta3)/8d0

      case(27)                  !-1111

         ris = cli4pt5 - basis8(x) 
     &+ basis3_6(x)
     &*ll1mx - (basis2_3(x)*ll1mx**2)/2d0 + (ll2
     &*ll1mx**3)/6d0 - (ll1mx**3*ll1px)/6d0

      case(28)                  !0-1-1-1

         ris = pi**4/90d0 - basis4(x) - basis3_4(x)
     &*ll1px - (pi**2*ll1px**2)/12d0 - (basis2_2(x)
     &*ll1px**2)/2d0 - (llx*ll1px**3)/3d0 
     &+ ll1px**4/8d0

      case(29)                  !0-1-10

         ris = -basis2_2(x)**2/2d0 - basis3_4(x)*llx 
     &- (pi**2*llx*ll1px)/6d0 - basis2_2(x)*llx*ll1px 
     &- (llx**2*ll1px**2)/2d0 + (llx*ll1px**3)/6d0 
     &+ llx*zeta3

      case(30)                  !0-1-11

         bcflag_save=bcflag
         ris = -basis18(x)
         bcflag=bcflag_save

      case(31)                  !0-10-1

         ris = pi**4/45d0 + basis2_2(x)**2/2d0 - 2*basis2(x) 
     &- 2*basis4(x) - 2*basis6(x) + 2*basis3_2(x)
     &*ll1px + (pi**2*ll1px**2)/6d0 + (llx
     &*ll1px**3)/3d0 - ll1px**4/6d0 - 2*ll1px*zeta3

      case(32)                  !0-100

         ris = -3*basis2(x) + 2*basis3_2(x)*llx 
     &- (basis2_2(x)*llx**2)/2d0

      case(33)                  !0-101

         bcflag_save=bcflag
         b16 = basis16(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = pi**4/180d0 - basis2_2(x)*basis2_1(x) - b16 
     &+ 2*basis3(x) + 2*basis2(x) + 2*basis1(x) - 2*basis5(x) 
     &- 2*basis4(x) - 2*basis6(x) - basis13(x)/2d0 
     &+ basis14(x)/2d0 - 2*basis3_2(x)*ll1mx 
     &- (pi**2*ll1mx**2)/8d0 - ll1mx**4/16d0 
     &+ (ll1mx**3*llx)/6d0 - 2*basis3_1(x)*ll1px 
     &+ (pi**2*ll1mx*ll1px)/12d0 + (ll1mx**3
     &*ll1px)/12d0 - (ll1mx**2*llx*ll1px)/2d0 
     &+ (5*pi**2*ll1px**2)/24d0 + (ll1mx**2*ll1px**2)
     &/8d0 - (ll1mx*llx*ll1px**2)/2d0 + (ll1mx
     &*ll1px**3)/12d0 + (llx*ll1px**3)/6d0 
     &- (7*ll1px**4)/48d0 - (3*ll1mx*zeta3)/2d0 
     &- (3*ll1px*zeta3)/2d0

      case(34)                  !0-11-1

         bcflag_save=bcflag
         b18 = basis18(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = pi**4/90d0 + (pi**2*basis2_2(x))/12d0 - basis2_3(x)
     &*basis2_2(x) + basis2_2(x)**2/2d0 
     &+ b18 
     &+ pi**4/480d0 - (pi**2*basis2_2(x))/12d0 
     &+ basis2_3(x)*basis2_2(x) - basis2_2(x)**2/2d0 
     &- 3*cli4pt5 + basis3(x) 
     &+ basis2(x) + basis1(x) + basis15(x)/2d0 
     &- basis9(x)/2d0 + basis10(x)/2d0 
     &+ 2*basis6(x) - 2*basis11(x) 
     &+ 3*basis7(x) - basis13(x)/4d0 + (basis2_2(x)
     &*ll2**2)/2d0 + basis3_4(x)*ll1mx - basis2_2(x)
     &*ll2*ll1mx - 2*basis3_1(x)*ll1px 
     &- 2*basis3_8(x)*ll1px - 2*basis3_7(x)
     &*ll1px + (pi**2*ll2*ll1px)/4d0 - (ll2**3
     &*ll1px)/2d0 + (5*pi**2*ll1mx*ll1px)/12d0 
     &+ basis2_2(x)*ll1mx*ll1px - ll2*ll1mx**2
     &*ll1px + (ll1mx**3*ll1px)/3d0 - ll1mx**2
     &*llx*ll1px - (5*pi**2*ll1px**2)/12d0 
     &+ (basis2_3(x)*ll1px**2)/2d0 - (basis2_2(x)
     &*ll1px**2)/2d0 + (basis2_1(x)*ll1px**2)/2d0 
     &+ ll2**2*ll1px**2 + (3*ll2*ll1mx
     &*ll1px**2)/2d0 + 2*ll1mx*llx*ll1px**2 
     &- (3*ll2*ll1px**3)/2d0 - (2*ll1mx
     &*ll1px**3)/3d0 - llx*ll1px**3 
     &+ (5*ll1px**4)/8d0 - (7*ll1mx*zeta3)/8d0 
     &- (5*ll1px*zeta3)/4d0
     &+ 6*cli4pt5 - basis2(x)/2d0 
     &- (3*basis1(x))/2d0 - (3*basis15(x))/4d0 
     &- basis4(x) - 2*basis6(x) + 3*basis11(x) 
     &- 6*basis7(x) - (basis2_2(x)*ll2**2)/2d0 + basis2_2(x)
     &*ll2*ll1mx - basis3_6(x)*ll1px 
     &+ basis3_3(x)*ll1px - basis3_2(x)*ll1px 
     &+ 3*basis3_1(x)
     &*ll1px + 3*basis3_8(x)*ll1px 
     &+ 2*basis3_7(x)*ll1px - (7*pi**2*ll2
     &*ll1px)/12d0 + (7*ll2**3*ll1px)/6d0 
     &- (5*pi**2*ll1mx*ll1px)/12d0 - (ll2**2*ll1mx
     &*ll1px)/2d0 + (3*ll2*ll1mx**2*ll1px)/2d0 
     &- (ll1mx**3*ll1px)/2d0 + (3*ll1mx**2*llx
     &*ll1px)/2d0 + (17*pi**2*ll1px**2)/24d0 
     &- (basis2_3(x)*ll1px**2)/2d0 + (basis2_2(x)
     &*ll1px**2)/2d0 - (basis2_1(x)*ll1px**2)/2d0 
     &- (7*ll2**2*ll1px**2)/4d0 - (3*ll2*ll1mx
     &*ll1px**2)/2d0 - 2*ll1mx*llx*ll1px**2 
     &+ 2*ll2*ll1px**3 + (ll1mx*ll1px**3)/2d0 
     &+ (7*llx*ll1px**3)/6d0 - (19*ll1px**4)/24d0 
     &+ (17*ll1px*zeta3)/8d0

      case(35)                  !0-110

         bcflag_save=bcflag
         b16 = basis16(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = pi**4/45d0 + basis2_2(x)*basis2_1(x) + b16 
     &+ 4*cli4pt5 + basis2(x) - 3*basis1(x) 
     &- basis15(x)/2d0 - 2*basis4(x) 
     &- 2*basis6(x) + 4*basis11(x) 
     &- 4*basis7(x) - basis3_6(x)*llx 
     &+ basis3_3(x)*llx - basis3_2(x)*llx 
     &+ basis3_1(x)*llx 
     &+ basis3_8(x)*llx - (pi**2*ll2*llx)/12d0 
     &+ (ll2**3*llx)/6d0 - (pi**2*ll1mx*llx)/12d0 
     &+ basis2_2(x)*ll1mx*llx 
     &- (ll2**2*ll1mx*llx)
     &/2d0 + (ll2*ll1mx**2*llx)/2d0 - (ll1mx**3
     &*llx)/6d0 + (ll1mx**2*llx**2)/2d0 + 2*basis3_1(x)
     &*ll1px - (pi**2*ll2*ll1px)/3d0 + (2*ll2**3
     &*ll1px)/3d0 + (pi**2*ll1px**2)/3d0 - ll2**2
     &*ll1px**2 + (2*ll2*ll1px**3)/3d0 + (llx
     &*ll1px**3)/3d0 - ll1px**4/3d0 - (llx*zeta3)/8d0 
     &+ (3*ll1px*zeta3)/2d0

      case(36)                  !0-111

         ris = -pi**4/72d0 - cli4pt5 - basis8(x) 
     &- basis3(x) - basis2(x)/2d0 + basis1(x)/2d0 
     &+ 2*basis5(x) - basis12(x) 
     &+ basis15(x)/4d0 + 2*basis4(x) 
     &+ 2*basis6(x) - 2*basis11(x) 
     &+ 2*basis7(x) + basis13(x)/4d0 
     &- basis14(x)/4d0 + basis3_6(x)*ll1mx 
     &- basis3_3(x)*ll1mx + basis3_2(x)*ll1mx 
     &- basis3_1(x)*ll1mx - basis3_8(x)*ll1mx 
     &+ (3*pi**2*ll1mx**2)/16d0 - (basis2_2(x)*ll1mx**2)/2d0 
     &+ (ll2**2*ll1mx**2)/4d0 - (ll2*ll1mx**3)/3d0 
     &+ (19*ll1mx**4)/96d0 - (7*ll1mx**3*llx)/12d0 
     &+ (pi**2*ll2*ll1px)/6d0 - (ll2**3*ll1px)/3d0 
     &- (pi**2*ll1mx*ll1px)/24d0 - (ll1mx**3
     &*ll1px)/24d0 + (ll1mx**2*llx*ll1px)/4d0 
     &- (13*pi**2*ll1px**2)/48d0 + (ll2**2*ll1px**2)/2d0 
     &- (ll1mx**2*ll1px**2)/16d0 + (ll1mx*llx
     &*ll1px**2)/4d0 - (ll2*ll1px**3)/3d0 - (ll1mx
     &*ll1px**3)/24d0 - (llx*ll1px**3)/4d0 
     &+ (23*ll1px**4)/96d0 + (7*ll1mx*zeta3)/4d0

      case(37)                  !00-1-1

         ris = -pi**4/90d0 + basis2(x) + basis4(x) 
     &+ basis6(x) - basis3_2(x)*ll1px 
     &- (pi**2*ll1px**2)/12d0 - (llx*ll1px**3)/6d0 
     &+ ll1px**4/12d0 + ll1px*zeta3

      case(38)                  !00-10

         ris = 3*basis2(x) - basis3_2(x)*llx

      case(39)                  !00-11

         ris = -pi**4/72d0 - 2*cli4pt5 - basis3(x) 
     &- (3*basis2(x))/2d0 + basis1(x)/2d0 + basis5(x) 
     &+ basis15(x)/4d0 + 2*basis4(x) 
     &+ 2*basis6(x) - 2*basis11(x) 
     &+ 2*basis7(x) + basis13(x)/4d0 
     &- basis14(x)/4d0 + basis3_2(x)*ll1mx 
     &+ (pi**2*ll1mx**2)/16d0 + ll1mx**4/32d0 
     &- (ll1mx**3*llx)/12d0 + (pi**2*ll2*ll1px)/6d0 
     &- (ll2**3*ll1px)/3d0 - (pi**2*ll1mx
     &*ll1px)/24d0 - (ll1mx**3*ll1px)/24d0 
     &+ (ll1mx**2*llx*ll1px)/4d0 - (13*pi**2
     &*ll1px**2)/48d0 + (ll2**2*ll1px**2)/2d0 
     &- (ll1mx**2*ll1px**2)/16d0 + (ll1mx*llx
     &*ll1px**2)/4d0 - (ll2*ll1px**3)/3d0 
     &- (ll1mx*ll1px**3)/24d0 - (llx*ll1px**3)/4d0 
     &+ (23*ll1px**4)/96d0 + (3*ll1mx*zeta3)/4d0

      case(40)                  !000-1

         ris = -basis2(x)

      case(41)                  !0000

         ris = llx**4/24d0

      case(42)                  !0001

         ris = basis1(x)

      case(43)                  !001-1

         ris = pi**4/90d0 + 2*cli4pt5 + basis2(x)/2d0 
     &- (3*basis1(x))/2d0 - basis15(x)/4d0 
     &- basis4(x) - basis6(x) + 2*basis11(x) 
     &- 2*basis7(x) + basis3_1(x)*ll1px - (pi**2*ll2
     &*ll1px)/6d0 + (ll2**3*ll1px)/3d0 
     &+ (pi**2*ll1px**2)/6d0 - (ll2**2*ll1px**2)/2d0 
     &+ (ll2*ll1px**3)/3d0 + (llx*ll1px**3)/6d0 
     &- ll1px**4/6d0 + (3*ll1px*zeta3)/4d0

      case(44)                  !0010

         ris = -3*basis1(x) + basis3_1(x)*llx

      case(45)                  !0011
         
         ris = pi**4/90d0 - basis3(x) + basis1(x) + basis5(x) 
     &- basis3_1(x)*ll1mx + (pi**2*ll1mx**2)/12d0 
     &+ ll1mx**4/24d0 - (ll1mx**3*llx)/6d0 
     &+ ll1mx*zeta3

      case(46)                  !01-1-1
         ris = -pi**4/90d0 - 3*cli4pt5 - basis2(x)/2d0 
     &+ basis1(x)/2d0 + basis15(x)/4d0 + basis4(x) 
     &- basis11(x) + 3*basis7(x) + basis3_2(x)
     &*ll1px - basis3_1(x)*ll1px - basis3_8(x)
     &*ll1px + basis3_4(x)*ll1px 
     &- basis3_7(x)*ll1px - basis3_5(x)
     &*ll1px + (pi**2*ll2*ll1px)/6d0 - (ll2**3
     &*ll1px)/3d0 + (pi**2*ll1mx*ll1px)/6d0 
     &- (ll2*ll1mx**2*ll1px)/2d0 + (ll1mx**3
     &*ll1px)/6d0 - (ll1mx**2*llx*ll1px)/2d0 
     &- (pi**2*ll1px**2)/8d0 + (basis2_1(x)*ll1px**2)/2d0 
     &+ (ll2**2*ll1px**2)/4d0 + ll2*ll1mx
     &*ll1px**2 + ll1mx*llx*ll1px**2 
     &- (ll2*ll1px**3)/2d0 - (ll1mx*ll1px**3)/2d0 
     &- (llx*ll1px**3)/6d0 + ll1px**4/6d0 
     &- (3*ll1px*zeta3)/4d0

      case(47)                  !01-10

         bcflag_save=bcflag         
         b16 = basis16(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = -pi**4/45d0 - b16 - 4*cli4pt5 
     &- basis2(x) 
     &+ 3*basis1(x) + basis15(x)/2d0 + 2*basis4(x) 
     &+ 2*basis6(x) - 4*basis11(x) 
     &+ 4*basis7(x) + basis3_2(x)*llx - basis3_1(x)*llx 
     &- basis3_8(x)*llx + basis3_4(x)*llx 
     &- basis3_7(x)*llx - basis3_5(x)*llx 
     &- (pi**2*ll2*llx)/12d0 + (ll2**3*llx)/6d0 
     &+ (pi**2*ll1mx*llx)/6d0 - (ll2*ll1mx**2
     &*llx)/2d0 + (ll1mx**3*llx)/6d0 - (ll1mx**2
     &*llx**2)/2d0 - 2*basis3_1(x)*ll1px + (pi**2*ll2
     &*ll1px)/3d0 - (2*ll2**3*ll1px)/3d0 + (pi**2
     &*llx*ll1px)/12d0 + basis2_1(x)*llx*ll1px 
     &- (ll2**2*llx*ll1px)/2d0 + ll2*ll1mx
     &*llx*ll1px + ll1mx*llx**2*ll1px 
     &- (pi**2*ll1px**2)/3d0 + ll2**2*ll1px**2 
     &- (ll1mx*llx*ll1px**2)/2d0 - (2*ll2
     &*ll1px**3)/3d0 - (llx*ll1px**3)/3d0 
     &+ ll1px**4/3d0 + (7*llx*zeta3)/8d0 
     &- (3*ll1px*zeta3)/2d0

      case(48)                  !01-11

         bcflag_save=bcflag
         b17= basis17(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)


         ris = pi**4/72d0 + (pi**2*basis2_1(x))/12d0 - basis2_3(x)
     &*basis2_1(x) + basis2_2(x)*basis2_1(x) -basis2_1(x)**2/2d0 
     &- (b17-( -pi**4/288d0 
     &- (pi**2*basis2_1(x))/12d0 + basis2_3(x)
     &*basis2_1(x) - basis2_2(x)*basis2_1(x) +basis2_1(x)**2/2d0 
     &- 5*cli4pt5 
     &+ basis8(x) - 2*basis3(x) + basis2(x) + basis1(x) 
     &- 2*basis12(x) + basis15(x)/2d0 
     &+ 3*basis4(x) - basis9(x)/2d0 
     &+ basis10(x)/2d0 + 2*basis6(x) 
     &- 4*basis11(x) + 4*basis7(x) 
     &+ basis13(x)/4d0 + (basis2_1(x)*ll2**2)/2d0 
     &- 2*basis3_8(x)*ll1mx + (pi**2*ll2
     &*ll1mx)/12d0 - basis2_1(x)*ll2*ll1mx 
     &- (ll2**3*ll1mx)/6d0 + (pi**2*ll1mx**2)/8d0 
     &+ (basis2_3(x)*ll1mx**2)/2d0 - (basis2_2(x)
     &*ll1mx**2)/2d0 + (basis2_1(x)*ll1mx**2)/2d0 
     &+ (ll2**2*ll1mx**2)/2d0-(2*ll2*ll1mx**3)/3d0 
     &+ (7*ll1mx**4)/24d0 - (ll1mx**3*llx)/3d0 
     &- basis3_3(x)*ll1px - 2*basis3_1(x)*ll1px 
     &+ (pi**2*ll2*ll1px)/3d0-(2*ll2**3*ll1px)/3d0 
     &- (3*pi**2*ll1px**2)/8d0 + ll2**2*ll1px**2 
     &- (2*ll2*ll1px**3)/3d0 - (llx*ll1px**3)/3d0 
     &+ (3*ll1px**4)/8d0 + (7*ll1mx*zeta3)/4d0 
     &- (5*ll1px*zeta3)/8d0))
     &+ 6*cli4pt5 + 3*basis3(x) 
     &- basis2(x)/2d0 - (3*basis1(x))/2d0 - 2*basis5(x) 
     &+ 3*basis12(x) - (3*basis15(x))/4d0 
     &- 4*basis4(x) - 4*basis6(x) + 6*basis11(x) 
     &- 6*basis7(x) - basis13(x)/4d0 
     &+ basis14(x)/4d0 - (basis2_1(x)*ll2**2)/2d0 
     &- basis3_2(x)*ll1mx + basis3_1(x)*ll1mx 
     &+ 3*basis3_8(x)*ll1mx - basis3_4(x)
     &*ll1mx + basis3_7(x)*ll1mx 
     &+ basis3_5(x)*ll1mx + (pi**2*ll2
     &*ll1mx)/12d0 + basis2_1(x)*ll2*ll1mx
     &-(ll2**3
     &*ll1mx)/6d0 - (17*pi**2*ll1mx**2)/48d0 
     &- (basis2_3(x)*ll1mx**2)/2d0 
     &+ (basis2_2(x)*ll1mx**2)/2d0 
     &-(basis2_1(x)*ll1mx**2)/2d0 
     &- (ll2**2*ll1mx**2)/4d0 + ll2*ll1mx**3 
     &- (47*ll1mx**4)/96d0 + (11*ll1mx**3*llx)/12d0 
     &+ 2*basis3_1(x)*ll1px - (pi**2*ll2*ll1px)/2d0 
     &+ ll2**3*ll1px - (pi**2*ll1mx*ll1px)/24d0 
     &- basis2_1(x)*ll1mx*ll1px + (ll2**2*ll1mx
     &*ll1px)/2d0 - ll2*ll1mx**2*ll1px 
     &+ (ll1mx**3*ll1px)/24d0 - (5*ll1mx**2*llx
     &*ll1px)/4d0 + (29*pi**2*ll1px**2)/48d0 
     &- (3*ll2**2*ll1px**2)/2d0 + (9*ll1mx**2
     &*ll1px**2)/16d0 - (ll1mx*llx*ll1px**2)/4d0 
     &+ ll2*ll1px**3 + (ll1mx*ll1px**3)/24d0 
     &+ (7*llx*ll1px**3)/12d0 - (55*ll1px**4)/96d0 
     &- (29*ll1mx*zeta3)/8d0 + (3*ll1px*zeta3)/2d0

      case(49)                  !010-1

         bcflag_save=bcflag
         ris = basis16(x)
         bcflag=bcflag_save

      case(50)                  !0100

         ris = 3*basis1(x) - 2*basis3_1(x)*llx 
     &+ (basis2_1(x)*llx**2)/2d0

      case(51)                  !0101

         ris = -pi**4/45d0 + basis2_1(x)**2/2d0 + 2*basis3(x) 
     &- 2*basis1(x) - 2*basis5(x) + 2*basis3_1(x)*ll1mx 
     &- (pi**2*ll1mx**2)/6d0 - ll1mx**4/12d0 
     &+ (ll1mx**3*llx)/3d0 - 2*ll1mx*zeta3

      case(52)                  !011-1

         bcflag_save=bcflag
         ris = basis17(x)
         bcflag=bcflag_save

      case(53)                  !0110

         ris = -basis2_1(x)**2/2d0 - basis3_3(x)*llx 
     &+ (pi**2*ll1mx*llx)/6d0 - basis2_1(x)*ll1mx*llx 
     &- (ll1mx**2*llx**2)/2d0 + llx*zeta3

      case(54)                  !0111

         ris = pi**4/90d0 - basis3(x) + basis3_3(x)*ll1mx 
     &- (pi**2*ll1mx**2)/12d0 + (basis2_1(x)*ll1mx**2)/2d0 
     &+ (ll1mx**3*llx)/3d0

      case(55)                  !1-1-1-1

         ris = cli4pt5 - basis7(x) 
     &+ basis3_5(x)
     &*ll1px - (pi**2*ll1px**2)/12d0 + (basis2_3(x)
     &*ll1px**2)/2d0 +(ll2**2*ll1px**2)/2d0 - (ll2
     &*ll1mx*ll1px**2)/2d0 - (ll2*ll1px**3)/3d0 
     &+ (ll1mx*ll1px**3)/3d0

      case(56)                  !1-1-10

         bcflag_save=bcflag
         b18 = basis18(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = -pi**4/480d0 + basis2_2(x)**2/2d0 
     &+ b18 
     &+ pi**4/480d0 - (pi**2*basis2_2(x))/12d0 
     &+ basis2_3(x)*basis2_2(x) - basis2_2(x)**2/2d0 
     &- 3*cli4pt5 + basis3(x) 
     &+ basis2(x) + basis1(x) + basis15(x)/2d0 
     &- basis9(x)/2d0 + basis10(x)/2d0 
     &+ 2*basis6(x) - 2*basis11(x) 
     &+ 3*basis7(x) - basis13(x)/4d0 + (basis2_2(x)
     &*ll2**2)/2d0 + basis3_4(x)*ll1mx - basis2_2(x)
     &*ll2*ll1mx - 2*basis3_1(x)*ll1px 
     &- 2*basis3_8(x)*ll1px - 2*basis3_7(x)
     &*ll1px + (pi**2*ll2*ll1px)/4d0 - (ll2**3
     &*ll1px)/2d0 + (5*pi**2*ll1mx*ll1px)/12d0 
     &+ basis2_2(x)*ll1mx*ll1px - ll2*ll1mx**2
     &*ll1px + (ll1mx**3*ll1px)/3d0 - ll1mx**2
     &*llx*ll1px - (5*pi**2*ll1px**2)/12d0 
     &+ (basis2_3(x)*ll1px**2)/2d0 - (basis2_2(x)
     &*ll1px**2)/2d0 + (basis2_1(x)*ll1px**2)/2d0 
     &+ ll2**2*ll1px**2 + (3*ll2*ll1mx
     &*ll1px**2)/2d0 + 2*ll1mx*llx*ll1px**2 
     &- (3*ll2*ll1px**3)/2d0 - (2*ll1mx
     &*ll1px**3)/3d0 - llx*ll1px**3 
     &+ (5*ll1px**4)/8d0 - (7*ll1mx*zeta3)/8d0 
     &- (5*ll1px*zeta3)/4d0
     &+ 3*cli4pt5 - basis3(x) 
     &- basis2(x) - basis1(x) - basis15(x)/2d0 
     &+ basis9(x)/2d0 - basis10(x)/2d0 
     &- 2*basis6(x) + 2*basis11(x) 
     &- 3*basis7(x) + basis13(x)/4d0 
     &+ basis3_5(x)*llx + (pi**2*ll2*llx)/12d0 
     &- (ll2**3*llx)/6d0 + 2*basis3_1(x)*ll1px 
     &+ 2*basis3_8(x)*ll1px + 2*basis3_7(x)
     &*ll1px - (pi**2*ll2*ll1px)/4d0 + (ll2**3
     &*ll1px)/2d0 - (pi**2*ll1mx*ll1px)/4d0 
     &+ ll2*ll1mx**2*ll1px - (ll1mx**3
     &*ll1px)/3d0 - (pi**2*llx*ll1px)/6d0 
     &+ basis2_3(x)*llx*ll1px + ll2**2*llx
     &*ll1px - ll2*ll1mx*llx*ll1px 
     &+ ll1mx**2*llx*ll1px + (5*pi**2*ll1px**2)
     &/12d0 - (basis2_3(x)*ll1px**2)/2d0 + (basis2_2(x)
     &*ll1px**2)/2d0 - (basis2_1(x)*ll1px**2)/2d0
     &-ll2**2
     &*ll1px**2 - (3*ll2*ll1mx*ll1px**2)/2d0 
     &- (ll2*llx*ll1px**2)/2d0 - ll1mx*llx
     &*ll1px**2 + (3*ll2*ll1px**3)/2d0 + (ll1mx
     &*ll1px**3)/2d0 + llx*ll1px**3 - (5*ll1px**4)
     &/8d0 - (ll1mx*zeta3)/8d0 - (7*llx*zeta3)/8d0 
     &+ (5*ll1px*zeta3)/4d0

      case(57)                  !1-1-11

         ris = -pi**4/288d0 + (pi**2*basis2_3(x))/12d0 
     &- basis2_3(x)**2/2d0 + (pi**2*ll2**2)/24d0 
     &- (basis2_3(x)*ll2**2)/2d0 - ll2**4/8d0 
     &- basis3_5(x)*ll1mx - (pi**2*ll2*ll1mx)/6d0 
     &+ basis2_3(x)*ll2*ll1mx + (2*ll2**3
     &*ll1mx)/3d0 - (ll2**2*ll1mx**2)/2d0 
     &+ (pi**2*ll1mx*ll1px)/6d0 - basis2_3(x)
     &*ll1mx*ll1px - ll2**2*ll1mx*ll1px 
     &+ ll2*ll1mx**2*ll1px + (ll2*ll1mx
     &*ll1px**2)/2d0 - (ll1mx**2*ll1px**2)/2d0 
     &+ (7*ll1mx*zeta3)/8d0

      case(58)                  !1-10-1

         bcflag_save=bcflag
         b18 = basis18(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = (11*pi**4)/720d0 - basis2_2(x)**2/2d0 
     &- b18 
     &- (pi**4/480d0 - (pi**2*basis2_2(x))/12d0 
     &+ basis2_3(x)*basis2_2(x) - basis2_2(x)**2/2d0 
     &- 3*cli4pt5 + basis3(x) 
     &+ basis2(x) + basis1(x) + basis15(x)/2d0 
     &- basis9(x)/2d0 + basis10(x)/2d0 
     &+ 2*basis6(x) - 2*basis11(x) 
     &+ 3*basis7(x) - basis13(x)/4d0 + (basis2_2(x)
     &*ll2**2)/2d0 + basis3_4(x)*ll1mx - basis2_2(x)
     &*ll2*ll1mx - 2*basis3_1(x)*ll1px 
     &- 2*basis3_8(x)*ll1px - 2*basis3_7(x)
     &*ll1px + (pi**2*ll2*ll1px)/4d0 - (ll2**3
     &*ll1px)/2d0 + (5*pi**2*ll1mx*ll1px)/12d0 
     &+ basis2_2(x)*ll1mx*ll1px - ll2*ll1mx**2
     &*ll1px + (ll1mx**3*ll1px)/3d0 - ll1mx**2
     &*llx*ll1px - (5*pi**2*ll1px**2)/12d0 
     &+ (basis2_3(x)*ll1px**2)/2d0 - (basis2_2(x)
     &*ll1px**2)/2d0 + (basis2_1(x)*ll1px**2)/2d0 
     &+ ll2**2*ll1px**2 + (3*ll2*ll1mx
     &*ll1px**2)/2d0 + 2*ll1mx*llx*ll1px**2 
     &- (3*ll2*ll1px**3)/2d0 - (2*ll1mx
     &*ll1px**3)/3d0 - llx*ll1px**3 
     &+ (5*ll1px**4)/8d0 - (7*ll1mx*zeta3)/8d0 
     &- (5*ll1px*zeta3)/4d0)
     &+ 2*basis3(x) + (3*basis2(x))/2d0 
     &+ basis1(x)/2d0 + basis15(x)/4d0 - basis4(x) 
     &- basis9(x) + basis10(x) 
     &+ 2*basis6(x) - basis11(x) - basis13(x)/2d0 
     &- basis3_6(x)*ll1px + basis3_3(x)*ll1px 
     &- basis3_2(x)*ll1px - basis3_1(x)*ll1px 
     &- basis3_8(x)*ll1px - 2*basis3_7(x)
     &*ll1px - (pi**2*ll2*ll1px)/12d0 + (ll2**3
     &*ll1px)/6d0 + (pi**2*ll1mx*ll1px)/12d0 
     &- (ll2**2*ll1mx*ll1px)/2d0 - (ll2
     &*ll1mx**2*ll1px)/2d0 + (ll1mx**3*ll1px)
     &/6d0 - (ll1mx**2*llx*ll1px)/2d0 - (pi**2
     &*ll1px**2)/8d0 + (basis2_3(x)*ll1px**2)/2d0 
     &- (basis2_2(x)*ll1px**2)/2d0 
     &+(basis2_1(x)*ll1px**2)/2d0 
     &+ (ll2**2*ll1px**2)/4d0 + (3*ll2*ll1mx
     &*ll1px**2)/2d0 + ll1mx*llx*ll1px**2 
     &- ll2*ll1px**3 - (ll1mx*ll1px**3)/2d0 
     &- (5*llx*ll1px**3)/6d0 + (11*ll1px**4)/24d0 
     &+ (ll1mx*zeta3)/4d0 - (3*ll1px*zeta3)/8d0

      case(59)                  !1-100

         ris = pi**4/72d0 + 2*cli4pt5 + basis3(x) 
     &+ (3*basis2(x))/2d0 - basis1(x)/2d0 - basis5(x) 
     &- basis15(x)/4d0 - 2*basis4(x) 
     &- 2*basis6(x) + 2*basis11(x) 
     &- 2*basis7(x) - basis13(x)/4d0 
     &+ basis14(x)/4d0 - (pi**2*ll1mx**2)/16d0 
     &- ll1mx**4/32d0 - basis3_6(x)*llx + basis3_3(x)
     &*llx - basis3_2(x)*llx + basis3_1(x)*llx 
     &+ basis3_8(x)*llx - (pi**2*ll2*llx)/12d0 
     &+ (ll2**3*llx)/6d0 - (pi**2*ll1mx*llx)/12d0 
     &- (ll2**2*ll1mx*llx)/2d0 + (ll2*ll1mx**2
     &*llx)/2d0 - (ll1mx**3*llx)/12d0 - (pi**2*llx**2)
     &/24d0 + (basis2_3(x)*llx**2)/2d0 + (ll2**2
     &*llx**2)/4d0 - (ll2*ll1mx*llx**2)/2d0 
     &+ (ll1mx**2*llx**2)/2d0 - (pi**2*ll2*ll1px)
     &/6d0 + (ll2**3*ll1px)/3d0 + (pi**2*ll1mx
     &*ll1px)/24d0 + (ll1mx**3*ll1px)/24d0 
     &- (ll1mx**2*llx*ll1px)/4d0 + (13*pi**2
     &*ll1px**2)/48d0 - (ll2**2*ll1px**2)/2d0 
     &+ (ll1mx**2*ll1px**2)/16d0 - (ll1mx*llx
     &*ll1px**2)/4d0 + (ll2*ll1px**3)/3d0 + (ll1mx
     &*ll1px**3)/24d0 + (llx*ll1px**3)/4d0 
     &- (23*ll1px**4)/96d0 - (3*ll1mx*zeta3)/4d0 
     &- (llx*zeta3)/8d0

      case(60)                  !1-101

         bcflag_save=bcflag
         b17= basis17(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)


         ris = -pi**4/72d0 + basis2_2(x)*basis2_1(x)
     &- basis2_1(x)**2/2d0 
     &-(b17-( -pi**4/288d0 
     &- (pi**2*basis2_1(x))/12d0 + basis2_3(x)
     &*basis2_1(x) - basis2_2(x)*basis2_1(x) 
     &+ basis2_1(x)**2/2d0 
     &- 5*cli4pt5 
     &+ basis8(x) - 2*basis3(x) + basis2(x) + basis1(x) 
     &- 2*basis12(x) + basis15(x)/2d0 
     &+ 3*basis4(x) - basis9(x)/2d0 
     &+ basis10(x)/2d0 + 2*basis6(x) 
     &- 4*basis11(x) + 4*basis7(x) 
     &+ basis13(x)/4d0 + (basis2_1(x)*ll2**2)/2d0 
     &- 2*basis3_8(x)*ll1mx + (pi**2*ll2
     &*ll1mx)/12d0 - basis2_1(x)*ll2*ll1mx 
     &- (ll2**3*ll1mx)/6d0 + (pi**2*ll1mx**2)/8d0 
     &+ (basis2_3(x)*ll1mx**2)/2d0 - (basis2_2(x)
     &*ll1mx**2)/2d0 + (basis2_1(x)*ll1mx**2)/2d0 
     &+ (ll2**2*ll1mx**2)/2d0-(2*ll2*ll1mx**3)/3d0 
     &+ (7*ll1mx**4)/24d0 - (ll1mx**3*llx)/3d0 
     &- basis3_3(x)*ll1px - 2*basis3_1(x)*ll1px 
     &+ (pi**2*ll2*ll1px)/3d0-(2*ll2**3*ll1px)/3d0 
     &- (3*pi**2*ll1px**2)/8d0 + ll2**2*ll1px**2 
     &- (2*ll2*ll1px**3)/3d0 - (llx*ll1px**3)/3d0 
     &+ (3*ll1px**4)/8d0 + (7*ll1mx*zeta3)/4d0 
     &- (5*ll1px*zeta3)/8d0))
     &+ 4*cli4pt5 
     &- 2*basis8(x) + basis3(x) - (3*basis2(x))/2d0 
     &- basis1(x)/2d0 + 2*basis5(x) + basis12(x) 
     &- basis15(x)/4d0 + 2*basis11(x) 
     &- 2*basis7(x) + basis13(x)/4d0 
     &- basis14(x)/4d0 + basis3_6(x)*ll1mx 
     &- basis3_3(x)*ll1mx + basis3_2(x)*ll1mx 
     &- basis3_1(x)
     &*ll1mx + basis3_8(x)*ll1mx - (pi**2*ll2
     &*ll1mx)/12d0 + (ll2**3*ll1mx)/6d0 + (5*pi**2
     &*ll1mx**2)/48d0 - (basis2_3(x)*ll1mx**2)/2d0 
     &+ (basis2_2(x)*ll1mx**2)/2d0 
     &-(basis2_1(x)*ll1mx**2)/2d0 
     &- (ll2**2*ll1mx**2)/4d0 + (ll2*ll1mx**3)/3d0 
     &- (3*ll1mx**4)/32d0 - (ll1mx**3*llx)/4d0 
     &+ 2*basis3_1(x)*ll1px - (pi**2*ll2*ll1px)/6d0 
     &+ (ll2**3*ll1px)/3d0 - (pi**2*ll1mx*ll1px)
     &/24d0 - (ll1mx**3*ll1px)/24d0 + (ll1mx**2
     &*llx*ll1px)/4d0 + (pi**2*ll1px**2)/16d0 
     &- (ll2**2*ll1px**2)/2d0 - (ll1mx**2
     &*ll1px**2)/16d0 + (ll1mx*llx*ll1px**2)/4d0 
     &+ (ll2*ll1px**3)/3d0 - (ll1mx*ll1px**3)/24d0 
     &+ (llx*ll1px**3)/12d0 - (3*ll1px**4)/32d0 
     &+ (5*ll1mx*zeta3)/8d0 + (3*ll1px*zeta3)/2d0

      case(61)                  !1-11-1

         ris = (-23*pi**4)/1440d0 - (pi**2*basis2_3(x))/12d0 
     &+ basis2_3(x)**2/2d0 - 2*basis8(x) 
     &- 2*basis10(x) + 2*basis7(x) 
     &- (pi**2*ll2**2)/24d0 + (basis2_3(x)*ll2**2)/2d0 
     &+ ll2**4/8d0 - (pi**2*ll2*ll1mx)/12d0 
     &- basis2_3(x)*ll2*ll1mx - (ll2**3
     &*ll1mx)/6d0 + 2*basis3_6(x)*ll1px 
     &+ (pi**2*ll2*ll1px)/3d0-(2*ll2**3*ll1px)/3d0 
     &+ ll2**2*ll1mx*ll1px - (pi**2*ll1px**2)/6d0 
     &+ (ll2**2*ll1px**2)/2d0 - ll2*ll1mx
     &*ll1px**2 + (ll1mx*ll1px**3)/3d0 
     &- ll1px**4/12d0 + (ll1mx*zeta3)/4d0 
     &- 2*ll1px*zeta3

      case(62)                  !1-110

         bcflag_save=bcflag
         b17= basis17(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = -pi**4/72d0 - basis2_2(x)*basis2_1(x)
     &+ basis2_1(x)**2/2d0 
     &+ b17-( -pi**4/288d0 
     &- (pi**2*basis2_1(x))/12d0 + basis2_3(x)
     &*basis2_1(x) - basis2_2(x)*basis2_1(x) +basis2_1(x)**2/2d0 
     &- 5*cli4pt5 
     &+ basis8(x) - 2*basis3(x) + basis2(x) + basis1(x) 
     &- 2*basis12(x) + basis15(x)/2d0 
     &+ 3*basis4(x) - basis9(x)/2d0 
     &+ basis10(x)/2d0 + 2*basis6(x) 
     &- 4*basis11(x) + 4*basis7(x) 
     &+ basis13(x)/4d0 + (basis2_1(x)*ll2**2)/2d0 
     &- 2*basis3_8(x)*ll1mx + (pi**2*ll2
     &*ll1mx)/12d0 - basis2_1(x)*ll2*ll1mx 
     &- (ll2**3*ll1mx)/6d0 + (pi**2*ll1mx**2)/8d0 
     &+ (basis2_3(x)*ll1mx**2)/2d0 - (basis2_2(x)
     &*ll1mx**2)/2d0 + (basis2_1(x)*ll1mx**2)/2d0 
     &+ (ll2**2*ll1mx**2)/2d0-(2*ll2*ll1mx**3)/3d0 
     &+ (7*ll1mx**4)/24d0 - (ll1mx**3*llx)/3d0 
     &- basis3_3(x)*ll1px - 2*basis3_1(x)*ll1px 
     &+ (pi**2*ll2*ll1px)/3d0-(2*ll2**3*ll1px)/3d0 
     &- (3*pi**2*ll1px**2)/8d0 + ll2**2*ll1px**2 
     &- (2*ll2*ll1px**3)/3d0 - (llx*ll1px**3)/3d0 
     &+ (3*ll1px**4)/8d0 + (7*ll1mx*zeta3)/4d0 
     &- (5*ll1px*zeta3)/8d0) - 6*cli4pt5 
     &- 3*basis3(x) + basis2(x)/2d0 + (3*basis1(x))/2d0 
     &+ 2*basis5(x) - 3*basis12(x) 
     &+ (3*basis15(x))/4d0 + 4*basis4(x) 
     &+ 4*basis6(x) - 6*basis11(x) 
     &+ 6*basis7(x) + basis13(x)/4d0 
     &- basis14(x)/4d0 - 2*basis3_8(x)
     &*ll1mx + (3*pi**2*ll1mx**2)/16d0 
     &+ (basis2_3(x)*ll1mx**2)/2d0 - (basis2_2(x)
     &*ll1mx**2)/2d0 + (basis2_1(x)*ll1mx**2)/2d0 
     &+ (ll2**2*ll1mx**2)/4d0 - (ll2*ll1mx**3)/2d0 
     &+ (31*ll1mx**4)/96d0 + 2*basis3_6(x)*llx 
     &+ (pi**2*ll2*llx)/6d0 - (ll2**3*llx)/3d0 
     &- (pi**2*ll1mx*llx)/12d0 - basis2_3(x)
     &*ll1mx*llx + (ll2**2*ll1mx*llx)/2d0 
     &- (5*ll1mx**3*llx)/12d0 - 2*basis3_1(x)*ll1px 
     &+ (pi**2*ll2*ll1px)/2d0 - ll2**3*ll1px 
     &- (pi**2*ll1mx*ll1px)/24d0 - (ll1mx**3
     &*ll1px)/24d0 + (ll1mx**2*llx*ll1px)/4d0 
     &- (29*pi**2*ll1px**2)/48d0 + (3*ll2**2
     &*ll1px**2)/2d0 - (ll1mx**2*ll1px**2)/16d0 
     &+ (ll1mx*llx*ll1px**2)/4d0 - ll2
     &*ll1px**3 - (ll1mx*ll1px**3)/24d0 
     &- (7*llx*ll1px**3)/12d0 + (55*ll1px**4)/96d0 
     &+ (11*ll1mx*zeta3)/4d0 - (7*llx*zeta3)/4d0 
     &- (3*ll1px*zeta3)/2d0

      case(63)                  !1-111

         ris = -3*cli4pt5 + 3*basis8(x) 
     &- 2*basis3_6(x)*ll1mx + (pi**2*ll2
     &*ll1mx)/12d0 - (ll2**3*ll1mx)/6d0 
     &+ (basis2_3(x)*ll1mx**2)/2d0 
     &- (7*ll1mx*zeta3)/8d0

      case(64)                  !10-1-1

         ris = -pi**4/480d0 - basis3(x) + basis9(x)/2d0 
     &- basis10(x)/2d0 + basis13(x)/4d0 
     &+ basis3_6(x)*ll1px - basis3_3(x)*ll1px 
     &- basis3_4(x)*ll1px + basis3_7(x)
     &*ll1px + basis3_5(x)*ll1px + (pi**2*ll2
     &*ll1px)/6d0 - (ll2**3*ll1px)/3d0 + (ll2**2
     &*ll1mx*ll1px)/2d0 - (pi**2*ll1px**2)/6d0 
     &- (basis2_1(x)*ll1px**2)/2d0
     &+(ll2**2*ll1px**2)/2d0 
     &- ll2*ll1mx*ll1px**2 - (ll1mx*llx
     &*ll1px**2)/2d0 + (ll1mx*ll1px**3)/2d0 
     &- (ll1mx*zeta3)/8d0 - (ll1px*zeta3)/8d0

      case(65)                  !10-10

         bcflag_save=bcflag
         b16 = basis16(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = -pi**4/180d0 + b16 - 2*basis3(x) 
     &- 2*basis2(x) - 2*basis1(x) + 2*basis5(x) 
     &+ 2*basis4(x) + 2*basis6(x) + basis13(x)/2d0 
     &- basis14(x)/2d0 + (pi**2*ll1mx**2)/8d0 
     &+ ll1mx**4/16d0 + basis3_6(x)*llx - basis3_3(x)
     &*llx - basis3_4(x)*llx + basis3_7(x)
     &*llx + basis3_5(x)*llx + (pi**2*ll2*llx)
     &/6d0 - (ll2**3*llx)/3d0 - (pi**2*ll1mx*llx)/12d0 
     &+ (ll2**2*ll1mx*llx)/2d0 - (ll1mx**3
     &*llx)/6d0 + 2*basis3_1(x)*ll1px - (pi**2*ll1mx
     &*ll1px)/12d0 - (ll1mx**3*ll1px)/12d0 
     &- (pi**2*llx*ll1px)/12d0 - basis2_1(x)*llx*ll1px 
     &+ (ll2**2*llx*ll1px)/2d0 - ll2*ll1mx
     &*llx*ll1px + (ll1mx**2*llx*ll1px)/2d0 
     &- ll1mx*llx**2*ll1px - (5*pi**2*ll1px**2)
     &/24d0 - (ll1mx**2*ll1px**2)/8d0 + ll1mx*llx
     &*ll1px**2 - (ll1mx*ll1px**3)/12d0 - (llx
     &*ll1px**3)/6d0 + (7*ll1px**4)/48d0 + (3*ll1mx
     &*zeta3)/2d0 - (3*llx*zeta3)/4d0 + (3*ll1px*zeta3)/2d0

      case(66)                  !10-11

         bcflag_save=bcflag
         b17= basis17(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = pi**4/72d0 - (pi**2*basis2_1(x))/12d0 + basis2_3(x)
     &*basis2_1(x) - basis2_2(x)*basis2_1(x) +basis2_1(x)**2/2d0 
     &+ b17-( -pi**4/288d0 
     &- (pi**2*basis2_1(x))/12d0 + basis2_3(x)
     &*basis2_1(x) - basis2_2(x)*basis2_1(x) +basis2_1(x)**2/2d0 
     &- 5*cli4pt5 
     &+ basis8(x) - 2*basis3(x) + basis2(x) + basis1(x) 
     &- 2*basis12(x) + basis15(x)/2d0 
     &+ 3*basis4(x) - basis9(x)/2d0 
     &+ basis10(x)/2d0 + 2*basis6(x) 
     &- 4*basis11(x) + 4*basis7(x) 
     &+ basis13(x)/4d0 + (basis2_1(x)*ll2**2)/2d0 
     &- 2*basis3_8(x)*ll1mx + (pi**2*ll2
     &*ll1mx)/12d0 - basis2_1(x)*ll2*ll1mx 
     &- (ll2**3*ll1mx)/6d0 + (pi**2*ll1mx**2)/8d0 
     &+ (basis2_3(x)*ll1mx**2)/2d0 - (basis2_2(x)
     &*ll1mx**2)/2d0 + (basis2_1(x)*ll1mx**2)/2d0 
     &+ (ll2**2*ll1mx**2)/2d0-(2*ll2*ll1mx**3)/3d0 
     &+ (7*ll1mx**4)/24d0 - (ll1mx**3*llx)/3d0 
     &- basis3_3(x)*ll1px - 2*basis3_1(x)*ll1px 
     &+ (pi**2*ll2*ll1px)/3d0-(2*ll2**3*ll1px)/3d0 
     &- (3*pi**2*ll1px**2)/8d0 + ll2**2*ll1px**2 
     &- (2*ll2*ll1px**3)/3d0 - (llx*ll1px**3)/3d0 
     &+ (3*ll1px**4)/8d0 + (7*ll1mx*zeta3)/4d0 
     &- (5*ll1px*zeta3)/8d0)
     &- 4*cli4pt5 + 2*basis8(x)
     & - basis3(x) + (3*basis2(x))/2d0 + basis1(x)/2d0 
     &- 2*basis5(x) - basis12(x) 
     &+ basis15(x)/4d0 - 2*basis11(x) 
     &+ 2*basis7(x) - basis13(x)/4d0 
     &+ basis14(x)/4d0 + (basis2_1(x)*ll2**2)/2d0 
     &- basis3_6(x)*ll1mx + basis3_3(x)*ll1mx 
     &- 2*basis3_8(x)*ll1mx + basis3_4(x)
     &*ll1mx - basis3_7(x)*ll1mx 
     &- basis3_5(x)*ll1mx - basis2_1(x)*ll2*ll1mx 
     &+ (pi**2*ll1mx**2)/16d0 + (basis2_3(x)
     &*ll1mx**2)/2d0 - (basis2_2(x)*ll1mx**2)/2d0 
     &+ (basis2_1(x)*ll1mx**2)/2d0
     &+(ll2**2*ll1mx**2)/4d0 
     &- (5*ll2*ll1mx**3)/6d0 + (25*ll1mx**4)/96d0 
     &- (ll1mx**3*llx)/4d0 - 2*basis3_1(x)*ll1px 
     &+ (pi**2*ll2*ll1px)/6d0 - (ll2**3*ll1px)/3d0 
     &+ (pi**2*ll1mx*ll1px)/8d0 + basis2_1(x)*ll1mx
     &*ll1px - (ll2**2*ll1mx*ll1px)/2d0 
     &+ ll2*ll1mx**2*ll1px + (ll1mx**3
     &*ll1px)/24d0 + (3*ll1mx**2*llx*ll1px)/4d0 
     &- (pi**2*ll1px**2)/16d0 + (ll2**2*ll1px**2)/2d0 
     &- (7*ll1mx**2*ll1px**2)/16d0 - (ll1mx*llx
     &*ll1px**2)/4d0 - (ll2*ll1px**3)/3d0 
     &+ (ll1mx*ll1px**3)/24d0 - (llx*ll1px**3)
     &/12d0 + (3*ll1px**4)/32d0 + (ll1mx*zeta3)/4d0 
     &- (3*ll1px*zeta3)/2d0

      case(67)                  !100-1

         bcflag_save=bcflag
         b16 = basis16(x)
         bcflag=bcflag_save         

         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = pi**4/360d0 - b16 + basis3(x) + basis2(x) 
     &+ basis1(x) - basis5(x) - basis4(x) 
     &- basis6(x) - basis13(x)/4d0 + basis14(x)
     &/4d0 - (pi**2*ll1mx**2)/16d0 - ll1mx**4/32d0 
     &+ (ll1mx**3*llx)/12d0 - basis3_1(x)*ll1px 
     &+ (pi**2*ll1mx*ll1px)/24d0 + (ll1mx**3
     &*ll1px)/24d0 - (ll1mx**2*llx*ll1px)/4d0 
     &+ (5*pi**2*ll1px**2)/48d0 + (ll1mx**2*ll1px**2)
     &/16d0 - (ll1mx*llx*ll1px**2)/4d0 + (ll1mx
     &*ll1px**3)/24d0 + (llx*ll1px**3)/12d0 
     &- (7*ll1px**4)/96d0 - (3*ll1mx*zeta3)/4d0 
     &- (3*ll1px*zeta3)/4d0

      case(68)                  !1000

         ris = -basis1(x)+basis3_1(x)*llx
     &-(basis2_1(x)*llx**2)/2d0 
     &- (ll1mx*llx**3)/6d0

      case(69)                  !1001

         ris = -basis2_1(x)**2/2d0 - basis3_1(x)*ll1mx

      case(70)                  !101-1

         bcflag_save=bcflag
         b17= basis17(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = -pi**4/144d0 + (pi**2*basis2_1(x))/12d0 - basis2_3(x)
     &*basis2_1(x) + basis2_2(x)*basis2_1(x) -basis2_1(x)**2/2d0 
     &+ 4*cli4pt5 
     &- (b17-( -pi**4/288d0 
     &- (pi**2*basis2_1(x))/12d0 + basis2_3(x)
     &*basis2_1(x) - basis2_2(x)*basis2_1(x) +basis2_1(x)**2/2d0 
     &- 5*cli4pt5 
     &+ basis8(x) - 2*basis3(x) + basis2(x) + basis1(x) 
     &- 2*basis12(x) + basis15(x)/2d0 
     &+ 3*basis4(x) - basis9(x)/2d0 
     &+ basis10(x)/2d0 + 2*basis6(x) 
     &- 4*basis11(x) + 4*basis7(x) 
     &+ basis13(x)/4d0 + (basis2_1(x)*ll2**2)/2d0 
     &- 2*basis3_8(x)*ll1mx + (pi**2*ll2
     &*ll1mx)/12d0 - basis2_1(x)*ll2*ll1mx 
     &- (ll2**3*ll1mx)/6d0 + (pi**2*ll1mx**2)/8d0 
     &+ (basis2_3(x)*ll1mx**2)/2d0 - (basis2_2(x)
     &*ll1mx**2)/2d0 + (basis2_1(x)*ll1mx**2)/2d0 
     &+ (ll2**2*ll1mx**2)/2d0-(2*ll2*ll1mx**3)/3d0 
     &+ (7*ll1mx**4)/24d0 - (ll1mx**3*llx)/3d0 
     &- basis3_3(x)*ll1px - 2*basis3_1(x)*ll1px 
     &+ (pi**2*ll2*ll1px)/3d0-(2*ll2**3*ll1px)/3d0 
     &- (3*pi**2*ll1px**2)/8d0 + ll2**2*ll1px**2 
     &- (2*ll2*ll1px**3)/3d0 - (llx*ll1px**3)/3d0 
     &+ (3*ll1px**4)/8d0 + (7*ll1mx*zeta3)/4d0 
     &- (5*ll1px*zeta3)/8d0))
     &- 2*basis8(x) + basis3(x) - (3*basis2(x))/2d0 
     &- basis1(x)/2d0 + 2*basis5(x) + basis12(x) 
     &- basis15(x)/4d0 - 2*basis4(x) 
     &+ basis9(x) - basis10(x) 
     &+ 2*basis11(x) - 2*basis7(x) 
     &- basis13(x)/4d0 - basis14(x)/4d0 
     &- (basis2_1(x)*ll2**2)/2d0 + 2*basis3_8(x)*ll1mx 
     &-(pi**2*ll2*ll1mx)/6d0
     &+basis2_1(x)*ll2*ll1mx 
     &+ (ll2**3*ll1mx)/3d0 - (pi**2*ll1mx**2)/16d0 
     &- (basis2_3(x)*ll1mx**2)/2d0 + (basis2_2(x)
     &*ll1mx**2)/2d0 - (basis2_1(x)*ll1mx**2)/2d0 
     &- (3*ll2**2*ll1mx**2)/4d0 + (5*ll2*ll1mx**3)
     &/6d0 - (25*ll1mx**4)/96d0 + (ll1mx**3*llx)/4d0 
     &+ 2*basis3_3(x)*ll1px + 2*basis3_1(x)*ll1px - (pi**2
     &*ll2*ll1px)/6d0 + (ll2**3*ll1px)/3d0 
     &- (pi**2*ll1mx*ll1px)/24d0 - (ll1mx**3
     &*ll1px)/24d0 + (ll1mx**2*llx*ll1px)/4d0 
     &+ (7*pi**2*ll1px**2)/48d0 - (ll2**2*ll1px**2)/2d0 
     &- (ll1mx**2*ll1px**2)/16d0 + (ll1mx*llx
     &*ll1px**2)/4d0 + (ll2*ll1px**3)/3d0 - (ll1mx
     &*ll1px**3)/24d0 + (llx*ll1px**3)/12d0 
     &- (17*ll1px**4)/96d0 - (3*ll1mx*zeta3)/4d0 
     &- (ll1px*zeta3)/4d0

      case(71)                  !1010

         ris = pi**4/45d0 + basis2_1(x)**2/2d0 - 2*basis3(x) 
     &+ 2*basis1(x) + 2*basis5(x) + (pi**2*ll1mx**2)/6d0 
     &+ ll1mx**4/12d0 + 2*basis3_3(x)*llx - (pi**2
     &*ll1mx*llx)/3d0 + basis2_1(x)*ll1mx*llx 
     &- (ll1mx**3*llx)/3d0 + ll1mx**2*llx**2 
     &+ 2*ll1mx*zeta3 - 2*llx*zeta3

      case(72)                  !1011

         ris = -pi**4/30d0 + 3*basis3(x) - 2*basis3_3(x)*ll1mx 
     &+ (pi**2*ll1mx**2)/12d0 - (basis2_1(x)*ll1mx**2)/2d0 
     &- (ll1mx**3*llx)/2d0 - ll1mx*zeta3

      case(73)                  !11-1-1

         ris = (7*pi**4)/720d0 + basis8(x) 
     &+ basis10(x) - basis7(x) + (pi**2*ll2
     &*ll1mx)/12d0 - (ll2**3*ll1mx)/6d0 + (ll2**2
     &*ll1mx**2)/4d0 - basis3_6(x)*ll1px 
     &- (pi**2*ll2*ll1px)/6d0 + (ll2**3*ll1px)/3d0 
     &- (ll2**2*ll1mx*ll1px)/2d0 + (pi**2
     &*ll1px**2)/12d0 - (ll2**2*ll1px**2)/4d0 
     &+ (ll2*ll1mx*ll1px**2)/2d0 - (ll1mx
     &*ll1px**3)/6d0 + ll1px**4/24d0 - (ll1mx*zeta3)
     &/8d0 + ll1px*zeta3

      case(74)                  !11-10

         ris = pi**4/72d0 + cli4pt5 + basis8(x) 
     &+ basis3(x) + basis2(x)/2d0 - basis1(x)/2d0 
     &- 2*basis5(x) + basis12(x) 
     &- basis15(x)/4d0 - 2*basis4(x) 
     &- 2*basis6(x) + 2*basis11(x) 
     &- 2*basis7(x) - basis13(x)/4d0 
     &+ basis14(x)/4d0 + (pi**2*ll2*ll1mx)/12d0 
     &- (ll2**3*ll1mx)/6d0 - (5*pi**2*ll1mx**2)/48d0 
     &+ (ll2**2*ll1mx**2)/4d0 - (ll2*ll1mx**3)/6d0 
     &- ll1mx**4/32d0 - basis3_6(x)*llx - (pi**2
     &*ll2*llx)/12d0 + (ll2**3*llx)/6d0 + (pi**2
     &*ll1mx*llx)/12d0 - (ll2**2*ll1mx*llx)/2d0 
     &+ (ll2*ll1mx**2*llx)/2d0 + (ll1mx**3*llx)
     &/12d0 - (pi**2*ll2*ll1px)/6d0+(ll2**3*ll1px)
     &/3d0 + (pi**2*ll1mx*ll1px)/24d0 + (ll1mx**3
     &*ll1px)/24d0 - (ll1mx**2*llx*ll1px)/4d0 
     &+ (13*pi**2*ll1px**2)/48d0 - (ll2**2*ll1px**2)/2d0 
     &+ (ll1mx**2*ll1px**2)/16d0 - (ll1mx*llx
     &*ll1px**2)/4d0 + (ll2*ll1px**3)/3d0 + (ll1mx
     &*ll1px**3)/24d0 + (llx*ll1px**3)/4d0 
     &- (23*ll1px**4)/96d0 - (13*ll1mx*zeta3)/8d0 
     &+ (7*llx*zeta3)/8d0

      case(75)                  !11-11
         
         ris = 3*cli4pt5 - 3*basis8(x) 
     &+ basis3_6(x)*ll1mx - (pi**2*ll2*ll1mx)/6d0 
     &+ (ll2**3*ll1mx)/3d0 + (pi**2*ll1mx**2)/24d0 
     &- (ll2**2*ll1mx**2)/4d0 + (7*ll1mx*zeta3)/4d0

      case(76)                  !110-1

         ris = -pi**4/288d0 + basis4(x) 
     &- basis9(x)/2d0 + basis10(x)/2d0 
     &+ basis13(x)/4d0 + (pi**2*ll1mx**2)/24d0 
     &- basis3_3(x)*ll1px - (pi**2*ll1px**2)/24d0 
     &+ ll1px**4/24d0 + (5*ll1mx*zeta3)/8d0 
     &+ (7*ll1px*zeta3)/8d0

      case(77)                  !1100

         ris = -pi**4/90d0 + basis3(x) - basis1(x) - basis5(x) 
     &- (pi**2*ll1mx**2)/12d0 - ll1mx**4/24d0 - basis3_3(x)
     &*llx + (pi**2*ll1mx*llx)/6d0 + (ll1mx**3
     &*llx)/6d0 - (ll1mx**2*llx**2)/4d0 
     &- ll1mx*zeta3 + llx*zeta3

      case(78)                  !1101

         ris = pi**4/30d0 - 3*basis3(x) + basis3_3(x)*ll1mx 
     &+ (pi**2*ll1mx**2)/12d0 + 2*ll1mx*zeta3

      case(79)                  !111-1

         ris = -cli4pt5 + basis8(x) 
     &+ (pi**2*ll2*ll1mx)/12d0 -(ll2**3*ll1mx)/6d0 
     &- (pi**2*ll1mx**2)/24d0 + (ll2**2*ll1mx**2)/4d0 
     &- (ll2*ll1mx**3)/6d0 - (7*ll1mx*zeta3)/8d0

      case(80)                  !1110

         ris = -pi**4/90d0 + basis3(x) - (pi**2*ll1mx**2)/12d0 
     &- ll1mx*zeta3

      case(81)                  !1111

         ris = ll1mx**4/24d0
         
      end select
      
c --- set the imaginary part back to zero if it has been modified to
c --- get the branch cuts right (and should be zero). Also, set imaginary
c --- part of result to zero if x is between 0 and 1.

      if (bcflag.eq.1) then
         x = x - dcmplx(0d0,1d-60)
         xre = dreal(x)
         if (xre.ge.0d0.and.xre.le.1d0) then
            ris = dcmplx(dreal(ris),0d0)
         endif
      endif

      HPL4else=ris
      return
      end function
c-source-file HPL4.f
C=============================================================================
C---  HPLs  of Rank 4
C=============================================================================
c--- main forking function
      double complex function HPL4(n1, n2, n3, n4, x)
      implicit none
      integer n1,n2,n3,n4
      double complex x,ris
      double complex HPL4at0,HPL4at1,HPL4atm1
      double complex HPL4ar1,HPL4arm1,HPL4ar0,HPL4else
      double precision rad1,radm1,rad0

      rad1 = 0.01d0
      radm1 = 0.025d0
      rad0 = 0.025d0

      if ((abs(n1).gt.1).or.(abs(n2).gt.1).or.(abs(n3).gt.1)
     &     .or.(abs(n4).gt.1)) then
         print*, ""
         print*, "****************"
         print*, "Error in HPL4:"
         print*, "Indices",n1,n2,n3,n4," out of range !"
         print*, "Aborting..."
         print*,"****************"
         stop
      endif

      ris = dcmplx(0d0,0d0)
      
      if (x.eq.dcmplx(0d0,0d0)) then
c         print*, "I'm in 1"
         ris = HPL4at0(n1,n2,n3,n4)
      elseif (x.eq.dcmplx(1d0,0d0)) then
c         print*, "I'm in 2"
         ris = HPL4at1(n1,n2,n3,n4)
      elseif (x.eq.dcmplx(-1d0,0d0)) then
c         print*, "I'm in 3"
         ris = HPL4atm1(n1,n2,n3,n4)
      elseif (abs(x-dcmplx(1d0,0d0)).lt.rad1) then
c         print*, "I'm in 4"
         ris = HPL4ar1(n1,n2,n3,n4,x)
      elseif (abs(x+dcmplx(1d0,0d0)).lt.radm1) then
c         print*, "I'm in 5"
         ris = HPL4arm1(n1,n2,n3,n4,x)
      elseif (abs(x-dcmplx(0d0,0d0)).lt.rad0) then
c         print*, "I'm in 6"
         ris = HPL4ar0(n1,n2,n3,n4,x)
      else
c         print*, "I'm in 7"
         ris = HPL4else(n1,n2,n3,n4,x)
      endif
      HPL4=ris

c      write(*,'(I1,I1,I1,I1,D16.16,F16.16,F16.16)') 
c     &     n1,n2,n3,n4,x,abs(x),abs(x-dcmplx(1d0,0d0))

c      print*, n1,n2,n3,n4,x,abs(x),ris

      return
      end function
c ------------------------------------------------
      double complex function HPL4at0(n1, n2, n3, n4)
      implicit none
      integer n1,n2,n3,n4,j
      double complex ris
    
      j=1+(n4+1)*1+(n3+1)*3+(n2+1)*9+(n1+1)*27
      ris = dcmplx(0d0,0d0)

      if (j.eq.41) then             
         print*, ""
         print*, "****************"
         print*, "ERROR in HPL4: "
         print*, "HPL4(",n1,",",n2,",",n3,",",n4
     &        ,",0) is divergent!"
         print*, "Aborting..."
         print*,"****************"
         stop
      endif
      HPL4at0=ris
      return
      end function
c ------------------------------------------------
c --- Real part of HPL4     
      double precision function HPL4real(n1,n2,n3,n4,xr,xi)
      implicit none
      double precision xr,xi
      integer n1,n2,n3,n4
      double complex x,HPL4
      x=dcmplx(xr,xi)
      HPL4real = dreal(HPL4(n1,n2,n3,n4,x))
      return
      end

c --- Imaginary part of HPL4 
      double precision function HPL4im(n1,n2,n3,n4,xr,xi)
      implicit none
      double precision xr,xi
      integer n1,n2,n3,n4
      double complex x,HPL4
      x=dcmplx(xr,xi)
      HPL4im = dimag(HPL4(n1,n2,n3,n4,x))
      return
      end
