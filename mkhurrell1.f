c  cd77 -ezget -lats mkhurrell1.f -o mkhurrell1

c To compile (once LD_LIBRARY_PATH includes lats lib)
c [durack1@oceanonly 150219_AMIPForcingData]$ cd77 -ezget -lats mkhurrell1.f -o mkhurrell1 -fcray-pointer -std=legacy

c For debugging
c [durack1@oceanonly 150219_AMIPForcingData]$ valgrind --log-file=mkhurrell1.valout2 --track-origins=yes mkhurrell1

c 3456789012345678901234567890123456789012345678901234567890123456789012

c    21 August 1997
c    program description and comments added 8 February 2002
c    Code tweaks and recompile on oceanonly 10 June 2015 - PJD

c    Karl E. Taylor
c    PCMDI
c    taylor13@llnl.gov
c
c
c
c *********************************************************************
c *********************************************************************

c  GENERAL PROGRAM DESCRIPTION

c *********************************************************************

c  Purpose:

c    to create on some specified "target" grid an artificial mid-month 
c    sea ice fraction or SST data set (referred to here as the 
c    "boundary condition data" set) that, when linearly interpolated 
c    (in time), produces the observed monthly means (referred to here
c    as the "observed data").  

c  REFERENCE:  
c       
c    Taylor, K.E., D. Williamson, and F. Zwiers (2000): The sea surface
c          temperature and sea-ice concentration boundary conditions for
c          AMIP II simulations. PCMDI Report No. 60 and UCRL-MI-125597, 
c          Lawrence Livermore National Laboratory, Livermore, CA, 25 pp.
c          
c          pdf file available at 
c              http://www-pcmdi.llnl.gov/publications/ab60.html

c  Input files:

c      monthly mean (observed) SSTs and/or sea ice fraction
    
c      optional input: 
c               a file containing a field on the target grid which
c                       can be used to define the target grid
c                     (Alternatively, the target grid can be defined 
c                     within this program.)
c               a land/sea mask for the input data and/or
c                        a land/sea mask for the output grid.

c  Output files:

c      text file identified by suffix ".out" (e.g., gisstbc_sst_1x1.out)
c          The output filename is generated based on user-specified   
c             input parameters.
c          The file contains a list of the input parameters specified
c             by the user and global and zonal statistical information 
c             concerning the data processed.

c      optional output:
c         (1) The boundary condition data set for the years requested.
c               (The filenames contain the string "bc" and don't 
c                include the string "clim" (e.g. amipbc_sst_1978.nc).)
c         (2) A boundary condition data set for the year preceeding
c               the years requested. (The filenames contain the string
c               "spinup" (e.g., spinup_sic_360x180_1955.nc).)
c         (3) A climatological boundary condition data set based on the
c               years specified for computing this climatology.
c               (The filenames contain the strings "bc" and "clim"
c               (e.g., amipbc_sic_360x180_1979_2001_clim.nc).)
c         (4) Observed monthly mean data on the target grid.
c               (The filenames contain the string "obs" and don't
c                include the string "clim" (e.g. amipobs_sst_1978.nc).)
c         (5) Climatological observed monthly means based on the years
c               specified for computing this climatology. (The 
c               filenames contain the strings "obs" and "clim" (e.g.,
c               amipobs_sic_360x180_1979_2001_clim.nc).)
c         (6) Statistics useful in quality control of the data.
c               (The filenames contain the string "bcinfo", (e.g., 
c                bcinfo_sic_360x180.nc).)



c  Libraries used:
c    The following two libraries are required unless the input data 
c         (observations) are in pp format.  (if pp format input files 
c         are used, then the user can remove the call to subroutine 
c         getobs and also remove the subroutine itself.

c      ezget library: high level interface to read, mask, and regrid 
c                     input data. (required unless input, observational
c                     data files are pp format)
c      cdms library:  lower level interface to read data in netcdf, 
c                     drs, hdf or grib/grads format.  (required by
c                     ezget)
c      netcdf library: interface to netCDF files (required if input
c                      or output is netcdf format; if lats or cdms 
c                      libraries are linked, unsatisfied external 
c                      reference error may result unless drs library 
c                      is also linked.)
c      drs library:   interface to drs files (required if output is
c                     drs format; if cdms library is linked, 
c                     unsatisfied external reference error may result
c                     unless drs library is also linked.) 
c      hdf library interface to hdf files (if cdms library is linked, 
c                     unsatisfied external reference error may result
c                    unless hdf library is also linked.) 
c      lats library: lats contain the output subroutines for writing
c                    netCDF and grib files. (required unless output 
c                    is pp format or ascii)
c                 .

c  Overview of algorithm:

c     The following assumes familiarity with the information contained
c     in the reference listed above.

c        Away from regions of near-freezing temperature, the SSTs are
c     generated by solving a set of N linear equations where N
c     is the number of months considered.  The mid-month (boundary
c     condition) temperatures for a given month depend on the observed
c     monthly-mean temperature of that month, and also the temperatures
c     of the preceeding month and following month.  Thus, 
c     the mid-month temperature for the first and last months 
c     formally require observed temperatures before and after 
c     the time period considered.  In practice the dependence 
c     on the temperatures outside the period considered is fairly weak. 
c     To close the problem mathematically, however, we must impose some
c     sort of condition on the temperatures before and following the
c     period of interest.  

c        If observations are available for several months before and 
c     after the period of interest, then these can be used along with 
c     imposition of a periodic boundary condition to fully constrain the
c     problem.  (By "periodic boundary condition," we mean  that the 
c     first month of the available observations is assumed to follow    
c     the last month).  The periodic boundary condition should only be 
c     imposed on an  observational period that comprises an integral 
c     number of years, so, for example, if the first month is March, the
c     last month must be February.

c        The periodic condition is a convenient way of closing the
c     problem mathematically, but is of course unrealistic.  In  
c     practice, however, it really affects only the mid-month 
c     temperatures of a few months: the temperatures for the 2 or 3 
c     months at the beginning and end of the observational period 
c     considered.  If the observational record extends well outside 
c     (both before and after) the period for which boundary condition 
c     data are needed, then this is not a problem.  If, however, we need
c     boundary condition data for every month for which observations are
c     available, then the months at the beginning and end will be  
c     sensitive to this assumption.

c        If you really want to simulate the full time period for which 
c     observed monthly means are available, the "end" effects
c     can be reduced by generating artifical "observed" temperatures
c     for a brief period preceeding and following the actual 
c     observational period.  In this code we do this as follows: The 
c     artificial "observed" temperature anomalies are assumed to decay  
c     to zero in the first few months prior to and following the  
c     observational period, so that the temperature approaches  
c     climatology.  The rate of decay is based on the the   
c     autocorrelation function for the monthly mean SST time-series    
c     (analyzed for the years 1979-1998).  The spatial mean 
c     (area-weighted over all ocean grid cells) of the autocorrelation 
c     function, at lags of 1, 2, 3, ..., and 8 months have been
c     computed and are specified within this code.  The correlations 
c     for lags greater than 8 months are small, and are assigned so 
c     that the correlation reaches zero (in a smooth way) at month 12. 
c     [Note: for sea ice fraction the correlations are based on 
c     observations for lags less than 6 months, and the correlation is 
c     assumed to be zero for lags greater than 7 months.]

c        In summary, to minimize dependence on the boundary end 
c     conditions (i.e., the values of "observations" specfied for the 
c     period before and after the interval where observations are 
c     actually available), it is best to use only the mid-month values 
c     that are generated for a subinterval away from the boundary.  In  
c     practice this means the mid-month values for the 2 or 3 months at 
c     the beginning and the end of the observational period should be  
c     excluded from use in forcing a GCM simulation.  

c        Consider the following example.  Suppose observations
c     are available for the months January, 1956 through August, 2002.
c     Suppose, further, we would like to prescribe the SST boundary
c     condition in an atmospheric GCM simulation of these same months
c     (1/1956 - 8/2002).  The recommended approach would be to generate
c     artificial observational data for the year preceeding 1/56 and the
c     year following 8/02.  We also require that we treat an integral 
c     number of years.  We could conservatively choose to extend the 
c     "observations" to 12/03 at the end, and to also tack on a year
c     at the beginning (i.e., 12 months starting at 1/55).  
c     Since the artifical "observational" data approaches climatology
c     outside the time period of interest, the user must also 
c     specify which years will contribute to the climatology. If there
c     is no significant temperature trend over the period considered,
c     the climatology could be justifiably based on any interval of
c     several (say, 10 or more) years.  

c        If, however, the data exhibit a noticable trend, then the
c     user should probably avoid trying to simulate the entire 
c     period because the current code is unable to handle this case
c     accurately; it assumes that the artificial data outside both ends
c     of the observational period "decay" toward (i.e., approach) the 
c     same climatology.  It would be better to assume that at the 
c     beginning of the period and at the end of the period, the values 
c     approached the different characteristic climatologies.

c        In any case, it is wise to avoid the last few months for which
c     observations are currently available.  This is because as 
c     observations become available for succeeding months, these
c     will replace the artificially generated "observations and will
c     have some effect on boundary condition data for the few months 
c     near the end of the observational period.

c     The user can specify the beginning and ending months of the
c        following:
c            
c        1) the interval for which observational monthly mean data
c               are available.
c        2) the entire interval over which the analysis will be
c               performed, including buffer periods prior to and
c               following the observed data, which are recommended to 
c               total a year or two.
c        3) the interval over which the climatology is computed.
c        4) the interval defining which months will be included
c               as output (in the form of SST mid-month boundary
c               conditions).

c        Away from regions of near-freezing temperature, the mid-month
c     SSTs [represented here by S] satisfy:

c               A S = B

c      where S is a column vector of dimension N representing the 
c           mid-month values that constitute the SST boundary condition,
c           B is a column vector of dimension N, representing the 
c           observed monthly mean time-series, and A (except for 2 
c           elements) is a tri-diagnal matrix of dimension NxN.  The
c           two unconforming elements account for the assumption of 
c           periodicity. (The following shows the matrix structure.)
c

c          |b(1) c(1)             ...    a(1) | 
c          |a(2) b(2) c(2)                    |
c          | .   a(3) b(3) c(3)   ...         |
c      A = |                   .              |
c          |                   .              |
c          |                   .              |
c          |       ...    a(N-1) b(N-1) c(N-1)|
c          |c(N)   ...       .    a(N)   b(N) |
 
c      where 
c          a(i) = aa(i)/8
c          b(i) = 1 - aa(i)/8 - cc(i)/8
c          c(i) = cc(i)/8

c          aa(i) = 2*n(i)/(n(i)+n(i-1))
c          cc(i) = 2*n(i)/(n(i)+n(i+1))

c      where n(i) = number of days in month i.

c     Note that if all the months were of equal length, then aa=cc=1 
c          and a=1/8, b=3/4, and c=1/8.

c        For sea ice, the above procedure has to be modified because
c     the sea ice fraction is constrained to be between 0 and 1.  
c     Similarly, the water temperature cannot fall below its freezing
c     point, so this also places a physical constraint that is not 
c     always consistent with the above procedure.  In these cases 
c     the equations that must be solved have a similar structure as 
c     shown above, but the coefficients (a's, b's, and c's) depend on 
c     temperature or the sea ice fraction.  Thus, the equations in this
c     case are no longer linear, but they can be solved using an
c     iterative Newton-Raphson approach.

c     The Jacobian that is required under this approach is 
c     not generated analytically, but is approximated numerically.

c     There is another constraint necessary to ensure a unique
c     solution in the nonlinear case.  To see why the constraint
c     is needed, consider a grid cell where the ocean is ice-free
c     year around, except for 1 month, when the ice fraction is
c     10%.   The mid-month value for this month is not uniquely
c     determined because the cell could be covered by little
c     ice over the entire month, or by lots of ice for only a short
c     time in the middle of the month.  The algorithm relied on 
c     in this code tries to minimize the absolute difference between
c     mid-month values that exceed the maximum physically allowed
c     value (S_max=100% for sea ice) and S_max, while still yielding 
c     the correct monthly means: i.e.,
 
c                 where S > S_max, minimize (S - S_max) 
         
c     Similarly, for values that are less than the minimum physically
c     allowed value (S_min = 0% for sea ice):

c                 where S < S_min, minimize (S_min - S)

c     In the example above this results in the following mid-month 
c     values for sea ice (assuming for simplicity that all months
c     are of the same length):    S(i-1) = S(i+1) = -20%  and 
c     S(i) = 20%  where i is the month with mean sea ice fraction of 
c     10%.  These mid-month values when linearly interpolated give a 
c     sea ice fraction for month i that starts at 0%, linearly grows to
c     20% at the middle of the month, and then linearly decreases to 0%
c     at the end of the month.  For the month preceeding and the month
c     following month i, the linear interpolation leads to negative
c     values, but recall that the model's alogrithm will "clip" these
c     values, setting them to 0%.

c  Further details:
c               
c    SPECIFY the first and last month and year for period in which 
c        observed monthly mean data will be read. (The first month
c        read must not preceed mon1, iyr1, and the last month must
c        not follow monn, iyrn).  In the code set, for example:
           
c      mon1rd = 1            !AMIP
c      iyr1rd = 1956         !AMIP

c      monnrd = 8            !AMIP
c      iyrnrd = 2002         !AMIP

c    SPECIFY first and last month and year for entire period that will
c        be treated (i.e., the period of interest plus buffers at ends).
c        Note that the entire period treated should be an integral   
c        number of years.  In the code set, for example:

c      mon1 = 1            !AMIP
c      iyr1 = 1955         !AMIP

c      monn = 12           !AMIP
c      iyrn = 2003         !AMIP

c     SPECIFY the first and last month and year that will be included in
c         climatological mean.  (This must be an interval within the 
c         observed period).  In the code set, for example:

c      mon1clm = 1            !AMIP
c      iyr1clm = 1979         !AMIP

c      monnclm = 12           !AMIP
c      iyrnclm = 2000         !AMIP

c     SPECIFY the first and last month and year written to the output 
c         file.  (Try to begin at east a few months after the mon1rd,
c         iyr1rd, and end a few months before monnrd, iyrnrd, to avoid
c         sensitivity to the artificial data outside the observed 
c         period.)  In the code, set, for example:
 
c      mon1out = 1            !AMIP
c      iyr1out = 1956         !AMIP
c      
c      monnout = 6            !AMIP
c      iyrnout = 2002         !AMIP


c The structure of the input files should be as follows.  The data can
c reside in 1 or more input files covering the time period of interst.
c The last month of data found in each file (except chronologically the
c last input file) must be December.  Thus, for example you could have 
c the data organized roughly by decade as follows:

c     file 1:  march-december 1959
c     file 2: January 1960 - December 1969
c     file 3: January 1970 - December 1979
c     file 4: January 1980 - December 1989
c     file 5: January 1990 - December 1999
c     file 6: January 2000 - August 2001

c This is just one example; you could have all the data in a single 
c file.  Note that you must make a list of the input files as input to 
c my code, so if you're going to treat lots of years, the list might be 
c long if you store each year in a separate file.




c *********************************************************************

c THINGS TO DO ???:
c
c   generalize to treat 30-day calendars (in addition to realistic 
c      calendars
c
c     check whether program works if nlat/mlat .ne. an integer.
c   checked 7/14/03 that nlat/mlat need not be integer
c
c
c    don't open files if they already exist (change from unknown to new)
c          -- abort
c
c ???  *** GISST *** are units correct for sea ice??? (K or %?)
c
c   mismatch between reading and writing temp.pp files (.pp not there?)

c
c    SPECIFY FOLLOWING PARAMETERS:

c     nmon = number of months in entire period treated (including
c              buffers preceeding and following actual months of 
c              interest.
c     nlon = number of longitudes in output grid
c     nlat = number of latitudes in output grid
c     mlat = number of latitudes at one time (i.e., latitudes per
c              chunk) 
c     nzon = number of zones for which various summary statistics will
c              be calculated.  (Ideally, nlat/nzon = an integer)
c
      implicit none
      integer nmon, nlon, nlat, mlat, nzon, nzonp, n1, n2, nmon12, 
     &    nlagm, nchunks, niofiles

c *****************************************************************

      parameter (nmon=147*12, nlon=360, nlat=180, mlat=180, ! PJD Oceanonly 1870-2014 - '147*12' requires changing in subroutines
c      parameter (nmon=147*12, nlon=360, nlat=180, mlat=90, ! PJD Oceanonly 1870-2014
c      parameter (nmon=4*12, nlon=12, nlat=6, mlat=4, !test
c      parameter (nmon=27*12, nlon=288, nlat=181, mlat=46, !bala
c      parameter (nmon=51*12, nlon=360, nlat=180,  mlat=60, !AMIP            
c      parameter (nmon=142*12, nlon=360, nlat=180, mlat=45, !obs           
c      parameter (nmon=142*12, nlon=32, nlat=64,  mlat=64, !    
c      parameter (nmon=142*12, nlon=96, nlat=73,  mlat=73, !    
c      parameter (nmon=142*12, nlon=288, nlat=217,  mlat=31, !    
c      parameter (nmon=142*12, nlon=192, nlat=145,  mlat=29, !    
     &       nzon=6, nzonp=nzon+1, n1=nzonp*21, n2=nzonp*41,
     &       nmon12=12, nlagm=13, nchunks=(nlat-1)/mlat+1,
     &       niofiles=150)
c *****************************************************************

      logical caseindp
      integer i, j, icntl, maxiter, ntlat, ntlon
      integer lagcalc, iie, ie, ii
      integer m, mp, n, mm, nn, i1, kk, isrc, isrc1, nfilesi, nfileso
      integer ismc, jsmc, ismf, jsmf, ik, isea, je, k, ig, j1, jchnk, jn
      INTEGER  mons, m1, m2, mon1spin, iyr1spin, nmonspin

      real conv, r, wtmin, bbmin, dcbias, fac, tmin, offset,
     &    tmax, histo, varmin, dt, amissin, amissout, amissasc,
     &    rlat0, rlon0
      real obsclim(nlon,nlat,nmon12), sst(nlon,mlat,nmon), alons(nlon),
     &    alats(nlat), centmon(nlon,nlat,nmon12), gauss(nlat)
      real omax(nlon,nlat), omin(nlon,nlat), cmax(nlon,nlat),
     &    cmin(nlon,nlat), ovarmon(nlon,nlat), cvarmon(nlon,nlat),
     &    sstwts(nlon,nlat), wtfrac(nlon,nlat), 
     &    space(nlon,nlat,nmon12), spinup(nlon,nlat,nmon12)
      real woseas(nzon), wcseas(nzon), wovar(nzonp), wcvar(nzonp), 
     &   wocorr(nzonp), wccorr(nzonp), 
     &   acvarmon(nzonp), aovarmon(nzonp), aobsclim(nmon12,nzon),
     &   acentmon(nmon12,nzon), bpp(19)
      real ocorrel(nlagm,nzonp), ccorrel(nlagm,nzonp), correl(nmon12)
      real elem1(2), elemn(2), vecout(nmon),
     &     vecin(nmon), a(nmon), c(nmon), ac(nmon12), cc(nmon12)
      integer icount(nzonp,-20:20), jcount(nzonp,-10:10),
     &     kcount(nzonp,-20:20), mcount(nzonp,-10:10), lagm(nlagm)
      integer monlen(nmon12,2), lpp(45), jyr1out(niofiles), kyr1out(3),
     &    iyr1sst(niofiles), iyr1sic(niofiles), iyr1in(niofiles)
      character*80 title, grid
      character*120 fbcinfo, fbc, sftfileo, sftfilei, dum,
     &    inputfil(niofiles), outfile, pathout, tempbc(nchunks),
     &    tempobs(nchunks), fclimbc, fclimobs, fspinup, fobs, 
     &    inputsic(niofiles), inputsst(niofiles), outcopy,
     &    sftfilis, sftfilic  ! PJD Oceanonly 1870-2014
      character*16 vname, varin, model, sftnameo, sftnamei
      character varout*6, varout1*3        
      character*40 units, units2
      character*120 sourceo, sourceb
      character*30 suff
      character*3 suffix, c3
      character*4 ayr
      character*16 abbrev, src, src1  ! PJD Oceanonly 1870-2014
      character*9 cclim
      character*15 fnclim, blnk
c      double precision sum
      character*16 outftype, calclim, calgreg, gtype, typfil, inftype
      character*80 center
      character*256 parmtabl, dfltparm 
      integer idvar, ibasemon, ibaseyr, idvara(6), lbfc, iend, iregrid
      integer iout(6)
      integer lats, mon1rd, iyr1rd, monnrd, iyrnrd, mon1clm, iyr1clm,
     &     monnclm, iyrnclm, mon1, iyr1, monn, iyrn, mon1out, iyr1out, 
     &     monnout, iyrnout, mon1in, msklndi, msklndo
   
      data lagm/0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/
      data ig/1/
      data jyr1out/niofiles*-999/, iyr1sst/niofiles*-999/,
     &      iyr1sic/niofiles*-999/, kyr1out/3*-999/

      data lpp/45*0/, bpp/19*0.0/, blnk/''/

      data monlen/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31,
     &            31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
c        *** GISST ***   might assume months of equal length
c      data monlen/30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
c     &            30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30/
c *****************************************************************

c    SPECIFY controlling parameters:

c     icntl = 1 for SST, 2 for sic, and 3 for SST and sic


c     SPECIFY output wanted: 1=output 0=no output
c           iout(1) = boundary conditions for period specified by  
c           cd77 -ezget -lats mkgisst16.f -o mkgisst16      mon1out, iyr1out, monnout, iyrnout
c           iout(2) = boundary conditions for year preceeding period
c                       specified in 1
c           iout(3) = boundary conditions based on climatology
c           iout(4) = observed monthly means (normally =0 if input in 
c                       pp format, unless you want to change the units)
c           iout(5) = climatology of monthly means
c           iout(6) = diagnostic statistics (normally =0 for GISST)

      iout(1) = 1
      iout(2) = 1
      iout(3) = 1
      iout(4) = 1  
c      iout(4) = 0     !GISST
      iout(5) = 1
      iout(6) = 1

c  1 for sst 2 for ice 3 for both
          icntl = 3

c     SPECIFY source for boundary condition data (used to control 
c        description stored in output file cd77 -ezget -lats mkgisst16.f -o mkgisst16) 
c        (<16 characters)  
c        usually set to 'AMIP' or 'GISST' 
c        isrc is the length of the significant part of the character
c             string 

      src = 'hurrell'            !AMIP.
      isrc = 7                !AMIP
c      src = 'AMIP'            !AMIP.
c      isrc = 4                !AMIP
c      src = 'GISST'           !GISST
c      isrc = 5                !GISST

c       SPECIFY src1, which is an abbreviated lower case version of src 
c           used in generating output filenames (must not contain 
c           blanks) (< 16 characters)

      src1 = 'amip'           !AMIP
c      src1 = 'smip'           !SMIP
c      src1 = 'gisst'           !GISST
      isrc1 = index(src1, ' ') - 1

c    SPECIFY either 'obs' (if not regridding input data and EzGet input 
c       has been specified) or a brief indicator of the target grid.
c          (< 17 characters).

      model = 'obs' ! PJD Oceanonly 1870-2014
c      model = 'obs'       !obs
c      model = 'CCM3-T159'
c      model = 'pcmdi'     !AMIP  & test & bala 

c     SPECIFY the research center that has requested the data (written
c        as a comment in ascii and lats output files).  For grib and 
c        grads files, the center must appear in the parameter table 
c        lists (see parmtabl and dfltparm).
c        (< 81 characters)
c      center = 'UKMO'    !GISST
      center =  'pcmdi'   !AMIP

c    SPECIFY grid or model i.d. for use in generating file names 
c        (< 17 characters)
c         usually indicating grid.
c          examples: ncep, 360x180, 96x73, T42, T108,  ....
c      abbrev = '360x180'   !normal for GISST
c
c     SPECIFY grid type, which controls weights and for lats, controls
c              output grid type:
c              'gaussian' or 'cosine' or standard model acronym
c              (gtype is used only to control weighting if if both 
c               input and output are in pp format)
c     
c      gtype = 'cosine'   !GISST

c     SPECIFY complete grid description (needed for ascii output only)
c    for ascii files following is output to identify grid type:
c      grid = '1 x 1 degree uniformly-spaced longitude/latitude grid'
c      grid = 'T62 Gaussian grid'

c     SPECIFY file type for output? ('pp', 'drs', 'coards', 'cdf',   
c                  'grads', 'grib', 'ascii', 'all' (i.e., pp, drs, 
c                  coards, grib, and as.cii), or 'coards&grib' (i.e.,
c                    coards and grib).
c ???       why coards and not cdf?  coards expresses mid-month values
c                in hours since, whereas cdf expresses as year, month
c                 and fraction of day (e.g. 19790116.5)
c                 also cdf leads to an error in writing "fixed" fields
c                  (i.e., without a time dimension; e.g. bcinfo files)
c ???       why grib and not grads?
c ???       explain here differences between different lats formats


c    SPECIFY path to all output from this program
c      pathout = '/pcmdi/tmp1/doutriau/bcs/ecmwf'

c     SPECIFY file type for input?  ('pp' or 'ezget')
      inftype = 'ezget'       !AMIP
c      inftype = 'pp'           !GISST


      if (inftype(1:2) .ne. 'pp') then
c       land/sea mask for output (not used when reading pp format
c       note: land/sea mask is used for target grid even if msklndo=0 
c          0 = no mask
c          1 = 0.0 to 100.0 (% land)
c          3 = 0.0, 1.0, 2.0 OR 0, 1, 2 indicating land, ocean,
c                 and sea ice, respectively
c          4 = detailed geography mask

c *********************************************************************
c *********************************************************************
c               AMIP II Boundary Conditions Changes to do (start)
c *********************************************************************
c *********************************************************************

c      original grid (iregrid=0) or regrid (iregrid=1)?

c       iregrid=1    ! AMIP
       iregrid=0    ! obs ! PJD Oceanonly 1870-2014

          abbrev = '360x180' ! PJD Oceanonly 1870-2014
c          abbrev = 'T21'
c          abbrev = '64x32'
c          abbrev = '288x181'   !bala
c          abbrev = '12x6'     !test
c          abbrev = '360x180'  !obs 

      gtype = 'cosine' ! PJD Oceanonly 1870-2014
c           gtype  = 'gaussian'
c         gtype = 'cosine'
          
C    the following will appear only on ascii files:
      grid = '1 x 1 degree uniformly-spaced longitude/latitude grid' ! PJD Oceanonly 1870-2014
c      grid = '0.5 x 0.5 degree uniformly-spaced longitude/latitude grid'
c       grid = 'Gaussian grid'
c      grid = 'uniformly spaced'

c     outftype options include 'drs', 'coards', 'grib', 'ascii', 
c                     'coards&grib' and 'all'
c          outftype = 'grib'        test
c          outftype = 'coards&grib'
c          outftype = 'coards&asc'
c          outftype = 'notdrs'       !obs
c          outftype = 'all'        !AMIP & bala
c          outftype = 'ascii'
          outftype = 'coards'

          pathout = '/work/durack1/Shared/150219_AMIPForcingData/'
     & // '360x180'
c     !PJD Oceanonly 1870-2014
c          pathout = '/pcmdi/tobala/288x181/ORIG/'    !bala
c          pathout = '/pcmdi/zooks1/SSTCICE/360x180/'     !obs
c          pathout = '/pcmdi/zooks1/SSTCICE/96x73/'     !AMIP
c          pathout = '/pcmdi/zooks1/SSTCICE/T21/'     !AMIP
c          pathout = '/pcmdi/zooks1/SSTCICE/12x6/'    !test
c          pathout = '/pcmdi/zooks1/SSTCICE/T255_NtoS/'   !AMIP
c          pathout = '/pcmdi/zooks1/SSTCICE/SMIP_T42/'   !AMIP


          msklndo=0
          sftfileo = 'none'
c          sftfileo = '/pcmdi/staff/longterm/doutriau/ldseamsk/amipII/' 
c     &            // 'pcmdi_sftlf_ccm_T159.nc'
c          sftfileo = '/pcmdi/roseland0/amip/bcs/codes/' 
c     &            // 'Uniform-grid-451x720.nc'
c     &            // 'Gaussian-grid-T340.nc'
          sftnameo = 'sftlf'
c
c
c   define output grid (ignored if sftfileo .ne. 'none' or if
c                  iregrid=0)

      if ((sftfileo .eq. 'none') .or. (iregrid .eq. 1)) then 
c         ntlat = number of latitudes on target grid
c         rlat0 = first latitude location 
c                 (for gaussian grids, the sign of this scalar is only
c                  meaningful; it determines whether the grid begins 
c                  near the south pole (negative sign) or north pole 
c                  (+ sign))
c    
c         ntlon = number of longitude grid cells
c         rlon0 = first longitude location
c      
        ntlat = nlat
        rlat0 = -89.5  !AMIP obs
c        rlat0 = 75.    !test
c        rlat0 = 90.0   !bala
c        rlat0 = -90.    !Hadley
        ntlon = nlon
        rlon0 = 0.5     !obs ! PJD Oceanonly 1870-2014
c         rlon0 = 0.5    !AMIP regular
c         rlon0 = 0.    ! gaussian ! PJD Oceanonly 1870-2014
c        rlon0 = -177.5  ! GISS
c         rlon0 = 0.0     ! Hadley

        if (msklndo .ne. 0) then 
          print*, ' Error in mkgisst -- if sftfileo = "none" then '
          print*, '      msklndo must be set to 0'
          stop
        endif

      endif

      



c ***************************************************************************
c ***************************************************************************
c               AMIP II Boundary Conditions Changes to do (done)
c ***************************************************************************
c ***************************************************************************
c       land/sea mask for input (see above for key) 

         msklndi = 0
c       msklndi = 1

        sftnamei = 'sftlf'

        sftfilis = 
     &     '/work/durack1/Shared/150219_AMIPForcingData/SSTCICE/OBS/'
     & // 'Hurrell_Shea/sftlf_360x180.nc' ! PJD Oceanonly 1870-2014
c     &     '/pcmdi/AMIP2/data_fixed/ldseamsk/obs/sftlf_180x360.nc'
        sftfilic = 
     &     '/work/durack1/Shared/150219_AMIPForcingData/SSTCICE/OBS/'
     & // 'Hurrell_Shea/sftlf_360x180.nc' ! PJD Oceanonly 1870-2014
c     &     '/pcmdi/AMIP2/data_fixed/ldseamsk/obs/sftlf_180x360.nc'

c    &  '/pcmdi/staff/longterm/doutriau/ldseamsk/bcamip2/mask.1deg.ctl'
c        sftfilis = 
c     &        '/pcmdi/staff/doutriau/ldseamsk/bcamip2/'//
c     &        'sftlf_hadisst10_sst.nc'  
c        sftfilic = 
c     &        '/pcmdi/staff/doutriau/ldseamsk/bcamip2/'//
c     &        'sftlf_hadisst10_sic.nc'  
c
c        sftnamei = 'sftl'
c        sftfilis = 
c     &      '/pcmdi/staff/longterm/gleckler/gleckler/ktaylor/sftl.nc'
c        sftfilic = 
c     &      '/pcmdi/staff/longterm/gleckler/gleckler/ktaylor/sftl.nc'
c        sftfilei = '/pcmdi/share1/obs/amip2/bcs/fnl/amip2.bcs.mask.nc'
c        sftfilei = '/pcmdi/doutriau/geog/sealand/data/sftl_1x1.nc'

      else

        sftnameo = ' '
        msklndo = 0
        sftfileo = ' '
        sftnamei = ' '
        msklndi = 0
        sftfilei = ' '

      endif


c    SPECIFY first month and year for entire period that will be treated 
c       (period of interest plus buffers at ends)
      mon1 = 1            !hurrell
      iyr1 = 1869         !hurrell
c      mon1 = 1            !AMIP
c      iyr1 = 1955         !AMIP
c      mon1 = 1            !GISST
c      iyr1 = 1869         !GISST
c      mon1 = 1            !bala
c      iyr1 = 1976         !bala
c      mon1 = 1            !test
c      iyr1 = 19August78         !test

c    last month and year for entire period that will be treated (period
c        of interest plus buffers at ends).  Note the entire period 
c        treated should be an integral number of years.
      monn = 12           !hurrell
      iyrn = 2015 ! PJD Oceanonly 1870-2014
c      iyrn = 2010         !hurrell
c      monn = 12           !AMIP & bala
c      iyrn = 2007         !AMIP & bala
c      monn = 12           !GISST
c      iyrn = 1997         !GISST
c      monn = 12           !test2
c      iyrn = 1873         !test2
c      monn = 12           !test
c      iyrn = 1981         !test

c    first month and year for period in which observed monthly mean data
c           will be read (must not preceed mon1, iyr1)            
      mon1rd = 1            !AMIP
      iyr1rd = 1870         !AMIP
c      mon1rd = 1            !GISST
c      iyr1rd = 1870         !GISST
c      mon1rd = 1            !bala
c      iyr1rd = 1977         !bala
c      mon1rd = 1            !test
c      iyr1rd = 1979         !test

c    last month and year for period in which observed monthly mean data
c           will be read (must not follow monn, iyrn) 
      monnrd = 5 ! PJD Oceanonly 1870-2014
      iyrnrd = 2015 ! PJD Oceanonly 1870-2014
c      monnrd = 6            !AMIP & bala
c      iyrnrd = 2010          !AMIP & bala
c      iyrnrd = 1981         !test3
c      monnrd = 12           !GISST
c      iyrnrd = 1996         !GISST
c      monnrd = 12           !test2
c      iyrnrd = 1872         !test2
c      monnrd = 2            !test
c      iyrnrd = 1981         !test

c     first month and year that will be included in climatological mean
      mon1clm = 1            !1/18/07
      iyr1clm = 1988         !1/18/07
c      mon1clm = 1            !SMIP
c      iyr1clm = 1979         !SMIP 
c      mon1clm = 1            !AMIP & bala
c      iyr1clm = 1979         !AMIP & bala
c      iyr1clm = 1971         !SMIP
c      mon1clm = 1            !GISST
c      iyr1clm = 1870         !GISST
c      mon1clm = 1            !test2
c      iyr1clm = 1870         !test2
c      mon1clm = 1            !test
c      iyr1clm = 1979         !test

c     last month and year that will be included in climatological mean
      monnclm = 12            !1/18/07
      iyrnclm = 2007          !1/18/07
c      monnclm = 12            !SMIP
c      iyrnclm = 2001          !SMIP
c      monnclm = 12            !AMIP & bala
c      iyrnclm = 2001          !AMIP & bala
c      iyrnclm = 2000          !SMIP
c      iyrnclm = 1981         !test3
c      monnclm = 12           !GISST
c      iyrnclm = 1996         !GISST
c      monnclm = 12           !test2
c      iyrnclm = 1872         !test2
c      monnclm = 12           !test
c      iyrnclm = 1980         !test

c     first month and year written to output file 
      mon1out = 1            ! hurrell
      iyr1out = 1870         ! hurrell
c      mon1out = 1            !AMIP
c      iyr1out = 1956         !AMIP
c      mon1out = 1            !GISST
c      iyr1out = 1870         !GISST
c      mon1out = 1            !bala
c      iyr1out = 1978         !bala
c      mon1out = 1            !test
c      iyr1out = 1979         !test

c     last month and year written to output file
      monnout = 3 ! PJD Oceanonly 1870-2014
      iyrnout = 2015 ! PJD Oceanonly 1870-2014
c      monnout = 3            ! 1/18/07
c      iyrnout = 2010          ! 1/18/07
c      monnout = 6            !AMIP & test
c      iyrnout = 2005         !AMIP & test
c      monnout = 12           !AMIP
c      iyrnout = 2001         !AMIP
c      iyrnout = 1981         !test3
c      monnout = 12           !GISST
c      iyrnout = 1996         !GISST
c      monnout = 1            !test2
c      iyrnout = 1872         !test2
c      monnout = 1            !test
c      iyrnout = 1980         !test

c     first month found in first input file (these values are 
c          ignored if the 'pp' input format has been selected)
c          data is assumed to be stored by month and if non-pp
c          format, data is extracted by time-index (not actual
c          time-coordinate value).
          
      mon1in = 1           !AMIP
c      mon1in = 0            !GISST



c     SPECIFY path to SST concentration data (not used if processing 
c              only sic data) Also if more than 1 SST input file, then
c              specify year of first month of data in each file.
c              (Note, last month in each file, except possibly last
c               file, should be December) 
        inputsst(1) = 
     &     '/work/durack1/Shared/150219_AMIPForcingData/'
     & //  'SST_NEW3/'
     & //  'MODEL.SST.HAD187001-198110.OI198111-201505.nc' ! PJD Oceanonly 1870-2014

c     &     '/pcmdi/zooks1/SSTCICE/OBS/Hurrell_Shea/' //
c     &     'MODEL.SST.HAD187001-198110.OI198111-201006.nc'
c     &     '/pcmdi/zooks1/SSTCICE/sst.ctl'
c     &     '/pcmdi/roseland0/amip2/bcs/data/latest/sst/sst.ctl'
c     &     '/pcmdi/zooks1/SSTCICE/sst.ctl'
c     &     '/zooks1/SSTCICE/sst_new.ctl'
c     &     '/pcmdi/doutriau/amip2/bcs/data/sst/amip2.bcs.sst.ctl'
c     &    '/pcmdi/staff/longterm/fiorino/sst/ac_corr/amip2.bcs.sst.ctl'
        iyr1sst(1) = 1870
c        iyr1sst(1) = 1956
c        inputsst(1) = 
c     &    '/pcmdi/staff/longterm/ktaylor/GISST/Gisst2.3a.pp_1870_1900'
c        iyr1sst(1) = 1870
c        inputsst(2) = 
c     &    '/pcmdi/staff/longterm/ktaylor/GISST/Gisst2.3a.pp_1901_1930'
c        iyr1sst(2) = 1901
c        inputsst(3) = 
c     &    '/pcmdi/staff/longterm/ktaylor/GISST/Gisst2.3a.pp_1931_1960'
c        iyr1sst(3) = 1931
c        inputsst(4) = 
c     &    '/pcmdi/staff/longterm/ktaylor/GISST/Gisst2.3a.pp_1961_1996'
c        iyr1sst(4) = 1961

c     SPECIFY path to sea ice concentration data (not used if processing 
c              only sst data) Also if more than 1 sic input file, then
c              specify year of first month of data in each file.
c              (Note, last month in each file, except possibly last
c               file, should be December)
        inputsic(1) = 
     &     '/work/durack1/Shared/150219_AMIPForcingData/'
     & //  'SST_NEW3/'
     & //  'MODEL.ICE.HAD187001-198110.OI198111-201505.nc' ! PJD Oceanonly 1870-2014

c     &     '/pcmdi/zooks1/SSTCICE/OBS/Hurrell_Shea/' //
c     &     'MODEL.ICE.HAD187001-198110.OI198111-201006.nc'
c     &     '/pcmdi/zooks1/SSTCICE/sic.ctl'
c     &     '/pcmdi/roseland0/amip2/bcs/data/latest/sic/sic.ctl'
c     &     '/pcmdi/zooks1/SSTCICE/sic.ctl'
c     &     '/zooks1/SSTCICE/sic_new.ctl'
c     &     '/pcmdi/doutriau/amip2/bcs/data/sic/amip2.bcs.sic.ctl'
c     &      '/pcmdi/staff/longterm/fiorino/sic/orig/amip2.bcs.sic.ctl'
        iyr1sic(1) = 1870
c        iyr1sic(1) = 1956
c        inputsic(1) = '/pcmdi/staff/longterm/ktaylor/GISST/' //
c     &               'Gisst2.3a_ice.pp_1870_1900'
c        iyr1sic(1) = 1870
c        inputsic(2) = '/pcmdi/staff/longterm/ktaylor/GISST/' //
c     &                'Gisst2.3a_ice.pp_1901_1930'
c        iyr1sic(2) = 1901
c        inputsic(3) = '/pcmdi/staff/longterm/ktaylor/GISST/' //
c     &                'Gisst2.3a_ice.pp_1931_1960'
c        iyr1sic(3) = 1931
c        inputsic(4) = '/pcmdi/staff/longterm/ktaylor/GISST/' //
c     &                'Gisst2.3a_ice.pp_1961_1996'
c        iyr1sic(4) = 1961


c    SPECIFY year of first month of boundary condition and monthly mean
c         data stored in each output file (i.e., fbc and fobs files)
c        It is not necessary to define jyr1out(1), because this must be
c           set equal to iyr1out (defined above).
c        If you want all data placed in a single file, it is not
c            necessary to define jyr1out
c        If you want each year to be written to separate files, then set
c           (jyr1out(i), i=2,number of years) equal to consecutive years
c           (e.g., 1975, 1976, 1977)
c        If you want each decade of data to be written to separate files
c           then set (jyr1out(i), i=2,number of decades) equal to years
c           separated by increments of decades (e.g., 1942, 1952, 1962)
c           note: jyr1out(2) - jyr1out(1) can differ from the other 
c              increments
c     jyr1out(1) = iyr1out      
c     do 26 i=2,44         !AMIP
c       jyr1out(i) = jyr1out(1) + i - 1
c  26 continue
c      do 26 i=2,13        !GISST
c        jyr1out(i) = 1870 + 10*(i - 1)
c   26 continue

      jyr1out(1) = iyr1out
      ii = iyrnout-iyr1out+1
      if (ii .gt. niofiles) then
         print*, 'Error found.  You must dimension niofiles .ge. ', ii
         stop
      endif

      if (ii .gt. 1) then
        do 26 i=2,ii         !AMIP
          jyr1out(i) = jyr1out(1) + i - 1
   26   continue
       endif            

      write(fnclim, "(i4, '_', i4, '_clim')") iyr1clm, iyrnclm
      write(cclim, "(i4, '-', i4)") iyr1clm, iyrnclm

c *****************************************************************


c      If ocean contribution to target grid cell is less than wtmin, 
c          then values from land will be included in calculating
c          value on target grid.

       wtmin = 0.5

c    parameters controlling nonlinear solver:
c       maximum number of iterations?

       maxiter = 200

c    smallest diagnal element allowed in Jacobian?

       bbmin = 0.001

c   calculate lag correlations?
c       0 = no
c       1 = calculate monthly lag correlations

      lagcalc = 0
    
c     base-time?  (used only for lats output)

      ibasemon = 1
      ibaseyr = 1979         !AMIP
c      ibaseyr = 1870         !GISST 

      if (src(1:isrc) .eq. 'GISST') then
c       ***GISST*** check source description of GISST data:
        sourceo = 'GISST 2.3a data: Parker et al., Hadley Centre '//
     &         'Clim. Res. Tech. Note 63, 1995 & Rayner et al.,' //
     &         ' HCCRTN 74, 1996'
        sourceb = 'Based on GISST 2.3a data: Parker et al., Hadley'//
     &         ' Centre Clim. Res. Tech. Note 63, 1995 & Rayner et '//
     &         'al., HCCRTN 74, 1996'

      elseif (src(1:isrc) .eq. 'AMIP') then
        sourceo = 'Observations on 1x1 deg grid (Fiorino, 1996): ' //
     &                  'http://www-pcmdi.llnl.gov/amip2/' //
     &                  'AMIP2EXPDSN/BCS_OBS/amip2_bcs.htm'
        sourceb = 'See: http://www-pcmdi.llnl.gov/amip2/' //
     &                  'AMIP2EXPDSN/BCS/amip2bcs.html'
      elseif (src(1:isrc) .eq. 'hurrell') then
         sourceo = 'Hurrell SST/sea ice consistency criteria applied '//
     &         'to merged HadISST (1870-01 1981-10) & NCEP-0I2 ' //
     &         '(1981-11 to now) data'
         sourceb = 'Based on Hurrell SST/sea ice consistency criteria'//
     &         ' applied to merged HadISST (1870-01 1981-10) & ' //
     &         'NCEP-0I2 (1981-11 to ?)'


      else
        print*, "Error: src should be either 'GISST' or 'AMIP' or ' //
     &      'hurrell')"
        call exit(1)
      endif

c     SPECIFY missing data indicators: 
c        amissin  :  missing data value in input files
c        amissout :  value assigned to missing data in output
c                      files (except ascii files)
c        amissasc :  value assigned to missing data in ascii
c                      output files

      amissin = 1.0e20       !AMIP
c      amissin = -1.0e30       !GISST

      amissout = 1.0e20      !AMIP
c      amissout = -1.0e30      !GISST

      amissasc = 0.0

c      calendar?  ('climatology', 'julian', 'gregorian' or 'no_leap')
c          note: 'gregorian' is only option for 'coards' file type
c          (parameter ignored if pp, drs, or ascii output selected)
c ???    generalize for 30-day months

      calgreg = 'gregorian'

        if ((outftype(1:6) .eq. 'coards') .or. 
     &        (outftype(1:3) .eq. 'all') .or.
     &        (outftype(1:11) .eq. 'coards&grib') .or.
     &        (outftype(1:10) .eq. 'coards&asc') .or.
     &        (outftype(1:6) .eq. 'notdrs')) then
          calclim = 'gregorian'
        else
          calclim = 'climatology'
        endif

        if (caseindp(calgreg , 'no_leap')) monlen(2,2) = 28
 
   22 if ((icntl .eq. 2) .or. (icntl .eq. 3)) then
        varout = 'sicbcs'
        varout1 = 'sic'
      else 
        varout = 'tosbcs'
        varout1 = 'tos'
      endif

      if (varout(1:3) .eq. 'sic') then

        sftfilei = sftfilic

c       convergence criteria
c         this is used to determine how closely the monthly means
c         determined by the boundary condition data must match
c         the observed before convergence is declared.
c         The units are the same as the output data.
c         the value should not be set too small, because
c         this is also used in calculating derivative with 
c         respect to time used in generating the jacobi matrix
c         required for Newton-Raphson type interative method
c         for solving nonlinear equations.  
c         In short, set it smaller than the measurement precision 
c         of the input data, but don't approach the precision of
c         a single precision floating point number (which is ~1e-6
c         times a typical value of the data).  We recommend 
c         .01 if sea ice concentration is expressed in percent;
c         .0001 if in fraction.
    
        conv = 1.e-2    !AMIP
c        conv = 1.e-4    !bala
c ???        conv = 0.5e-2

c       constant to subtract to reduce truncation errors
        dcbias = 0.0

c       factor to multiply input data by (before any cropping).
c          Otput will be scaled by this factor:
         fac = 1.0          ! Hurrell
c        fac = 100.0        !AMIP
c        fac = 1.0          !bala
c  ???        *** GISST *** check
c        fac = 10.0            !GISST
c        number added to data (after scaling by fac)
        offset = 0.0

c        output units for sea ice (also variance units)
c ???    check following for GISST
        units = '%'            !AMIP
        units2 = '%**2'        !AMIP
c        units = 'fraction'     !bala
c        units2 = ''            !bala

c       prescribed limits (in output units) used (if necessary) to crop 
c         data and accounted for in creating boundary condition data
        tmin = 0.0            !AMIP
        tmax = 100.0          !AMIP
c        tmin = 0.0            !bala
c        tmax = 1.0            !bala
c          *** GISST ***  check whether tenths are stored as character 
c          data?
c         tmin = 0.0            !GISST
c         tmax = 100.0           !GISST

c       smoothing when values go from mininum to maximum
c          maximum jump in monthly means allowed is
c          tmax-tmin-dt  (in output units)
c        dt = 2.0
c        dt = (tmax-tmin)/50.
        dt = (tmax-tmin)/100.   !AMIP 1/14/2003
c
        if (src(1:isrc) .eq. 'hurrell') then
           varin = 'SEAICE'
        else
           varin = 'sic'
        endif

        suffix = 'sic'
        suff = 'Sea Ice Concentration (%)'       !AMIP
c        suff = 'Sea Ice Fraction'                !bala
c ???      suff = 'Sea Ice Concentration (%)'         !GISST
c          length of above string
        iie = 25

c       histogram fine increment?
        histo = (tmax-tmin)/100.0   !AMIP
c ???    histo = (tmax-tmin)/10.0    !GISST

c       minimum variance for grid cell to be included in area means
c         (expressed in output units)
        varmin = 1.0

c       field code for pp-format output
        lbfc = 37

c       lag correlations for monthly mean anomalies to be used in
c          generating artificial data

        if (iregrid .eq. 0) then
c          following are from AMIP II data on original 1x1 degree grid
           correl(1) = 0.53
           correl(2) = 0.23
           correl(3) = 0.13
           correl(4) = 0.08
           correl(5) = 0.05
c          following are NOT taken from observations
           correl(6) = 0.02
           correl(7) = 0.00
           correl(8) = 0.00
           correl(9) = 0.00
           correl(10) = 0.00
           correl(11) = 0.00
           correl(12) = 0.00
c          observed values:
c           correl(6) = 0.05
c           correl(7) = 0.05
c           correl(8) = 0.06
c           correl(9) = 0.06
c           correl(10) = 0.08
c           correl(11) = 0.09
c           correl(12) = 0.09
        else
c           following are from AMIP II observations regridded to 3x3 
c                 degree resolution
           correl(1) = 0.58
           correl(2) = 0.27
           correl(3) = 0.15
           correl(4) = 0.09
           correl(5) = 0.06
c          following are NOT taken from observations
           correl(6) = 0.03
           correl(7) = 0.00
           correl(8) = 0.00
           correl(9) = 0.00
           correl(10) = 0.00
           correl(11) = 0.00
           correl(12) = 0.00
c          observed values:
c           correl(6) = 0.05
c           correl(7) = 0.05
c           correl(8) = 0.06
c           correl(9) = 0.07
c           correl(10) = 0.08
c           correl(11) = 0.09
c           correl(12) = 0.08
         endif

         do 31 i=1,niofiles
           inputfil(i) = inputsic(i)
           iyr1in(i) = iyr1sic(i)
   31    continue

c         if output is pp format, but input is ezget, 
c           then generate pp header
          if ((outftype(1:2) .eq. 'pp') .and. 
     &         (inftype(1:2) .ne. 'pp')) then

            print*, 'at present input file type must be in'
            print*, 'in pp format if output is in pp format'

            stop
c           call genpphd
          endif



      else

        sftfilei = sftfilis

c       convergence criteria
c         this is used to determine how closely the monthly means
c         determined by the boundary condition data must match
c         the observed before convergence is declared.
c         The units are the same as the output data.
c         the value should not be set too small, because
c         this is also used in calculating derivative with 
c         respect to time used in generating the jacobi matrix
c         required for Newton-Raphson type interative method
c         for solving nonlinear equations.  
c         In short, set it smaller than the measurement precision 
c         of the input data, but don't approach the precision of
c         a single precision floating point number (which is ~1e-6
c         times a typical value of the data).  We recommend 
c         .001 for SST data in Cesius or kelvin.
    
        conv = 1.e-3

c       freezing point in output units 
c        ***GISST*** check whether data read are in C or K.
        dcbias = 273.16      !AMIP
c        dcbias = 0.0         !bala
c        dcbias = 0.0         !GISST

c       factor to multiply input data by (before any cropping).
c          Output will be scaled by this factor:
        fac = 1.0          !AMIP & GISST
c        number added to data (after scaling by fac)
c        to convert from Kelvin to Celsius, set offset = -273.16
        offset = 273.16    ! hurrell to convert input from C to K
c        offset = 0.0       !AMIP
c        offset = -273.16   !bala
c        offset = 0.0       !GISST
c ???    offset = +-271.35   !GISST

c        output units for sea ice (also variance units)
c        ***GISST**** check units
c        units for sst (also variance units)
        units = 'K'          !AMIP
        units2 = 'K**2'      !AMIP
c        units = 'C'          !GISST & bala
c        units2 = 'C**2'      !GISST & bala


c       prescribed limits (in output units) used (if necessary) to crop 
c         data and accounted for in creating boundary condition data
c   
          tmin = 271.36           !ERA40 & AMIP
          tmax = 1000.0           !AMIP
c          tmin = -1.8             !bala
c          tmax = 1000.0           !bala
c          tmin = -1.8             !GISST  this is appropriate for C
c          tmax =  1000.0          !GISST

c       smoothing when values go from mininum to maximum
c          maximum jump in monthly means allowed is
c          tmax-tmin-dt  (in output units)
        dt = 0.001

        if (src(1:isrc) .eq. 'hurrell') then
           varin = 'SST'
        else
           varin = 'sst'
        endif

        suffix = 'sst'
        suff = 'SST (K)'    !AMIP
c        suff = 'SST (C)'    !bala
c         suff = 'SST (C)'    !GISST  check units

        iie = 7

c       histogram fine increment?
        histo = 0.1

c       minimum variance for grid cell to be included in area means
c         (expressed in output units)
        varmin = 0.0

c       field code for pp-format output
        lbfc = 16

c       lag correlations for monthly mean anomalies to be used in
c          generating artificial data

        if (iregrid .eq. 0) then
c          following are from AMIP II data on original 1x1 degree grid
           correl(1) = 0.68
           correl(2) = 0.46
           correl(3) = 0.33
           correl(4) = 0.26
           correl(5) = 0.20
           correl(6) = 0.17
           correl(7) = 0.14
           correl(8) = 0.11
c          following are NOT taken from observations
           correl(9) = 0.08
           correl(10) = 0.05
           correl(11) = 0.02
           correl(12) = 0.00
c          observed values:
c           correl(9) = 0.09
c           correl(10) = 0.07
c           correl(11) = 0.06
c           correl(12) = 0.03
            
        else
c           following are from AMIP II observations regridded to 3x3 
c                 degree resolution
           correl(1) = 0.69
           correl(2) = 0.47
           correl(3) = 0.34
           correl(4) = 0.27
           correl(5) = 0.21
           correl(6) = 0.17
           correl(7) = 0.14
           correl(8) = 0.12
c          following are NOT taken from observations
           correl(9) = 0.09
           correl(10) = 0.06
           correl(11) = 0.03
           correl(12) = 0.00
c          observed values:
c           correl(9) = 0.10
c           correl(10) = 0.08
c           correl(11) = 0.06
c           correl(12) = 0.03
        endif

         do 32 i=1,niofiles
           inputfil(i) = inputsst(i)
           iyr1in(i) = iyr1sst(i)
   32    continue

c         if output is pp format, but input is ezget, 
c           then generate pp header
          if ((outftype(1:2) .eq. 'pp') .and. 
     &         (inftype(1:2) .ne. 'pp')) then

            print*, 'at present input file type must be in'
            print*, 'in pp format if output is in pp format'

            stop
c           call genpphd
          endif

      endif

      tmin = tmin - dcbias
      tmax = tmax - dcbias

      do 33 i=1,niofiles
        if (iyr1in(i) .ne. -999) nfilesi = i
        if (jyr1out(i) .ne. -999) nfileso = i
   33 continue          

      je = index(pathout, ' ') - 1

      ie = index(abbrev, ' ')
      if (ie .eq. 0) then
        ie = 16
      else
        ie = ie-1
      endif

      outfile = pathout(1:je)//
     &     '/'//src1(1:isrc1)//'bc_'//suffix//'_'//abbrev(1:ie)//'.out'
c   example:   outfile = gisstbc_sst_1x1.out

c ***      fclimbc =  pathout(1:je)//'/climbc_'//suffix//'_'//abbrev(1:ie)
      fclimbc =  pathout(1:je)//
     &         '/'//src1(1:isrc1)//'bc_'//suffix//'_'//abbrev(1:ie)
      fbcinfo =  pathout(1:je)//'/bcinfo_'//suffix//'_'//abbrev(1:ie)
      fbc = pathout(1:je)//
     &       '/'//src1(1:isrc1)//'bc_'//suffix//'_'//abbrev(1:ie)
c ***     if (src1(1:4) .eq. 'amip') then
c         fspinup =  pathout(1:je)//'/bc1978_'//suffix//'_'//abbrev(1:ie)
c      else
         fspinup =  pathout(1:je)//'/spinup_'//suffix//'_'//abbrev(1:ie)
c      endif
      fobs = pathout(1:je)//
     &          '/'//src1(1:isrc1)//'obs_'//suffix//'_'//abbrev(1:ie)
c ***      fclimobs = pathout(1:je)//'/climobs_'//suffix//'_'//abbrev(1:ie)
      fclimobs = pathout(1:je)//
     &    '/'//src1(1:isrc1)//'obs_'//suffix//'_'//abbrev(1:ie)

c      the following only used for lats output
c      parmtabl = '/home/taylor13/pcmdi/util/ketgrib.parms'
c      dfltparm = '/usr/local/lats/table/amip2.lats.table' ! sunOS
c      dfltparm = '/usr/local/lib/lats/amip2.lats.table' ! stargate
c      dfltparm = '/usr/local/lib/lats/amip2.parms' ! linux
c      dfltparm = '/work/durack1/Shared/150219_AMIPForcingData/src/lats/'
c     & // 'amip2.parms' ! PJD Oceanonly 1870-2014
      parmtabl = '/work/durack1/Shared/150219_AMIPForcingData/'
     & // 'ketgrib.parms' ! PJD Oceanonly 1870-2014
      dfltparm = '/work/durack1/Shared/150219_AMIPForcingData/'
     & // 'ketgrib.parms' ! PJD Oceanonly 1870-2014

      open(9, file=outfile, status='new')

      print*, ' '
      print*, 'Output from monthly boundary condition generator '
      print*, 'varin = ', varin
      print*, 'varout = ', varout
      print*, 'varout1 = ', varout1
      print*, 'src = ', src
      print*, 'src1 = ', src1
      print*, 'model = ', model
      print*, 'grid = ', grid
      print*, 'gtype = ', gtype
      print*, 'inftype = ', inftype
      print*, 'outftype = ', outftype
      print*, 'calclim = ', calclim
      print*, 'center = ', center
      print*, 'abbrev = ', abbrev
      print*, 'suffix = ', suffix
      print*, 'suff = ', suff

      print*, ' '
      print*, 'icntl = ', icntl
      print*, 'lagcalc = ', lagcalc
      print*, 'iout = ', iout

      print*, ' '
      print*, 'nlon = ', nlon
      print*, 'nlat = ', nlat
      print*, 'mlat = ', mlat
      print*, 'nchunks = ', nchunks

      print*, ' '
      print*, 'nmon = ', nmon
      print*, 'mon1rd = ', mon1rd
      print*, 'iyr1rd = ', iyr1rd
      print*, 'monnrd = ', monnrd
      print*, 'iyrnrd = ', iyrnrd

      print*, 'mon1clm = ', mon1clm
      print*, 'iyr1clm = ', iyr1clm
      print*, 'monnclm = ', monnclm
      print*, 'iyrnclm = ', iyrnclm

      print*, 'mon1   = ', mon1
      print*, 'iyr1   = ', iyr1
      print*, 'monn   = ', monn
      print*, 'iyrn   = ', iyrn

      print*, 'mon1out = ', mon1out
      print*, 'iyr1out = ', iyr1out
      print*, 'monnout = ', monnout
      print*, 'iyrnout = ', iyrnout

      print*, 'mon1in  = ', mon1in

      print*, 'ibasemon = ', ibasemon
      print*, 'ibaseyr = ', ibaseyr

      print*, ' '
      print*, 'conv = ', conv
      print*, 'wtmin = ', wtmin
      print*, 'maxiter = ', maxiter
      print*, 'bbmin = ', bbmin
      print*, 'histo = ', histo
      print*, 'varmin = ', varmin
      print*, 'amissin = ', amissin
      print*, 'amissout = ', amissout
      print*, 'amissasc = ', amissasc

      print*, ' '
      print*, 'correl = ', correl

      print*, ' '
      print*, 'monlen = ', monlen

      print*, ' '
      print*, 'tmin = ', tmin
      print*, 'tmax = ', tmax
      print*, 'dt = ', dt
      print*, 'dcbias = ', dcbias
      print*, 'fac =  ', fac
      print*, 'offset =  ', offset

      print*, ' '
      do 41 i=1,nfilesi
        print*, 'iyr1in = ', iyr1in(i)
        print*, 'inputfil = ', inputfil(i)
   41 continue

      print*, ' '
      print*, 'msklndi = ', msklndi
      print*, 'sftnamei = ', sftnamei
      print*, 'sftfilei = ', sftfilei
      print*, 'msklndo = ', msklndo
      print*, 'sftnameo = ', sftnameo
      print*, 'sftfileo = ', sftfileo

      print*, ' '
      print*, 'fclimbc = ', fclimbc
      print*, 'fbcinfo = ', fbcinfo
      print*, 'fspinup = ', fspinup
      print*, 'fclimobs = ', fclimobs
      print*, 'outfile = ', outfile
      print*, ' '
      print*, 'jyr1out = ', (jyr1out(i), i=1,nfileso)
      print*, 'fbc = ', fbc
      print*, 'fobs = ', fobs

      write(9,*) ' '
      write(9,*) 'Output from monthly boundary condition generator '
      write(9,*) 'varin = ', varin
      write(9,*) 'varout = ', varout
      write(9,*) 'varout1 = ', varout1
      write(9,*) 'src = ', src
      write(9,*) 'src1 = ', src1
      write(9,*) 'model = ', model
      write(9,*) 'grid = ', grid
      write(9,*) 'gtype = ', gtype
      write(9,*) 'inftype = ', inftype
      write(9,*) 'outftype = ', outftype
      write(9,*) 'calclim = ', calclim
      write(9,*) 'center = ', center
      write(9,*) 'abbrev = ', abbrev
      write(9,*) 'suffix = ', suffix
      write(9,*) 'suff = ', suff

      write(9,*) ' '
      write(9,*) 'icntl = ', icntl
      write(9,*) 'lagcalc = ', lagcalc
      write(9,*) 'iout = ', iout

      write(9,*) ' '
      write(9,*) 'nlon = ', nlon
      write(9,*) 'nlat = ', nlat
      write(9,*) 'mlat = ', mlat
      write(9,*) 'nchunks = ', nchunks

      write(9,*) ' '
      write(9,*) 'nmon = ', nmon
      write(9,*) 'mon1rd = ', mon1rd
      write(9,*) 'iyr1rd = ', iyr1rd
      write(9,*) 'monnrd = ', monnrd
      write(9,*) 'iyrnrd = ', iyrnrd

      write(9,*) 'mon1clm = ', mon1clm
      write(9,*) 'iyr1clm = ', iyr1clm
      write(9,*) 'monnclm = ', monnclm
      write(9,*) 'iyrnclm = ', iyrnclm

      write(9,*) 'mon1   = ', mon1
      write(9,*) 'iyr1   = ', iyr1
      write(9,*) 'monn   = ', monn
      write(9,*) 'iyrn   = ', iyrn

      write(9,*) 'mon1out = ', mon1out
      write(9,*) 'iyr1out = ', iyr1out
      write(9,*) 'monnout = ', monnout
      write(9,*) 'iyrnout = ', iyrnout

      write(9,*) 'mon1in  = ', mon1in

      write(9,*) 'ibasemon = ', ibasemon
      write(9,*) 'ibaseyr = ', ibaseyr

      write(9,*) ' '
      write(9,*) 'conv = ', conv
      write(9,*) 'wtmin = ', wtmin
      write(9,*) 'maxiter = ', maxiter
      write(9,*) 'bbmin = ', bbmin
      write(9,*) 'histo = ', histo
      write(9,*) 'varmin = ', varmin
      write(9,*) 'amissin = ', amissin
      write(9,*) 'amissout = ', amissout
      write(9,*) 'amissasc = ', amissasc

      write(9,*) ' '
      write(9,*) 'correl = ', correl

      write(9,*) ' '
      write(9,*) 'monlen = ', monlen

      write(9,*) ' '
      write(9,*) 'tmin = ', tmin
      write(9,*) 'tmax = ', tmax
      write(9,*) 'dt = ', dt
      write(9,*) 'dcbias = ', dcbias
      write(9,*) 'fac =  ', fac
      write(9,*) 'offset =  ', offset

      write(9,*) ' '
      do 42 i=1,nfilesi
        write(9,*) 'iyr1in = ', iyr1in(i)
        write(9,*) 'inputfil = ', inputfil(i)
   42 continue

      write(9,*) ' '
      write(9,*) 'msklndi = ', msklndi
      write(9,*) 'sftnamei = ', sftnamei
      write(9,*) 'sftfilei = ', sftfilei
      write(9,*) 'msklndo = ', msklndo
      write(9,*) 'sftnameo = ', sftnameo
      write(9,*) 'sftfileo = ', sftfileo

      write(9,*) ' '
      write(9,*) 'fclimbc = ', fclimbc
      write(9,*) 'fbcinfo = ', fbcinfo
      write(9,*) 'fspinup = ', fspinup
      write(9,*) 'fclimobs = ', fclimobs
      write(9,*) 'outfile = ', outfile

      write(9,*) ' '
      write(9,*) 'jyr1out = ', (jyr1out(i), i=1,nfileso)
      write(9,*) 'fbc = ', fbc
      write(9,*) 'fobs = ', fobs
      write(9,*) ' '

c   checked 7/14/03 that nlat/mlat need not be integer
c      if (mod(nlat,mlat) .ne. 0) then
c        print*, 'error in specification of mlat'
c        print*, 'nlat/mlat must be an integer'
c        stop
c      endif

      if (inftype .ne. 'pp') then
        if ((mon1rd+nmon12*iyr1rd) .lt. (mon1in+nmon12*iyr1in(1))) then
          print*, 'error in time specifications:  first month read'
          print*, 'must not preceed first month in input file.'
          print*, 'check mon1rd, iyr1rd, mon1in, iyr1in(1)'
          stop
        endif
      endif

      if ((mon1rd+nmon12*iyr1rd) .gt. (monnrd+nmon12*(iyrnrd-1))) then
        print*, 'error in time specifications:  last month read must'
        print*, 'follow first month read by at least 1 year'
        print*, 'check mon1rd, iyr1rd, monnrd, iyrnrd'
        stop
      endif

      if ((mon1+nmon12*iyr1) .gt. (mon1rd+nmon12*iyr1rd)) then
        print*, 'error in time specifications:  first month read must '
        print*, 'not precede first month to be processed'
        print*, 'check mon1, iyr1, mon1rd, iyr1rd'
        stop
      endif

      if ((monn+nmon12*iyrn) .lt. (monnrd+nmon12*iyrnrd)) then
        print*, 'error in time specifications: last month read must '
        print*, 'not follow last month processed'
        print*, 'check monn, iyrn, monnrd, iyrnrd'
        stop
      endif

      if (mod((monn+nmon12-mon1+1),nmon12) .ne. 0) then
        print*, 'error in time specifications:  only integral number'
        print*, 'of years allowed'
        print*, 'check monn, mon1'
        stop
      endif

      if ((monn-mon1+nmon12*(iyrn-iyr1)+1) .ne. nmon) then
        print*, 'error in time specifications:  parameter nmon '
        print*, 'must be consistent with first and last months'
        print*, 'specified'
        print*, 'check nmon, mon1, iyr1, monn, iyrn'
        stop
      endif

      if ((mon1clm+nmon12*iyr1clm) .lt. (mon1rd+nmon12*iyr1rd)) then
        print*, 'error in time specifications:  first month '
        print*, 'contributing to climatology must not precede first'
        print*, 'month read'
        print*, 'check mon1clm, iyr1clm, mon1rd, iyr1rd'
        stop
      endif

      if ((monnclm+nmon12*iyrnclm) .gt. (monnrd+nmon12*iyrnrd)) then
        print*, 'error in time specifications:  last month '
        print*, 'contributing to climatology must not follow last'
        print*, 'month read'
        print*, 'check monnclm, iyrnclm, monnrd, iyrnrd'
        stop
      endif


c   calculate jacobian elements for climatological years

      m = nmon12
      mp = 1
      do 63 n=1,nmon12
        mm = m
        m = mp
        mp = mp + 1
        if (mp .eq. (nmon12+1)) mp = 1

        ac(n) = 2.0*monlen(m,1)/float(monlen(m,1)+monlen(mm,1))
        cc(n) = 2.0*monlen(m,1)/float(monlen(m,1)+monlen(mp,1))

   63 continue

c   calculate jacobian elements for all years

      i = iyr1
      j = 1
      if (((mod(i,4) .eq. 0) .and. (mod(i,100) .ne. 0))
     &    .or. (mod(i,400) .eq. 0)) j = 2
      m = mod((mon1+nmon12-2), nmon12) + 1
      mp = mod((mon1-1), nmon12) + 1

      do 65 n=1,nmon

        mm = m 
        m  = mp
        mp = mp + 1

        if (mp .eq. (nmon12+1)) then
          mp = 1
          i = i + 1
          j = 1
          if (((mod(i,4) .eq. 0) .and. (mod(i,100) .ne. 0)) .or.
     &        (mod(i,400) .eq. 0)) j = 2
        endif

        a(n) = 2.0*monlen(m,j)/float(monlen(m,j)+monlen(mm,j))
        c(n) = 2.0*monlen(m,j)/float(monlen(m,j)+monlen(mp,j))

   65 continue

      do 100 i=1,nzonp
        aovarmon(i) = 0.0
        acvarmon(i) = 0.0
        wovar(i) = 0.0
        wcvar(i) = 0.0
        wocorr(i) = 0.0
        wccorr(i) = 0.0
        if (i .lt. nzonp) then
          woseas(i) = 0.0
          wcseas(i) = 0.0
          do 70 m=1,nmon12
            acentmon(m,i) = 0.0
            aobsclim(m,i) = 0.0
   70     continue
        endif
        do 75 m=1,nlagm
          ccorrel(m,i) = 0.0
          ocorrel(m,i) = 0.0
   75   continue

        do 120 m=-10,10
          jcount(i,m) = 0
          mcount(i,m) = 0
  120   continue

        do 130 m=-20,20
          icount(i,m) = 0
          kcount(i,m) = 0
  130   continue

  100 continue

c      test
c      do 101 m=1,nlon
c        alons(m) = -180. + (m-0.5)*360./nlon
c  101 continue

c      do 102 m=1,nlat
c        alats(m) = 90. - (m-0.5)*180./nlat
c  102 continue

c   loop through latitude chunks

      ismc = 0
      jsmc = 0
      ismf = 0
      jsmf = 0

      isea = 0

      jn = 0

      do 1000 jchnk=1,nchunks

        j1 = jn + 1
        jn = min0((j1+mlat-1), nlat)


c    obtain monthly means


      lats = jn - j1 + 1

      m1 = (iyr1rd-iyr1)*nmon12 + mon1rd - mon1 + 1
      m2 = (iyrnrd-iyr1)*nmon12 + monnrd - mon1 + 1
      mons = m2 - m1 + 1

      print*, 'Reading input for chunk ', jchnk

c  test 
c        goto 2345

      if (inftype .eq. 'pp') then

        call readpp(inputfil, nlat, nlon, mlat, j1, lats, iyr1rd,
     &           mon1rd, iyr1in, mons, amissin, amissout, alons, alats, 
     &           sst(1,1,m1), sstwts, wtfrac, lpp, bpp, space)

        lpp(29) = 0
        bpp(5) = 0.0
        bpp(18) = amissout
        bpp(19) = 1.0

      else

        call getobs(iregrid, gtype, msklndi, sftfilei, sftnamei, 
     &        msklndo, sftfileo, sftnameo, inputfil, varin,
     &        nlon, nlat, mlat, j1, jn, iyr1in(1), mon1in, iyr1rd, 
     &        mon1rd, iyrnrd, monnrd, nmon12, amissin, amissout, wtmin,
     &        sst(1,1,m1), sstwts(1,j1), wtfrac(1,j1), alons, alats(j1),
     &        ntlat, rlat0, ntlon, rlon0)

        lpp(15) = nlat*nlon
c         the following indicates a regular lat/long grid boxes 
c             (grid points are box centres)
c ???      change for Gaussian grid?
        lpp(16) = 2
        lpp(18) = nlat
        lpp(19) = nlon
        lpp(22) = 2
        lpp(23) = lbfc

        bpp(11) = 90.0
c ???  following incorrect for Gaussian grid
        bpp(14) = alats(1) + 2.*alats(1)/(nlat-1)
        bpp(15) = -2.*alats(1)/(nlat-1)
        bpp(16) = alons(1)
        bpp(17) = alons(2) - alons(1)
        bpp(18) = amissout
        bpp(19) = 1.0

      endif

c  test
c 2345 continue

      print*, 'Finished reading input for chunk ', jchnk

      print*, 'latitudes ', alats(j1), ' through ', alats(jn)
     &       , ' read successfully.'
      print*, ' '

      write(9,*) 'latitudes ', alats(j1), ' through ', alats(jn)
     &       , ' read successfully.'
      write(9,*) ' '

      do 19 j=j1,jn
        do 18 i=1,nlon
          if (wtfrac(i,j) .gt. wtmin) then
            do 17 m=m1,m2
              sst(i,j-j1+1,m) = fac*sst(i,j-j1+1,m) + offset
   17       continue
          endif
   18   continue
   19 continue

      if (iout(4) .gt. 0) then

c         generate file name for chunk of data

          write(c3, '(i3.3)') jchnk
          tempobs(jchnk) = pathout(1:je)//'/tempobs_'//c3

          nn = (iyrnout-iyr1out)*nmon12 + monnout-mon1out + 1
          kk = (iyr1out-iyr1)*nmon12 + mon1out-mon1 + 1
          call wrtpp(tempobs(jchnk), 1, nlon, mlat, alats(j1), lats, 
     &        monlen, iyr1out, mon1out, nn, lpp, bpp, sst(1,1,kk))

          tempobs(jchnk) = pathout(1:je)//'/tempobs_'//c3//'.pp'
      endif

      do 29 j=j1,jn
        do 28 i=1,nlon
          if (wtfrac(i,j) .gt. wtmin) then
            do 27 m=m1,m2
              sst(i,j-j1+1,m) = sst(i,j-j1+1,m) - dcbias
   27       continue
          endif
   28   continue
   29 continue
  
c      compute climatological means

      write(9,*) 
      write(9,*) 'Calculating climatological values'
c test  
c      go to 3456

      call calcclim(j1, jn, nlon, nlat, mlat, nmon, 
     &    ismc, jsmc, nzon, nzonp, nlagm, nmon12, mon1clm, iyr1clm,
     &    iyr1, mon1, iyr1rd, mon1rd, iyrnrd, monnrd,
     &    monnclm, iyrnclm, lagcalc, maxiter, amissout, histo, dcbias, 
     &    tmin, tmax, dt, wtmin, conv, bbmin, varmin, lagm, 
     &    isea, icount, jcount, kcount, mcount, wtfrac, sstwts, ovarmon,
     &    aovarmon, obsclim, aobsclim, sst, woseas, wovar, 
     &    ocorrel, centmon, acentmon, cmax, cmin, omax, omin,
     &    vecin, vecout, alons, alats, ac, cc, wcseas, wocorr)

      print*, ' '
      print*, ' '
      print*, ' '
      print*, '***********  finished climatology *************'
      print*, ' '
      print*, ' '
      print*, ' '

c     solve for full mid-month values

      write(9,*) 
      write(9,*) 'Calculating actual monthly values'

      call calcfull(j1, jn, nlon, nlat, mlat, nmon, 
     &    iyr1rd, mon1rd, iyr1, mon1, iyrn, monn, iyrnrd, monnrd,
     &    mon1clm, monnclm, iyr1clm, iyrnclm, ismf, jsmf,
     &    nzon, nzonp, nmon12, lagcalc, nlagm,
     &    maxiter, wtmin, amissout, tmin, tmax,  
     &    dt, conv, bbmin, varmin, lagm, sstwts, wtfrac, 
     &    wcseas, cvarmon, acvarmon, correl, ccorrel, 
     &    obsclim, vecin, vecout, centmon, sst, alons, alats, a, c, 
     &    ovarmon, wcvar, wccorr)


      call finish(j1, jn, nlon, nlat, mlat, mon1, nmon, nmon12, wtmin, 
     &     dcbias, amissout, wtfrac, centmon, obsclim, sst)

c test 
c 3456  continue

        if (iout(1) .gt. 0) then

c         generate file name for chunk of data

          write(c3, '(i3.3)') jchnk
          tempbc(jchnk) = pathout(1:je)//'/tempbc_'//c3

          nn = (iyrnout-iyr1out)*nmon12 + monnout-mon1out + 1
          kk = (iyr1out-iyr1)*nmon12 + mon1out-mon1 + 1

          call wrtpp(tempbc(jchnk), 0, nlon, mlat, alats(j1), lats, 
     &       monlen, iyr1out, mon1out, nn, lpp, bpp, sst(1,1,kk))
          tempbc(jchnk) = pathout(1:je)//'/tempbc_'//c3//'.pp'

        endif

        if (iout(2) .gt. 0) then
c         store spin-up data:  starts in first January preceeding
c           (and nearest to) first output month (or if this month is  
c           unavailable, starts at imon1, iyr1); ends in December of 
c           same year

c          m1 = (iyr1out-iyr1)*nmon12 - mon1 + 2
c          m2 = m1 + nmon12 - 1
c          mon1spin = 1
c          if (m1 .gt. 0) then
c            mon1spin = 1
c          else
c            m1 = 1
c            mon1spin = mon1
c          endif
c          m2 = max0(m2, 1)
c          nmonspin = m2-m1+1
c          iyr1spin = iyr1 + (mon1-2+m1)/nmon12

          mon1spin = 1
          if (mon1out .eq. 1) then
            iyr1spin = iyr1out -1
          else
            iyr1spin = iyr1out
          endif

          if (iyr1spin .le. iyr1) then
            iyr1spin = iyr1
            mon1spin = mon1
          endif

          kyr1out(1) = iyr1spin

          m1 = (iyr1spin - iyr1)*nmon12 + (mon1spin - mon1) + 1
          m2 = m1 + (nmon12-mon1spin)
          nmonspin = m2-m1+1

          do 325 m=m1,m2
            do 320 j=j1,jn
              do 310 i=1,nlon
                spinup(i,j,m-m1+1) = sst(i,j-j1+1,m)
  310         continue
  320       continue
  325     continue

        endif

1000  continue

      print*, ' '
      print*, ismc, ' climatological monthly means were smoothed'
      print*, jsmc, ' grid cells were affected'
      print*, ' '
      write(9,*) ' '
      write(9,*) ismc, ' climatological monthly means were smoothed'
      write(9,*) jsmc, ' grid cells were affected'
      write(9,*) ' '

      write(9,*) ' '
      write(9,*) 'ocean grid cells = ', isea
      write(9,*) ' '
c
      print*, ' '
      print*, ismf, ' monthly means were smoothed'
      print*, jsmf, ' grid cells were affected'
      print*, ' '
      write(9,*) ' '
      write(9,*) ismf, ' monthly means were smoothed'
      write(9,*) jsmf, ' grid cells were affected'
      write(9,*) ' '

        do 191 j=1,nzonp
          if (wovar(j) .gt. 0.0) aovarmon(j) = aovarmon(j)/wovar(j)
  191   continue

        do 195 j=1,nzon

          if (woseas(j) .gt. 0.0) then
            do 194 n=1,nmon12
              aobsclim(n,j) = aobsclim(n,j)/woseas(j)
  194       continue
          else
            aovarmon(j) = amissout
            do 196 n=1,nmon12
              aobsclim(n,j) = amissout
  196       continue
          endif

  195   continue

      if (lagcalc .gt. 0) then
        do 260 ik=1,nzonp

          if (wocorr(ik) .gt. 0.0) then
            do 255 k=1,nlagm
              ocorrel(k,ik) = ocorrel(k,ik)/wocorr(ik)
  255       continue
c          else
c            do 261 k=1,nlagm
c              ocorrel(k,ik) = amissout
c  261       continue
          endif

  260   continue

        do 560 ik=1,nzonp
          if (wccorr(ik) .gt. 0.0) then
            do 555 k=1,nlagm
                ccorrel(k,ik) = ccorrel(k,ik)/wccorr(ik)
  555       continue
c          else
c            do 556 k=1,nlagm
c              ccorrel(k,ik) = amissout
c  556       continue
          endif
  560   continue

      endif

      do 505 ik=1,nzon
        do 502 n=-20,20
          icount(nzonp,n) = icount(nzonp,n)+icount(ik,n)
          kcount(nzonp,n) = kcount(nzonp,n)+kcount(ik,n)
  502   continue
        do 503 n=-10,10
          jcount(nzonp,n) = jcount(nzonp,n)+jcount(ik,n)
          mcount(nzonp,n) = mcount(nzonp,n)+mcount(ik,n)
  503   continue

        do 504 m=1,nmon12
          if (wcseas(ik) .gt. 0.0) then
            acentmon(m,ik) = acentmon(m,ik)/wcseas(ik)
          else
            acentmon(m,ik) = amissout
          endif
  504   continue
  505 continue

      do 906 ik=1,nzonp
        if (wcvar(ik) .gt. 0.0) then
           acvarmon(ik) = acvarmon(ik)/wcvar(ik)
        elseif (ik .lt. nzonp) then
          if (wcseas(ik) .eq. 0.0) acvarmon(ik) = amissout
        endif
  906 continue

      write(9,*) 
     &'Number of cells where difference in climatological max is '
      write(9,*) 'within +- ', histo,  ' of 0.0: '
      write(9,*) (icount(ik,0), ik=1,nzonp)
      write(9,*) ' '
      write(9,*) 'number of cells where difference between constructed'
      write(9,*) 'mid-point climatological monthly maximum and observed'
      write(9,*) 'climatological monthly mean exceeds different values:'
      write(9,*) ' '

      if (alats(1) .lt. alats(2)) then
        write(9,*) 'value       ',
     &   '-90 & -60 -60 & -30  -30 & 0     0 & 30   30 & 60   60 & 90',
     &   '  Global'
      else
        write(9,*) 'value       ',
     &   ' 90 & 60   60 & 30    30 & 0    0 & -30  -30 & -60 -60 & -90',
     &   ' Global'
      endif
      write(9,*) 'exceeded'
      do 201 i=-20,20
        if (i .ne. 0) then
          r = i*histo
          write(9,'(f10.3, 7i10)') r, (icount(ik,i), ik=1,nzonp)
        endif
  201 continue
      write(9,*) ' '
      do 301 i=-10,10
        if (i .ne. 0) then
          r = i*10*histo
          write(9,'(f10.3, 7i10)') r, (jcount(ik,i), ik=1,nzonp)
        endif
  301 continue    
      write(9,*) ' '
      write(9,*) 
     &'Number of cells where difference in climatological min is '
      write(9,*) 'within  +- ', histo,  ' of 0.0: '
      write(9,*) (kcount(ik,0), ik=1,nzonp)
      write(9,*) ' '
      write(9,*)'number of cells where -(difference between constructed'
      write(9,*) 'mid-point climatological monthly minimum and observed'
      write(9,*)
     &  'climatological monthly mean minimum) exceeds different values:'
      write(9,*) ' '
      if (alats(1) .lt. alats(2)) then
        write(9,*) 'value       ',
     &   '-90 & -60 -60 & -30  -30 & 0     0 & 30   30 & 60   60 & 90',
     &   '  Global'
      else
        write(9,*) 'value       ',
     &   ' 90 & 60   60 & 30    30 & 0    0 & -30  -30 & -60 -60 & -90',
     &   ' Global'
      endif
      write(9,*) 'exceeded'
      do 202 i=-20,20
        if (i .ne. 0) then
          r = i*histo
          write(9,'(f10.3, 7i10)') r, (kcount(ik,i), ik=1,nzonp)
        endif
  202 continue
      write(9,*) ' '
      do 302 i=-10,10
        if (i .ne. 0) then
          r = i*10*histo
          write(9,'(f10.3, 7i10)') r, (mcount(ik,i), ik=1,nzonp)
        endif
  302 continue 

      write(9,*) ' '
      write(9,*) ' '
      write(9,*) 'Climatological Seasonal Variation by Zone'
      write(9,*) ' '

      if (alats(2) .lt. alats(1)) then
        write(9,*) 'month  ',
     & '    90 & 60     60 & 30       30 & 0      0 & -30   -30 & -60',
     &   '  -60 & -90'
      else
        write(9,*) 'month  ',
     & '   -90 & -60   -60 & -30     -30 & 0      0 & 30     30 & 60 ',
     &   '    60 & 90'
      endif
      write(9,*)'         avg   mid   avg   mid   avg   mid   avg',
     &   '   mid   avg   mid   avg   mid'
      write(9,'((i2, 6x, 12f6.2))')
     &     (n, (aobsclim(n,ik), acentmon(n,ik), ik=1,nzon), n=1,nmon12)
      write(9,'(/ "area    ", 12f6.3)')
     &             (woseas(ik), wcseas(ik), ik=1,nzon)
c
      write(9,*) ' '
      write(9,*) ' '
      write(9,*) 'Variance of Monthly Anomalies by Zone'
      write(9,*) ' '

      if (alats(2) .lt. alats(1)) then
        write(9,*) 
     &  '          90 & 60    60 & 30     30 & 0    0 & -30 -30 & -60',
     &   ' -60 & -90  global'
      else
        write(9,*) 
     &  '         -90 & -60  -60 & -30   -30 & 0    0 & 30   30 & 60 ',
     &   '  60 & 90   global'
      endif
      write(9,'("mon avg.", 7f10.3)') (aovarmon(ik), ik=1,nzonp)
      write(9,'("mid mon ", 7f10.3)') (acvarmon(ik), ik=1,nzonp)
      write(9,*) ' '
      write(9,'("mon area", 7f10.5)') (wovar(ik), ik=1,nzonp)
      write(9,'("mid area", 7f10.5)') (wcvar(ik), ik=1,nzonp)
      write(9,*) ' '

      if (lagcalc .gt. 0) then

      write(9,*) ' '
      write(9,*) 'Variance and Auto Correlation by Zone'
      write(9,*) ' '
      write(9,*) '       Monthly Mean Anomalies '
      if (alats(2) .lt. alats(1)) then
        write(9,*) 
     &  '          90 & 60    60 & 30     30 & 0    0 & -30 -30 & -60',
     &   ' -60 & -90  global'
      else
        write(9,*) 
     &  '         -90 & -60  -60 & -30   -30 & 0    0 & 30   30 & 60 ',
     &   '  60 & 90   global'
      endif
      if (alats(2) .lt. alats(1)) then
        write(9,*) 
     &  '          90 & 60    60 & 30     30 & 0    0 & -30 -30 & -60',
     &   ' -60 & -90  global'
      else
        write(9,*) 
     &  '         -90 & -60  -60 & -30   -30 & 0    0 & 30   30 & 60 ',
     &   '  60 & 90   global'
      endif
      write(9,'("variance:", 7f10.3 )') (ocorrel(1,ik), ik=1,nzonp)
      write(9,*) 'lag          Auto Correlation'
      write(9,'(i8, 7f10.3)') (lagm(k), (ocorrel(k,ik), ik=1,nzonp), 
     &       k=2,nlagm)
      write(9,*) ' '
      write(9,'("area    ", 7f10.5)') (wocorr(ik), ik=1,nzonp)

      write(9,*) ' '
      write(9,*) '       Mid-Month Anomalies '
      if (alats(2) .lt. alats(1)) then
        write(9,*) 
     &  '          90 & 60    60 & 30     30 & 0    0 & -30 -30 & -60',
     &   ' -60 & -90  global'
      else
        write(9,*) 
     &  '         -90 & -60  -60 & -30   -30 & 0    0 & 30   30 & 60 ',
     &   '  60 & 90   global'
      endif
      write(9,'("variance:", 7f10.3 )') (ccorrel(1,ik), ik=1,nzonp)
      write(9,*) 'lag          Auto Correlation'
      write(9,'(i8, 7f10.3)') (lagm(k), (ccorrel(k,ik), ik=1,nzonp), 
     &       k=2,nlagm)
      write(9,*) ' '
      write(9,'("area    ", 7f10.5)') (wccorr(ik), ik=1,nzonp)

      write(9,*) ' '

      endif

      if ((outftype(1:2) .eq. 'pp') .or. (outftype(1:3) .eq. 'all'))
     &        then

        if (iout(1) .gt. 0) then
          call ppconcat (tempbc, fbc, fnclim, nmon12, nlon, nlat, 
     &         nchunks, iyr1out, mon1out, iyrnout, monnout, jyr1out,
     &         space)
        endif

        if (iout(2) .gt. 0) then
c         write spin-up data

          write(ayr, '(i4)') iyr1spin
          iend = index(fspinup, ' ') - 1
          outcopy = fspinup(1:iend) // '_' //ayr//' '
            
          call wrtpp(outcopy, 0, nlon, nlat, alats(1), nlat, 
     &        monlen, iyr1spin, mon1spin, nmonspin, lpp, bpp, spinup)
        endif

        if (iout(3) .gt. 0) then
c         write climatology b.c.
          iend = index(fclimbc, ' ') - 1
          outcopy = fclimbc(1:iend) // '_' // fnclim // ' '
          call wrtpp(outcopy, 0, nlon, nlat, alats(1),
     &        nlat, monlen, 0, 1, nmon12, lpp, bpp, centmon)
        endif

        if (iout(4) .gt. 0) then
          call ppconcat (tempobs, fobs, fnclim, nmon12, nlon, nlat, 
     &         nchunks, iyr1out, mon1out, iyrnout, monnout, jyr1out,
     &         space)
        endif

        if (iout(5) .gt. 0) then
c         write climatology of obs
          iend = index(fclimobs, ' ') - 1
          outcopy = fclimobs(1:iend) // '_' // fnclim // ' '
          call wrtpp(outcopy, 1, nlon, nlat, alats(1), nlat, 
     &        monlen, 0, 1, nmon12, lpp, bpp, obsclim)
        endif

      endif



      if ((outftype(1:3) .eq. 'cdf') .or. 
     &        (outftype(1:6) .eq. 'coards') .or.
     &        (outftype(1:5) .eq. 'grads') .or.
     &        (outftype(1:4) .eq. 'grib') .or.
     &        (outftype(1:11) .eq. 'coards&grib') .or.
     &        (outftype(1:10) .eq. 'coards&asc') .or.
     &        (outftype(1:6) .eq. 'notdrs') .or.
     &        (outftype(1:3) .eq. 'all')) then

        if ((outftype(1:3) .eq. 'all') .or. 
     &     (outftype(1:6) .eq. 'notdrs') .or.
     &     (outftype(1:11) .eq. 'coards&grib') .or. 
     &     (outftype(1:10) .eq. 'coards&asc')) then
          typfil = 'coards'
          if (outftype(1:10) .eq. 'coards&asc') then
             i1 = 1
          else
             i1 = 0
          endif
        else 
          typfil = outftype
          i1 = 1
        endif

  330   continue

        if (typfil(1:6) .eq. 'coards') then
          calclim = 'gregorian'
        else
          calclim = 'climatology'
        endif
 
        if (iout(1) .gt. 0) then

          title = src(1:isrc)//
     &      ' Boundary Condition Data: Mid-Month '//suff(1:iie)

          nn = (iyrnout-iyr1out)*nmon12 + monnout-mon1out + 1
          call wrtlats(ig, 1, 1, 1, 1, 1, 1,  
     &      'mid', 'instantaneous', varout, fbc, fnclim, typfil, 
     &      dfltparm, gtype, calgreg, title, sourceb, model, center, 
     &      amissout, nlon, nlat, nchunks, alons, alats, nn,  
     &      mon1out, iyr1out, ibasemon, ibaseyr, idvar, space,
     &      tempbc, jyr1out)
          ig = 0

        endif

        if (iout(2) .gt. 0) then
      
          title = src(1:isrc)//
     &          ' Spin-up Boundary Condition Data: Mid-Month '
     &          //suff(1:iie)
          call wrtlats(ig, 1, 1, 1, 1, 1, 0, 
     &      'mid', 'instantaneous', varout, fspinup, fnclim, typfil, 
     &      dfltparm, gtype, calgreg, title, sourceb, model, center,  
     &      amissout, nlon, nlat, nchunks, alons, alats, nmonspin, 
     &      mon1spin, iyr1spin, ibasemon, ibaseyr, idvar, spinup,
     &      dum, kyr1out)
          ig = 0

        endif

        if (iout(3) .gt. 0) then


          title = cclim //
     &       ' Clim. Bdry. Condition Data: Mid-Month '
     &       //suff(1:iie)

          call wrtlats(ig, 1, 1, 1, 1, 1, 0,  
     &      'mid', 'instantaneous', varout, fclimbc, fnclim, typfil, 
     &      dfltparm, gtype, calclim, title, sourceb, model, center, 
     &      amissout, nlon, nlat, nchunks, alons, alats, nmon12, 
     &      1, 2, 1, 2, idvar, centmon, dum, jyr1out)
          ig = 0

        endif
        
        if (iout(4) .gt. 0) then

          title = ' Observed Monthly Mean '//suff(1:iie)
          nn = (iyrnout-iyr1out)*nmon12 + monnout-mon1out + 1

          call wrtlats(ig, 1, 1, 1, 1, 1, 1, 
     &      'mean', 'average', varout1, fobs, fnclim, typfil, 
     &      dfltparm, gtype, calgreg, title, sourceo, model, center, 
     &      amissout, nlon, nlat, nchunks, alons, alats, nn, 
     &      mon1out, iyr1out, ibasemon, ibaseyr, idvar, space, 
     &      tempobs, jyr1out)
          ig = 0
        endif
       
        if (iout(5) .gt. 0) then

          title =  cclim //
     &            ' Observed Climatology: Monthly Mean '//suff(1:iie)
          call wrtlats(ig, 1, 1, 1, 1, 1, 0, 
     &      'mean','average', varout1, fclimobs, fnclim, typfil, 
     &      dfltparm, gtype, calclim, title, sourceo, model, center, 
     &      amissout, nlon, nlat, nchunks, alons, alats, nmon12, 
     &      1, 2, 1, 2, idvar, obsclim, dum, jyr1out)
          ig = 0

        endif

        if (iout(6) .gt. 0) then
cjfp no calls of wrtlats, effectively ignoring iout(6):
c        if (iout(6) .gt. 0 .and. iout(6).lt.0 ) then
          vname = varout1//'max'
          title = 
     &        ' Maximum (Monthly) Mean Climatological '// suff(1:iie)

          call wrtlats(ig, 1, 1, 1, 0, 0, 0, 
     &      'none', 'average', vname, fbcinfo, blnk, typfil, 
     &      parmtabl, gtype, calclim, title, sourceo, model, center, 
     &      amissout, nlon, nlat, nchunks, alons, alats, 0, 
     &      1, 0, 1, 0, idvara(1), omax, dum, jyr1out)
          ig  = 0


          vname = varout//'max'
          title = 
     &        ' Maximum (Mid-Month) Climatological '// suff(1:iie)

          call wrtlats(0, 0, 1, 1, 0, 0, 0, 
     &      'none', 'average', vname, fbcinfo, blnk, typfil, 
     &      'none', gtype, calclim, title, sourceb, model, center, 
     &      amissout, nlon, nlat, nchunks, alons, alats, 0, 
     &      1, 0, 1, 0, idvara(2), cmax, dum, jyr1out)

          vname = varout1//'min'
          title = 
     &        ' Minimum (Monthly) Mean Climatological '// suff(1:iie)

          call wrtlats(0, 0, 1, 1, 0, 0, 0,
     &      'none', 'average', vname, fbcinfo, blnk, typfil, 
     &      'none', gtype, calclim, title, sourceo, model, center, 
     &      amissout, nlon, nlat, nchunks, alons, alats, 0, 
     &      1, 0, 1, 0, idvara(3), omin, dum, jyr1out)

          vname = varout//'min'
          title = 
     &       ' Minimum (Mid-Month) Climatological '// suff(1:iie)

          call wrtlats(0, 0, 1, 1, 0, 0, 0, 
     &      'none', 'average', vname, fbcinfo, blnk, typfil, 
     &      'none', gtype, calclim, title, sourceb, model, center, 
     &      amissout, nlon, nlat, nchunks, alons, alats, 0, 
     &      1, 0, 1, 0, idvara(4), cmin, dum, jyr1out)

          vname = varout1//'var'
          m1 = index(suff, '(')
          if (m1 .gt. 0) then
            m1 = m1 - 1
          else
            m1 = iie
          endif
          
          title = 
     &        ' Variance of Monthly Mean Anomalies: '// suff(1:m1) 

          call wrtlats(0, 0, 1, 1, 0, 0, 0, 
     &      'none', 'average', vname, fbcinfo, blnk, typfil, 
     &      'none', gtype, calclim, title, sourceo, model, center, 
     &      amissout, nlon, nlat, nchunks, alons, alats, 0, 
     &      1, 0, 1, 0, idvara(5), ovarmon, dum, jyr1out)

          vname = varout//'var'
          title = 
     &        ' Variance of Mid-Month Anomalies: '// suff(1:m1) 

          call wrtlats(0, 0, 1, 1, 0, 0, 0, 
     &      'none', 'average', vname, fbcinfo, blnk, typfil, 
     &      'none', gtype, calclim, title, sourceb, model, center, 
     &      amissout, nlon, nlat, nchunks, alons, alats, 0, 
     &      1, 0, 1, 0, idvara(6), cvarmon, dum, jyr1out)


          vname = varout1//'max'
          call wrtlats(0, 0, 0, 1, 1, 0, 0, 
     &      'none', 'average', vname, fbcinfo, blnk, typfil, 
     &      'none', gtype, calclim, title, sourceo, model, center, 
     &      amissout, nlon, nlat, nchunks, alons, alats, 0, 
     &      1, 0, 1, 0, idvara(1), omax, dum, jyr1out)

          vname = varout//'max'
          call wrtlats(0, 0, 0, 1, 1, 0, 0, 
     &      'none', 'average', vname, fbcinfo, blnk, typfil, 
     &      'none', gtype, calclim, title, sourceb, model, center, 
     &      amissout, nlon, nlat, nchunks, alons, alats, 0, 
     &      1, 0, 1, 0, idvara(2), cmax, dum, jyr1out)

          vname = varout1//'min'
          call wrtlats(0, 0, 0, 1, 1, 0, 0, 
     &      'none', 'average', vname, fbcinfo, blnk, typfil, 
     &      'none', gtype, calclim, title, sourceo, model, center, 
     &      amissout, nlon, nlat, nchunks, alons, alats, 0, 
     &      1, 0, 1, 0, idvara(3), omin, dum, jyr1out)

          vname = varout//'min'
          call wrtlats(0, 0, 0, 1, 1, 0, 0, 
     &      'none', 'average', vname, fbcinfo, blnk, typfil, 
     &      'none', gtype, calclim, title, sourceb, model, center, 
     &      amissout, nlon, nlat, nchunks, alons, alats, 0, 
     &      1, 0, 1, 0, idvara(4), cmin, dum, jyr1out)

          vname = varout1//'var'
          call wrtlats(0, 0, 0, 1, 1, 0, 0, 
     &      'none', 'average', vname, fbcinfo, blnk, typfil, 
     &      'none', gtype, calclim, title, sourceo, model, center, 
     &      amissout, nlon, nlat, nchunks, alons, alats, 0, 
     &      1, 0, 1, 0, idvara(5), ovarmon, dum, jyr1out)

          vname = varout//'var'
          call wrtlats(0, 0, 0, 1, 1, 1, 0, 
     &      'none', 'average', vname, fbcinfo, blnk, typfil, 
     &      'none', gtype, calclim, title, sourceb, model, center, 
     &      amissout, nlon, nlat, nchunks, alons, alats, 0, 
     &      1, 0, 1, 0, idvara(6), cvarmon, dum, jyr1out)
    
        endif

        if (i1 .eq. 0) then
          typfil = 'grib'
          i1 = 1
          go to 330
        endif

      endif

      if ((outftype(1:3) .eq. 'drs') .or. (outftype(1:3).eq.'all'))
     &          then

        elem1(1) = alons(1)
        elemn(1) = alons(nlon)
        elem1(2) = alats(1)
        elemn(2) = alats(nlat)

        if (caseindp(gtype, 'gau') 
     +     .or. caseindp(gtype, 'bmr') .or. caseindp(gtype, 'ccc') 
     +     .or. caseindp(gtype, 'cnr') .or. caseindp(gtype, 'col') 
     +     .or. caseindp(gtype, 'csi') .or. caseindp(gtype, 'der') 
     +     .or. caseindp(gtype, 'ecm') .or. caseindp(gtype, 'gfd') 
     +     .or. caseindp(gtype, 'mgo') .or. caseindp(gtype, 'mpi')
     +     .or. caseindp(gtype, 'nca') .or. caseindp(gtype, 'nmc') 
     +     .or. caseindp(gtype, 'rpn') .or. caseindp(gtype, 'sun') 
     +     .or. caseindp(gtype, 'uga') .or. caseindp(gtype, 'sng')
     +     .or. caseindp(gtype, 'gen')
     +     .or. caseindp(gtype, 'ccm') .or. caseindp(gtype, 'ccs')
     +     .or. caseindp(gtype, 'ech') .or. caseindp(gtype, 'nce')) then

           do 108 i=1,nlat
             gauss(i) = alats(i)
  108      continue

        elseif (caseindp(gtype, 'cos') 
     +     .or. caseindp(gtype, 'dnm') .or. caseindp(gtype, 'gis') 
     +     .or. caseindp(gtype, 'csu') .or. caseindp(gtype, 'ucl')
     +     .or. caseindp(gtype, 'gla') .or. caseindp(gtype, 'gsf') 
     +     .or. caseindp(gtype, 'iap') .or. caseindp(gtype, 'mri') 
     +     .or. caseindp(gtype, 'uiu') 
     +     .or. caseindp(gtype, 'jma') .or. caseindp(gtype, 'nrl') 
     +     .or. caseindp(gtype, 'ukm') .or. caseindp(gtype, 'yon')) 
     +      then

          gauss(1) = 0.0
          gauss(2) = 0.0

        else
          print*, 'Grid type not recognized '
          print*, 'Grid type = ', gtype
          write(9,*) 'Grid type not recognized '
          write(9,*) 'Grid type = ', gtype
          stop
        endif


        if (iout(1) .gt. 0) then

          vname = varout
          title = src(1:isrc)//
     &      ' Boundary Condition Data: Mid-Month '//suff(1:iie)
     
          call wrtdrs(1, 1, 1, fbc, fnclim, 
     &         nlon, nlat, nmon12, nchunks, iyr1out, mon1out,
     &         iyrnout, monnout, tempbc, jyr1out, elem1, elemn,
     &         sourceb, vname, title, units, gauss, space)

        endif

        if (iout(2) .gt. 0) then

          vname = varout
          title = src(1:isrc)//
     &          ' Spin-up Boundary Condition Data: Mid-Month '
     &          //suff(1:iie)

           m2 = mon1spin+nmonspin-1

          call wrtdrs(1, 1, 0, fspinup, fnclim,
     &         nlon, nlat, nmon12, nchunks, iyr1spin, mon1spin,
     &         iyr1spin, m2, tempbc, kyr1out, elem1, elemn,
     &         sourceb, vname, title, units, gauss, spinup)

        endif

        if (iout(3) .gt. 0) then

          vname = varout
          title = cclim //
     &       ' Clim. Bdry. Condition Data: Mid-Month '
     &       //suff(1:iie)

          call wrtdrs(1, 1, 0, fclimbc, fnclim,
     &         nlon, nlat, nmon12, nchunks, 0, 1,
     &         0, nmon12, dum, jyr1out, elem1, elemn,
     &         sourceb, vname, title, units, gauss, centmon)

        endif

        if (iout(4) .gt. 0) then

          vname = varout1
          title = ' Observed Monthly Mean '//suff(1:iie)

          call wrtdrs(1, 1, 1, fobs, fnclim,
     &         nlon, nlat, nmon12, nchunks, iyr1out, mon1out,
     &         iyrnout, monnout, tempobs, jyr1out, elem1, elemn,
     &         sourceo, vname, title, units, gauss, space)

        endif

        if (iout(5) .gt. 0) then

          vname = varout1
          title =  cclim //
     &             ' Observed Climatology: Monthly Mean '//suff(1:iie)

          call wrtdrs(1, 1, 0, fclimobs, fnclim,
     &         nlon, nlat, nmon12, nchunks, 0, 1,
     &         0, nmon12, dum, jyr1out, elem1, elemn,
     &         sourceo, vname, title, units, gauss, obsclim)

        endif
   
      endif
  
       if ((outftype(1:3) .eq. 'drs') 
     &     .or. (outftype(1:3) .eq. 'all')) then

        if (iout(6) .gt. 0) then

          vname = varout1//'max'
          title = 
     &        ' Maximum (Monthly) Mean Climatological '// suff(1:iie)

          call wrtdrs(1, 0, 0, fbcinfo, blnk,
     &         nlon, nlat, nmon12, nchunks, 0, 1,
     &         0, 1, dum, jyr1out, elem1, elemn,
     &         sourceo, vname, title, units, gauss, omax)

          vname = varout//'max'
          title = 
     &        ' Maximum (Mid-Month) Climatological '// suff(1:iie)

          call wrtdrs(0, 0, 0, fbcinfo, blnk,
     &         nlon, nlat, nmon12, nchunks, 0, 1,
     &         0, 1, dum, jyr1out, elem1, elemn,
     &         sourceb, vname, title, units, gauss, cmax)

          vname = varout1//'min'
          title = 
     &        ' Minimum (Monthly) Mean Climatological '// suff(1:iie)

          call wrtdrs(0, 0, 0, fbcinfo, blnk,
     &         nlon, nlat, nmon12, nchunks, 0, 1,
     &         0, 1, dum, jyr1out, elem1, elemn,
     &         sourceo, vname, title, units, gauss, omin)

          vname = varout//'min'
          title = 
     &       ' Minimum (Mid-Month) Climatological '// suff(1:iie)

          call wrtdrs(0, 0, 0, fbcinfo, blnk,
     &         nlon, nlat, nmon12, nchunks, 0, 1,
     &         0, 1, dum, jyr1out, elem1, elemn,
     &         sourceb, vname, title, units, gauss, cmin)

          vname = varout1//'var'
          m1 = index(suff, '(')
          if (m1 .gt. 0) then
            m1 = m1 - 1
          else
            m1 = iie
          endif

          title = 
     &        ' Variance of Monthly Mean Anomalies: '// suff(1:m1) 

          call wrtdrs(0, 0, 0, fbcinfo, blnk,
     &         nlon, nlat, nmon12, nchunks, 0, 1,
     &         0, 1, dum, jyr1out, elem1, elemn,
     &         sourceo, vname, title, units2, gauss, ovarmon)

          vname = varout//'var'
          title = 
     &        ' Variance of Mid-Month Anomalies: '// suff(1:m1) 

          call wrtdrs(0, 1, 0, fbcinfo, blnk,
     &         nlon, nlat, nmon12, nchunks, 0, 1,
     &         0, 1, dum, jyr1out, elem1, elemn,
     &         sourceb, vname, title, units2, gauss, cvarmon)

        endif

      endif


      if ((outftype(1:5) .eq. 'ascii') .or. (outftype(1:3) .eq. 'all')
     &      .or. (outftype(1:10) .eq. 'coards&asc') .or. 
     &           (outftype(1:6) .eq. 'notdrs')) then

c         should call wrtasc last (if 'all') because missing data value
c          is different from other formats and sub wrtasc alters
c           the output arrays.

         elem1(1) = alons(1)
         elemn(1) = alons(nlon)
         elem1(2) = alats(1)
         elemn(2) = alats(nlat)

         if (iout(1) .gt. 0) then
           vname = varout
           title = src(1:isrc)//
     &      ' Boundary Condition Data: Mid-Month '//suff(1:iie)

           call wrtasc(fbc, fnclim, nlon, nlat, nmon12, iyr1out, 
     &        mon1out,
     &        iyrnout, monnout, elem1, elemn, model, center, sourceb, 
     &        vname, title, units, grid, amissout, amissasc, 1, 
     &        nchunks, tempbc, jyr1out, space)
         endif

         if (iout(2) .gt. 0) then
           vname = varout
           title = src(1:isrc)//
     &       ' Spin-up Boundary Condition Data: Mid-Month '
     &       //suff(1:iie)

           m2 = mon1spin+nmonspin-1
           call wrtasc(fspinup, fnclim, nlon, nlat, nmon12, iyr1spin, 
     &         mon1spin, 
     &         iyr1spin, m2, elem1, elemn, model, center, sourceb, 
     &         vname, title, units, grid, amissout, amissasc, 0,
     &         nchunks, dum, kyr1out, spinup)
         endif

         if (iout(3) .gt. 0) then
           vname = varout
           title = cclim //
     &         ' Clim. Bdry. Condition Data: Mid-Month '
     &         //suff(1:iie)

           call wrtasc(fclimbc, fnclim, nlon, nlat, nmon12, 0, 1, 
     &         0, nmon12, elem1, elemn, model, center, sourceb, 
     &         vname, title, units, grid, amissout, amissasc, 0,
     &         nchunks, dum, jyr1out, centmon)
         endif

         if (iout(4) .gt. 0) then
           vname = varout1
           title = ' Observed Monthly Mean '//suff(1:iie)

           call wrtasc(fobs, fnclim, nlon, nlat, nmon12, iyr1out, 
     &        mon1out,
     &        iyrnout, monnout, elem1, elemn, model, center, sourceo, 
     &        vname, title, units, grid, amissout, amissasc, 1, 
     &        nchunks, tempobs, jyr1out, space)

         endif

         if (iout(5) .gt. 0) then
           vname = varout1
           title =  cclim //
     &             ' Observed Climatology: Monthly Mean '//suff(1:iie)

           call wrtasc(fclimobs, fnclim, nlon, nlat, nmon12, 0, 1,
     &        0, nmon12, elem1, elemn, model, center, sourceo,
     &        vname, title, units, grid, amissout, amissasc, 0,
     &        nchunks, dum, jyr1out, obsclim)
         endif

      endif

      close(9)

      if (icntl .eq. 3) then
        icntl = 1
        go to 22
      endif

c      syscomm = 'rm '//pathout(1:je)//'/temp*.pp'
c      call system(syscomm)

      print*, ' '
      print*, ' '
      print*, 'Program mkgisst finished execution.'
      print*, ' '
      print*, ' '
      call exit(1)
      end




      subroutine finish(j1, jn, nlon, nlat, mlat, mon1, nmon, nmon12,
     &     wtmin, dcbias, amissout, wtfrac, centmon, obsclim, sst)

      implicit none
      integer nlon, nlat, mon1, nmon, nmon12, mlat
      real wtmin, dcbias, amissout
      real wtfrac(nlon,nlat), centmon(nlon,nlat,nmon12), 
     &    obsclim(nlon,nlat,nmon12), sst(nlon,mlat,nmon)

      integer i, j, m, n, j1, jn

c    calculate full mid-month time-series

      do 990 j=j1,jn
        do 980 i=1,nlon

          if (wtfrac(i,j) .le. wtmin) then

            do 920 n=1,nmon
              sst(i,j-j1+1,n) = amissout
  920       continue

          else

            do 930 m=1,nmon12
              centmon(i,j,m) = centmon(i,j,m) + dcbias
              obsclim(i,j,m) = obsclim(i,j,m) + dcbias
  930       continue   

            do 970 n=1,nmon
              m = mod((n+mon1-2), nmon12) + 1
              sst(i,j-j1+1,n) = sst(i,j-j1+1,n) + centmon(i,j,m)
  970       continue

          endif

  980   continue
  990 continue

      return
      end



      subroutine calcclim(j1, jn, nlon, nlat, mlat, nmon, 
     &    ismc, jsmc, nzon, nzonp, nlagm, nmon12, mon1clm, iyr1clm,  
     &    iyr1, mon1, iyr1rd, mon1rd, iyrnrd, monnrd,
     &    monnclm, iyrnclm, lagcalc, maxiter, amiss, histo, dcbias, 
     &    tmin, tmax, dt, wtmin, conv, bbmin, varmin, lagm, 
     &    isea, icount, jcount, kcount, mcount, wtfrac, sstwts, ovarmon,
     &    aovarmon, obsclim, aobsclim, sst, woseas, wovar, ocorrel,
     &    centmon, acentmon, cmax, cmin, omax, omin, vecin,
     &    vecout, alons, alats, ac, cc, wcseas, wocorr)
  
      implicit none
      integer nlon, nlat, nmon, nzon, j1, jn, mon1clm, iyr1clm, monnclm,
     &      iyrnclm, nzonp, nlagm, nmon12, 
     &      lagcalc, maxiter, isea, mlat,
     &      iyr1, mon1, iyr1rd, mon1rd, iyrnrd, monnrd

      real amiss, histo, dcbias, tmin, tmax, dt, wtmin, conv,
     &      bbmin, varmin

      integer icount(nzonp,-20:20), jcount(nzonp,-10:10),
     &      kcount(nzonp,-20:20), mcount(nzonp, -10:10), lagm(nlagm)

      real wtfrac(nlon,nlat), sstwts(nlon,nlat), ovarmon(nlon,nlat),
     &      aovarmon(nzonp), obsclim(nlon,nlat,nmon12), 
     &      aobsclim(nmon12,nzon), sst(nlon,mlat,nmon), 
     &      woseas(nzon), wovar(nzonp), ocorrel(nlagm,nzonp), 
     &      centmon(nlon,nlat,nmon12), 
     &      acentmon(nmon12,nzon), cmax(nlon,nlat), cmin(nlon,nlat),
     &      omax(nlon,nlat), omin(nlon,nlat), vecin(nmon), 
     &      vecout(nmon), alons(nlon), alats(nlat), ac(nmon12),
     &      cc(nmon12), wcseas(nzon), wocorr(nzonp)

      integer ismc, jsmc, i, j, ik, n, m, nn, mn, k, mm, momax, 
     &      momin, mmax, mmin, ii, jj, kk, mmnn, ji, m1, m2, m3, m4

      real var, ocentmax, ocentmin, centmax, centmin, diffmax, diffmin,
     &      rmin, rmax, rrmin, rrmax

      double precision sum, avg1, avg2, covar

        mmnn = 0
        rrmin = 0.0
        rrmax = 0.0

      m1 = nmon12*(iyr1rd-iyr1) + mon1rd-mon1 + 1
      m2 = nmon12*(iyrnrd-iyr1) + monnrd-mon1 + 1

      m3 = nmon12*(iyr1clm-iyr1) + mon1clm-mon1 + 1
      m4 = nmon12*(iyrnclm-iyr1) + monnclm-mon1 + 1

        do 190 j=j1,jn
          ji = j-j1+1
          ik = (j-1)/(nlat/nzon) + 1
          ik = min0(ik, nzon)

          print*, 'computing climatology for latitude ', alats(j)
          do 180 i=1,nlon

            if (wtfrac(i,j) .le. wtmin) then

              ovarmon(i,j) = amiss
              do 148 n=1,nmon12
                obsclim(i,j,n) = amiss
  148         continue

            else

              do 150 n=1,nmon12
                obsclim(i,j,n) = 0.0
  150         continue

              mn = 0
              rmin = 0.0
              rmax = 0.0
              do 155 n=m1,m2
                if (sst(i,ji,n) .gt. tmax) then
                  if (sst(i,ji,n) .gt. tmax+1.0e-2) then
                    print*, 'tmax exceeded ', i, j, n, sst(i,ji,n)
                    rmax = amax1(rmax,  (sst(i,ji,n)-tmax))
                    mn = mn + 1
                  endif
                  sst(i,ji,n) = tmax
                elseif (sst(i,ji,n) .lt. tmin) then
                  if (sst(i,ji,n) .lt. tmin-dt) then
c                    print*, 'tmin exceeded ', i, j, n, sst(i,ji,n)
                    if ((tmin-sst(i,ji,n)) .lt. 900.) then
                      rmin = amin1(rmin,  (sst(i,ji,n)-tmin))
                      mn = mn + 1
                    endif
                  endif
                  sst(i,ji,n) = tmin                  
                endif
  155         continue

              if (mn .gt. 0) then

                 mmnn = mmnn + 1
                 rrmin = amin1(rrmin, rmin)
                 rrmax = amax1(rrmax, rmax)
                 write(9,*) ' '
                 write(9,*) 
     &              ' WARNING -- observed value exceeds limits at'
                 write(9,*) mn, ' time points at latitude ', alats(j)
                 write(9,*) 'and longitude ', alons(i)
                 write(9,*) 'max error = ', rmax,
     &                      '  min error = ', rmin

              endif

c              compute climatology  

              woseas(ik) = woseas(ik) + sstwts(i,j)

              do 170 n=1,nmon12

                m = mod((mon1clm+n-2), nmon12) + 1
                obsclim(i,j,m)=0.0
                nn = 0

                do 160 k = m3+n-1, m4, nmon12 
                  nn = nn + 1
                  obsclim(i,j,m) = obsclim(i,j,m)+sst(i,ji,k)
  160           continue

                obsclim(i,j,m) = obsclim(i,j,m)/nn
                aobsclim(m,ik) = aobsclim(m,ik) + 
     &                               obsclim(i,j,m)*sstwts(i,j)

  170         continue

c            remove climatology to generate anomalies

              do 175 n=m1,m2
                mm = mod((mon1rd+n-m1-1), nmon12) + 1
                sst(i,ji,n) = sst(i,ji,n) - obsclim(i,j,mm)
  175         continue
              
              sum = 0.0
              avg1 = 0.0
              do 177 n=m3,m4
                sum = sum + sst(i,ji,n)*sst(i,ji,n)
                avg1 = avg1 + sst(i,ji,n)
  177         continue

              ovarmon(i,j) =
     &         (sum-(avg1*avg1/(m4-m3+1)))/float(m4-m3)
              if (ovarmon(i,j) .ge. varmin) then
                aovarmon(ik) = aovarmon(ik) + ovarmon(i,j)*sstwts(i,j)
                aovarmon(nzonp) = aovarmon(nzonp) + 
     &                              ovarmon(i,j)*sstwts(i,j)
                wovar(ik) = wovar(ik) + sstwts(i,j)
                wovar(nzonp) = wovar(nzonp) + sstwts(i,j)
              endif

            endif

  180     continue

  190   continue

        if (mmnn .gt. 0) then 
          write(9,*) ' '
          write(9,*) ' WARNING -- observed value exceeds limits at'
          write(9,*) mmnn, ' grid cells'
          write(9,*) 'max error = ', rrmax, '  min error = ', rrmin
        endif

      if (lagcalc .gt. 0) then

        do 250 j=j1,jn

          print*, 'computing lag correlations for latitude ', alats(j)
        ik = (j-1)/(nlat/nzon) + 1
        ik = min0(ik, nzon)
        ji = j-j1+1

        do 240 i=1,nlon

          if ((wtfrac(i,j) .gt. wtmin) .and. (ovarmon(i,j) .ge. varmin))
     &            then

            wocorr(ik) = wocorr(ik) + sstwts(i,j)
            wocorr(nzonp) = wocorr(nzonp) + sstwts(i,j)

            do 230 k=1,nlagm

              covar = 0.0
              avg1 = 0.0
              avg2 = 0.0
              
              nn = m4 - lagm(k)

              do 220 n=m3,nn
                mm = n+lagm(k)
                covar = covar + sst(i,ji,n)*sst(i,ji,mm)
                avg1 = avg1 + sst(i,ji,n)
                avg2 = avg2 + sst(i,ji,mm)
  220         continue

              if (k .eq. 1) then
c               calculate variance
                var = (covar - (avg1*avg2/(nn-m3+1)))/(nn-m3)
                ocorrel(k,ik) = ocorrel(k,ik) + sstwts(i,j)*var
                ocorrel(k,nzonp) = ocorrel(k,nzonp) + sstwts(i,j)*var
              elseif (var .gt. 0.) then
c               calculate correlations
                ocorrel(k,ik) = ocorrel(k,ik) + sstwts(i,j)*
     &            ((covar - (avg1*avg2/(nn-m3+1)))/(nn-m3))/var 
                ocorrel(k,nzonp) = ocorrel(k,nzonp) + sstwts(i,j)*
     &            ((covar - (avg1*avg2/(nn-m3+1)))/(nn-m3))/var 
              endif

  230       continue

          endif

  240   continue
  250   continue

      endif

c    solve for climatological mid-month values

      do 500 j=j1,jn

        print*, 
     &      'Computing climatological mid-month values for latitude ',
     &                alats(j)

        ik = (j-1)/(nlat/nzon) + 1
        ik = min0(ik, nzon)

        do 400 i=1,nlon

          if (wtfrac(i,j) .le. wtmin)then

            do 308 m=1,nmon12
              centmon(i,j,m) = amiss
  308       continue

            cmax(i,j) = amiss
            cmin(i,j) = amiss
            omax(i,j) = amiss
            omin(i,j) = amiss

          else

            do 310 m=1,nmon12
              vecin(m) = obsclim(i,j,m)
              vecout(m) = obsclim(i,j,m)
  310       continue

            call solvmid(alons(i), alats(j), nmon12, conv, dt, tmin, 
     &          tmax, bbmin, maxiter, ac, cc, vecin, vecout, ismc, jsmc)

            do 320 m=1,nmon12
              centmon(i,j,m) = vecout(m)
  320       continue

            ocentmax = -1.e20
            ocentmin = 1.e20
            centmax = -1.e20
            centmin = 1.e20

            do 390 m=1,nmon12

c             find max and min values and months of max and min sst

              if (obsclim(i,j,m) .gt. ocentmax) then
                ocentmax = obsclim(i,j,m)
                momax = m
              endif
              if (obsclim(i,j,m) .lt. ocentmin) then
                ocentmin = obsclim(i,j,m)
                momin = m
              endif  

              if (centmon(i,j,m) .gt. centmax) then
                centmax = centmon(i,j,m)
                mmax = m
              endif
              if (centmon(i,j,m) .lt. centmin) then
                centmin = centmon(i,j,m)
                mmin = m
              endif  

  390       continue

            cmax(i,j) = centmax + dcbias
            cmin(i,j) = centmin + dcbias
            omax(i,j) = ocentmax + dcbias
            omin(i,j) = ocentmin + dcbias

            wcseas(ik) = wcseas(ik) + sstwts(i,j)
            do 395 m=1,nmon12
              acentmon(m,ik) = acentmon(m,ik) + 
     &                      centmon(i,j,m)*sstwts(i,j)
  395       continue

c           count ocean grid cells
            isea = isea + 1

            diffmax = centmax - ocentmax
            diffmin = ocentmin - centmin

           if ((diffmax .gt. 700.0) .or. (diffmin .gt. 700.0)) then
             print*, 'exceeds 700 at lat = ', alats(j), 
     &              ' alon = ', alons(i)
             print*, 'kk = ', kk, ' obs = ', (obsclim(i,j,m), m=1,12) 
             print*, 'kk = ', kk, ' mid = ', (centmon(i,j,m), m=1,12)
           endif


            ii = diffmax/histo
            if (ii .gt. 20) ii=20
            if (ii .lt. -20) ii=-20
            jj = 0.1*diffmax/histo
            if (jj .gt. 10) jj=10
            if (jj .lt. -10) jj=-10

            if (ii .gt. 0) then
              do 84 m=1,ii
                icount(ik,m) = icount(ik,m) + 1
   84         continue
            elseif (ii .lt. 0) then
              do 85 m=-1,ii,-1
                icount(ik,m) = icount(ik,m) + 1
   85         continue
            else
              icount(ik,0) = icount(ik,0) + 1
            endif
            if (jj .gt. 0) then
              do 86 m=1,jj
                jcount(ik,m) = jcount(ik,m) + 1
   86         continue
            elseif (jj .lt. 0) then
              do 87 m=-1,jj,-1
                jcount(ik,m) = jcount(ik,m) + 1
   87         continue
            else
              jcount(ik,0) = jcount(ik,0) + 1             
            endif
               
            kk = diffmin/histo
            if (kk .gt. 20) kk=20
            if (kk .lt. -20) kk=-20
            mm = 0.1*diffmin/histo
            if (mm .gt. 10) mm=10
            if (mm .lt. -10) mm=-10

            if (kk .gt. 0) then
              do 94 m=1,kk
                kcount(ik,m) = kcount(ik,m) + 1
   94         continue
            elseif (kk .lt. 0) then
              do 95 m=-1,kk,-1
                kcount(ik,m) = kcount(ik,m) + 1
   95         continue
            else
              kcount(ik,0) = kcount(ik,0) + 1
            endif
            if (mm .gt. 0) then
              do 96 m=1,mm
                mcount(ik,m) = mcount(ik,m) + 1
   96         continue
            elseif (mm .lt. 0) then
              do 97 m=-1,mm,-1
                mcount(ik,m) = mcount(ik,m) + 1
   97         continue
            else
              mcount(ik,0) = mcount(ik,0) + 1
            endif

          endif

  400   continue
  500 continue

      return
      end


      subroutine calcfull(j1, jn, nlon, nlat, mlat, nmon,
     &    iyr1rd, mon1rd, iyr1, mon1, iyrn, monn, iyrnrd, monnrd,
     &    mon1clm, monnclm, iyr1clm, iyrnclm, ismf, jsmf,
     &    nzon, nzonp, nmon12, lagcalc, nlagm,
     &    maxiter, wtmin, amiss, tmin, tmax,  
     &    dt, conv, bbmin, varmin, lagm, sstwts, wtfrac, 
     &    wcseas, cvarmon, acvarmon, correl, ccorrel,  
     &    obsclim, vecin, vecout, centmon, sst, alons, alats, a, c, 
     &    ovarmon, wcvar, wccorr)

      implicit none
      integer nlon, nlat, nmon, iyr1rd, mon1rd, iyr1, mon1, iyrn, monn,
     &    mon1clm, monnclm, iyr1clm, iyrnclm, ismf, jsmf, 
     &   iyrnrd, monnrd, nzon, nzonp, nmon12, lagcalc, nlagm,
     &    maxiter, mlat, ji, nprior, nafter, j1, jn

      real wtmin, amiss, tmin, tmax, dt, conv, bbmin, varmin

      integer lagm(nlagm)

      real sstwts(nlon,nlat), wtfrac(nlon,nlat), wcseas(nzon),
     &    cvarmon(nlon,nlat), acvarmon(nzonp), 
     &    correl(nmon12), ccorrel(nlagm,nzonp),
     &    obsclim(nlon,nlat,nmon12), vecin(nmon), vecout(nmon), 
     &    centmon(nlon,nlat,nmon12), sst(nlon,mlat,nmon), alons(nlon), 
     &    alats(nlat), a(nmon), c(nmon), ovarmon(nlon,nlat), 
     &    wcvar(nzonp), wccorr(nzonp)

      integer i, j, k, m, mm, ik, nn, n, m3, m4

      real tt, s, var, cc

      double precision avg1, avg2, covar, sum

      m3 = nmon12*(iyr1clm-iyr1) + mon1clm-mon1 + 1
      m4 = nmon12*(iyrnclm-iyr1) + monnclm-mon1 + 1

      do 900 j=j1,jn

        print*, 'Calculating mid-month values for latitude ', alats(j)
        ik = (j-1)/(nlat/nzon) + 1
        ik = min0(ik, nzon)
        ji = j-j1+1

        do 800 i=1,nlon 

          if (wtfrac(i,j) .le. wtmin) then

            cvarmon(i,j) = amiss

          else


c           fill in some months prior to the 1st observed month and
c              after the last observed month

            nprior = nmon12*(iyr1rd-iyr1) + mon1rd - mon1 

            if (nprior .gt. 0) then
              do 708 m=1,nprior

                if ((nprior+1-m) .gt. nmon12) then
                  cc = 0.0
                else
                  cc = correl(nprior+1-m)
                endif
                mm = mod((mon1-2+m), nmon12) + 1  
                tt = sst(i,ji,nprior+1)*cc
                if ((obsclim(i,j,mm)+tt) .gt. tmax) 
     &                tt = tmax-obsclim(i,j,mm)
                if ((obsclim(i,j,mm)+tt) .lt. tmin) 
     &                tt = tmin-obsclim(i,j,mm)
                sst(i,ji,m) = tt

  708         continue
            endif

            nafter = nmon12*(iyrn-iyrnrd) + monn - monnrd 

            if (nafter .gt. 0) then
              do 709 m=1,nafter

                if (m .gt. nmon12) then
                  cc = 0.0
                else
                  cc = correl(m)
                endif 
                mm = mod((monnrd-1+m), nmon12) + 1
                tt = sst(i,ji,nmon-nafter)*cc
                if ((obsclim(i,j,mm)+tt) .gt. tmax) 
     &                tt = tmax-obsclim(i,j,mm)
                if ((obsclim(i,j,mm)+tt) .lt. tmin) 
     &                tt = tmin-obsclim(i,j,mm)

                sst(i,ji,nmon-nafter+m) = tt

  709         continue

            endif
          
c            copy data into vecin, vecout
            do 710 n=1,nmon
              mm = mod((mon1+n-2), nmon12) + 1
              vecin(n) = sst(i,ji,n) + obsclim(i,j,mm)
              vecout(n) = vecin(n)
  710       continue


            call solvmid(alons(i), alats(j), nmon, conv, dt, tmin, tmax,
     &            bbmin, maxiter, a, c, vecin, vecout, ismf, jsmf)

            sum = 0.0
            avg1 = 0.0

            do 718 n=m3, m4
              mm = mod((mon1clm+n-m3-1), nmon12) + 1
              s = amax1(tmin, amin1(tmax, vecout(n))) - 
     &            amax1(tmin, amin1(tmax, centmon(i,j,mm)))
              sum = sum + s*s
              avg1 = avg1 + s
  718       continue

            cvarmon(i,j) = (sum-avg1*avg1/(m4-m3+1))/(m4-m3)

            do 720 n=1,nmon
              mm = mod((mon1+n-2), nmon12) + 1
              sst(i,ji,n) = vecout(n) - centmon(i,j,mm)
  720       continue


            if (ovarmon(i,j) .ge. varmin) then

              wcvar(ik) = wcvar(ik) + sstwts(i,j)
              acvarmon(ik) = acvarmon(ik) + cvarmon(i,j)*sstwts(i,j)
              wcvar(nzonp) = wcvar(nzonp) + sstwts(i,j)
              acvarmon(nzonp) = acvarmon(nzonp) + 
     &                                     cvarmon(i,j)*sstwts(i,j)

            endif

          endif

  800   continue
  900 continue

      if (lagcalc .gt. 0) then

        do 550 j=j1,jn
          print*,
     &      'calculating mid-month lag correlations for latitude ', 
     &      alats(j)

        ji = j-j1+1
        ik = (j-1)/(nlat/nzon) + 1
        ik = min0(ik, nzon)
        do 540 i=1,nlon

          if ((wtfrac(i,j) .gt. wtmin) .and. (ovarmon(i,j) .ge. varmin))
     &          then

            wccorr(ik) = wccorr(ik) + sstwts(i,j)
            wccorr(nzonp) = wccorr(nzonp) + sstwts(i,j)

            do 530 k=1,nlagm

              covar = 0.0
              avg1 = 0.0
              avg2 = 0.0
              
              nn = m4 - lagm(k)

              do 520 n=m3,nn
                mm = n+lagm(k)
                covar = covar + sst(i,ji,n)*sst(i,ji,mm)
                avg1 = avg1 + sst(i,ji,n)
                avg2 = avg2 + sst(i,ji,mm)
  520         continue

              if (k .eq. 1) then
c               calculate variance
                var = (covar - (avg1*avg2/(nn-m3+1)))/(nn-m3)
                ccorrel(k,ik) = ccorrel(k,ik) + sstwts(i,j)*var
                ccorrel(k,nzonp) = ccorrel(k,nzonp) + sstwts(i,j)*var
              elseif (var .gt. 0.) then
c               calculate correlations
                ccorrel(k,ik) = ccorrel(k,ik) + sstwts(i,j)*
     &            ((covar - (avg1*avg2/(nn-m3+1)))/(nn-m3))/var 
                ccorrel(k,nzonp) = ccorrel(k,nzonp) + sstwts(i,j)*
     &            ((covar - (avg1*avg2/(nn-m3+1)))/(nn-m3))/var 
              endif

  530       continue

          endif

  540   continue
  550   continue

      endif

      return
      end




      subroutine solvmid(alon, alat, nmon, conv, dt, tmin, tmax, 
     &     bbmin, maxiter, a, c, obsmean, ss, icnt, jcnt)

      implicit none
      integer nmont, nmon12

      parameter (nmont=147*12, nmon12=12) ! PJD Oceanonly 1870-2014

      integer nmon, icnt, jcnt, maxiter
      real conv, tmin, tmax, dt, bbmin, alon, alat
      real obsmean(nmon), a(nmon), c(nmon), ss(nmon)

      integer i, n, imethod, n1, n2, nn, jj, i1, i2, i3, nnn, jend,
     &      kk, j, k, kkk, nm, np
      integer jbeg(nmont)
      real relax, residmax, resid, dxm, dxp, s1, s2, addmax, addmin
      real r(nmont), avg(nmont), aa(nmont), bb(nmont), cc(nmont), 
     &     add(nmont)
      double precision s(nmont), sum

      imethod = 1
c ???  check following value
      relax = 1.0

      if (nmon .gt. nmont) then
        print*, 'error-- nmont not declared large enough in '
        print*, 'subroutine solvmid'
        stop
      endif

c    check for occurance where obs monthly means are consecutively
c      at upper and lower limits. If so, smooth data, being careful
c      to preserve annual mean.

      do 48 n=1,nmon
        add(n) = 0.0
   48 continue

      n2 = nmon
      do 50 n=1,nmon

        n1 = n2
        n2 = n

        if ((obsmean(n2)-obsmean(n1)) .gt. (tmax-tmin-2.*dt)) then
          add(n1) = add(n1) + 
     &           (0.5*(obsmean(n2)-obsmean(n1)-tmax+tmin) + dt)/c(n1)
          add(n2) = add(n2) - 
     &           (0.5*(obsmean(n2)-obsmean(n1)-tmax+tmin) + dt)/a(n2)
        elseif ((obsmean(n1)-obsmean(n2)) .gt. (tmax-tmin-2.*dt)) then
          add(n1) = add(n1) - 
     &           (0.5*(obsmean(n1)-obsmean(n2)-tmax+tmin) + dt)/c(n1)
          add(n2) = add(n2) + 
     &           (0.5*(obsmean(n1)-obsmean(n2)-tmax+tmin) + dt)/a(n2)
        endif         

   50 continue
  
      nn = 0
      addmax = 0.0
      addmin = 0.0
      do 51 n=1,nmon
         
        if (add(n) .ne. 0.0) then

           if (jcnt .lt. 5000) then
              print*, 'jump from extreme to extreme at time ', n
              write(9,*) 'jump from extreme to extreme at time ', n
           endif

           if (jcnt .eq. 5000) then
              print*, ' **********************************'
              print*, ' No more jump messages will be written'
              print*, ' **********************************'
           endif

           if (jcnt .lt. 5) then
          print*,
     &     'monthly means go from one limit to the other in 1 month'
          print*, 'alat = ', alat, ' alon = ', alon, ' n = ', n,
     &         ' add = ', add(n)
             n1 = mod((n+nmon-2), nmon) + 1
             n2 = mod(n, nmon) + 1

          print*, 'observed mean for 3 cells: ', obsmean(n1), 
     &          obsmean(n), obsmean(n2)
           endif
          if (jcnt .le. 2) 
     &        print*, (obsmean(n1), n1=1,nmon)

          addmax = amax1(addmax, add(n))
          addmin = amin1(addmin, add(n))
          obsmean(n) = obsmean(n) + add(n)
          nn = nn + 1
        endif
   51 continue

      icnt = icnt + nn
      if (nn .gt. 0) then
        jcnt = jcnt + 1
c        if (nmon .eq. nmon12) then
c          print*, 'Climatology: '
c          write(9,*) 'Climatology: '
c        endif
        if (jcnt .le. 1000) then
          print*,  nn, 
     &      ' monthly values smoothed at lat,lon', alat, alon
          print*, 'max added = ', addmax, 
     &             '  max subtracted = ', addmin
          write(9,*)  nn, 
     &      ' monthly values smoothed at lat,lon', alat, alon
          write(9,*) 'max added = ', addmax, 
     &             '  max subtracted = ', addmin
        endif
        if (jcnt .eq. 50) then
          print*, ' '
          if (nmon .eq. nmon12) then
            print*, 
     &      'No more warnings will be printed concerning smoothing '//
     &       'of climatological data'
          else
            print*, 
     &      'No more warnings will be printed concerning smoothing '//
     &       'of monthly data'
          endif
          print*, ' '
          write(9,*) ' '
          if (nmon .eq. nmon12) then
            write(9,*)
     &      'No more warnings will be printed concerning smoothing '//
     &       'of climatological data'
          else
            write(9,*)
     &      'No more warnings will be printed concerning smoothing '//
     &       'of monthly data'
          endif
          write(9,*) ' '
        endif
      endif

c    check if all are le tmin or all are ge tmax

      if (obsmean(1) .le. (tmin+0.01*dt)) then

        do 80 i=2,nmon
          if (obsmean(i) .gt. (tmin+0.01*dt)) go to 99
   80   continue
        do 85 i=1,nmon
          ss(i) = tmin
   85   continue
c        if (nmon .eq. nmon12) print*, 'Climatology: '
c        print*, 'all values were at minimum at this grid cell:'
c        print*, 'latitude = ', alat, ' longitude = ', alon 
        return

      elseif (obsmean(1) .ge. (tmax-0.01*dt)) then

        do 90 i=2,nmon
          if (obsmean(i) .lt. (tmax-0.01*dt)) go to 99
   90   continue
        do 95 i=1,nmon
          ss(i) = tmax
   95   continue
c        if (nmon .eq. nmon12) print*, 'Climatology: '
c        print*, 'all values were at maximum at this grid cell:'
c        print*, 'latitude = ', alat, ' longitude = ', alon 
        return

      endif

   99 jj = 0
      do 100 i=1,nmon

        i1 = i
        i2 = mod(i,nmon) + 1
        i3 = mod((i+1), nmon) + 1
     
        if (((obsmean(i1) .le. tmin+0.01*dt) .and. 
     &       (obsmean(i2) .le. tmin+0.01*dt) .and.
     &       (obsmean(i3) .gt. tmin+0.01*dt)) .or.
     &      ((obsmean(i1) .ge. tmax-0.01*dt) .and. 
     &       (obsmean(i2) .ge. tmax-0.01*dt) .and.
     &       (obsmean(i3) .lt. tmax-0.01*dt))) then
          jj = jj + 1
          jbeg(jj) = i2
        endif

  100 continue

      if (jj .eq. 0) then

c       simple cyclic treatment

c         latest approximation of means (given mid-month values)

         nnn = 0
  105    nnn = nnn + 1
         sum = 0.0
         residmax = 0.0

         do 110 n = 1, nmon

           nm = mod((n+nmon-2), nmon) + 1
           np = mod(n, nmon) + 1

           bb(n) = 0.0
           avg(n) = 0.0

           if (nnn .lt. imethod) then

             call approx(tmin, tmax, a(n), c(n), ss(nm), ss(n),
     &                ss(np), aa(n), bb(n), cc(n), avg(n))

           else

             call numer(conv, tmin, tmax, bbmin, a(n), c(n), ss(nm), 
     &                ss(n), ss(np), aa(n), bb(n), cc(n), avg(n))

           endif


           r(n) = obsmean(n) - avg(n) 
           sum = sum + r(n)**2
           residmax = amax1(residmax, abs(r(n)))

110      continue

         resid = dsqrt(sum)/nmon


         if (residmax .gt. conv) then

           if (nnn .gt. maxiter*0.5) then
             print*, 'iteration = ', nnn, ' residual = ', resid,
     &          ' maximum residual = ', residmax
           endif

           if (nnn .gt. maxiter*0.9) then
  
             print*, ' '
             print*, 'latitude = ', alat, ' longitude = ', alon 

             do 1234 n=1,nmon
               write(*,'(8(1pe10.2))') obsmean(n), avg(n), r(n), 
     &              s(n), ss(n), aa(n), bb(n), cc(n)
 1234        continue

           endif

           if (nnn .gt. maxiter) then

             print*, 'latitude = ', alat, ' longitude = ', alon 
             print*, 'does not converge'
             write(9,*) 'latitude = ', alat, ' longitude = ', alon 
             write(9,*) 'does not converge'
c             call exit(1)

           else

c            solve for new estimate of mid-month values

             call cyclic(alon, alat, aa, bb, cc, cc(nmon), aa(1), r, s,
     &                nmon)

             do 120 n=1,nmon
               ss(n) = ss(n) + relax*s(n)
  120        continue

c            if ss exceeds tmax or tmin, then it should exceed it no
c                more than absolutely necessary:

             do 130 n=1,nmon

               nm = mod((n+nmon-2), nmon) + 1
               np = mod(n, nmon) + 1

               if (ss(n) .gt. tmax) then

                 if (ss(nm) .le. tmax) then
                   dxm = (ss(n)-tmax)/((ss(n)-ss(nm))*a(n))
                 else 
                   dxm = 0.0
                 endif

                 if (ss(np) .le. tmax) then 
                   dxp = (ss(n)-tmax)/((ss(n)-ss(np))*c(n))
                 else
                   dxp = 0.0
                 endif

                 if ((dxm .gt. 0.5) .and. (dxp .gt. 0.5)) then
                   s1 = tmax + (tmax-ss(nm))*a(n)/(2.-a(n))
                   s2 = tmax + (tmax-ss(np))*c(n)/(2.-c(n))
                   ss(n) = amin1(s1, s2)
                 elseif ((dxm .eq. 0.0) .and. (dxp .eq. 0.0)) then
                   ss(n) = tmax
                 elseif ((dxp .eq. 0.0) .and. (dxm .gt. 0.5)) then
                   ss(n) = tmax + (tmax-ss(nm))*a(n)/(2.-a(n))
                 elseif ((dxm .eq. 0.0) .and. (dxp .gt. 0.5)) then
                   ss(n) = tmax + (tmax-ss(np))*c(n)/(2.-c(n))
                 endif

               elseif (ss(n) .lt. tmin) then

                 if (ss(nm) .ge. tmin) then
                   dxm = (ss(n)-tmin)/((ss(n)-ss(nm))*a(n))
                 else 
                   dxm = 0.0
                 endif

                 if (ss(np) .ge. tmin) then 
                   dxp = (ss(n)-tmin)/((ss(n)-ss(np))*c(n))
                 else
                   dxp = 0.0
                 endif

                 if ((dxm .gt. 0.5) .and. (dxp .gt. 0.5)) then
                   s1 = tmin + (tmin-ss(nm))*a(n)/(2.-a(n))
                   s2 = tmin + (tmin-ss(np))*c(n)/(2.-c(n))
                   ss(n) = amin1(s1, s2)
                 elseif ((dxm .eq. 0.0) .and. (dxp .eq. 0.0)) then
                   ss(n) = tmin
                 elseif ((dxp .eq. 0.0) .and. (dxm .gt. 0.5)) then
                   ss(n) = tmin + (tmin-ss(nm))*a(n)/(2.-a(n))
                 elseif ((dxm .eq. 0.0) .and. (dxp .gt. 0.5)) then
                   ss(n) = tmin + (tmin-ss(np))*c(n)/(2.-c(n))
                 endif
               endif

  130        continue

             go to 105

           endif

         endif

      else

c        treat independent segments

         do 300 j=1,jj
           jend = jbeg(j)
  150      jend = jend + 1
           i1 = mod((jend-1), nmon) + 1
           i2 = mod(jend, nmon) + 1
           if (((obsmean(i1) .le. tmin+0.01*dt) .and.
     &          (obsmean(i2) .le. tmin+0.01*dt)) .or.
     &         ((obsmean(i1) .ge. tmax-0.01*dt) .and.
     &          (obsmean(i2) .ge. tmax-0.01*dt))) then

c            calculate values for interval jbeg(j) to jend
c
c            latest approximation of means (given mid-month values)

             nnn = 0
  205        nnn = nnn + 1
             kk = jend - jbeg(j) + 1

             n = jbeg(j)
             avg(1) = obsmean(n)
             r(1) = 0.0
             n = mod((jend-1), nmon) + 1
             avg(kk) = obsmean(n)
             r(kk) = 0.0

             sum = 0.0
             residmax = 0.0

             do 210 k = 2, kk-1

               nm = mod((k+jbeg(j)-3), nmon) + 1
               n  = mod((k+jbeg(j)-2), nmon) + 1
               np = mod((k+jbeg(j)-1), nmon) + 1

               bb(k) = 0.0
               avg(k) = 0.0

               if (nnn .lt. imethod) then

                 call approx(tmin, tmax, a(n), c(n), ss(nm), ss(n),
     &                ss(np), aa(k), bb(k), cc(k), avg(k))

               else

                 call numer(conv, tmin, tmax, bbmin, a(n), c(n), ss(nm),
     &                ss(n), ss(np), aa(k), bb(k), cc(k), avg(k))

               endif

               r(k) = obsmean(n) - avg(k) 
               sum = sum + r(k)**2
               residmax = amax1(residmax, abs(r(k)))

210          continue

             resid = dsqrt(sum)/(kk-2)

             if (residmax .gt. conv) then

               if (nnn .gt. maxiter*0.9) then
                  print*, 'iter = ', nnn, ' kk = ', kk, ' residual = ',
     &               resid, ' maximum residual = ', residmax
c                  print*, ss(nm), ss(n), ss(np)
               endif

               if (nnn .gt. maxiter*0.9) then
  
                 print*, ' '
                 print*, 'latitude = ', alat, ' longitude = ', alon 

                 do 2234 k=1,kk
                   n  = mod((k+jbeg(j)-2), nmon) + 1
                   write(*,'(8(1pe10.2))') obsmean(n), avg(k), r(k), 
     &                  s(k), ss(n), aa(k), bb(k), cc(k)
 2234            continue

               endif

               if (nnn .gt. maxiter) then

                 print*, 'latitude = ', alat, ' longitude = ', alon 
                 print*, 'does not converge'
                 write(9,*) 'latitude = ', alat, ' longitude = ', alon 
                 write(9,*) 'does not converge'
c                 call exit(1)

               else

c                solve for new estimate of mid-month values

                 kkk = kk - 2
                 call tridag(alon, alat, aa(2), bb(2), cc(2), r(2),
     &                 s(2), kkk)

                 do 220 k=2,kk-1
                   n  = mod((k+jbeg(j)-2), nmon) + 1
                   ss(n) = ss(n) + relax*s(k)
  220            continue

c               if ss exceeds tmax or tmin, then it should exceed it no
c                    more than absolutely necessary:

                 n  = mod((jbeg(j)-1), nmon) + 1
                 np  = mod(jbeg(j), nmon) + 1

                 if (obsmean(n) .ge. (tmax-0.01*dt)) then 
                   ss(n) = 
     &                amax1(tmax, (tmax + (tmax-ss(np))*c(n)/(2.-c(n))))
                 else
                   ss(n) = 
     &                amin1(tmin, (tmin + (tmin-ss(np))*c(n)/(2.-c(n))))
                 endif

                 nm  = mod((jend+nmon-2), nmon) + 1
                 n  = mod((jend-1), nmon) + 1

                 if (obsmean(n) .ge. (tmax-0.01*dt)) then 
                   ss(n) = 
     &                amax1(tmax, (tmax + (tmax-ss(nm))*a(n)/(2.-a(n))))
                 else
                   ss(n) = 
     &                amin1(tmin, (tmin + (tmin-ss(nm))*a(n)/(2.-a(n))))
                 endif

                 do 230 k=2,kk-1

                   nm = mod((k+jbeg(j)+nmon-3), nmon) + 1
                   n  = mod((k+jbeg(j)-2), nmon) + 1
                   np = mod((k+jbeg(j)-1), nmon) + 1

                   if (ss(n) .gt. tmax) then

                     if (ss(nm) .le. tmax) then
                       dxm = (ss(n)-tmax)/((ss(n)-ss(nm))*a(n))
                     else 
                       dxm = 0.0
                     endif

                     if (ss(np) .le. tmax) then 
                       dxp = (ss(n)-tmax)/((ss(n)-ss(np))*c(n))
                     else
                       dxp = 0.0
                     endif

                     if ((dxm .gt. 0.5) .and. (dxp .gt. 0.5)) then
                       s1 = tmax + (tmax-ss(nm))*a(n)/(2.-a(n))
                       s2 = tmax + (tmax-ss(np))*c(n)/(2.-c(n))
                       ss(n) = amin1(s1, s2)
                     elseif ((dxm .eq. 0.0) .and. (dxp .eq. 0.0)) then
                       ss(n) = tmax
                     elseif ((dxp .eq. 0.0) .and. (dxm .gt. 0.5)) then
                       ss(n) = tmax + (tmax-ss(nm))*a(n)/(2.-a(n))
                     elseif ((dxm .eq. 0.0) .and. (dxp .gt. 0.5)) then
                       ss(n) = tmax + (tmax-ss(np))*c(n)/(2.-c(n))
                     endif

                   elseif (ss(n) .lt. tmin) then

                     if (ss(nm) .ge. tmin) then
                       dxm = (ss(n)-tmin)/((ss(n)-ss(nm))*a(n))
                     else 
                       dxm = 0.0
                     endif

                     if (ss(np) .ge. tmin) then 
                       dxp = (ss(n)-tmin)/((ss(n)-ss(np))*c(n))
                     else
                       dxp = 0.0
                     endif

                     if ((dxm .gt. 0.5) .and. (dxp .gt. 0.5)) then
                       s1 = tmin + (tmin-ss(nm))*a(n)/(2.-a(n))
                       s2 = tmin + (tmin-ss(np))*c(n)/(2.-c(n))
                       ss(n) = amax1(s1, s2)
                     elseif ((dxm .eq. 0.0) .and. (dxp .eq. 0.0)) then
                       ss(n) = tmin
                     elseif ((dxp .eq. 0.0) .and. (dxm .gt. 0.5)) then
                       ss(n) = tmin + (tmin-ss(nm))*a(n)/(2.-a(n))
                     elseif ((dxm .eq. 0.0) .and. (dxp .gt. 0.5)) then
                       ss(n) = tmin + (tmin-ss(np))*c(n)/(2.-c(n))
                     endif
                   endif

  230            continue

                 go to 205

               endif

             endif

             go to 300
           else
             go to 150
           endif

           if ((ss(kk) .gt. 700.0) .or. (ss(kk) .lt. -700.0)) then
             print*, 'exceeds 700 at lat = ', alat, 
     &              ' alon = ', alon
             print*, 'kk = ', kk, ' obs = ', obsmean(kk-1), obsmean(kk),
     &               obsmean(kk+1)
             print*, 'kk = ', kk, '  ss = ', ss(kk-1), ss(kk), 
     &              ss(kk+1)
           endif

  300    continue

c        fill in values where consecutive means are outside limits

         do 250 i=1,nmon
           i1 = mod((i-2+nmon), nmon) + 1
           i2 = mod((i-1), nmon) + 1
           i3 = mod(i, nmon) + 1
           if ((obsmean(i1) .le. (tmin+0.01*dt)) .and.
     &         (obsmean(i2) .le. (tmin+0.01*dt)) .and. 
     &         (obsmean(i3) .le. (tmin+0.01*dt))) then
             ss(i2) = tmin
           elseif ((obsmean(i1) .ge. (tmax-0.01*dt)) .and.
     &             (obsmean(i2) .ge. (tmax-0.01*dt)) .and.
     &             (obsmean(i3) .ge. (tmax-0.01*dt))) then
             ss(i2) = tmax
           endif
 
  250    continue

      endif


           return

      end


      subroutine numer(conv, tmin, tmax, bbmin, a, c, ssm, ss, ssp, aa,
     &                 bb, cc, avg)
c *********************************************************************

      implicit none
 
      real conv, tmin, tmax, bbmin, a, c, ssm, ss, ssp, aa, bb, cc, avg
      real ssmm, ssmp, sssm, sssp, sspm, sspp, r
      real amean

      avg = amean(tmin,tmax,a,c,ssm,ss,ssp)

      ssmm = ssm - conv
      ssmp = ssm + conv
      sssm = ss  - conv
      sssp = ss  + conv
      sspm = ssp - conv
      sspp = ssp + conv

      aa = (amean(tmin,tmax,a,c,ssmp,ss,ssp) - 
     &      amean(tmin,tmax,a,c,ssmm,ss,ssp)) / (2.*conv)

      bb = (amean(tmin,tmax,a,c,ssm,sssp,ssp) - 
     &      amean(tmin,tmax,a,c,ssm,sssm,ssp)) / (2.*conv)

      cc = (amean(tmin,tmax,a,c,ssm,ss,sspp) - 
     &      amean(tmin,tmax,a,c,ssm,ss,sspm)) / (2.*conv)

      aa = amin1(aa, bb)
      cc = amin1(cc, bb)

      if (bb .lt. bbmin) then
        bb = bbmin
        r = 0.2*bbmin
        aa = amax1(r, aa)
        cc = amax1(r, cc)
      endif

      return
      end




      function amean(tmin, tmax, a, c, ssm, ss, ssp)
c *********************************************************************

      implicit none

      real amean
      real tmin, tmax, a, c, ssm, ss, ssp
      real dx, dy, avg

      avg = 0.0

      if (ss .le.  tmin) then
      
        if (ssm .le. tmin) then

          avg = avg + tmin*0.5
                   
        elseif (ssm .ge. tmax) then

          dx = (ss-tmin)/((ss-ssm)*a)
          dy = (ss-tmax)/((ss-ssm)*a)
          if (dx .ge. 0.5) then
            avg = avg + tmin*0.5
          elseif (dy .le. 0.5) then
            avg = avg + tmin*dx + tmax*(.5-dy) + (dy-dx)*.5*(tmin+tmax)
          else
            avg = avg + 
     &              tmin*dx + (0.5-dx)*0.5*(tmin + ss - 0.5*a*(ss-ssm))
          endif

        else

          dx = (ss-tmin)/((ss-ssm)*a)
          if (dx .ge. 0.5) then
            avg = avg + tmin*0.5
          else
            avg = avg +
     &              tmin*dx + (0.5-dx)*0.5*(tmin + ss - 0.5*a*(ss-ssm))
          endif

        endif

        if (ssp .le. tmin) then

          avg = avg + tmin*0.5

        elseif (ssp .ge. tmax) then

          dx = (ss-tmin)/((ss-ssp)*c)
          dy = (ss-tmax)/((ss-ssp)*c)
          if (dx .ge. 0.5) then
            avg = avg + tmin*0.5
          elseif (dy .le. 0.5) then
            avg = avg + tmin*dx + tmax*(.5-dy) + (dy-dx)*.5*(tmin+tmax)
          else
            avg = avg + 
     &              tmin*dx + (0.5-dx)*0.5*(tmin + ss - 0.5*c*(ss-ssp))
          endif

        else

          dx = (ss-tmin)/((ss-ssp)*c)
          if (dx .ge. 0.5) then
            avg = avg + tmin*0.5
          else
            avg = avg + 
     &              tmin*dx + (0.5-dx)*0.5*(tmin + ss - 0.5*c*(ss-ssp))
          endif

        endif

      elseif (ss .ge. tmax) then

        if (ssm .ge. tmax) then

          avg = avg + tmax*0.5

        elseif (ssm .le. tmin) then

          dx = (ss-tmax)/((ss-ssm)*a)
          dy = (ss-tmin)/((ss-ssm)*a)
          if (dx .ge. 0.5) then
            avg = avg + tmax*0.5
          elseif (dy .le. 0.5) then
            avg = avg + tmax*dx + tmin*(.5-dy) + (dy-dx)*.5*(tmin+tmax)
          else
            avg = avg + 
     &              tmax*dx + (0.5-dx)*0.5*(tmax + ss - 0.5*a*(ss-ssm))
          endif

        else

          dx = (ss-tmax)/((ss-ssm)*a)
          if (dx .ge. 0.5) then
            avg = avg + tmax*0.5
          else
            avg = avg + 
     &        tmax*dx + (0.5-dx)*0.5*(tmax + ss - 0.5*a*(ss-ssm))
          endif

        endif
        
        if (ssp .ge. tmax) then

          avg = avg + tmax*0.5

        elseif (ssp .le. tmin) then

          dx = (ss-tmax)/((ss-ssp)*c)
          dy = (ss-tmin)/((ss-ssp)*c)
          if (dx .ge. 0.5) then
            avg = avg + tmax*0.5
          elseif (dy .le. 0.5) then
            avg = avg + tmax*dx + tmin*(.5-dy) + (dy-dx)*.5*(tmin+tmax)
          else
            avg = avg + 
     &              tmax*dx + (0.5-dx)*0.5*(tmax + ss - 0.5*c*(ss-ssp))
          endif

        else

          dx = (ss-tmax)/((ss-ssp)*c)
          if (dx .ge. 0.5) then
            avg = avg + tmax*0.5
          else
            avg = avg + 
     &              tmax*dx + (0.5-dx)*0.5*(tmax + ss - 0.5*c*(ss-ssp))
          endif

        endif

      else

        if (ssm .le. tmin) then

          dx = (ss-tmin)/((ss-ssm)*a)
          if (dx .ge. 0.5) then
            avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssm)*a)
          else
            avg = avg + tmin*(.5-dx) + dx*0.5*(tmin + ss)
          endif

        elseif (ssm .ge. tmax) then

          dx = (ss-tmax)/((ss-ssm)*a)
          if (dx .ge. 0.5) then
            avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssm)*a)
          else
            avg = avg + tmax*(.5-dx) + dx*0.5*(tmax + ss)
          endif

        else

          avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssm)*a)

        endif
        
        if (ssp .le. tmin) then

          dx = (ss-tmin)/((ss-ssp)*c)
          if (dx .ge. 0.5) then
            avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssp)*c)
          else
            avg = avg + tmin*(.5-dx) + dx*0.5*(tmin + ss)
          endif

        elseif (ssp .ge. tmax) then

          dx = (ss-tmax)/((ss-ssp)*c)
          if (dx .ge. 0.5) then
            avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssp)*c)
          else
            avg = avg + tmax*(.5-dx) + dx*0.5*(tmax + ss)
          endif

        else

          avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssp)*c)

        endif

      endif

      amean = avg

      return
      end



      subroutine approx(tmin, tmax, a, c, ssm, ss, ssp, aa, bb, cc, avg)

c *********************************************************************

      implicit none

      real tmin, tmax, a, c, ssm, ss, ssp, aa, bb, cc, avg
      real dx, dy



      if (ss .le.  tmin) then
      
        if (ssm .le. tmin) then

          avg = avg + tmin*0.5
          aa = a/32.
          bb = bb + 0.125 - a/32.
                   
        elseif (ssm .ge. tmax) then

          dx = (ss-tmin)/((ss-ssm)*a)
          dy = (ss-tmax)/((ss-ssm)*a)
          if (dx .ge. 0.5) then
            avg = avg + tmin*0.5
            aa = a/32.
            bb = bb + 0.125 - a/32.
          elseif (dy .le. 0.5) then
            avg = avg + tmin*dx + tmax*(.5-dy) + (dy-dx)*.5*(tmin+tmax)
            aa = a/16.
            bb = bb + 0.25 - a/16. 
          else
            avg = avg + 
     &              tmin*dx + (0.5-dx)*0.5*(tmin + ss - 0.5*a*(ss-ssm))
            aa = a/16.
            bb = bb + 0.25 - a/16.
          endif

        else

          dx = (ss-tmin)/((ss-ssm)*a)
          if (dx .ge. 0.5) then
            avg = avg + tmin*0.5
            aa = a/32.
            bb = bb + 0.125 - a/32.
          else
            avg = avg +
     &              tmin*dx + (0.5-dx)*0.5*(tmin + ss - 0.5*a*(ss-ssm))
            aa = a/16.
            bb = bb + 0.25 - a/16.
          endif

        endif

        if (ssp .le. tmin) then

          avg = avg + tmin*0.5
          cc = c/32.
          bb = bb + 0.125 - c/32.

        elseif (ssp .ge. tmax) then

          dx = (ss-tmin)/((ss-ssp)*c)
          dy = (ss-tmax)/((ss-ssp)*c)
          if (dx .ge. 0.5) then
            avg = avg + tmin*0.5
            cc = c/32.
            bb = bb + 0.125 - c/32.
          elseif (dy .le. 0.5) then
            avg = avg + tmin*dx + tmax*(.5-dy) + (dy-dx)*.5*(tmin+tmax)
            cc = c/16.
            bb = bb + 0.25 - c/16.
          else
            avg = avg + 
     &              tmin*dx + (0.5-dx)*0.5*(tmin + ss - 0.5*c*(ss-ssp))
            cc = c/16.
            bb = bb + 0.25 - c/16.
          endif

        else

          dx = (ss-tmin)/((ss-ssp)*c)
          if (dx .ge. 0.5) then
            avg = avg + tmin*0.5
            cc = c/32.
            bb = bb + 0.125 - c/32.
          else
            avg = avg + 
     &              tmin*dx + (0.5-dx)*0.5*(tmin + ss - 0.5*c*(ss-ssp))
            cc = c/16.
            bb = bb + 0.25 - c/16.
          endif

        endif

      elseif (ss .ge. tmax) then

        if (ssm .ge. tmax) then

          avg = avg + tmax*0.5
          aa = a/32.
          bb = bb + 0.125 - a/32.

        elseif (ssm .le. tmin) then

          dx = (ss-tmax)/((ss-ssm)*a)
          dy = (ss-tmin)/((ss-ssm)*a)
          if (dx .ge. 0.5) then
            avg = avg + tmax*0.5
            aa = a/32.
            bb = bb + 0.125 - a/32.
          elseif (dy .le. 0.5) then
            avg = avg + tmax*dx + tmin*(.5-dy) + (dy-dx)*.5*(tmin+tmax)
            aa = a/16.
            bb = bb + 0.25 - a/16.
          else
            avg = avg + 
     &              tmax*dx + (0.5-dx)*0.5*(tmax + ss - 0.5*a*(ss-ssm))
            aa = a/16.
            bb = bb + 0.25 - a/16.
          endif

        else

          dx = (ss-tmax)/((ss-ssm)*a)
          if (dx .ge. 0.5) then
            avg = avg + tmax*0.5
            aa = a/32.
            bb = bb + 0.125 - a/32.
          else
            avg = avg + 
     &        tmax*dx + (0.5-dx)*0.5*(tmax + ss - 0.5*a*(ss-ssm))
            aa = a/16.
            bb = bb + 0.25 - a/16.
          endif

        endif
        
        if (ssp .ge. tmax) then

          avg = avg + tmax*0.5
          cc = c/32.
          bb = bb + 0.125 - c/32.

        elseif (ssp .le. tmin) then

          dx = (ss-tmax)/((ss-ssp)*c)
          dy = (ss-tmin)/((ss-ssp)*c)
          if (dx .ge. 0.5) then
            cc = c/32.
            bb = bb + 0.125 - c/32.
            avg = avg + tmax*0.5
          elseif (dy .le. 0.5) then
            avg = avg + tmax*dx + tmin*(.5-dy) + (dy-dx)*.5*(tmin+tmax)
            cc = c/16.
            bb = bb + 0.25 - c/16.
          else
            avg = avg + 
     &              tmax*dx + (0.5-dx)*0.5*(tmax + ss - 0.5*c*(ss-ssp))
            cc = c/16.
            bb = bb + 0.25 - c/16.
          endif

        else

          dx = (ss-tmax)/((ss-ssp)*c)
          if (dx .ge. 0.5) then
            avg = avg + tmax*0.5
            cc = c/32.
            bb = bb + 0.125 - c/32.
          else
            avg = avg + 
     &              tmax*dx + (0.5-dx)*0.5*(tmax + ss - 0.5*c*(ss-ssp))
            cc = c/16.
            bb = bb + 0.25 - c/16.
          endif

        endif

      else

        if (ssm .le. tmin) then

          dx = (ss-tmin)/((ss-ssm)*a)
          if (dx .ge. 0.5) then
            avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssm)*a)
            aa = a/16.
            bb = bb + 0.25 - a/16.
          else
            avg = avg + tmin*(.5-dx) + dx*0.5*(tmin + ss)
            aa = a/16.
            bb = bb + 0.25 - a/16.
          endif

        elseif (ssm .ge. tmax) then

          dx = (ss-tmax)/((ss-ssm)*a)
          if (dx .ge. 0.5) then
            avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssm)*a)
            aa = a/16.
            bb = bb + 0.25 - a/16.
          else
            avg = avg + tmax*(.5-dx) + dx*0.5*(tmax + ss)
            aa = a/16.
            bb = bb + 0.25 - a/16.
          endif

        else

          avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssm)*a)
          aa = a/8.
          bb = bb + 0.5 - a/8.

        endif
        
        if (ssp .le. tmin) then

          dx = (ss-tmin)/((ss-ssp)*c)
          if (dx .ge. 0.5) then
            avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssp)*c)
            cc = c/16.
            bb = bb + 0.25 - c/16.
          else
            avg = avg + tmin*(.5-dx) + dx*0.5*(tmin + ss)
            cc = c/16.
            bb = bb + 0.25 - c/16.
          endif

        elseif (ssp .ge. tmax) then

          dx = (ss-tmax)/((ss-ssp)*c)
          if (dx .ge. 0.5) then
            avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssp)*c)
            cc = c/16.
            bb = bb + 0.25 - c/16.
          else
            avg = avg + tmax*(.5-dx) + dx*0.5*(tmax + ss)
            cc = c/16.
            bb = bb + 0.25 - c/16.
          endif

        else

          avg = avg + 0.5*0.5*(2.*ss - 0.5*(ss-ssp)*c)
          cc = c/8.
          bb = bb + 0.5 - c/8.

        endif

      endif

      return
      end


                

      SUBROUTINE tridag(alon,alat,a,b,c,r,u,n)
      implicit none
      INTEGER n,nmax
      REAL alon,alat,a(n),b(n),c(n),r(n)
      double precision u(n)
      PARAMETER (nmax=147*12) ! PJD Oceanonly 1870-2014
      INTEGER j
      REAL bet, gam(nmax)
      if (nmax .lt. n) then
           print*, 'Error nmax not declared large enough'
           print*, 'in tridag'
           stop
      endif
      if(b(1).eq.0.) then
          print*, 'longitude = ', alon, '  latitude = ', alat
          pause 'tridag: rewrite equations'
      endif
      bet=b(1)
      u(1)=r(1)/bet
      if (n .gt. 1) then
        do 11 j=2,n
          gam(j)=c(j-1)/bet
          bet=b(j)-a(j)*gam(j)
          if(bet.eq.0.) then
            print*, 'longitude = ', alon, '  latitude = ', alat
            pause 'tridag failed'
          endif
          u(j)=(r(j)-a(j)*u(j-1))/bet
11      continue
        do 12 j=n-1,1,-1
          u(j)=u(j)-gam(j+1)*u(j+1)
12      continue
      endif
      return
      END


      SUBROUTINE cyclic(alon,alat,a,b,c,alpha,beta,r,x,n)
      implicit none
      INTEGER n,nmax
      real alon,alat,alpha,beta,a(n),b(n),c(n),r(n)
      double precision x(n)
      PARAMETER (nmax=147*12) ! PJD Oceanonly 1870-2014
CU    USES tridag
      INTEGER i
      REAL fact,gamma,bb(nmax),u(nmax)
      double precision z(nmax)
      if(n.le.2)pause 'n too small in cyclic'
      if(n.gt.nmax)pause 'nmax too small in cyclic'
      gamma=-b(1)
      bb(1)=b(1)-gamma
      bb(n)=b(n)-alpha*beta/gamma
      do 11 i=2,n-1
        bb(i)=b(i)
11    continue
      call tridag(alon,alat,a,bb,c,r,x,n)
      u(1)=gamma
      u(n)=alpha
      do 12 i=2,n-1
        u(i)=0.
12    continue
      call tridag(alon,alat,a,bb,c,u,z,n)
      fact=(x(1)+beta*x(n)/gamma)/(1.+z(1)+beta*z(n)/gamma)
      do 13 i=1,n
        x(i)=x(i)-fact*z(i)
13    continue
      return
      END


      subroutine wrtdrs1(lcreate, lclose, lconcat, outfile,
     &         nlon, nlat, nmon12, nchunks, iyr1, mon1,
     &         iyrn, monn, chnkname, jyr1out, elem1, elemn,
     &         source, vname, title, units, gauss, array)

      implicit none
      integer lcreate, lclose, lconcat, nlon, nlat, nmon12,
     &          nchunks, iyr1, mon1, iyrn, monn, jyr1out(*)
      real array(*),elem1(2), elemn(2), gauss(*)
      character*16 vname
      character*120 source
      character*80 title
      character*120 outfile, chnkname(nchunks) ! PJD Oceanonly - 1870-2014
      character*40 units

      print*, 'Your version of the code is unable to write '
      print*, 'DRS-formatted files.'

      call exit(1)

      return
      end
      


      subroutine wrtdrs(lcreate, lclose, lconcat, outfile, fnclim,
     &         nlon, nlat, nmon12, nchunks, iyr1, mon1,
     &         iyrn, monn, chnkname, jyr1out, elem1, elemn,
     &         source, vname, title, units, gauss, array)
     
c #include "drsdef.h"
c      include '/usr/local/include/drsdef.h'
c      include '/work/durack1/Shared/150219_AMIPForcingData/src/'
c     & // 'libdrs/lib/drsdef.h'
      include '/work/durack1/Shared/150219_AMIPForcingData/src/drsdef.h'
      
c      implicit none

      integer lcreate, lclose, lconcat, nlon, nlat, nmon12,
     &          nchunks, iyr1, mon1, iyrn, monn, jyr1out(*)
      real array(*),elem1(2), elemn(2), gauss(*)
      character*16 vname, ayr
      character*15 fnclim
      character*120 source
      character*80 title
      character*120 outfile, chnkname(nchunks), fullfile ! PJD Oceanonly - 1870-2014
      character*40 units
      integer seterr, aslun, cluvdb, setname, setdim, 
     +    putdat, cllun, setvdim, putvdim
      integer ierr, ichunk, ix, iy, nnn, mm, iyr, iendyr, imon,
     &    jmon, iend, m, n2, n1, n, lpp(45), nn, nmon, mmm
      real bpp(19)
      real amonths(150)

      write(36,*) 'lcreate = ', lcreate
      write(36,*) 'lclose = ', lclose
      write(36,*) 'lconcat = ', lconcat
      write(36,*) 'outfile = ', outfile
      write(36,*) 'nlon = ', nlon
      write(36,*) 'nlat = ', nlat
      write(36,*) 'nmon12 = ', nmon12
      write(36,*) 'nchunks = ', nchunks
      write(36,*) 'iyr1 = ', iyr1 
      write(36,*) 'mon1 = ', mon1
      write(36,*) 'iyrn = ', iyrn
      write(36,*) 'monn = ', monn
      write(36,*) 'chnkname = ', chnkname
      write(36,*) 'jyr1out = ', jyr1out(1)
      write(36,*) 'elem1 = ', elem1
      write(36,*) 'elemn = ', elemn
      write(36,*) 'gauss = ', gauss(1), gauss(2), gauss(47), gauss(48)

      mm = 0
      do 195 m=mon1,monn
        mm = mm + 1
        amonths(mm) = float(m)
  195 continue

      nmon = monn - mon1 + 1

      if (lconcat .gt. 0) then
        do 200 ichunk=1,nchunks
          OPEN (12+ichunk,FILE=chnkname(ichunk),STATUS='OLD',
     &           FORM='UNFORMATTED')
  200   continue
      endif

      nnn = 0
      mm = 0
      mmm = 0

c     LOOP OVER YEARS

      iyr = iyr1 - 1
 1000 iyr = iyr + 1

        if ((lcreate .gt. 0) .and. 
     &      ((iyr .eq. iyr1) .or. (iyr .eq. jyr1out(mm+1)))) then
          mm = mm + 1
          if (iyr1 .eq. 0) then
            ayr = ''
            if (index(fnclim, 'clim') .gt. 0) ayr = '_' // fnclim
            iendyr = 0
          elseif (iyr .eq. iyrn) then
            write(ayr, '("_", i4)') iyr
            iendyr = iyr
          elseif (jyr1out(mm+1) .eq. -999) then
            write(ayr, '("_", i4, "_", i4)') iyr, iyrn
            iendyr = iyrn
          elseif (jyr1out(mm) .eq. (jyr1out(mm+1)-1)) then  
            write(ayr, '("_", i4)') iyr
            iendyr = iyr
          else
            iendyr = jyr1out(mm+1) - 1
            iendyr = min0(iendyr, iyrn)
            write(ayr, '("_", i4, "_", i4)') iyr, iendyr
          endif

          iend = index(outfile, ' ') - 1
          fullfile = outfile(1:iend)//ayr

          print*, 'Opening drs file:'
          print*, fullfile
          ierr=seterr(6,IDRS_WARNING)
          ierr=seterr(6,IDRS_FATAL)
          ierr=aslun(10,fullfile,11,' ',IDRS_CREATE)
          if  (gauss(1) .ne. gauss(2)) then
            ierr=cluvdb()
            ierr = setname(' ','Latitude','Gaussian Latitudes',
     &       'Degrees',' ')
            ierr=putvdim(10,nlat,gauss,ix,iy)
          endif

          if (iyr .eq. iyr1) then
            imon = mon1
          else
            imon = 1
          endif
   
          if (iendyr .eq. iyrn) then
            jmon = monn
          else
            jmon = nmon12
          endif

          nmon = (iendyr-iyr)*nmon12 + jmon-imon + 1

          if (nmon .gt. 150) then
            print*, 'amonths not dimensioned large enough in'
            print*, 'subroutine wrtdrs.  It should be dimensioned'
            print*, 'no smaller than ', nmon
            call exit(1)
          endif

          do 305 m=1,nmon
            amonths(m) = float(mmm+m)
  305     continue
          mmm = mmm + nmon

          if (lconcat .gt. 0) then
            ierr=cluvdb()
            ierr = setname(' ','time','Month', 'Month',' ')
            ierr=putvdim(10,nmon,amonths,ix,iy)
          endif

        endif      
   
        iyr = iendyr

        print*, title
        write(9,*) title

        if (lconcat .eq. 0) then

          ierr=cluvdb()

          ierr=setname(source,vname,title,units,' ')
          ierr=setdim(1,'longitude','Degrees',nlon,elem1(1),elemn(1))
          if (gauss(1) .ne. gauss(2)) then
             ierr=setvdim(2,' ','Latitude',' ',' ', elem1(2),elemn(2))
          else
             ierr=setdim(2,'latitude','Degrees',nlat,elem1(2),elemn(2))
          endif

          nn = nnn + 1
          nnn = nnn + nlon*nlat*(jmon-imon+1)

          if (jmon .gt. imon) ierr=setdim(3,'time','Month',nmon,
     &            amonths(1), amonths(nmon))

          ierr=putdat(10,array(nn))

        else

          ierr=cluvdb()

          ierr=setname(source,vname,title,units,' ')
          ierr=setdim(1,'longitude','Degrees',nlon,elem1(1),elemn(1))
          if (gauss(1) .ne. gauss(2)) then
             ierr=setvdim(2,' ','Latitude',' ',' ', elem1(2),elemn(2))
          else
             ierr=setdim(2,'latitude','Degrees',nlat,elem1(2),elemn(2))
          endif
          
          imon = 1
  250     jmon = min0(nmon, (imon+nmon12-1))

          ierr = setvdim(3,' ','time', ' ',' ',
     &              amonths(imon), amonths(jmon))

c          loop over months; read pp format temporary files

c         Process each chunk for this month
          n2 = 0
          do 400 m=imon,jmon
             do 300 ichunk = 1,nchunks

c              read chunk in pp-format
               read(12+ichunk)lpp,bpp
               n1 = n2+1
               n2 = n1+lpp(15)-1
               read(12+ichunk) (array(n), n=n1,n2)

 300         continue
 400      continue

          ierr=putdat(10,array)

          imon = jmon + 1
          if (imon .le. nmon) go to 250

        endif
 
        if ((lclose .gt. 0) .and. ((iyr .eq. iyrn) .or. 
     &     ((iyr+1) .eq. jyr1out(mm+1)))) ierr=cllun(10)

      if (iyr .lt. iyrn) go to 1000

      if (lconcat .gt. 0) then
        do 600 ichunk=1,nchunks
          close(12+ichunk)
  600   continue
      endif
c       
      return
      end


      subroutine wrtasc(outfile,fnclim,nlon,nlat,nmon12,iyr1,mon1,iyrn,
     &      monn,elem1,elemn,model,center,source,vname,title,units,
     &      grid,amissout,amissasc,lconcat,nchunks,chnkname,
     &      jyr1out, array)
        
      implicit none

      integer nlon, nlat, iyr1, iyrn, mon1, monn, nmon12, 
     &      lconcat, nchunks, mm, iendyr, jyr1out(*)
      real array(*),elem1(2), elemn(2), amissout, amissasc
      character*16 vname, model
      character*15 fnclim
      character*120 source, fullfile
      character*80 title, grid, center
      character*120 outfile, chnkname(nchunks)
      character*20 ayr
      character*40 units
      integer n, nn, nnn, imon, jmon, iyr, iend, n1, n2, ichunk,
     &    m, lpp(45)
      real bpp(19), fmax, fmin
      
      print*, 'Preparing to write ascii data ...'
      print*, title
      write(9,*) title

      nnn = 0
      if (index(units, '%') .eq. 0)  then
         fmax=9999.9
         fmin=-999.99
      else
         fmax=99999.0
         fmin=-9999.9
      endif


      if (lconcat .gt. 0) then
        do 200 ichunk=1,nchunks
          OPEN (12+ichunk,FILE=chnkname(ichunk),STATUS='OLD',
     &           FORM='UNFORMATTED')
  200   continue
      endif
    
      mm = 0

      do 1000 iyr = iyr1, iyrn

        if (iyr .eq. iyr1) then
          imon = mon1
        else
          imon = 1
        endif
 
        if (iyr .eq. iyrn) then
          jmon = monn
        else
          jmon = nmon12
        endif

        if ((iyr .eq. iyr1) .or. (iyr .eq. jyr1out(mm+1))) then
          mm = mm + 1
          if (iyr1 .eq. 0) then
            ayr = ''
            if (index(fnclim, 'clim') .gt. 0) ayr = "_" // 
     &                   fnclim(1:index(fnclim,' ')-1) // '.asc '
          elseif (iyr .eq. iyrn) then
            write(ayr, '("_", i4, ".asc")') iyr
            iendyr = iyrn
          elseif (jyr1out(mm+1) .eq. -999) then
            write(ayr, '("_", i4, "_", i4,".asc")') iyr, iyrn
            iendyr = iyrn
          elseif (jyr1out(mm) .eq. (jyr1out(mm+1)-1)) then  
            write(ayr, '("_", i4,".asc")') iyr
            iendyr = iyr
          else
            iendyr = jyr1out(mm+1) - 1
            iendyr = min0(iendyr, iyrn)
            write(ayr, '("_", i4, "_", i4,".asc")') iyr, iendyr
          endif

          iend = index(outfile, ' ') - 1
          fullfile = outfile(1:iend)//ayr
          print*, 'Opening ', fullfile
          open(12, file=fullfile, status='new')

          write(12,'(a16 / a80 / a80 / a40 / a40 / a60, 2x, a16 / a80)') 
     &        vname, title, source(1:80), source(81:120), units, center,
     &        model, grid
          write(12,'(i5, " = number of longitudes")') nlon
          write(12,'(i5, " = number of latitudes")')  nlat
          write(12,'(f10.4, " = first longitude")')      elem1(1)
          write(12,'(f10.4, " = last longitude")')       elemn(1)
          write(12,'(f10.4, " = first latitude")')       elem1(2)
          write(12,'(f10.4, " = last latitude")')        elemn(2)

          if (iyr1 .eq. 0) then
            write(12,'("    0 -     0 = climatology")')   
          else
            write(12,'(i5, " - ", i5, " = years written")') iyr, iendyr
          endif
          write(12,'(i5, 
     &    " = first month (jan, feb, mar, ... = 1, 2, 3, ...)")') imon
          write(12,'(i5, 
     &    " = last month (jan, feb, mar, ... = 1, 2, 3, ...)")') jmon
          write(12,'(1x)')

        endif

        if (lconcat .eq. 0) then

          nn = nnn + 1 
          nnn = nnn + nlon*nlat*(jmon-imon+1)

        else

c          loop over months; read pp format temporary files

c          Process each chunk for this month
           n2 = 0
           do 400 m=imon,jmon
             do 300 ichunk = 1,nchunks

c              read chunk in pp-format
               read(12+ichunk)lpp,bpp
               n1 = n2+1
               n2 = n1+lpp(15)-1
               read(12+ichunk) (array(n), n=n1,n2)

 300         continue
 400       continue

           nn = 1
           nnn = n2

        endif

        do 100 n=nn,nnn
          if (array(n) .eq. amissout) array(n) = amissasc
          if ((array(n) .gt. fmax) .or. 
     &        (array(n) .lt. fmin)) then
            print*, 'Error:  output data contains numbers like ', 
     &              array(n)
            print*, 'which are too large to express in the specified ',
     &              'format'
            call exit(1)
          endif
  100   continue

        if (index(units, '%') .eq. 0) then
           write(12,'(10f8.3)') (array(n), n=nn,nnn)
        else
           write(12,'(10f8.2)') (array(n), n=nn,nnn)
        endif
 
        if ((iyr .eq. iyrn) .or. ((iyr+1) .eq. jyr1out(mm+1)))
     &        close(12)

 1000 continue

      if (lconcat .gt. 0) then
        do 600 ichunk=1,nchunks
          close(12+ichunk)
  600   continue
      endif

c       
      return
      end


      SUBROUTINE readpp(filename, nlat, nlon, mlat, lat1, lats, iyr1rd, 
     &           mon1rd, iyr1in, mons, amissin, amissout, alons, alats, 
     &           array, amask, wtfrac, lpp, bpp, space)

c       Read in every month to month value from pp format file within 
c       the given chunk
c
c       assumes data stored from north to south

c       NB Uses unit number 10
c         IN:
c           filename : Filename of pp format file to be read
c           nlat     : total number of latitudes spanning globe
c           nlon     : Size of input data longitude dimension
c           mlat     : Size of input data latitude dimension
c           lat1     : Index of first latitude to be read
c           lats     : Number of latitudes to be read
c           iyr1in(1): Year at beginning of period to be read
c           iyr1in(i): for i>1, 1st year of succeeding input files
c           iyr1rd   : Year at beginning of period to be read
c           mon1rd   : Month at beginning of period to be read
c           mons     : Number of months to be read
c           amissin  : Value that indicates missing data in file read
c           amissout : value assigned to missing data before returned
c           space    : scratch space needed to read in data
c         OUT:
c           alons    : Vector containing the values of each longitude 
c                          in array
c           alats    : Vector containing the values of each latitude in 
c                          array
c           array    : Chunk of data to be returned
c           amask    : proportional to grid-cell area; 0. if missing
c           wtfrac   : 0. or 1. if missing or not.
c           lpp      : integer header
c           bpp      : real header

      IMPLICIT NONE
      CHARACTER*120 filename(*)
      INTEGER      nlat,lat1,lats,nlon,mlat,iyr1rd,mon1rd,iyr1in(*),mons
      REAL         array(nlon,mlat,mons),amask(nlon,nlat),amissin,
     &               amissout, wtfrac(nlon,nlat), space(nlon,nlat)
      REAL         alons(nlon),alats(nlat)
      INTEGER      lyear,lmon
      INTEGER      i,j,cnt,ncnt,nn, ij
c        Part one of pp-header
      INTEGER      lpp(45)                
c        Part two of pp-header
      REAL         bpp(19)  
      REAL         radlat1,radlat2,sindiff,pi

      pi=4.0*atan(1.)
c        3.1415926

      print*, 'Opening input file: ', filename(1)
      OPEN (10,FILE=filename(1),STATUS='OLD',FORM='UNFORMATTED')

      lyear = 0
      lmon = 0
      ij = 1

c     input file counter
      nn = 1

c     number of data points read (excluding missing data)
      ncnt = 0
      
c       This maintains an index of the current month
      cnt = 1    

c        Repeat until mons months have been processed
      DO WHILE (cnt.LE.mons) 

        if ((lyear .eq. (iyr1in(nn+1)-1)) .and. (lmon .eq. 12)) then
          ij = cnt - 1
          print*, 'found cnt= ', ij,' year= ',lyear,' month=', lmon
          close(10)
          nn = nn + 1
          print*, 'Opening input file: ', filename(nn)
          OPEN (10,FILE=filename(nn),STATUS='OLD',FORM='UNFORMATTED')
          ij = 0
        endif

c         Read in next field from data file
        READ (10,END=999)lpp,bpp

c       Extract the month and year for this input data array
        lyear = lpp(1)
        lmon = lpp(2)

c         If the data is for the period that we want
        IF ((lyear.GT.iyr1rd).OR.
     &     ((lyear.EQ.iyr1rd).AND.(lmon.GE.mon1rd))) THEN

          if (nlat .ne. lpp(18)) then
          print*, 'ERROR:  For input data in pp format, the number of' 
          print*, '   latitudes stored in input file should equal'
          print*, '   the number declared in the parameter statement' 
          print*, '   number stored = ', lpp(18)
          print*, '   nlat          = ', nlat
          stop
          endif

          if (nlon .ne. lpp(19)) then
          print*, 'ERROR:  For input data in pp format, the number of' 
          print*, '   longitudes stored in input file should equal'
          print*, '   the number declared in the parameter statement' 
          print*, '   number stored = ', lpp(19)  
          print*, '   nlon          = ', nlon
          stop
          endif


          READ (10,END=999) space

          if ((ij .eq. 1) .or. (cnt .eq. mons))
     &     print*, 'found cnt= ', cnt,' year= ',lyear,' month=', lmon
          ij=0
c          
          if (cnt .eq. 1) then

            DO j = lat1, lat1+lats-1

c               Set vector of latitude values
              alats(j) = bpp(14) + j*bpp(15)

c               Find the difference in sines of north/south boundaries  
c                 of square (assumes that if data stored from north to 
c                 south then bpp(15)<0.)

              radlat1 = (2.0*pi/360.0) * (bpp(14) + (j-0.5)*bpp(15))
              radlat2 = (2.0*pi/360.0) * (bpp(14) + (j+0.5)*bpp(15))

              sindiff = SIN(radlat1) - SIN(radlat2)

                do i = 1, nlon
                  amask(i,j) = abs(sindiff) / (2.*nlon)
                  wtfrac(i,j) = 1.0
                enddo
            enddo

c               Assign values to alons
            do i=1,nlon
              alons(i) = bpp(16) + i*bpp(17)
            enddo

          endif

          do j = lat1, lat1+lats-1

            DO i = 1,nlon

              IF (space(i,j) .NE. amissin) THEN
                ncnt = ncnt + 1
                array(i,j - lat1 + 1,cnt) = space(i,j)
             ELSE
                array(i,j - lat1 + 1,cnt) = amissout
                amask(i,j) = 0.0
                wtfrac(i,j) = 0.0
              ENDIF

            ENDDO
          ENDDO

c         Increment count of months used
          cnt = cnt + 1  

        else

          read(10) 

        ENDIF
  
      ENDDO

      go to 1000

  999 CONTINUE

      print*, 'Error in subroutine readpp:'
      print*, 'reached end of file before all data were read'
      print*, 'Check specification of months and years you '
      print*, 'want to read and input file name.'
      print*, 'Also check specification of nlon and nlat'

      call exit(1)


 1000 continue

      PRINT*, 'number of longitude x latitude x month data read'
      print*, '   (excluding missing data) = ', ncnt
      ncnt = ncnt/(cnt-1)
      PRINT*, 'number of months read = ', cnt-1

      CLOSE (10)
      print*, 'returning to main program'

      RETURN
      END

      subroutine wrtpp(filename, iobs, nlon, mlat, alat1, lats, 
     &        monlen, iyear1, mon1, nmon, lpp, bpp, array)

c         Write the adjusted data to a series of intermediate files,
c       each containing a chunk of data, which can later be adjusted
c       to pp-format
c         Uses logical unit number 11
c       IN:
c         Filename  : The name of the file to which the chunk will be 
c                       written
c         iobs      : 1 if observed data; 0 otherwise
c         nlon,mlat : The dimensions of the data array
c         alat1     : First latitude of chunk
c         lats      : Number of latitudes in each chunk
c         iyr1      : Year at which chunk begins
c         mon1      : Month at which chunk begins
c         nmon      : Number of months in chunk
c         array     : Chunk itself

c       OUT         :No value is returned, but the chunk stored in array
c                      is written to unit number 11, in pp format to
c                      retain info about nlat/nlon coverage
c            

      IMPLICIT NONE

      INTEGER       nlon,lats,iyear1,mon1,nmon,mlat,iobs,mon
      REAL          array(nlon,mlat,nmon),alat1
      CHARACTER*120  filename, outfile
      INTEGER       imon, ijk
      REAL          bpp(*), bbpp(19), day
      INTEGER       lpp(*), llpp(45), i,j,iyr, iyr1, imon1, iday1,
     &            ihour1, iyr2, imon2, iday2, ihour2, lastyr, ij, iend
      integer monlen(12,2)

      lastyr = (iyear1*12+mon1+nmon-2)/12

c       Define pp-header
      DO i=1,19
        llpp(i)=lpp(i)
        bbpp(i)=bpp(i)
      ENDDO
      DO i=20,45
        llpp(i)=lpp(i)
      ENDDO

      llpp(15) = nlon*lats   ! Size of chunk
      llpp(18) = lats        ! Height of chunk

c ???    check the following
      if (iobs .gt. 0) then
        llpp(25) = 128
        if (iyear1 .eq. 0) llpp(13) = 31
        if (iyear1 .ne. 0) llpp(13) = 21
      else
        llpp(25) = 0
        if (iyear1 .eq. 0) llpp(13) = 31
        if (iyear1 .ne. 0) llpp(13) = 1
      endif
      
c       Top of chunk
      bbpp(14) = alat1 - bbpp(15) 

      iend = index(filename, ' ') - 1
      outfile = filename(1:iend)//'.pp'

      print*, 'opening pp file:  '
      print*,  outfile 
      ijk = index(outfile, 'temp')
      if (ijk .eq. 0) then
        OPEN (11,FILE=outfile,STATUS='NEW',FORM='UNFORMATTED')
      else
        OPEN (11,FILE=outfile,STATUS='unknown',FORM='UNFORMATTED')
      endif

c       Initialise the month and year to the beginning of chunk
      iyr = iyear1
      imon = mon1 - 1

      DO mon = 1,nmon

        imon = imon + 1
        IF (imon.EQ.13) THEN
          iyr = iyr + 1
          imon = 1
        ENDIF

        if (iyear1 .eq. 0) then
c          we have climatological data

          if (iobs .gt. 0) then
c           we have monthly mean data

            ihour1 = 0
            iday1 = 1
            imon1 = imon
c ???    check the following
            iyr1 = iyear1

            ihour2 = 0
            iday2 = 1 
            imon2 = imon + 1
c ???    check the following
            iyr2 = lastyr

            if (imon2 .eq. 13) imon2 = 1

          else
c           we have mid-month data

            day = monlen(imon,1)/2.0 + 1.01
            iday1 = day
            if ((day - iday1) .gt. 0.1) then
                ihour1 = 12
            else
                ihour1 = 0
            endif

c ???    check the following
            iyr1 = iyear1
            imon1 = imon

            ihour2 = ihour1
            iday2 = iday1
            imon2 = imon
c ???    check the following
            iyr2 = lastyr

          endif

        else
c         we don't have climatology

          if (iobs .gt. 0) then
c           we have monthly mean data

            ihour1 = 0
            iday1 = 1
            imon1 = imon
            iyr1 = iyr

            ihour2 = 0
            iday2 = 1 
            imon2 = imon + 1
            iyr2 = iyr 

            if (imon2 .eq. 13) then
              imon2 = 1
              iyr2 = iyr2 + 1
            endif

          else

c           we have mid-month data

c ???     generalize for Julian, 360-day, 365-day years
            ij = 1
            if (((mod(iyr,4) .eq. 0) .and. (mod(iyr,100) .ne. 0)) .or. 
     &          (mod(iyr,400) .eq. 0)) ij = 2

            day = monlen(imon,ij)/2.0 + 1.01
            iday1 = day
            if ((day - iday1) .gt. 0.1) then
              ihour1 = 12
            else
                ihour1 = 0
            endif

            imon1 = imon
            iyr1 = iyr

            ihour2 = ihour1
            iday2 = iday1
            imon2 = imon1
            iyr2 = iyr

          endif

        endif

        llpp(1)  = iyr1        
        llpp(2)  = imon1       
        llpp(3)  = iday1       
        llpp(4)  = ihour1      
        llpp(7)  = iyr2        
        llpp(8)  = imon2       
        llpp(9)  = iday2      
        llpp(10) = ihour2      

c          Write each month of chunk in pp-format
        WRITE(11) llpp, bbpp
        WRITE(11) ((array(i,j,mon), i=1,nlon), j=1,lats)

      ENDDO
 
      CLOSE (11)

      RETURN
      END

      SUBROUTINE ppconcat(chnkname, outfile, fnclim, nmon12, nlon, nlat,
     &         nchunks, iyr1out, mon1out, iyrnout, monnout, jyr1out,
     &         aa)
c       Concatenate all the intermediary files back into a monthly 
c          pp file
c       Note, it uses logical units 12 through 12+number of chunks
c       IN:
c          chnkname  : The names of the files that each of the chunks 
c                        are stored in, in N-S order
c          outfile     : Name of output file
c          nmon12      : Number of months in a year
c          nlon,nlat   : Dimensions of output file
c          nchunks     : The number of chunks
c          iyr1out     : The first year for which output is written
c          mon1out     : The first month for which output is written
c          iyrnout     : The last year for which output is written
c          monnout     : The last month for which output is written
c          jyr1out(*)  : First year of each output file
c          aa(nlon,nlat): scratch space for reading and writing
c       OUT:
c                      : No output is returned, but the chunks stored in 
c                        the files given by chnkname have been 
c                        concatenated together

      IMPLICIT NONE
      integer         nmon12, nlon, nlat, nchunks, iyr1out, mon1out,
     &                iyrnout, monnout, lpp(45), jyr1out(*)
      CHARACTER*120    chnkname(nchunks),outfile, fullfile
      real            bpp(19), aa(nlon,nlat)
      character*20    ayr
      character*15    fnclim

      INTEGER         imon, mm, iendyr, iyr, m, jmon, i, j,
     &                i2, j1, j2, ichunk, iend
      real            alat0

      DO ichunk = 1,nchunks
        OPEN (12+ichunk,FILE=chnkname(ichunk),STATUS='OLD',
     &           FORM='UNFORMATTED')
      ENDDO

      mm = 0

      do 1000 iyr = iyr1out, iyrnout

        if (iyr .eq. iyr1out) then
          imon = mon1out
        else
          imon = 1
        endif
 
        if (iyr .eq. iyrnout) then
          jmon = monnout
        else
          jmon = nmon12
        endif

        if ((iyr .eq. iyr1out) .or. (iyr .eq. jyr1out(mm+1))) then
          mm = mm + 1
          if (iyr1out .eq. 0) then
            ayr = fnclim(1:index(fnclim,' ')-1) // '.pp'
          elseif (iyr .eq. iyrnout) then
            write(ayr, '("_", i4, ".pp")') iyr
          elseif (jyr1out(mm+1) .eq. -999) then
            write(ayr, '("_", i4, "_", i4,".pp")') iyr, iyrnout
          elseif (jyr1out(mm) .eq. (jyr1out(mm+1)-1)) then  
            write(ayr, '("_", i4,".pp")') iyr
          else
            iendyr = jyr1out(mm+1) - 1
            iendyr = min0(iendyr, iyrnout)
            write(ayr, '("_", i4, "_", i4,".pp")') iyr, iendyr
          endif

          iend = index(outfile, ' ') - 1
          fullfile = outfile(1:iend)//ayr

          print*, 'Opening pp file:'
          print*, fullfile

          OPEN (12,FILE=fullfile,STATUS='unknown',FORM='UNFORMATTED')

        endif

c       Concatenate each month consecutively

        do 500 m=imon,jmon
  
c        Process each chunk for this year
          j2 = 0
          DO 300 ichunk = 1,nchunks
 
c           Read chunk in pp-format
            READ(12+ichunk)lpp,bpp
            if (ichunk .eq. 1) alat0 = bpp(14)
            i2 = lpp(19)
            j1 = j2 + 1
            j2 = j1 + lpp(18) - 1
            read(12+ichunk) ((aa(i,j), i=1,i2), j=j1,j2)

  300     continue 

c         Assign the correct pp-header

          lpp(18)=    j2
          lpp(15) =   i2*j2
          bpp(14)=    alat0
        
c         Write the regenerated monthly field to unit 12
          WRITE(12)lpp,bpp
          WRITE(12) ((aa(i,j), i=1,i2), j=1,j2)
      
  500   continue

        if ((iyr .eq. iyrnout) .or. ((iyr+1) .eq. jyr1out(mm+1)))
     &        close(12)

 1000 continue

      DO ichunk = 1,nchunks
        close(12+ichunk)
      ENDDO

      RETURN
      END

      subroutine getobs1(iregrid, gtype, msklndi, sftfilei, sftnamei, 
     &        msklndo, sftfileo, sftnameo, inputfil, varin,
     &        nlon, nlat, mlat, lat1, latn, iyr1in, mon1in, iyr1rd, 
     &        mon1rd, iyrnrd, monnrd, nmon12, amissin, amissout, wtmin,
     &        sst, sstwts, wtfrac, alons, alats,
     &        ntlat, rlat0, ntlon, rlon0)

      implicit none

      integer msklndi, msklndo, nlon, nlat, mlat, lat1, latn, iyr1in,
     &         mon1in, iyr1rd, mon1rd, iyrnrd, monnrd, nmon12,
     &         ntlat, ntlon, iregrid
      real sst(nlon,mlat,*), alons(nlon), alats(*), rlat0, rlon0,
     &      sstwts(nlon,*), wtfrac(nlon,*), amissin, amissout, wtmin
      character*120 sftfilei, sftfileo, inputfil
      character*16 sftnamei, sftnameo, varin, gtype

      print*, 'Your version of the code cannot read files with'
      print*, 'the EzGet software (i.e., files other than pp-format).'

      call exit(1)

      return
      end



      subroutine getobs(iregrid, gtype, msklndi, sftfilei, sftnamei, 
     &        msklndo, sftfileo, sftnameo, inputfil, varin,
     &        nlon, nlat, mlat, lat1, latn, iyr1in, mon1in, iyr1rd, 
     &        mon1rd, iyrnrd, monnrd, nmon12, amissin, amissout, wtmin,
     &        sst, sstwts, wtfrac, alons, alats,
     &        ntlat, rlat0, ntlon, rlon0)

      implicit none

      integer nlon, nlat, mlat, lat1, latn, iyr1in, mon1in, iyr1rd,
     &         mon1rd, iyrnrd, monnrd, nmon12, msklndi, msklndo,
     &         ntlat, ntlon, iregrid

      real sst(nlon,mlat,*), alons(nlon), alats(*), rlat0, rlon0,
     &      sstwts(nlon,*), wtfrac(nlon,*), amissin, amissout, wtmin
      character*120 sftfilei, sftfileo, inputfil
      character*16 sftnamei, sftnameo, varin, gtype

      integer len, n, ierr, i3, i4, ierr1, i, j, ij, nn, ilon, ilat,
     &     lats, mon1, monn, mons, k, m, ii, jj, kk, isize, nnlon,
     &     nnlat

      real dlat, dlon

      integer miallocw, mifree

      pointer (ptwt, wts), (ptgeog, ageog), (pbdlat, bdlat),
     &     (ptland, sstland)
      real wts(*), ageog(*), bdlat(*), sstland(*)

      lats = latn - lat1 + 1
c   for AMIP input, the first
      mon1 = (iyr1rd-iyr1in)*nmon12 + mon1rd-mon1in + 1
      monn = (iyrnrd-iyr1in)*nmon12 + monnrd-mon1in + 1
      mons = monn - mon1 + 1

      call initget
      call defmisc('input missing value', 'real', amissin)
      call defmisc('output missing value', 'real', amissout)

      call defvar(2, varin, inputfil)

      if (iregrid .eq. 0) then

c       *** extract monthly mean data, possibly mask land, but
c             don't regrid

c        check latitude dimension
        call defdim(2, 1, 'longitude', 'unit', 'as saved', 0., 0., 0.)
        call defdim(2, 2, 'latitude', 'unit', 'as saved', 0., 0., 0.)
        call defdimi(2, 3, 'time', 'unit', 1, 1)
        call shape(2, ilon, ilat, i3, i4, isize)

        if (ilat .ne. nlat) then 
          print*, 'ERROR:  When reading observed data (without '
          print*, '   regridding), the number of latitudes stored in' 
          print*, '   input file should equal the number declared in '
          print*, '   the parameter statement.  You should '
          print*, '   set nlat = ', ilat  
          stop
        endif
        if (ilon .ne. nlon) then 
          print*, 'ERROR:  When reading observed data (without '
          print*, '   regridding), the number of longitudes stored in' 
          print*, '   input file should equal the number declared in '
          print*, '   the parameter statement.  You should '
          print*, '   set nlon = ', ilon  
          stop
        endif

c       **get full weights

        call defdim(2, 1, 'longitude', 'width', 'as saved', 
     &                        0.0, 0.0, 360.0)
        call defdimi(2, 2, 'latitude', 'cosine', lat1, latn)
        call defdimi(2, 3, 'time', 'unit', 1, 1)

        ilon = nlon
        ilat = lats
        i3 = 1
        i4 = 0
        call getdata(2, nlon,mlat,1,0, ilon,ilat,i3,i4, wtfrac, sst)

        if (msklndi .gt. 0) then

c         * get weights associated with ocean and sea-ice only

          call defvar(1, sftnamei, sftfilei)

          call defmisc('mask type in', 'integer', msklndi)
          call defgeog(2, 'in', 1, 'ocean, seaice')
          call getdata(2, nlon,mlat,1,0, ilon,ilat,i3,i4, sstwts, sst)

        else

c         * copy full weights to sstwts
           
          do 100 j=1,ilat
            do 90 i=1,ilon
              sstwts(i,j) = wtfrac(i,j)
   90       continue
  100     continue

        endif

c        ** extract all monthly mean data, but don't mask or regrid

        call defdimi(2, 3, 'time', 'unit', mon1, monn)
        len = nlon*lats*mons
        call defmisc('data size', 'integer', len)
        call getfield(2,sst)


        if (mlat .ne. lats) then
c         rearrange in memory:

          do 140 k=mons,1,-1
            do 130 j=lats,1,-1
              do 120 i=nlon,1,-1
                m = (k-1)*nlon*lats + (j-1)*nlon + i - 1
                kk = m/(nlon*mlat) + 1
                jj = (m - (kk-1)*nlon*mlat)/nlon + 1
                ii = m - (kk-1)*nlon*mlat - (jj-1)*nlon + 1
                sst(i,j,k) = sst(ii,jj,kk)        
  120         continue
  130       continue
  140     continue

        endif

      else

c       *** extract monthly mean data, possibly mask land using 
c           geography from original and/or target grid, and regrid to a 
c           target grid.

       if (sftfileo .eq. 'none') then

           len = ntlon*lats
           ilon = ntlon
           ilat = ntlat
           dlon = 360./ntlon

           call defdim(2, 1, 'longitude', 'width', 'range', 
     &                        0.0, 0.0, 360.0)
           call defdimi(2, 3, 'time', 'unit', 1, 1)

           if (gtype .eq. 'cosine') then 

               dlat = -((2.*rlat0)/(ntlat-1))

c              IF (abs(rlat0) .eq. 90.) dlat=dlat*(1.-1.e-7)

               call defdim(2, 2, 'latitude', 'cosine', 'range',
     &                        0.0, 0.0, 0.0)
               call defregrd(2, 'uniform', 0, 'area-weighted', ntlat,
     &               rlat0, dlat, ntlon, rlon0, dlon)

           elseif (gtype .eq. 'gaussian') then 

               if (rlat0 .gt. 0.) then
                  call defdim(2, 2, 'latitude', 'cosine', 'range',
     &                        90., -90., 0.0)
               else
                  call defdim(2, 2, 'latitude', 'cosine', 'range',
     &                        -90., 90., 0.0)
               endif

               call defregrd(2, 'gaussian', 0, 'area-weighted', ntlat,
     &               0., 0., ntlon, rlon0, dlon)

           else

              print*, ' Error in mkgisst -- gtype must be either ',
     &             '"cosine" or "gaussian", not ', gtype 
              stop

           endif

           call shape(2, nnlon, nnlat, i3, i4, isize)

           i3 = nnlat+1
           ierr1 = miallocw(i3, pbdlat)
           call getedges(2, 2, bdlat)

           ilat = lats

           call defdim(2, 2, 'latitude', 'cosine', 'range',  
     &                  bdlat(lat1), bdlat(latn+1), 0.0)

           ierr = mifree(pbdlat)

           if (gtype .eq. 'cosine') then 

               i3 = 0

               call defregrd(2, 'uniform', 0, 'area-weighted', i3,
     &               rlat0, dlat, ntlon, rlon0, dlon)

           endif

        else

        call defvar(3, sftnameo, sftfileo)

c        check latitude dimension
        call defdim(3, 1, 'longitude', 'unit', 'as saved', 0., 0., 0.)
        call defdim(3, 2, 'latitude', 'unit', 'as saved', 0., 0., 0.)
        call defdimi(3, 3, 'time', 'unit', 1, 1)
        call shape(3, nnlon, nnlat, i3, i4, isize)

c       ** extract full area weights

        call defdim(2, 1, 'longitude', 'width', 'range', 
     &                        0.0, 0.0, 360.0)
        call defdim(2, 2, 'latitude', 'cosine', 'range',
     &                        0.0, 0.0, 0.0)
        call defdimi(2, 3, 'time', 'unit', 1, 1)
        call defdim(3, 1, 'longitude', 'width', 'as saved',
     &                        0.0, 0.0, 360.0)
        call defdimi(3, 2, 'latitude', gtype, lat1, latn)
c        call defdim(3, 2, 'latitude', gtype, 'as saved',
c     &                        0.0, 0.0, 0.0)

        call defregrd(2, 'to', 3, 'area-weighted', 
     &                        0, 0.0, 0.0, 0, 0.0, 0.0)

        len = nlon*lats

       endif

        if (nnlat .ne. nlat) then 
          print*, 'ERROR:  When regridding input data to a target grid,'
          print*, '   the number of latitudes on the target grid' 
          print*, '   should equal the number declared in '
          print*, '   the parameter statement.  You should '
          print*, '   set nlat = ', nnlat  
          stop
        endif
        if (nnlon .ne. nlon) then 
          print*, 'ERROR:  When regridding input data to a target grid,'
          print*, '   the number of longitudes on the target grid' 
          print*, '   should equal the number declared in '
          print*, '   the parameter statement.  You should '
          print*, '   set nlon = ', nnlon  
          stop
        endif

           call shape(2, nnlon, nnlat, i3, i4, isize)

        call getnogap(2, len, wtfrac, sst)

c       ** extract monthly data

        call defdimi(2, 3, 'time', 'unit', mon1, monn)

        if (msklndi .gt. 0) then
          call defvar(1, sftnamei, sftfilei)
          call defmisc('mask type in', 'integer', msklndi)
          call defgeog(2, 'in', 1, 'ocean, sea-ice')
        endif

        len = nlon*lats
        ierr1 = miallocw(len, ptgeog)

        if (msklndo .gt. 0) then

          call defmisc('data size', 'integer', len)
          call defmisc('mask type out', 'integer', msklndo)
          call defgeog(3, 'out', 3, 'ocean, sea-ice')
          call getvdata(3, ageog, sstwts)
          call defgeog(2, 'out', 3, 'ocean, sea-ice')

        else

          do 280 n=1,len
             ageog(n) = 1.
 280      continue
           
        endif

        len = nlon*mlat*mons
        ierr1 = miallocw(len, ptwt)
        ilon=nlon
        ilat=lats
        i3=mons
        i4=0

        call getdata(2, nlon,mlat,mons,0, ilon,ilat,i3,i4, wts, sst)
 
c       copy weights from first month to sstwts
        do 300 j=1,ilat
          do 290 i=1,ilon
            n = (j-1)*ilon + i
            sstwts(i,j) = wts(n)
  290     continue
  300   continue
  
      endif

      do 360 j=1,ilat
        do 350 i=1,ilon
          if (wtfrac(i,j) .gt. 0.0) wtfrac(i,j)=sstwts(i,j)/wtfrac(i,j)
  350   continue
  360 continue

      call getcoord(2, 1, alons)
      call getcoord(2, 2, alats)

      do 390 j=1,ilat ! PJD Oceanonly - 1870-2014
         if (abs(alats(j)) .lt. 0.00001) alats(j) = 0.0
         if (alats(j) .gt. 90.) alats(j) = 90.
         if (alats(j) .lt. -90.) alats(j) = -90.
 390  continue

      len = 0

      if ((iregrid .ne. 0) .and. (msklndi .ne. 0)) then

         if (msklndo .ne. 0) then

            do 430 j=1,ilat
               do 420 i=1,ilon
                  ij = (j-1)*ilon+i
                  if (wtfrac(i,j) .le. wtmin) then
                     if (ageog(ij) .gt. 0.0) then
                        CALL clrtable
c     print*,
c     &           'WARNING -- target ocean grid cell is incompatible'
c     print*, '     with Fiorino ocean'
c     print*, 'alats = ', alats(j), ' alons = ', alons(i), 
c     &             '  wtfrac = ', wtfrac(i,j)
                        write(9,*)
     &              'WARNING -- target ocean grid cell is incompatible'
                        write(9,*) '     with Fiorino ocean'
                        write(9,*) 'alats = ', alats(j), ' alons = ', 
     &                       alons(i), '  wtfrac = ', wtfrac(i,j)
                        len = len + 1 
                        
                        call defdim(3, 1, 'longitude', 'width', 
     &                       'nearest', alons(i), alons(i), 360.0)
                        call defdim(3, 2, 'latitude', gtype, 'nearest',
     &                       alats(j), alats(j), 0.0)
                        
                        call defgeog(2, 'in', 0, ' ')
                        call defgeog(2, 'out', 0, ' ')
                        
                        nn = mons+1
                        call getnogap(2, mons, wts(nn), wts)
                        
                        do 405 n=1,mons
                           sst(i,j,n) = wts(n)
 405                    continue
                        
                        wtfrac(i,j) = wtmin*1.01
                        
                     else
                        
                        do 410 n=1,mons
                           sst(i,j,n) = amissout
 410                    continue
                        
                     endif
                     
                  endif
 420           continue
 430        continue
            
            ierr = mifree(ptgeog)
            
         else

c           there *is* a land mask for the original grid, but none
c              for the target grid.  Need to fill in the all
c              land points on the target grid with data

            i3 = lats*nlon*mons
            ierr1 = miallocw(i3, ptland)

            call defgeog(2, 'in', 0, ' ')
            ilon = nlon
            ilat = lats
            i3 = mons
            i4 = 0
            call getdata(2, nlon,lats,mons,0, ilon,ilat,i3,i4, wts, 
     &           sstland)

            do 530 j=1,ilat
               do 520 i=1,ilon
                  ij = (j-1)*ilon+i
                  if (wtfrac(i,j) .le. wtmin) then
                     wtfrac(i,j) = wtmin*1.01
                     len = len+1
                     do 505 n=1,mons
                        i4 = (n-1)*ilat*ilon + ij
                        sst(i,j,n) = sstland(i4)
 505                 continue
                  endif
 520           continue
 530        continue

            ierr = mifree(ptland)

         endif
      endif

      if (iregrid .ne. 0) ierr = mifree(ptwt)

      call closeget

      if ((msklndo .ne. 0) .and. (len .gt. 0)) then
         print*, ' '
         print*, 'Model land/sea mask is incompatible with 1x1 '
         print*, ' land/sea mask at ', len, ' grid cells.'
         print*, ' For these grid cells, land values will be used.' 
         print*, ' '
         write(9,*) ' '
         write(9,*) 
     &        'Model land/sea mask is incompatible with 1x1 Fiorino '
         write(9,*) ' land/sea mask at ', len, ' grid cells.'
         write(9,*) ' For these grid cells, land values will be used.' 
         write(9,*) ' '
c         call exit(1)
      endif

      if ((msklndi .ne. 0) .and. (len .gt. 0) .and. (msklndo .eq. 0)) 
     &       then
         print*, ' '
         print*, 'Land values will be used for ', len, ' grid cells.'
         print*, ' '
         write(9,*) ' '
         write(9,*) 
     &         'Land values will be used for ', len, ' grid cells.'
         write(9,*) ' '
c         call exit(1)
      endif

      return
      end

      subroutine wrtlats1(ldefgrid, lcreate, ldefvar, ldefmisc, lwrite,
     &      lclose, lconcat, date, stat, varname, outfile, outftype, 
     &      fnclim, parmtabl, gtype, calendar, varcomm, filecomm, 
     &      model, center, amiss, nlon, nlat, nchunks, alons, alats, 
     &      nmon, mon1, iyr1, ibasemon, ibaseyr, idvar, array,
     &      chnkname, jyr1out)

      implicit none

      integer ldefgrid, ldefvar, lcreate, ldefmisc, lwrite, lclose, 
     &     lconcat, nlon, nlat, nchunks, nmon, mon1, iyr1, ibasemon,
     &     ibaseyr, idvar, jyr1out(*)

      real alons(*), alats(*), amiss, array(*)
      character*80 center, varcomm
      character*120 outfile, chnkname(nchunks), filecomm
c      character*256 parmtabl ! PJD Oceanonly - 1870-2014
      character*(*) parmtabl ! PJD Oceanonly - 1870-2014
      character*16 outftype, calendar, gtype, model
      character*15 fnclim
      character*(*) varname, stat, date

      print*, 'Your version of this code cannot write with the "lats"'
      print*, 'subroutine.  Only pp format or ascii is allowed.'

      call exit(1)

      return
      end


      subroutine wrtlats(ldefgrid, lcreate, ldefvar, ldefmisc, lwrite,
     &      lclose, lconcat, date, stat, varname, outfile, fnclim,
     &      outftype, parmtabl, gtype, calendar, varcomm, filecomm, 
     &      model, center, amiss, nlon, nlat, nchunks, alons, alats, 
     &      nmon, mon1, iyr1, ibasemon, ibaseyr, idvar, array,
     &      chnkname, jyr1out)

      implicit none
c      include '/usr/local/include/lats.inc' ! stargate
c      include '/pcmdi/ktaylor/rosinski/pcmdisst/lats.inc' ! zooks
c      include '/usr/local/lats/include/lats.inc' ! sunOS
c      include '/work/durack1/Shared/150219_AMIPForcingData/src/lats/lats.inc' ! PJD Oceanonly - 1870-2014
      include '/work/durack1/Shared/150219_AMIPForcingData/src/lats.inc' ! PJD Oceanonly - 1870-2014
      save idfile, igrid

      integer maxlon, maxlat
      parameter(maxlon=1440, maxlat=901)
      integer idvar, nlon, nlat, iyr1, nmon, mon1, ldefgrid, nchunks,
     &   lcreate, ldefvar, lwrite, lclose, ldefmisc, lconcat, 
     &   iendyr, jyr1out(*)
      real alons(*), alats(*), amiss, array(*)
      character*80 center, varcomm
      character*120 outfile, chnkname(nchunks), filecomm ! PJD Oceanonly - 1870-2014
      character*16 outftype, calendar, gtype, model
      character*15 fnclim
c      character*256 parmtabl ! PJD Oceanonly - 1870-2014
      character*(*) parmtabl ! PJD Oceanonly - 1870-2014
      character*(*) varname, stat, date

      integer latsgrid, latsvar, latscreate, latswrite,
     &    latsmissreal, latsclose, latsbasetime, latsparmtab
      logical caseindp

      real delta, day, bpp(19)
      integer igrid, igridtyp, i, ifiletyp, icalend, lastyr,
     &    inctime, istat, ilevel, n, nn, ierr1, nnn, iday, iend,
     &    ihour, idfile, imon, iyr, iday1, ihr1, m, mm,
     &    ibasemon, ibaseyr, ifreq, ichunk, n1, n2, lpp(45),
     &    llcre
      double precision alev, dlons(maxlon), dlats(maxlat)
      character gridtype*8, fullfile*120, ayr*16
      integer monlen(12,2)
      data monlen/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31,
     &            31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

      lastyr = (iyr1*12+mon1+nmon-2)/12

      if ((parmtabl(1:4) .ne. 'none') .and. 
     &    (parmtabl(1:7) .ne. 'default')) then

          ierr1 = latsparmtab(parmtabl)

          if (ierr1 .eq. 0) then
            print*, 'Error -11 in wrtlats '
            idvar = -11
            return
          endif
   
      endif

      if (ldefgrid .gt. 0) then

        print*, '    Defining grid with LATS ....'
        if ((nlon .gt. maxlon) .or. (nlat .gt. maxlat)) then

          print*, ' Error in subroutine wrtlats ....'
          print*, ' longitude and/or latitude dimension limits exceeded'
          print*, ' nlon = ', nlon, '  maxlon = ', maxlon
          print*, ' nlat = ', nlat, '  maxlat = ', maxlat

          idvar = -1
          return

        endif

        do 100 i=1,nlon
          dlons(i) = alons(i)
  100   continue
        do 150 i=1,nlat
          dlats(i) = alats(i)
  150   continue
  

        if (caseindp(gtype, 'gau')) then
          igridtyp = LATS_GAUSSIAN
          gridtype = 'gaussian'
        elseif (caseindp(gtype, 'cos')) then
          igridtyp = LATS_LINEAR
          gridtype = 'linear'
        elseif (caseindp(gtype, 'lmd') .or. caseindp(gtype, 'lmc')) then
          igridtyp = LATS_GENERIC
          gridtype = 'generic'
        elseif (caseindp(gtype, 'bmr') .or. caseindp(gtype, 'ccc') 
     +     .or. caseindp(gtype, 'cnr') .or. caseindp(gtype, 'col') 
     +     .or. caseindp(gtype, 'csi') .or. caseindp(gtype, 'der') 
     +     .or. caseindp(gtype, 'ecm') .or. caseindp(gtype, 'gfd') 
     +     .or. caseindp(gtype, 'mgo') .or. caseindp(gtype, 'mpi')
     +     .or. caseindp(gtype, 'nca') .or. caseindp(gtype, 'nmc') 
     +     .or. caseindp(gtype, 'rpn') .or. caseindp(gtype, 'sun') 
     +     .or. caseindp(gtype, 'uga') .or. caseindp(gtype, 'sng')
     +     .or. caseindp(gtype, 'gen')
     +     .or. caseindp(gtype, 'ccm') .or. caseindp(gtype, 'ccs')
     +     .or. caseindp(gtype, 'ech') .or. caseindp(gtype, 'nce')) then
          igridtyp = LATS_GAUSSIAN
          gridtype = 'gaussian'
        elseif (caseindp(gtype, 'dnm') .or. caseindp(gtype, 'gis') 
     +     .or. caseindp(gtype, 'csu') .or. caseindp(gtype, 'ucl')
     +     .or. caseindp(gtype, 'gla') .or. caseindp(gtype, 'gsf') 
     +     .or. caseindp(gtype, 'iap') .or. caseindp(gtype, 'mri') 
     +     .or. caseindp(gtype, 'uiu') 
     +     .or. caseindp(gtype, 'jma') .or. caseindp(gtype, 'nrl') 
     +     .or. caseindp(gtype, 'ukm') .or. caseindp(gtype, 'yon')) 
     +      then
          igridtyp = LATS_LINEAR
          gridtype = 'linear'
        else
          print*, 'Grid type not recognized in wrtlats '
          print*, 'Grid type = ', gtype
          write(9,*) 'Grid type not recognized in wrtlats '
          write(9,*) 'Grid type = ', gtype
          stop
        endif

        igrid = latsgrid(gridtype, igridtyp, nlon, dlons, nlat, dlats)

        if (igrid .eq. 0) then
          print*, 'Error in defining grid with latsgrid'
          idvar = -3
          return
        else
          print*, 'Grid defined successfully by LATS'
        endif

      endif

      if ((lwrite .gt. 0) .and. (lconcat .gt. 0)) then
        do 200 ichunk=1,nchunks
          OPEN (12+ichunk,FILE=chnkname(ichunk),STATUS='OLD',
     &           FORM='UNFORMATTED')
  200   continue
      endif

      if (nmon .gt. 0) then
        nn = nmon
      else
        nn = 1
      endif

      iyr = iyr1
      imon = mon1-1
      mm = 0

      do 1000 n=1,nn

        if (imon .eq. 12) then
          imon = 1
          iyr = iyr + 1
        else
          imon = imon+1
        endif

        llcre = 0
           
        if ((lcreate .gt. 0) .and. ((n .eq. 1) .or.  
     &      ((iyr .eq. jyr1out(mm+1)) .and. (imon .eq. 1)))) then

          mm = mm + 1

          if (caseindp(outftype, 'grads')) then 
            ifiletyp = LATS_GRADS_GRIB
          elseif (caseindp(outftype, 'grib')) then
            ifiletyp = LATS_GRIB_ONLY
          elseif (caseindp(outftype, 'cdf')) then
            ifiletyp = LATS_PCMDI
          elseif (caseindp(outftype, 'coards')) then
            ifiletyp = LATS_COARDS
          else
            idvar = -1
            print*, 'Error -1 in wrtlats '
            return
          endif

          if (caseindp(calendar, 'clim')) then
            icalend = LATS_CLIM
          elseif (caseindp(calendar, 'jul')) then
            icalend = LATS_JULIAN
          elseif (caseindp(calendar, 'no_')) then
            icalend = LATS_NOLEAP
            monlen(2,2) = 28
          else 
            icalend = LATS_STANDARD
          endif
  
          if ((caseindp(outftype, 'coards')) .and.
     &        (icalend .ne. LATS_STANDARD)) then
            print*, 'Error -2 in wrtlats '
            idvar = -2
            return
          endif
              
          if (caseindp(stat, 'average')) then
            istat = LATS_AVERAGE
          elseif (caseindp(stat, 'instant')) then
            istat = LATS_INSTANT
          endif

          if (date(1:4) .eq. 'none') then
            ifreq = LATS_FIXED
c              the following corrects a problem with lats.
            if (ifiletyp .eq. LATS_PCMDI) ifiletyp = LATS_COARDS
          else
            if (.not. caseindp(outftype, 'grads')) then
              ifreq = LATS_DAILY
            else
              ifreq = LATS_MONTHLY
            endif
          endif

          if (caseindp(outftype, 'grads') .or. (nmon .eq. 0)) then
            inctime = 1
          else
            inctime = 0
          endif

          alev = 0.0

c          iyr = 0 if no time-dimension; iyr=2 if climatology 
          if ((iyr .eq. 0) .or. (iyr .eq. 2)) then
            ayr = ''
            if (index(fnclim, 'clim') .gt. 0) ayr = '_' // fnclim 
          elseif (iyr .eq. lastyr) then
            write(ayr, '("_", i4)') iyr
          elseif (jyr1out(mm+1) .eq. -999) then
            write(ayr, '("_", i4, "_", i4)') iyr, lastyr
          elseif (jyr1out(mm) .eq. (jyr1out(mm+1)-1)) then  
            write(ayr, '("_", i4)') iyr
          else
            iendyr = jyr1out(mm+1) - 1
            iendyr = min0(iendyr, lastyr)
            write(ayr, '("_", i4, "_", i4)') iyr, iendyr
          endif
              
          iend = index(outfile, ' ') - 1
          fullfile = outfile(1:iend)//ayr

          print*, ' '
          print*, '     Creating ', outftype, ' file with LATS:'
          print*, '          ', fullfile        

          idfile = latscreate(fullfile, ifiletyp, icalend, 
     &         ifreq, inctime, center, model, filecomm)
          if (idfile .eq. 0) then
            print*, 'Error in creating file in LATS '
            idvar = -4
            return
          endif

          if (date(1:4) .ne. 'none') then
            iday1 = 1
            ihr1 = 0
            ierr1 = latsbasetime(idfile, ibaseyr, ibasemon, iday1, ihr1)
            if (ierr1 .eq. 0) then
              print*, 'Error in defining basetime in LATS'
              idvar = -9
              return
            endif
          endif
        
          llcre = 1

        endif

        if ((ldefvar .gt. 0) .and. 
     &           ((llcre .gt. 0) .or. (lcreate .eq. 0))) then

            print*, '    Defining a variable with LATS: ', varname

            ilevel = 0
            idvar = latsvar(idfile, varname, LATS_FLOAT, istat, igrid,
     &           ilevel,  varcomm)

            if (idvar .eq. 0) then
              print*, 'Error in defining variable in LATS'
              idvar = -5
              return
            endif
  
            if (ldefmisc .gt. 0) then
              delta = 1.0e-5*abs(amiss)
              ierr1 = latsmissreal(idfile, idvar, amiss, delta)

              if (ierr1 .eq. 0) then
                print*, 'Error in defining missing value in LATS '
                idvar = -6
                return
              endif
            endif

        endif

        if (lwrite .gt. 0) then

          if (lconcat .eq. 0) then 

            nnn = (n-1)*nlon*nlat + 1

          else

            nnn = 1

c           Process each chunk for this month
            n2 = 0
            do 300 ichunk = 1,nchunks

c             read chunk in pp-format
              read(12+ichunk)lpp,bpp
              n1 = n2+1
              n2 = n1+lpp(15)-1
              read(12+ichunk) (array(m), m=n1,n2)

  300       continue

          endif

          if ((nmon .eq. 0) .or. (date(1:4) .eq. 'none')) then

            iday = 0
            ihour = 0

          elseif (caseindp(outftype, 'grads')) then

            iday = 1
            ihour = 0

          elseif ((date(1:3) .eq. 'mid') .or. 
     &            (date(1:4) .eq. 'mean')) then

            if ((iyr .ne. 0) .and. (((mod(iyr,4) .eq. 0) .and. 
     &         (mod(iyr,100) .ne. 0)) .or. (mod(iyr,400) .eq. 0))) then
 

              day = monlen(imon,2)/2.0 + 1.01
              iday = day
              if ((day - iday) .gt. 0.1) then
                ihour = 12
              else
                ihour = 0
              endif

            else
              
              day = monlen(imon,1)/2.0 + 1.01
              iday = day
              if ((day - iday) .gt. 0.1) then
                ihour = 12
              else
                ihour = 0
              endif

            endif

          else

            print*, 'Unrecognized date indicator in wrtlats.'
            print*, 'You should pass either "mean" or "mid"'
            print*, 'to indicate whether data represent a'
            print*, 'monthly mean or a mid-month value.'
            idvar = -10
            return

          endif

            print*,'time out: ',iyr,imon,iday,ihour ! PJD Oceanonly 1870-2014

          ierr1 = latswrite(idfile, idvar, alev, iyr,
     &            imon, iday, ihour, array(nnn))

          if (ierr1 .eq. 0) then
            print*, 'Error in writing field for year ', iyr,
     &          ', month ', imon
            idvar = -7
            return
          endif

        endif
        
        if ((lclose .gt. 0) .and.
     &      ((n .ge. nmon) .or. 
     &       ((imon .eq. 12) .and. ((iyr+1) .eq. jyr1out(mm+1))))) then
          ierr1 = latsclose(idfile)
          if (ierr1 .eq. 0) then
            print*, 'Error in closing file in LATS'
            idvar = -8
            return
          else
            print*, 'File closed successfully by LATS: ', fullfile
            print*, ' '
          endif
        endif

 1000 continue

      if ((lwrite .gt. 0) .and. (lconcat .gt. 0)) then
        do 1100 ichunk=1,nchunks
          close(12+ichunk)
 1100   continue
      endif

      if (lwrite .gt. 0) then
        print*, 'Field written successfully by LATS: ', varname
      endif

      return
      end


c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     logical function caseindp
c     -------
c
c     Description:
c     -----------
c    This function tests for character string equivalence
c
c     Usage:
c     ------
c
c      if (caseindp(a, b)) then
c
c     ------
c     Date: 9/17/94
c     ----
c
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      logical function caseindp1(a, b)

c    This function tests for character string equivalence
c
c    compare what follows any leading blanks
c    check for equivalence of strings ignoring case.
c    only check through the last character of the shorter string 
c    (where "shorter" is computed after removing leading blanks)
c    The only exception to this is if 1 string is null or completely
c    blank; then equivalence requires both strings to be null or
c    completely blank.   

      character*(*) a, b
      character*1 c, d
      integer i, j, k, m, n, j1, k1, ii
      logical caseindp1 ! 150610 PJD added to suppress warnings

c     look for leading blanks and neglect

      j = len(a)
      k = len(b)

c     check whether both strings are length zero

      if ((j .eq. 0) .or. (k .eq. 0)) then
        if (j .eq. k) then
          caseindp1 = .true.
        else
          caseindp1 = .false.
        endif
        return
      endif


      j1 = 0
      do 10 i=1,j
        if (a(i:i) .ne. ' ') go to 20
        j1 = j1 + 1
   10 continue
   20 continue

      k1 = 0
      do 30 i=1,k
        if (b(i:i) .ne. ' ') go to 40
        k1 = k1 + 1
   30 continue
   40 continue

c   check whether either string is completely blank

      if ((j1 .eq. j) .or. (k1 .eq. k)) then
        if ((j1 .eq. j) .and. (k1 .eq. k)) then
          caseindp1 = .true.
        else
          caseindp1 = .false.
        endif
        return
      endif


c   look for trailing blanks and neglect
      
      j2 = j
      do 50 i=1,j
        if (a(j-i+1:j-i+1) .ne. ' ') go to 60
        j2 = j2 - 1
   50 continue
   60 continue

      k2 = k
      do 70 i=1,k
        if (b(k-i+1:k-i+1) .ne. ' ') go to 80
        k2 = k2 - 1
   70 continue
   80 continue

      j = j2 - j1
      k = k2 - k1
      ii = min0(j, k)
 
      do 100 i=1,ii
        m = i + j1
        n = i + k1
        c = a(m:m)
        if (lge(c,'A') .and. lle(c,'Z')) c = char(ichar(c)+32)
        d = b(n:n)
        if (lge(d,'A') .and. lle(d,'Z')) d = char(ichar(d)+32)
        if (c .ne. d) then
           caseindp1 = .false.
           return
        endif
  100 continue

      caseindp1 = .true. 
      return
      end
          


