c -*-Mode: Fortran;-*-
c  Module:      LATS FORTRAN include file
c 
c  Copyright:	1996, Regents of the University of California
c 		This software may not be distributed to others without
c 		permission of the author.
c 
c  Author:      Bob Drach, Lawrence Livermore National Laboratory
c               drach@llnl.gov
c 
c  Version:     $Id: lats.inc,v 1.10 1996/10/22 19:05:04 fiorino Exp $
c 
c  Revision History:
c 
c  $Log: lats.inc,v $
c Revision 1.10  1996/10/22  19:05:04  fiorino
c latsgrib bug in .ctl creator
c
c  Revision 1.9  1996/10/21 17:50:24  drach
c  - Declared LATS_MONTHLY_TABLE_COMP
c
c  Revision 1.8  1996/10/21 17:20:15  drach
c  - Fixed LATS_MONTHLY_TABLE_COMP option
c
c  Revision 1.7  1996/10/10 23:15:44  drach
c  - lats_create filetype changed to convention, with options LATS_PCMDI,
c    LATS_GRADS_GRIB, and LATS_COARDS.
c  - monthly data defaults to 16-bit compression
c  - LATS_MONTHLY_TABLE_COMP option added to override 16-bit compression
c  - AMIP II standard parameter file
c  - parameter file incorporates GRIB center and subcenter
c  - if time delta is positive, check that (new_time - old_time)=integer*delta
c
c  Revision 1.6  1996/08/27 19:41:05  drach
c  - Explicit definition of parameters as integers
c
c  Revision 1.5  1996/08/20 18:34:07  drach
c  - lats_create has a new argument: calendar
c  - lats_grid: longitude, latitude dimension vectors are now double
c    precision (double, C).
c  - lats_vert_dim: vector of levels is now double precision (double,
c    C). lats_vert_dim need not be called for single-value/surface
c    dimensions, if defined in the parameter table. Multi-valued vertical
c    dimensions, such as pressure levels, must be defined with
c    lats_vert_dim.
c  - lats_var: set level ID to 0 for implicitly defined surface
c    dimension.
c  - lats_write: level value is double precision (double, C).
c  - lats_parmtab: replaces routine lats_vartab.
c  - FORTRAN latserropt added: allows program to set error handling
c    options.
c  - The parameter file format changed.
c
c  Revision 1.4  1996/05/25 00:27:48  drach
c  - Added tables for vertical dimension types, time statistics, originating
c    centers, and quality control marks
c  - Modified signatures of lats_create and lats_vert_dim
c
c  Revision 1.3  1996/05/10 22:44:39  drach
c  - Initial version before GRIB driver added:
c  - Made grids, vertical dimensions file-independent
c
c  Revision 1.2  1996/05/04 01:11:10  drach
c  - Added name, units to lats_vert_dim
c  - Added missing data attribute (latsnc.c)
c
c  Revision 1.1  1996/05/03 18:48:08  drach
c  - Initial repository version
c
c 
c
      integer LATS_MAX_COMMENTS, LATS_MAX_LEVEL_DIMS, LATS_MAX_GRIDS 
      integer LATS_MAX_NAME, LATS_MAX_PARMS, LATS_MAX_PARM_LINE 
      integer LATS_MAX_PATH, LATS_MAX_RELUNITS, LATS_MAX_VAR_DIMS 
      integer LATS_LOW_LEVEL, LATS_MEDIUM_LEVEL, LATS_HIGH_LEVEL 
      integer LATS_NETCDF, LATS_GRIB, LATS_GAUSSIAN 
      integer LATS_LINEAR, LATS_GENERIC, LATS_SINGLE 
      integer LATS_INCREASING, LATS_DECREASING, LATS_HOURLY 
      integer LATS_DAILY, LATS_WEEKLY, LATS_MONTHLY 
      integer LATS_YEARLY, LATS_FIXED, LATS_AVERAGE 
      integer LATS_INSTANT, LATS_ACCUM, LATS_OTHER_TIME_STAT 
      integer LATS_FLOAT, LATS_INT, LATS_SINGLE_LEVEL 
      integer LATS_MULTI_LEVEL, LATS_STANDARD, LATS_JULIAN 
      integer LATS_NOLEAP, LATS_360, LATS_CLIM 
      integer LATS_CLIMLEAP, LATS_CLIM360, LATS_FATAL 
      integer LATS_VERBOSE, LATS_QC, LATS_PCMDI
      integer LATS_GRADS_GRIB, LATS_GRIB_ONLY
      integer LATS_COARDS, LATS_MONTHLY_TABLE_COMP
      integer LATS_FORECAST_HOURLY

c			     Max length of comment arguments
      parameter(LATS_MAX_COMMENTS = 256)
c			     Max number of level dimensions
      parameter(LATS_MAX_LEVEL_DIMS = 32)
c                            Max number of grids
      parameter(LATS_MAX_GRIDS = 32)
c			     Max name/units length
      parameter(LATS_MAX_NAME = 128)
c			     Max number of parameters
      parameter(LATS_MAX_PARMS = 500)
c			     Max characters in parameter file entry (one line)
      parameter(LATS_MAX_PARM_LINE = 512)
c			     Max file pathname length 
      parameter(LATS_MAX_PATH = 256)
c			     Max relative time units length
      parameter(LATS_MAX_RELUNITS = 64)
c			     Max number of dimensions in a variable
      parameter(LATS_MAX_VAR_DIMS = 4)

      parameter(LATS_LOW_LEVEL = 1, LATS_MEDIUM_LEVEL = 2,
     $     LATS_HIGH_LEVEL = 3)
      parameter(LATS_PCMDI = 1, LATS_GRIB_ONLY = 2,
     $     LATS_GRADS_GRIB = 3, LATS_COARDS = 4)
      parameter(LATS_NETCDF = 1, LATS_GRIB = 2)
      parameter(LATS_GAUSSIAN = 1, LATS_LINEAR = 2, LATS_GENERIC = 3)
      parameter(LATS_SINGLE = 1, LATS_INCREASING = 2,
     $     LATS_DECREASING = 3)
      parameter(LATS_HOURLY = 1, LATS_DAILY = 2, LATS_WEEKLY = 3,
     $     LATS_MONTHLY = 4, LATS_YEARLY = 5, LATS_FIXED = 6,
     $     LATS_MONTHLY_TABLE_COMP = 7, LATS_FORECAST_HOURLY = 8)
      parameter(LATS_AVERAGE = 1, LATS_INSTANT = 2, LATS_ACCUM = 3,
     $     LATS_OTHER_TIME_STAT = 4)
      parameter(LATS_FLOAT = 1, LATS_INT = 2)
      parameter(LATS_SINGLE_LEVEL = 1, LATS_MULTI_LEVEL = 2)

      parameter(LATS_STANDARD = X'00001111')
      parameter(LATS_JULIAN = X'00011111')
      parameter(LATS_NOLEAP = X'00001011')
      parameter(LATS_360 = X'00000011')
      parameter(LATS_CLIM = X'00001000')
      parameter(LATS_CLIMLEAP = X'00001100')
      parameter(LATS_CLIM360 = X'00000000')
      
      parameter(LATS_FATAL = X'00000001')
      parameter(LATS_VERBOSE = X'0000002')
      parameter(LATS_QC = X'00000001')
 
