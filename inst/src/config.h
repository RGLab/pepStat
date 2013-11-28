/* src/config.h.  Generated from config.h.in by configure.  */
#ifndef CONFIG_H
#define CONFIG_H

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "http://www.rglab.org"

/* Define to the full name of this package. */
#define PACKAGE_NAME "pepStat"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "pepStat 0.0.1"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "pepStat"

/* Define to the home page for this package. */
/* #undef PACKAGE_URL */

/* Define to the version of this package. */
#define PACKAGE_VERSION "0.0.1"

/* Define to 1 if you have the <dispatch/dispatch.h> header file. */
/* #undef HAVE_DISPATCH_DISPATCH_H */

/* Define to 1 if you have the <omp.h> header file. */
#define HAVE_OMP_H 1

#ifdef HAVE_DISPATCH_DISPATCH_H
  #include <dispatch/dispatch.h>
  #define DO_NORMALIZE(task, counter_name, n_times) \
    dispatch_apply(n_times, \
                   dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0),	\
                   ^(size_t counter_name) { task; });
  
#else // not HAVE_DISPATCH_DISPATCH_H
  #ifdef HAVE_OMP_H
    #include <omp.h>
    #define DO_NORMALIZE(task, counter_name, n_times) \
      _Pragma("omp parallel for")                     \
      for (int counter_name = 0; counter_name < n_times; ++counter_name) { \
        task; \
      }
  #else // not HAVE_OMP_H
	//    # include <foobar.h>
    #define DO_NORMALIZE(task, counter_name, n_times) \
      for (int counter_name = 0; counter_name < n_times; ++counter_name) { \
        task; \
      }
  #endif // HAVE_OMP_H
#endif // HAVE_DISPATCH_DISPATCH_H

#endif
