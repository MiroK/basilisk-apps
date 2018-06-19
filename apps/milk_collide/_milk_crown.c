#line 1 "milk_crown-cpp.c"
#line 1 "<built-in>"
#line 1 "<command-line>"
#line 1 "/usr/include/stdc-predef.h"
#line 1 "<command-line>"
#line 1 "milk_crown-cpp.c"
#if _XOPEN_SOURCE < 700
#undef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif
#if _GNU_SOURCE
#include <stdint.h>
#include <string.h>
#include <fenv.h>
#endif



#line 1 "common.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/common.h"
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#if _CADNA
# include <cadna.h>
#endif

#if __cplusplus
# define delete delete_qcc
# define right right_qcc
# define left left_qcc
# define norm norm_qcc
# define new new_qcc
#endif

#define pi 3.14159265358979
#undef HUGE
#define HUGE ((double)1e30)
#define nodata HUGE
#define _NVARMAX 65536
#define is_constant(v) ((v).i >= _NVARMAX)
#define constant(v) (is_constant(v) ? _constant[(v).i - _NVARMAX] : nodata)

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))
#define sq(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))
#define sign(x) ((x) > 0 ? 1 : -1)
#define noise() (1. - 2.*rand()/(double)RAND_MAX)
#define clamp(x,a,b) ((x) < (a) ? (a) : (x) > (b) ? (b) : (x))
#define swap(type,a,b) { type tmp = a; a = b; b = tmp; }
#define unmap(x,y)

#define trash(x)


#define systderr stderr
#define systdout stdout

#if 1
static FILE * qstderr (void);
static FILE * qstdout (void);
FILE * ferr, * fout;
#define not_mpi_compatible()\
do {\
  if (npe() > 1) {\
    fprintf (ferr, "%s() is not compatible with MPI (yet)\n", __func__);\
    exit (1);\
  }\
} while(0)\

#line 62

#else
# define qstderr() stderr
# define qstdout() stdout
# define ferr stderr
# define fout stdout
# define not_mpi_compatible()
#endif



#if MTRACE

struct {
  FILE * fp;
  size_t total, max;
  size_t overhead, maxoverhead;
  size_t nr;
  size_t startrss, maxrss;
  char * fname;
} pmtrace;

typedef struct {
  char * func, * file;
  size_t max, total;
  int line, id;
} pmfunc;

typedef struct {
  size_t id, size;
} pmdata;

static pmfunc * pmfuncs = NULL;
static int pmfuncn = 0;

#define sysmalloc malloc
#define syscalloc calloc
#define sysrealloc realloc
#define sysfree free
#define systrdup strdup

static int pmfunc_index (const char * func, const char * file, int line)
{
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn; i++, p++)
    if (p->line == line && !strcmp(func, p->func) && !strcmp(file, p->file))
      return p->id;
  pmfuncn++;
  pmfuncs = (pmfunc *) sysrealloc (pmfuncs, pmfuncn*sizeof(pmfunc));
  p = &pmfuncs[pmfuncn - 1];
  memset (p, 0, sizeof(pmfunc));
  p->func = systrdup(func);
  p->file = systrdup(file);
  p->line = line;
  p->id = pmfuncn;
  if (pmtrace.fp)
    fprintf (pmtrace.fp, "@ %d %s %s %d\n", pmfuncn, func, file, line);
  return pmfuncn;
}

static void pmfunc_trace (pmfunc * f, char c)
{
  if (pmtrace.fp)
    fprintf (pmtrace.fp, "%c %d %ld %ld %ld",
      c, f->id, pmtrace.nr, pmtrace.total, f->total);
#if _GNU_SOURCE
  if (pmtrace.nr % 1 == 0) {
    struct rusage usage;
    getrusage (RUSAGE_SELF, &usage);
    if (pmtrace.fp)
      fprintf (pmtrace.fp, " %ld", usage.ru_maxrss*1024);
    if (!pmtrace.nr)
      pmtrace.startrss = usage.ru_maxrss;
    if (usage.ru_maxrss > pmtrace.maxrss)
      pmtrace.maxrss = usage.ru_maxrss;
  }
#endif
  if (pmtrace.fp)
    fputc ('\n', pmtrace.fp);
  pmtrace.nr++;
}

static void * pmfunc_alloc (pmdata * d, size_t size,
       const char * func, const char * file, int line,
       char c)
{
  assert (d != NULL);
  d->id = pmfunc_index(func, file, line);
  d->size = size;
  pmfunc * f = &pmfuncs[d->id - 1];
  f->total += size;
  if (f->total > f->max)
    f->max = f->total;
  pmtrace.total += size;
  pmtrace.overhead += sizeof(pmdata);
  if (pmtrace.total > pmtrace.max) {
    pmtrace.max = pmtrace.total;
    pmtrace.maxoverhead = pmtrace.overhead;
  }
  pmfunc_trace (f, c);
  return ((char *)d) + sizeof(pmdata);
}

static void * pmfunc_free (void * ptr, char c)
{
  if (!ptr)
    return ptr;
  pmdata * d = (pmdata *) (((char *)ptr) - sizeof(pmdata));
  if (d->id < 1 || d->id > pmfuncn) {
    fputs ("*** MTRACE: ERROR!: corrupted free()", qstderr());
    if (d->size == 0)
      fputs (", possible double free()", qstderr());
    else
      fputs (", not traced?", qstderr());
    fputs (", aborting...\n", qstderr());
    abort();
    return ptr;
  }
  else {
    pmfunc * f = &pmfuncs[d->id - 1];
    if (f->total < d->size) {
      fprintf (qstderr(), "*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\n",
        f->total, d->size);
      abort();
    }
    else
      f->total -= d->size;
    if (pmtrace.total < d->size) {
      fprintf (qstderr(), "*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\n",
        pmtrace.total, d->size);
      abort();
    }
    else {
      pmtrace.total -= d->size;
      pmtrace.overhead -= sizeof(pmdata);
    }
    d->id = 0;
    d->size = 0;
    pmfunc_trace (f, c);
    return d;
  }
}

static void * pmalloc (size_t size,
         const char * func, const char * file, int line)
{
  return pmfunc_alloc ((pmdata *) sysmalloc (sizeof(pmdata) + size),
         size, func, file, line, '+');
}

static void * pcalloc (size_t nmemb, size_t size,
         const char * func, const char * file, int line)
{
  void * p = pmalloc (nmemb*size, func, file, line);
  return memset (p, 0, nmemb*size);
}

static void * prealloc (void * ptr, size_t size,
   const char * func, const char * file, int line)
{
  return pmfunc_alloc ((pmdata *) sysrealloc (pmfunc_free(ptr, '<'),
           sizeof(pmdata) + size),
         size, func, file, line, '>');
}

static void pfree (void * ptr,
     const char * func, const char * file, int line)
{
  sysfree (pmfunc_free (ptr, '-'));
}

static char * pstrdup (const char * s,
         const char * func, const char * file, int line)
{
  char * d = (char *) pmalloc (strlen(s) + 1, func, file, line);
  return strcpy (d, s);
}

#if MTRACE < 3
static int pmaxsort (const void * a, const void * b) {
  const pmfunc * p1 = a, * p2 = b;
  return p1->max < p2->max;
}
#endif

static int ptotalsort (const void * a, const void * b) {
  const pmfunc * p1 = (const pmfunc *) a, * p2 = (const pmfunc *) b;
  return p1->total < p2->total;
}

static void pmfuncs_free()
{
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn; i++, p++) {
    sysfree (p->func);
    sysfree (p->file);
  }
  sysfree (pmfuncs);
}

void pmuntrace (void)
{
#if MTRACE < 3
  fprintf (qstderr(),
    "*** MTRACE: max resident  set size: %10ld bytes\n"
    "*** MTRACE: max traced memory size: %10ld bytes"
    " (tracing overhead %.1g%%)\n"
    "%10s    %20s   %s\n",
    pmtrace.maxrss*1024,
    pmtrace.max, pmtrace.maxoverhead*100./pmtrace.max,
    "max bytes", "function", "file");
  qsort (pmfuncs, pmfuncn, sizeof(pmfunc), pmaxsort);
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn && p->max > 0; i++, p++)
    fprintf (qstderr(), "%10ld    %20s   %s:%d\n",
      p->max, p->func, p->file, p->line);

  if (pmtrace.fp) {
    char * fname = pmtrace.fname, * s;
    while ((s = strchr(fname,'/')))
      fname = s + 1;

    fputs ("load(\"`echo $BASILISK`/mtrace.plot\")\n", pmtrace.fp);
    fprintf (pmtrace.fp,
      "plot '%s' u 3:($6-%g) w l t 'ru_maxrss - %.3g',"
      "total(\"%s\") w l t 'total'",
      fname,
      pmtrace.startrss*1024.,
      pmtrace.startrss*1024.,
      fname);
    pmfunc * p = pmfuncs;
    for (int i = 0; i < pmfuncn && p->max > 0.01*pmtrace.max; i++, p++)
      fprintf (pmtrace.fp,
        ",func(\"%s\",%d) w l t '%s'",
        fname, p->id, p->func);
    fputc ('\n', pmtrace.fp);
    fprintf (qstderr(),
      "*** MTRACE: To get a graph use: tail -n 2 %s | gnuplot -persist\n",
      fname);
    fclose (pmtrace.fp);
    pmtrace.fp = NULL;
    sysfree (pmtrace.fname);
  }
#endif

  if (pmtrace.total > 0) {
    qsort (pmfuncs, pmfuncn, sizeof(pmfunc), ptotalsort);
    pmfunc * p = pmfuncs;
    for (int i = 0; i < pmfuncn && p->total > 0; i++, p++)
      fprintf (qstderr(), "%s:%d: error: %ld bytes leaked here\n",
        p->file, p->line, p->total);
    pmfuncs_free();
    exit(1);
  }
  else {
#if MTRACE < 3
    fputs ("*** MTRACE: No memory leaks\n", qstderr());
#endif
    pmfuncs_free();
  }
}

#else
# define pmalloc(s,func,file,line) malloc(s)
# define pcalloc(n,s,func,file,line) calloc(n,s)
# define prealloc(p,s,func,file,line) realloc(p,s)
# define pfree(p,func,file,line) free(p)
# define pstrdup(s,func,file,line) strdup(s)
#endif







typedef struct {
  void * p;
  long max, len;
} Array;

Array * array_new()
{
  Array * a = ((Array *) pmalloc ((1)*sizeof(Array),__func__,__FILE__,__LINE__));
  a->p = NULL;
  a->max = a->len = 0;
  return a;
}

void array_free (Array * a)
{
  if (a->max > 0)
    pfree (a->p,__func__,__FILE__,__LINE__);
  pfree (a,__func__,__FILE__,__LINE__);
}

void array_append (Array * a, void * elem, size_t size)
{
  if (a->len + size >= a->max) {
    a->max += max (size, 4096);
    a->p = prealloc (a->p, a->max,__func__,__FILE__,__LINE__);
  }
  memcpy (((char *)a->p) + a->len, elem, size);
  a->len += size;
}

void * array_shrink (Array * a)
{
  void * p = prealloc (a->p, a->len,__func__,__FILE__,__LINE__);
  pfree (a,__func__,__FILE__,__LINE__);
  return p;
}



#if 2 == 1
#include <extrae_user_events.h>

typedef struct {
  Array index, stack;
  extrae_type_t type;
} Trace;

Trace trace_func = {
  {NULL, 0, 0}, {NULL, 0, 0},
  60000010,
};

Trace trace_mpi_func = {
  {NULL, 0, 0}, {NULL, 0, 0},
  60000011,
};

static int lookup_func (Array * a, const char * func)
{
  for (int i = 0; i < a->len/sizeof(char *); i++) {
    char * s = ((char **)a->p)[i];
    if (!strcmp (func, s))
      return i + 1;
  }
  char * s = pstrdup (func,__func__,__FILE__,__LINE__);
  array_append (a, &s, sizeof(char *));
  return a->len;
}

static void trace_push (Trace * t, const char * func)
{
  int value = lookup_func (&t->index, func);
  Extrae_eventandcounters (t->type, value);
  array_append (&t->stack, &value, sizeof(int));
}

static void trace_pop (Trace * t, const char * func)
{
  assert (t->stack.len > 0);
  t->stack.len -= sizeof(int);
  int value = t->stack.len > 0 ?
    ((int *)t->stack.p)[t->stack.len/sizeof(int) - 1] : 0;
  Extrae_eventandcounters (t->type, value);
}

static void trace_define (Trace * t, char * description)
{
  if (t->index.len > 0) {
    extrae_value_t values[t->index.len/sizeof(char *) + 1];
    char * names[t->index.len/sizeof(char *) + 1],
      ** func = (char **) t->index.p;
    names[0] = "OTHER";
    values[0] = 0;
    unsigned len = 1;
    for (int i = 0; i < t->index.len/sizeof(char *); i++, func++) {
      names[len] = *func;
      values[len++] = i + 1;
    }
    Extrae_define_event_type (&t->type, description, &len, values, names);
  }
}

static void trace_free (Trace * t)
{
  char ** func = (char **) t->index.p;
  for (int i = 0; i < t->index.len/sizeof(char *); i++, func++)
    pfree (*func,__func__,__FILE__,__LINE__);
  pfree (t->index.p,__func__,__FILE__,__LINE__);
  pfree (t->stack.p,__func__,__FILE__,__LINE__);
}

static void trace_off()
{
  trace_define (&trace_func, "Basilisk functions");
  trace_define (&trace_mpi_func, "Basilisk functions (MPI-related)");
  trace_free (&trace_func);
  trace_free (&trace_mpi_func);
}






# define trace(func, file, line) trace_push (&trace_func, func)
# define end_trace(func, file, line) trace_pop (&trace_func, func)

#elif 2

typedef struct {
  char * func, * file;
  int line, calls;
  double total, self;
} TraceIndex;

struct {
  Array stack, index;
  double t0;
} Trace = {
  {NULL, 0, 0}, {NULL, 0, 0},
  -1
};

static void trace_add (const char * func, const char * file, int line,
         double total, double self)
{
  TraceIndex * t = (TraceIndex *) Trace.index.p;
  int i, len = Trace.index.len/sizeof(TraceIndex);
  for (i = 0; i < len; i++, t++)
    if (t->line == line && !strcmp (func, t->func) && !strcmp (file, t->file))
      break;
  if (i == len) {
    TraceIndex t = {pstrdup(func,__func__,__FILE__,__LINE__), pstrdup(file,__func__,__FILE__,__LINE__), line, 1, total, self};
    array_append (&Trace.index, &t, sizeof(TraceIndex));
  }
  else
    t->calls++, t->total += total, t->self += self;
}

static void trace (const char * func, const char * file, int line)
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  if (Trace.t0 < 0)
    Trace.t0 = tv.tv_sec + tv.tv_usec/1e6;
  double t[2] = {(tv.tv_sec - Trace.t0) + tv.tv_usec/1e6, 0.};
  array_append (&Trace.stack, t, 2*sizeof(double));




}

static void end_trace (const char * func, const char * file, int line)
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  double te = (tv.tv_sec - Trace.t0) + tv.tv_usec/1e6;
  double * t = (double *) Trace.stack.p;
  assert (Trace.stack.len >= 2*sizeof(double));
  t += Trace.stack.len/sizeof(double) - 2;
  Trace.stack.len -= 2*sizeof(double);
  double dt = te - t[0];




  trace_add (func, file, line, dt, dt - t[1]);
  if (Trace.stack.len >= 2*sizeof(double)) {
    t -= 2;
    t[1] += dt;
  }
}

static int compar_self (const void * p1, const void * p2)
{
  const TraceIndex * t1 = p1, * t2 = p2;
  return t1->self < t2->self;
}

static void trace_off()
{
  int i, len = Trace.index.len/sizeof(TraceIndex);
  double total = 0.;
  TraceIndex * t;
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    total += t->self;
  qsort (Trace.index.p, len, sizeof(TraceIndex), compar_self);
  fprintf (qstdout(), "   calls    total     self   %% total   function\n");
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++) {
    fprintf (qstdout(), "%8d   %6.2f   %6.2f     %4.1f%%   %s():%s:%d\n",
      t->calls, t->total, t->self, t->self*100./total,
      t->func, t->file, t->line);
    pfree (t->func,__func__,__FILE__,__LINE__); pfree (t->file,__func__,__FILE__,__LINE__);
  }

  pfree (Trace.index.p,__func__,__FILE__,__LINE__);
  Trace.index.p = NULL;
  Trace.index.len = Trace.index.max = 0;

  pfree (Trace.stack.p,__func__,__FILE__,__LINE__);
  Trace.stack.p = NULL;
  Trace.stack.len = Trace.stack.max = 0;
}

#else
# define trace(...)
# define end_trace(...)
#endif



#if _OPENMP

#include <omp.h>
#define OMP(x) _Pragma(#x)
#define tid() omp_get_thread_num()
#define pid() 0
#define npe() omp_get_num_threads()
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_double(v,op)

#elif 1

#include <mpi.h>
#define OMP(x)

static bool in_prof = false;
static double prof_start, _prof;
#define prof_start(name)\
  assert (!in_prof); in_prof = true;\
  prof_start = MPI_Wtime();\

#line 590

#define prof_stop()\
  assert (in_prof); in_prof = false;\
  _prof = MPI_Wtime();\
  mpi_time += _prof - prof_start;\

#line 595


#if FAKE_MPI
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_double(v,op)
#else

int mpi_all_reduce0 (void *sendbuf, void *recvbuf, int count,
       MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{ trace ("mpi_all_reduce0", "/home/miro3/Documents/Programming/basilisk/src/common.h", 604);
  { int _ret =  MPI_Allreduce (sendbuf, recvbuf, count, datatype, op, comm); end_trace("mpi_all_reduce0", "/home/miro3/Documents/Programming/basilisk/src/common.h", 605);  return _ret; }
 end_trace("mpi_all_reduce0", "/home/miro3/Documents/Programming/basilisk/src/common.h", 606); }
#define mpi_all_reduce(v,type,op) {\
  prof_start ("mpi_all_reduce");\
  union { int a; float b; double c;} global;\
  mpi_all_reduce0 (&(v), &global, 1, type, op, MPI_COMM_WORLD);\
  memcpy (&(v), &global, sizeof (v));\
  prof_stop();\
}\

#line 614

#define mpi_all_reduce_double(v,op) {\
  prof_start ("mpi_all_reduce");\
  double global, tmp = v;\
  mpi_all_reduce0 (&tmp, &global, 1, MPI_DOUBLE, op, MPI_COMM_WORLD);\
  v = global;\
  prof_stop();\
}\

#line 622


#endif

static int mpi_rank, mpi_npe;
#define tid() mpi_rank
#define pid() mpi_rank
#define npe() mpi_npe

#define QFILE FILE

static FILE * qstderr (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "log-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systderr;
  }
  return fp;
}

static FILE * qstdout (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "out-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systdout;
  }
  return fp;
}

static void finalize (void)
{
  MPI_Finalize();
}

void mpi_init()
{
  int initialized;
  MPI_Initialized (&initialized);
  if (!initialized) {
    MPI_Init (NULL, NULL);
    MPI_Comm_set_errhandler (MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
    atexit (finalize);
    MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size (MPI_COMM_WORLD, &mpi_npe);
    srand (mpi_rank + 1);
    if (mpi_rank > 0) {
      ferr = fopen ("/dev/null", "w");
      fout = fopen ("/dev/null", "w");
    }
    else {
      ferr = qstderr();
      fout = qstdout();
    }
    char * etrace = getenv ("MALLOC_TRACE"), name[80];
    if (etrace && mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      setenv ("MALLOC_TRACE", name, 1);
    }
#if MTRACE == 1
    etrace = getenv ("MTRACE");
    if (!etrace)
      etrace = "mtrace";
    if (mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      pmtrace.fp = fopen (name, "w");
      pmtrace.fname = systrdup(name);
    }
    else {
      pmtrace.fp = fopen (etrace, "w");
      pmtrace.fname = systrdup(etrace);
    }
#endif
  }
}

#else

#define OMP(x)
#define tid() 0
#define pid() 0
#define npe() 1
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_double(v,op)

#endif

void init_solver()
{
#if _CADNA
  cadna_init (-1);
#endif
#if 1
  mpi_init();
#elif MTRACE == 1
  char * etrace = getenv ("MTRACE");
  pmtrace.fp = fopen (etrace ? etrace : "mtrace", "w");
  pmtrace.fname = systrdup (etrace ? etrace : "mtrace");
#endif
}

#define OMP_PARALLEL() OMP(omp parallel)

#define NOT_UNUSED(x) (void)(x)

#define VARIABLES ;
#define val(a,k,l,m) data(k,l,m)[a.i]

double _val_higher_dimension = 0.;
#define _val_higher_dimension(x,a,b,c) _val_higher_dimension
#line 751 "/home/miro3/Documents/Programming/basilisk/src/common.h"
#if (_GNU_SOURCE || __APPLE__) && !_OPENMP && !_CADNA
double undefined;
# if __APPLE__
# include <stdint.h>
# include "fp_osx.h"
# endif
# define enable_fpe(flags) feenableexcept (flags)
# define disable_fpe(flags) fedisableexcept (flags)
static void set_fpe (void) {
  int64_t lnan = 0x7ff0000000000001;
  assert (sizeof (int64_t) == sizeof (double));
  memcpy (&undefined, &lnan, sizeof (double));
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}
#else
# define undefined ((double) DBL_MAX)
# define enable_fpe(flags)
# define disable_fpe(flags)
static void set_fpe (void) {}
#endif


typedef struct {
  long n;
  long tn;
  int depth;
  int maxdepth;
} Grid;
Grid * grid = NULL;

double X0 = 0., Y0 = 0., Z0 = 0.;

double L0 = 1.;




int N = 16;


typedef struct { int i; } scalar;

typedef struct {
  scalar x;

  scalar y;


  scalar z;

} vector;

typedef struct {
  vector x;

  vector y;


  vector z;

} tensor;

typedef struct {
  double x, y, z;
} coord;
#line 828 "/home/miro3/Documents/Programming/basilisk/src/common.h"
void normalize (coord * n)
{
  double norm = 0.;
  {
#line 831

    norm += sq(n->x);
#line 831

    norm += sq(n->y);
#line 831

    norm += sq(n->z);}
  norm = sqrt(norm);
  {
#line 834

    n->x /= norm;
#line 834

    n->y /= norm;
#line 834

    n->z /= norm;}
}

struct _origin { double x, y, z; };

void origin (struct _origin p) {
  X0 = p.x; Y0 = p.y; Z0 = p.z;
}

void size (double L) {
  L0 = L;
}

double zero (double s0, double s1, double s2) { return 0.; }
#line 857 "/home/miro3/Documents/Programming/basilisk/src/common.h"
  enum { right, left, top, bottom, front, back };

int nboundary = 2*3;



#define dirichlet(x) (2.*(x) - val(_s,0,0,0))
#define dirichlet_homogeneous() (- val(_s,0,0,0))
#define neumann(x) (Delta*(x) + val(_s,0,0,0))
#define neumann_homogeneous() (val(_s,0,0,0))

double * _constant = NULL;
extern size_t datasize;
typedef struct _Point Point;

#line 1 "grid/boundaries.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/grid/boundaries.h"


typedef struct _Boundary Boundary;

struct _Boundary {
  void (* destroy) (Boundary * b);
  void (* level) (const Boundary * b, scalar * list, int l);

  void (* restriction) (const Boundary * b, scalar * list, int l);
};

static Boundary ** boundaries = NULL;

void add_boundary (Boundary * b) {
  int len = 0;
  if (boundaries) {
    Boundary ** i = boundaries;
    while (*i++) len++;
  }
  boundaries = (Boundary * *) prealloc (boundaries, (len + 2)*sizeof(Boundary *),__func__,__FILE__,__LINE__);
  boundaries[len] = b;
  boundaries[len+1] = NULL;
}

void free_boundaries() {
  if (!boundaries)
    return;
  Boundary ** i = boundaries, * b;
  while ((b = *i++))
    if (b->destroy)
      b->destroy (b);
    else
      pfree (b,__func__,__FILE__,__LINE__);
  pfree (boundaries,__func__,__FILE__,__LINE__);
  boundaries = NULL;
}
#line 47 "/home/miro3/Documents/Programming/basilisk/src/grid/boundaries.h"
typedef struct {
  Boundary parent;
  int d;
} BoxBoundary;
#line 873 "/home/miro3/Documents/Programming/basilisk/src/common.h"



typedef struct {

#line 878 "/home/miro3/Documents/Programming/basilisk/src/common.h"

  double (** boundary) (Point, Point, scalar);
  double (** boundary_homogeneous) (Point, Point, scalar);
  double (* gradient) (double, double, double);
  void (* delete) (scalar);
  char * name;
  struct {
    int x;

    int y;


    int z;

  } d;
  vector v;
  bool face, nodump;

#line 17 "/home/miro3/Documents/Programming/basilisk/src/grid/multigrid-common.h"

  void (* prolongation) (Point, scalar);
  void (* restriction) (Point, scalar);

#line 8 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h"

  void (* refine) (Point, scalar);

#line 87 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h"

  void (* coarsen) (Point, scalar);

#line 76 "./fractions.h"

  vector n;

#line 26 "/home/miro3/Documents/Programming/basilisk/src/vof.h"

  scalar * tracers;
  bool inverse;

#line 20 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"

  scalar phi;

#line 21 "/home/miro3/Documents/Programming/basilisk/src/tension.h"

  double sigma;

} _Attributes;
_Attributes * _attribute;
#line 876 "/home/miro3/Documents/Programming/basilisk/src/common.h"























int list_len (scalar * list)
{
  if (!list) return 0;
  int ns = 0;
  if (list) for (scalar s = *list, *_i0 = list; ((scalar *)&s)->i >= 0; s = *++_i0) ns++;
  return ns;
}

scalar * list_append (scalar * list, scalar s)
{
  int len = list_len (list);
  list = (scalar *) prealloc (list, (len + 2)*sizeof(scalar),__func__,__FILE__,__LINE__);
  list[len] = s;
  list[len + 1].i = -1;
  return list;
}

scalar * list_add (scalar * list, scalar s)
{
  if (list) for (scalar t = *list, *_i1 = list; ((scalar *)&t)->i >= 0; t = *++_i1)
    if (t.i == s.i)
      return list;
  return list_append (list, s);
}

int list_lookup (scalar * l, scalar s)
{
  if (l != NULL)
    if (l) for (scalar s1 = *l, *_i2 = l; ((scalar *)&s1)->i >= 0; s1 = *++_i2)
      if (s1.i == s.i)
 return true;
  return false;
}

scalar * list_copy (scalar * l)
{
  scalar * list = NULL;
  if (l != NULL)
    if (l) for (scalar s = *l, *_i3 = l; ((scalar *)&s)->i >= 0; s = *++_i3)
      list = list_append (list, s);
  return list;
}

scalar * list_concat (scalar * l1, scalar * l2)
{
  scalar * l3 = list_copy (l1);
  if (l2) for (scalar s = *l2, *_i4 = l2; ((scalar *)&s)->i >= 0; s = *++_i4)
    l3 = list_append (l3, s);
  return l3;
}

void list_print (scalar * l, FILE * fp)
{
  int i = 0;
  if (l) for (scalar s = *l, *_i5 = l; ((scalar *)&s)->i >= 0; s = *++_i5)
    fprintf (fp, "%s%s", i++ == 0 ? "{" : ",", _attribute[s.i].name);
  fputs (i > 0 ? "}\n" : "{}\n", fp);
}

int vectors_len (vector * list)
{
  if (!list) return 0;
  int nv = 0;
  if (list) for (vector v = *list, *_i6 = list; ((scalar *)&v)->i >= 0; v = *++_i6) nv++;
  return nv;
}

vector * vectors_append (vector * list, vector v)
{
  int len = vectors_len (list);
  list = (vector *) prealloc (list, (len + 2)*sizeof(vector),__func__,__FILE__,__LINE__);
  list[len] = v;
  list[len + 1] = (vector){{-1}};
  return list;
}

vector * vectors_add (vector * list, vector v)
{
  if (list) for (vector w = *list, *_i7 = list; ((scalar *)&w)->i >= 0; w = *++_i7) {
    bool id = true;
    {
#line 979

      if (w.x.i != v.x.i)
 id = false;
#line 979

      if (w.y.i != v.y.i)
 id = false;
#line 979

      if (w.z.i != v.z.i)
 id = false;}
    if (id)
      return list;
  }
  return vectors_append (list, v);
}

vector * vectors_copy (vector * l)
{
  vector * list = NULL;
  if (l != NULL)
    if (l) for (vector v = *l, *_i8 = l; ((scalar *)&v)->i >= 0; v = *++_i8)
      list = vectors_append (list, v);
  return list;
}

vector * vectors_from_scalars (scalar * s)
{
  vector * list = NULL;
  while (s->i >= 0) {
    vector v;
    {
#line 1002
 {
      assert (s->i >= 0);
      v.x = *s++;
    }
#line 1002
 {
      assert (s->i >= 0);
      v.y = *s++;
    }
#line 1002
 {
      assert (s->i >= 0);
      v.z = *s++;
    }}
    list = vectors_append (list, v);
  }
  return list;
}

int tensors_len (tensor * list)
{
  if (!list) return 0;
  int nt = 0;
  if (list) for (tensor t = *list, *_i9 = list; ((scalar *)&t)->i >= 0; t = *++_i9) nt++;
  return nt;
}

tensor * tensors_append (tensor * list, tensor t)
{
  int len = tensors_len (list);
  list = (tensor *) prealloc (list, (len + 2)*sizeof(tensor),__func__,__FILE__,__LINE__);
  list[len] = t;
  list[len + 1] = (tensor){{{-1}}};
  return list;
}

tensor * tensors_from_vectors (vector * v)
{
  tensor * list = NULL;
  while (v->x.i >= 0) {
    tensor t;
    {
#line 1033
 {
      assert (v->x.i >= 0);
      t.x = *v++;
    }
#line 1033
 {
      assert (v->y.i >= 0);
      t.y = *v++;
    }
#line 1033
 {
      assert (v->z.i >= 0);
      t.z = *v++;
    }}
    list = tensors_append (list, t);
  }
  return list;
}

scalar * all = NULL;



scalar (* init_scalar) (scalar, const char *);
vector (* init_vector) (vector, const char *);
tensor (* init_tensor) (tensor, const char *);
vector (* init_face_vector) (vector, const char *);





typedef struct _Event Event;
typedef int (* Expr) (int *, double *, Event *);

struct _Event {
  int last, nexpr;
  int (* action) (const int, const double, Event *);
  Expr expr[3];
  int * arrayi;
  double * arrayt;
  char * file;
  int line;
  char * name;
  double t;
  int i, a;
  void * data;
  Event * next;
};

static Event * Events = NULL;

int iter = 0, inext = 0;
double t = 0, tnext = 0;
void init_events (void);
void event_register (Event event);
void _init_solver (void);



#if 1
static double mpi_time = 0.;
#endif

typedef struct {
  clock_t c;
  struct timeval tv;
  double tm;
} timer;

timer timer_start (void)
{
  timer t;
  t.c = clock();
  gettimeofday (&t.tv, NULL);
#if 1
  t.tm = mpi_time;
#endif
  return t;
}

double timer_elapsed (timer t)
{
  struct timeval tvend;
  gettimeofday (&tvend, NULL);
  return ((tvend.tv_sec - t.tv.tv_sec) +
   (tvend.tv_usec - t.tv.tv_usec)/1e6);
}



vector zerof= {{_NVARMAX + 0},{_NVARMAX + 1},{_NVARMAX + 2}};
vector unityf= {{_NVARMAX + 3},{_NVARMAX + 4},{_NVARMAX + 5}};
scalar unity= {_NVARMAX + 6};
scalar zeroc= {_NVARMAX + 7};



 vector fm = {{_NVARMAX + 3},{_NVARMAX + 4},{_NVARMAX + 5}};
 scalar cm = {(_NVARMAX + 6)};



static FILE ** qpopen_pipes = NULL;

FILE * qpopen (const char * command, const char * type)
{
  FILE * fp = popen (command, type);
  if (fp) {
    FILE ** i = qpopen_pipes;
    int n = 0;
    while (i && *i) { n++; i++; }
    qpopen_pipes = (FILE * *) prealloc (qpopen_pipes, (n + 2)*sizeof(FILE *),__func__,__FILE__,__LINE__);
    qpopen_pipes[n] = fp;
    qpopen_pipes[n+1] = NULL;
  }
  return fp;
}

int qpclose (FILE * fp)
{
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i == fp)
      *i = (FILE *) 1;
    i++;
  }
  return pclose (fp);
}

static void qpclose_all()
{
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i != (FILE *) 1)
      pclose (*i);
    i++;
  }
  pfree (qpopen_pipes,__func__,__FILE__,__LINE__);
  qpopen_pipes = NULL;
}






FILE * lfopen (const char * name, const char * mode)
{
  char fname[80];
  sprintf (fname, "%s-%d", name, pid());
  return fopen (fname, mode);
}



void * matrix_new (int n, int p, size_t size)
{
  void ** m = ((void * *) pmalloc ((n)*sizeof(void *),__func__,__FILE__,__LINE__));
  char * a = ((char *) pmalloc ((n*p*size)*sizeof(char),__func__,__FILE__,__LINE__));
  for (int i = 0; i < n; i++)
    m[i] = a + i*p*size;
  return m;
}

double matrix_inverse (double ** m, int n, double pivmin)
{
  int indxc[n], indxr[n], ipiv[n];
  int i, icol = 0, irow = 0, j, k, l, ll;
  double big, dum, pivinv, minpiv = HUGE;

  for (j = 0; j < n; j++)
    ipiv[j] = -1;

  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++)
      if (ipiv[j] != 0)
 for (k = 0; k < n; k++) {
   if (ipiv[k] == -1) {
     if (fabs (m[j][k]) >= big) {
       big = fabs (m[j][k]);
       irow = j;
       icol = k;
     }
   }
 }
    ipiv[icol]++;
    if (irow != icol)
      for (l = 0; l < n; l++)
 swap (double, m[irow][l], m[icol][l]);
    indxr[i] = irow;
    indxc[i] = icol;
    if (fabs (m[icol][icol]) <= pivmin)
      return 0.;
    if (fabs (m[icol][icol]) < minpiv)
      minpiv = fabs (m[icol][icol]);
    pivinv = 1.0/m[icol][icol];
    m[icol][icol] = 1.0;
    for (l = 0; l < n; l++) m[icol][l] *= pivinv;
    for (ll = 0; ll < n; ll++)
      if (ll != icol) {
 dum = m[ll][icol];
 m[ll][icol] = 0.0;
 for (l = 0; l < n; l++)
   m[ll][l] -= m[icol][l]*dum;
      }
  }
  for (l = n - 1; l >= 0; l--) {
    if (indxr[l] != indxc[l])
      for (k = 0; k < n; k++)
 swap (double, m[k][indxr[l]], m[k][indxc[l]]);
  }
  return minpiv;
}

void matrix_free (void * m)
{
  pfree (((void **) m)[0],__func__,__FILE__,__LINE__);
  pfree (m,__func__,__FILE__,__LINE__);
}
#line 14 "milk_crown-cpp.c"
static double _boundary0 (Point point, Point neighbor, scalar _s);
static double _boundary0_homogeneous (Point point, Point neighbor, scalar _s);
static double _boundary1 (Point point, Point neighbor, scalar _s);
static double _boundary1_homogeneous (Point point, Point neighbor, scalar _s);
static double _boundary2 (Point point, Point neighbor, scalar _s);
static double _boundary2_homogeneous (Point point, Point neighbor, scalar _s);
static double _boundary3 (Point point, Point neighbor, scalar _s);
static double _boundary3_homogeneous (Point point, Point neighbor, scalar _s);
static double _boundary4 (Point point, Point neighbor, scalar _s);
static double _boundary4_homogeneous (Point point, Point neighbor, scalar _s);
static double _boundary5 (Point point, Point neighbor, scalar _s);
static double _boundary5_homogeneous (Point point, Point neighbor, scalar _s);
#line 1 "milk_crown.c"
#line 1 "grid/octree.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/grid/octree.h"


#line 1 "grid/tree.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"
#line 1 "grid/mempool.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/grid/mempool.h"





typedef struct _Pool Pool;

struct _Pool {
  Pool * next;
};

typedef struct {
  char * first, * lastb;
  size_t size;
  size_t poolsize;
  Pool * pool, * last;
} Mempool;

typedef struct {
  char * next;
} FreeBlock;

Mempool * mempool_new (size_t poolsize, size_t size)
{

  assert (poolsize % 8 == 0);
  assert (size >= sizeof(FreeBlock));


  poolsize = min(1 << 20, poolsize + sizeof(Pool));
  Mempool * m = ((Mempool *) pcalloc (1, sizeof(Mempool),__func__,__FILE__,__LINE__));
  m->poolsize = poolsize;
  m->size = size;
  return m;
}

void mempool_destroy (Mempool * m)
{
  Pool * p = m->pool;
  while (p) {
    Pool * next = p->next;
    pfree (p,__func__,__FILE__,__LINE__);
    p = next;
  }
  pfree (m,__func__,__FILE__,__LINE__);
}

void * mempool_alloc (Mempool * m)
{
  if (!m->first) {

    Pool * p = (Pool *) pmalloc (m->poolsize,__func__,__FILE__,__LINE__);
    p->next = NULL;
    if (m->last)
      m->last->next = p;
    else
      m->pool = p;
    m->last = p;
    m->first = m->lastb = ((char *)m->last) + sizeof(Pool);
    FreeBlock * b = (FreeBlock *) m->first;
    b->next = NULL;
  }
  void * ret = m->first;
  FreeBlock * b = (FreeBlock *) ret;
  char * next = b->next;
  if (!next) {
    m->lastb += m->size;
    next = m->lastb;
    if (next + m->size > ((char *) m->last) + m->poolsize)
      next = NULL;
    else {
      FreeBlock * b = (FreeBlock *) next;
      b->next = NULL;
    }
  }
  m->first = next;
#if TRASH
  double * v = (double *) ret;
  for (int i = 0; i < m->size/sizeof(double); i++)
    v[i] = undefined;
#endif
  return ret;
}

void * mempool_alloc0 (Mempool * m)
{
  void * ret = mempool_alloc (m);
  memset (ret, 0, m->size);
  return ret;
}

void mempool_free (Mempool * m, void * p)
{
#if TRASH
  double * v = (double *) p;
  for (int i = 0; i < m->size/sizeof(double); i++)
    v[i] = undefined;
#endif
  FreeBlock * b = (FreeBlock *) p;
  b->next = m->first;
  m->first = (char *) p;
}
#line 2 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"
#line 22 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"
typedef struct {
  unsigned short flags;

  unsigned short neighbors;
  int pid;
} Cell;

enum {
  active = 1 << 0,
  leaf = 1 << 1,
  border = 1 << 2,
  vertex = 1 << 3,
  user = 4,

  face_x = 1 << 0

  , face_y = 1 << 1


  , face_z = 1 << 2

};

#define is_active(cell) ((cell).flags & active)
#define is_leaf(cell) ((cell).flags & leaf)
#define is_coarse() ((cell).neighbors > 0)
#define is_border(cell) ((cell).flags & border)
#define is_local(cell) ((cell).pid == pid())



typedef struct {
  int i;

  int j;


  int k;

} IndexLevel;

typedef struct {
  IndexLevel * p;
  int n, nm;
} CacheLevel;

typedef struct {
  int i;

  int j;


  int k;

  int level, flags;
} Index;

typedef struct {
  Index * p;
  int n, nm;
} Cache;




static void * new_refarray (size_t len, size_t size) {
  return pcalloc (len + 1, size,__func__,__FILE__,__LINE__);
}

static void refarray (void * p, size_t len, size_t size) {
  int * refcount = (int *)(((char *)p) + len*size);
  (*refcount)++;
}

static bool unrefarray (void * p, size_t len, size_t size) {
  int * refcount = (int *)(((char *)p) + len*size);
  (*refcount)--;
  if (*refcount == 0) {
    pfree (p,__func__,__FILE__,__LINE__);
    return true;
  }
  return false;
}




typedef struct {





  char **** m;

  Mempool * pool;
  int nc;
  int len;
} Layer;

static size_t _size (size_t depth)
{
  return (1 << depth) + 2*2;
}

static size_t poolsize (size_t depth, size_t size)
{






  return cube(_size(depth))*size;

}

static Layer * new_layer (int depth)
{
  Layer * l = ((Layer *) pmalloc ((1)*sizeof(Layer),__func__,__FILE__,__LINE__));
  l->len = _size (depth);
  if (depth == 0)
    l->pool = NULL;
  else {
    size_t size = sizeof(Cell) + datasize;


    l->pool = mempool_new (poolsize (depth, size), (1 << 3)*size);
  }





  l->m = ((char *** *) pcalloc (l->len, sizeof(char ***),__func__,__FILE__,__LINE__));

  l->nc = 0;
  return l;
}

static void destroy_layer (Layer * l)
{
  if (l->pool)
    mempool_destroy (l->pool);
  pfree (l->m,__func__,__FILE__,__LINE__);
  pfree (l,__func__,__FILE__,__LINE__);
}
#line 186 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"
static void layer_add_row (Layer * l, int i, int j)
{







  if (!l->m[i]) {
    l->m[i] = (char ***) new_refarray (l->len, sizeof (char **));
    l->nc++;
  }
  refarray (l->m[i], l->len, sizeof(char *));
  if (!l->m[i][j])
    l->m[i][j] = (char **) new_refarray (l->len, sizeof (char *));
  refarray (l->m[i][j], l->len, sizeof(char *));

}

static bool layer_remove_row (Layer * l, int i, int j)
{

  if (unrefarray (l->m[i][j], l->len, sizeof (char *)))
    l->m[i][j] = NULL;

  if (unrefarray (l->m[i], l->len, sizeof (char *))) {
    l->m[i] = NULL;
    if (--l->nc == 0) {
      destroy_layer (l);
      return true;
    }
    assert (l->nc >= 0);
  }
  return false;
}




typedef struct {
  Grid g;
  Layer ** L;

  Cache leaves;
  Cache faces;
  Cache vertices;
  Cache refined;
  CacheLevel * active;
  CacheLevel * prolongation;
  CacheLevel * boundary;

  CacheLevel * restriction;

  bool dirty;
} Tree;



struct _Point {

  int i;

  int j;


  int k;

  int level;
};
static Point last_point;



static void cache_level_append (CacheLevel * c, Point p)
{
  if (c->n >= c->nm) {
    c->nm += 128;
    c->p = (IndexLevel *) prealloc (c->p, (c->nm)*sizeof(IndexLevel),__func__,__FILE__,__LINE__);
  }
  c->p[c->n].i = p.i;

  c->p[c->n].j = p.j;


  c->p[c->n].k = p.k;

  c->n++;
}

static void cache_level_shrink (CacheLevel * c)
{
  if (c->nm > (c->n/128 + 1)*128) {
    c->nm = (c->n/128 + 1)*128;
    assert (c->nm > c->n);
    c->p = (IndexLevel *) prealloc (c->p, sizeof (Index)*c->nm,__func__,__FILE__,__LINE__);
  }
}

static void cache_append (Cache * c, Point p, unsigned short flags)
{
  if (c->n >= c->nm) {
    c->nm += 128;
    c->p = (Index *) prealloc (c->p, (c->nm)*sizeof(Index),__func__,__FILE__,__LINE__);
  }
  c->p[c->n].i = p.i;

  c->p[c->n].j = p.j;


  c->p[c->n].k = p.k;

  c->p[c->n].level = p.level;
  c->p[c->n].flags = flags;
  c->n++;
}

void cache_shrink (Cache * c)
{
  cache_level_shrink ((CacheLevel *)c);
}
#line 362 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"
#define allocated(a,l,n) (point.i+a >= 0 &&\
         point.i+a < (1 << point.level) + 2*2 &&\
         ((Tree *)grid)->L[point.level]->m[point.i+a] &&\
         point.j+l >= 0 &&\
         point.j+l < (1 << point.level) + 2*2 &&\
         ((Tree *)grid)->L[point.level]->m[point.i+a][point.j+l] &&\
         point.k+n >= 0 &&\
         point.k+n < (1 << point.level) + 2*2 &&\
         ((Tree *)grid)->L[point.level]->m[point.i+a][point.j+l]\
         [point.k+n])\

#line 372


#define NEIGHBOR(a,l,n) (((Tree *)grid)->L[point.level]->m[point.i+a][point.j+l]\
                                          [point.k+n])\

#line 376

#define PARENT(a,l,n) (((Tree *)grid)->L[point.level-1]->m[(point.i+2)/2+a]\
    [(point.j+2)/2+l][(point.k+2)/2+n])\

#line 379

#define allocated_child(a,l,n) (level < depth() &&\
         point.i > 0 && point.i <= (1 << level) + 2 &&\
         point.j > 0 && point.j <= (1 << level) + 2 &&\
         point.k > 0 && point.k <= (1 << level) + 2 &&\
         ((Tree *)grid)->L[point.level+1]->m[2*point.i-2 +a]\
   && ((Tree *)grid)->L[point.level+1]->m[2*point.i-2 +a][2*point.j-2 +l]\
   && ((Tree *)grid)->L[point.level+1]->m[2*point.i-2 +a][2*point.j-2 +l]\
         [2*point.k-2 +n])\

#line 388

#define CHILD(a,l,n) (((Tree *)grid)->L[point.level+1]->m[2*point.i-2 +a]\
      [2*point.j-2 +l][2*point.k-2 +n])\

#line 391


#define CELL(m) (*((Cell *)(m)))


#define depth() (grid->depth)
#define aparent(k,l,n) CELL(PARENT(k,l,n))
#define child(k,l,n) CELL(CHILD(k,l,n))


#define cell CELL(NEIGHBOR(0,0,0))
#define neighbor(k,l,n) CELL(NEIGHBOR(k,l,n))
#define neighborp(l,m,n) (Point) {\
    point.i + l,\
\
    point.j + m,\
\
\
    point.k + n,\
\
    point.level }\

#line 412



#define data(k,l,n) ((double *) (NEIGHBOR(k,l,n) + sizeof(Cell)))
#define fine(a,k,l,n) ((double *) (CHILD(k,l,n) + sizeof(Cell)))[a.i]
#define coarse(a,k,l,n) ((double *) (PARENT(k,l,n) + sizeof(Cell)))[a.i]

#define POINT_VARIABLES\
  VARIABLES\
  int level = point.level; NOT_UNUSED(level);\
\
\
\
\
\
\
\
  struct { int x, y, z; } child = {\
    2*((point.i+2)%2)-1, 2*((point.j+2)%2)-1, 2*((point.k+2)%2)-1\
  };\
\
  NOT_UNUSED(child);\
  Point parent = point; NOT_UNUSED(parent);\
  parent.level--;\
  parent.i = (point.i + 2)/2;\
\
  parent.j = (point.j + 2)/2;\
\
\
  parent.k = (point.k + 2)/2;\
\

#line 443


#line 1 "grid/foreach_cell.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/grid/foreach_cell.h"
#line 66 "/home/miro3/Documents/Programming/basilisk/src/grid/foreach_cell.h"
#define foreach_cell_root(root)\
  {\
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);\
    Point point;\
\
\
\
\
\
    int kg = 0; NOT_UNUSED(kg);\
    struct { int l, i, j, k, stage; } stack[20];\
\
    int _s = -1;\
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].k = root.k; stack[_s].stage = 0; };\
    while (_s >= 0) {\
      int stage;\
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; point.k = stack[_s].k; stage = stack[_s].stage; _s--; };\
      if (!allocated (0,0,0))\
 continue;\
      switch (stage) {\
      case 0: {\
 POINT_VARIABLES;\
\

#line 89

#define end_foreach_cell_root()\
        if (point.level < grid->depth) {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 1; };\
          { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };\
        }\
        break;\
      }\
\
      case 1: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 2; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;\
      case 2: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 3; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;\
      case 3: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 4; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;\
      case 4: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 5; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;\
      case 5: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 6; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;\
      case 6: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 7; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;\
      case 7: { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;\
\
      }\
    }\
  }\

#line 123


#define foreach_cell() {\
\
\
\
\
\
  Point root = {2,2,2,0};\
\
  foreach_cell_root (root)\

#line 134

#define end_foreach_cell() end_foreach_cell_root() }

#define foreach_cell_all() {\
  Point root = { .level = 0 };\
  for (root.i = 0; root.i <= 2*2; root.i++)\
\
    for (root.j = 0; root.j <= 2*2; root.j++)\
\
\
      for (root.k = 0; root.k <= 2*2; root.k++)\
\
 foreach_cell_root (root)\

#line 147

#define end_foreach_cell_all() end_foreach_cell_root() }

#define foreach_cell_post_root(condition, root)\
  {\
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);\
    Point point;\
\
\
\
\
\
    int kg = 0; NOT_UNUSED(kg);\
    struct { int l, i, j, k, stage; } stack[20];\
\
    int _s = -1;\
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].k = root.k; stack[_s].stage = 0; };\
    while (_s >= 0) {\
      int stage;\
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; point.k = stack[_s].k; stage = stack[_s].stage; _s--; };\
      if (!allocated (0,0,0))\
 continue;\
      switch (stage) {\
      case 0: {\
        POINT_VARIABLES;\
 if (point.level == grid->depth) {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 8; };\
 }\
 else {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 1; };\
   if (condition)\
     { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };\
 }\
 break;\
      }\
\
      case 1:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 2; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };\
 break;\
      case 2:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 3; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };\
 break;\
      case 3:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 4; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };\
 break;\
      case 4:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 5; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };\
 break;\
      case 5:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 6; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };\
 break;\
      case 6:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 7; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };\
 break;\
      case 7:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 8; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };\
 break;\
\
      default: {\
        POINT_VARIABLES;\
\

#line 244

#define end_foreach_cell_post_root()\
      }\
      }\
    }\
  }\

#line 250


#define foreach_cell_post(condition)\
  {\
\
\
\
\
\
    Point root = {2,2,2,0};\
\
    foreach_cell_post_root(condition, root)\

#line 262

#define end_foreach_cell_post() end_foreach_cell_post_root() }

#define foreach_cell_post_all(condition) {\
  Point root = { .level = 0 };\
  for (root.i = 0; root.i <= 2*2; root.i++)\
\
    for (root.j = 0; root.j <= 2*2; root.j++)\
\
\
      for (root.k = 0; root.k <= 2*2; root.k++)\
\
 foreach_cell_post_root (condition, root)\

#line 275

#define end_foreach_cell_post_all() end_foreach_cell_post_root() }

#define foreach_leaf() foreach_cell()\
  if (is_leaf (cell)) {\
    if (is_active(cell) && is_local(cell)) {\

#line 281

#define end_foreach_leaf() } continue; } end_foreach_cell()
#line 446 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"
#line 481 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"
#define foreach_child() {\
  int _i = 2*point.i - 2, _j = 2*point.j - 2, _k = 2*point.k - 2;\
  point.level++;\
  for (int _l = 0; _l < 2; _l++) {\
    point.i = _i + _l;\
    for (int _m = 0; _m < 2; _m++) {\
      point.j = _j + _m;\
      for (int _n = 0; _n < 2; _n++) {\
 point.k = _k + _n;\
 POINT_VARIABLES;\

#line 491

#define end_foreach_child()\
      }\
    }\
  }\
  point.i = (_i + 2)/2;point.j = (_j + 2)/2;point.k = (_k + 2)/2;\
  point.level--;\
}\

#line 499

#define foreach_child_break() _l = _m = _n = 2
#line 509 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"
#define foreach_cache(_cache) {\
  OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
\
\
\
\
\
  Point point = {2,2,2,0};\
\
  int _k; unsigned short _flags; NOT_UNUSED(_flags);\
  OMP(omp for schedule(static))\
  for (_k = 0; _k < _cache.n; _k++) {\
    point.i = _cache.p[_k].i;\
\
    point.j = _cache.p[_k].j;\
\
\
    point.k = _cache.p[_k].k;\
\
    point.level = _cache.p[_k].level;\
    _flags = _cache.p[_k].flags;\
    POINT_VARIABLES;\

#line 532

#define end_foreach_cache() } } }

#define foreach_cache_level(_cache,_l) {\
  OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
\
\
\
\
\
  Point point = {2,2,2,0};\
\
  point.level = _l;\
  int _k;\
  OMP(omp for schedule(static))\
  for (_k = 0; _k < _cache.n; _k++) {\
    point.i = _cache.p[_k].i;\
\
    point.j = _cache.p[_k].j;\
\
\
    point.k = _cache.p[_k].k;\
\
    POINT_VARIABLES;\

#line 557

#define end_foreach_cache_level() } } }

#define foreach_boundary_level(_l) {\
  if (_l <= depth()) {\
    { if (((Tree *)grid)->dirty) update_cache_f(); };\
    CacheLevel _boundary = ((Tree *)grid)->boundary[_l];\
    foreach_cache_level (_boundary,_l)\

#line 565

#define end_foreach_boundary_level() end_foreach_cache_level(); }}



#define foreach_boundary(_b) {\
  for (int _l = depth(); _l >= 0; _l--)\
    foreach_boundary_level(_l) {\
      if ((- cell.pid - 1) == _b)\
 for (int _d = 0; _d < 3; _d++) {\
   for (int _i = -1; _i <= 1; _i += 2) {\
     if (_d == 0) ig = _i; else if (_d == 1) jg = _i; else kg = _i;\
     if (allocated(-ig,-jg,-kg) &&\
  is_leaf (neighbor(-ig,-jg,-kg)) &&\
  !(neighbor(-ig,-jg,-kg).pid < 0) &&\
  is_local(neighbor(-ig,-jg,-kg))) {\
       point.i -= ig; x -= ig*Delta/2.;\
\
       point.j -= jg; y -= jg*Delta/2.;\
\
\
       point.k -= kg; z -= kg*Delta/2.;\
\

#line 588

#define end_foreach_boundary()\
       point.i += ig; x += ig*Delta/2.;\
\
       point.j += jg; y += jg*Delta/2.;\
\
\
       point.k += kg; z += kg*Delta/2.;\
\
            }\
   }\
   ig = jg = kg = 0;\
 }\
    } end_foreach_boundary_level(); }\

#line 602


#define foreach_halo(_name,_l) {\
  if (_l <= depth()) {\
    { if (((Tree *)grid)->dirty) update_cache_f(); };\
    CacheLevel _cache = ((Tree *)grid)->_name[_l];\
    foreach_cache_level (_cache, _l)\

#line 609

#define end_foreach_halo() end_foreach_cache_level(); }}

#line 1 "grid/neighbors.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/grid/neighbors.h"
#line 35 "/home/miro3/Documents/Programming/basilisk/src/grid/neighbors.h"
#define foreach_neighbor(_s) {\
  int _nn = _s + 0 ? _s + 0 : 2;\
  int _i = point.i, _j = point.j, _k = point.k;\
  for (int _l = - _nn; _l <= _nn; _l++) {\
    point.i = _i + _l;\
    for (int _m = - _nn; _m <= _nn; _m++) {\
      point.j = _j + _m;\
      for (int _n = - _nn; _n <= _nn; _n++) {\
 point.k = _k + _n;\
 POINT_VARIABLES;\

#line 45

#define end_foreach_neighbor()\
      }\
    }\
  }\
  point.i = _i; point.j = _j; point.k = _k;\
}\

#line 52

#define foreach_neighbor_break() _l = _m = _n = _nn + 1
#line 613 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"

static inline bool has_local_children (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 615 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"

   { foreach_child()
    if (is_local(cell))
      return true; end_foreach_child(); }
  return false;
}

static inline void cache_append_face (Point point, unsigned short flags)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 623 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"

  Tree * q = ((Tree *)grid);
  cache_append (&q->faces, point, flags);
#line 637 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"
  {
#line 637

    if (flags & face_x)
      for (int i = 0; i <= 1; i++)
 for (int j = 0; j <= 1; j++)
   if (!(neighbor(0,i,j).flags & vertex)) {
     cache_append (&q->vertices, neighborp(0,i,j), 0);
     neighbor(0,i,j).flags |= vertex;
   }
#line 637

    if (flags & face_y)
      for (int i = 0; i <= 1; i++)
 for (int j = 0; j <= 1; j++)
   if (!(neighbor(j,0,i).flags & vertex)) {
     cache_append (&q->vertices, neighborp(j,0,i), 0);
     neighbor(j,0,i).flags |= vertex;
   }
#line 637

    if (flags & face_z)
      for (int i = 0; i <= 1; i++)
 for (int j = 0; j <= 1; j++)
   if (!(neighbor(i,j,0).flags & vertex)) {
     cache_append (&q->vertices, neighborp(i,j,0), 0);
     neighbor(i,j,0).flags |= vertex;
   }}

}



static void update_cache_f (void)
{
  Tree * q = ((Tree *)grid);

   { foreach_cache(q->vertices){

#line 654 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"

    if (level <= depth() && allocated(0,0,0))
      cell.flags &= ~vertex; } end_foreach_cache(); }


  q->leaves.n = q->faces.n = q->vertices.n = 0;
  for (int l = 0; l <= depth(); l++)
    q->active[l].n = q->prolongation[l].n =
      q->boundary[l].n = q->restriction[l].n = 0;

  const unsigned short fboundary = 1 << user;
   { foreach_cell(){

#line 665 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"
 {



    if (is_local(cell) && is_active(cell)) {


      cache_level_append (&q->active[level], point);
    }
#line 689 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"
    if (!(cell.pid < 0)) {

       { foreach_neighbor (2)
 if (allocated(0,0,0) && (cell.pid < 0) && !(cell.flags & fboundary)) {
   cache_level_append (&q->boundary[level], point);
   cell.flags |= fboundary;
 } end_foreach_neighbor(); }
    }

    else if (level > 0 && is_local(aparent(0,0,0)))
      cache_level_append (&q->restriction[level], point);

    if (is_leaf (cell)) {
      if (is_local(cell)) {
 cache_append (&q->leaves, point, 0);

 unsigned short flags = 0;
 {
#line 706

   if ((neighbor(-1,0,0).pid < 0) || (!is_leaf(neighbor(-1,0,0)) && !neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) ||
       is_leaf(neighbor(-1,0,0)))
     flags |= face_x;
#line 706

   if ((neighbor(0,-1,0).pid < 0) || (!is_leaf(neighbor(0,-1,0)) && !neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) ||
       is_leaf(neighbor(0,-1,0)))
     flags |= face_y;
#line 706

   if ((neighbor(0,0,-1).pid < 0) || (!is_leaf(neighbor(0,0,-1)) && !neighbor(0,0,-1).neighbors && neighbor(0,0,-1).pid >= 0) ||
       is_leaf(neighbor(0,0,-1)))
     flags |= face_z;}
 if (flags)
   cache_append (&q->faces, point, flags);
 {
#line 712

   if ((neighbor(1,0,0).pid < 0) || (!is_leaf(neighbor(1,0,0)) && !neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0) ||
       (!is_local(neighbor(1,0,0)) && is_leaf(neighbor(1,0,0))))
     cache_append (&q->faces, neighborp(1,0,0), face_x);
#line 712

   if ((neighbor(0,1,0).pid < 0) || (!is_leaf(neighbor(0,1,0)) && !neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0) ||
       (!is_local(neighbor(0,1,0)) && is_leaf(neighbor(0,1,0))))
     cache_append (&q->faces, neighborp(0,1,0), face_y);
#line 712

   if ((neighbor(0,0,1).pid < 0) || (!is_leaf(neighbor(0,0,1)) && !neighbor(0,0,1).neighbors && neighbor(0,0,1).pid >= 0) ||
       (!is_local(neighbor(0,0,1)) && is_leaf(neighbor(0,0,1))))
     cache_append (&q->faces, neighborp(0,0,1), face_z);}

 for (int i = 0; i <= 1; i++)

   for (int j = 0; j <= 1; j++)


     for (int k = 0; k <= 1; k++)

       if (!(neighbor(i,j,k).flags & vertex)) {
  cache_append (&q->vertices, neighborp(i,j,k), 0);
  neighbor(i,j,k).flags |= vertex;
       }

        if (cell.neighbors > 0)
   cache_level_append (&q->prolongation[level], point);
      }
      else if (!(cell.pid < 0) || is_local(aparent(0,0,0))) {

 unsigned short flags = 0;
 {
#line 735

   if (allocated(-1,0,0) &&
       is_local(neighbor(-1,0,0)) && (!is_leaf(neighbor(-1,0,0)) && !neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0))
     flags |= face_x;
#line 735

   if (allocated(0,-1,0) &&
       is_local(neighbor(0,-1,0)) && (!is_leaf(neighbor(0,-1,0)) && !neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0))
     flags |= face_y;
#line 735

   if (allocated(0,0,-1) &&
       is_local(neighbor(0,0,-1)) && (!is_leaf(neighbor(0,0,-1)) && !neighbor(0,0,-1).neighbors && neighbor(0,0,-1).pid >= 0))
     flags |= face_z;}
 if (flags)
   cache_append_face (point, flags);
 {
#line 741

   if (allocated(1,0,0) && is_local(neighbor(1,0,0)) &&
       (!is_leaf(neighbor(1,0,0)) && !neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0))
     cache_append_face (neighborp(1,0,0), face_x);
#line 741

   if (allocated(0,1,0) && is_local(neighbor(0,1,0)) &&
       (!is_leaf(neighbor(0,1,0)) && !neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0))
     cache_append_face (neighborp(0,1,0), face_y);
#line 741

   if (allocated(0,0,1) && is_local(neighbor(0,0,1)) &&
       (!is_leaf(neighbor(0,0,1)) && !neighbor(0,0,1).neighbors && neighbor(0,0,1).pid >= 0))
     cache_append_face (neighborp(0,0,1), face_z);}
      }

      continue;

    }
  } } end_foreach_cell(); }


  cache_shrink (&q->leaves);
  cache_shrink (&q->faces);
  cache_shrink (&q->vertices);
  for (int l = 0; l <= depth(); l++) {
    cache_level_shrink (&q->active[l]);
    cache_level_shrink (&q->prolongation[l]);
    cache_level_shrink (&q->boundary[l]);
    cache_level_shrink (&q->restriction[l]);
}

  q->dirty = false;


  for (int l = depth(); l >= 0; l--)
     { foreach_boundary_level (l){

#line 767 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"

      cell.flags &= ~fboundary; } end_foreach_boundary_level(); }



  grid->n = q->leaves.n;

#if !1
  grid->tn = grid->n;
  grid->maxdepth = grid->depth;
#endif
}

#define foreach() { if (((Tree *)grid)->dirty) update_cache_f(); }; foreach_cache(((Tree *)grid)->leaves)
#define end_foreach() end_foreach_cache()

#define foreach_face_generic()\
  { if (((Tree *)grid)->dirty) update_cache_f(); };\
  foreach_cache(((Tree *)grid)->faces) 
#line 784

#define end_foreach_face_generic() end_foreach_cache()

#define is_face_x() (_flags & face_x)

#define is_face_y() (_flags & face_y)


#define is_face_z() (_flags & face_z)


#define foreach_vertex()\
  { if (((Tree *)grid)->dirty) update_cache_f(); };\
  foreach_cache(((Tree *)grid)->vertices) {\
    x -= Delta/2.;\
\
    y -= Delta/2.;\
\
\
    z -= Delta/2.;\
\

#line 806

#define end_foreach_vertex() } end_foreach_cache()
#line 818 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"
#define foreach_level(l) {\
  if (l <= depth()) {\
    { if (((Tree *)grid)->dirty) update_cache_f(); };\
    CacheLevel _active = ((Tree *)grid)->active[l];\
    foreach_cache_level (_active,l)\

#line 823

#define end_foreach_level() end_foreach_cache_level(); }}

#define foreach_coarse_level(l) foreach_level(l) if (!is_leaf(cell)) {
#define end_foreach_coarse_level() } end_foreach_level()

#define foreach_level_or_leaf(l) {\
  for (int _l1 = l; _l1 >= 0; _l1--)\
    foreach_level(_l1)\
      if (_l1 == l || is_leaf (cell)) {\

#line 833

#define end_foreach_level_or_leaf() } end_foreach_level(); }

#if TRASH
# undef trash
# define trash(list) reset(list, undefined)
#endif

void reset (void * alist, double val)
{
  scalar * list = (scalar *) alist;
  Tree * q = ((Tree *)grid);

  for (int l = 0; l <= depth(); l++) {
    Layer * L = q->L[l];
    for (int i = 0; i < L->len; i++)
      if (L->m[i])





 for (int j = 0; j < L->len; j++)
   if (L->m[i][j])





          for (int k = 0; k < L->len; k++)
     if (L->m[i][j][k])
       if (list) for (scalar s = *list, *_i10 = list; ((scalar *)&s)->i >= 0; s = *++_i10)
           if (!is_constant(s))
    ((double *)(L->m[i][j][k] + sizeof(Cell)))[s.i] = val;


  }
}

#define cache_level_resize(name, a)\
{\
  for (int i = 0; i <= depth() - a; i++)\
    pfree (q->name[i].p,__func__,__FILE__,__LINE__);\
  pfree (q->name,__func__,__FILE__,__LINE__);\
  q->name = ((CacheLevel *) pcalloc (depth() + 1, sizeof(CacheLevel),__func__,__FILE__,__LINE__));\
}\

#line 879


static void update_depth (int inc)
{
  Tree * q = ((Tree *)grid);
  grid->depth += inc;
  q->L = &(q->L[-1]);
  q->L = (Layer * *) prealloc (q->L, (grid->depth + 2)*sizeof(Layer *),__func__,__FILE__,__LINE__);
  q->L = &(q->L[1]);
  if (inc > 0)
    q->L[grid->depth] = new_layer (grid->depth);
  cache_level_resize (active, inc);
  cache_level_resize (prolongation, inc);
  cache_level_resize (boundary, inc);
  cache_level_resize (restriction, inc);
}

static void alloc_children (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 897 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"

  if (point.level == grid->depth)
    update_depth (+1);
  else if (allocated_child(0,0,0))
    return;
#line 929 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"
  Layer * L = ((Tree *)grid)->L[point.level + 1];
  size_t len = sizeof(Cell) + datasize;
  char * b = (char *) mempool_alloc0 (L->pool);
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++) {
      layer_add_row (L, 2*point.i - 2 + k, 2*point.j - 2 + l);
      for (int n = 0; n < 2; n++) {
 assert (!CHILD(k,l,n));
 CHILD(k,l,n) = b;
 b += len;
      }
    }


  int pid = cell.pid;
   { foreach_child() {
    cell.pid = pid;
#if TRASH
    if (all) for (scalar s = *all, *_i11 = all; ((scalar *)&s)->i >= 0; s = *++_i11)
      val(s,0,0,0) = undefined;
#endif
  } end_foreach_child(); }
}
#line 985 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"
static void free_children (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 986 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"


  Layer * L = ((Tree *)grid)->L[point.level + 1];
  assert (CHILD(0,0,0));
  mempool_free (L->pool, CHILD(0,0,0));
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++) {
      for (int n = 0; n < 2; n++)
 CHILD(k,l,n) = NULL;
      if (layer_remove_row (L, 2*point.i - 2 + k, 2*point.j - 2 + l)) {
 assert (point.level + 1 == grid->depth);
 update_depth (-1);
      }
    }
}


void increment_neighbors (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1004 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"

  ((Tree *)grid)->dirty = true;
   { foreach_neighbor (2/2)
    if (cell.neighbors++ == 0)
      alloc_children (point); end_foreach_neighbor(); }
}

void decrement_neighbors (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1012 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"

  ((Tree *)grid)->dirty = true;
   { foreach_neighbor (2/2)
    if (allocated(0,0,0)) {
      cell.neighbors--;
      if (cell.neighbors == 0)
 free_children (point);
    } end_foreach_neighbor(); }
  if (cell.neighbors) {
    int pid = cell.pid;
     { foreach_child() {
      cell.flags = 0;
      cell.pid = pid;
    } end_foreach_child(); }
  }
}

void realloc_scalar (void)
{

  Tree * q = ((Tree *)grid);
  size_t newlen = sizeof(Cell) + datasize;
  size_t oldlen = newlen - sizeof(double);

  size_t len = _size(0);
  for (int i = 0; i < len; i++)



    for (int j = 0; j < len; j++)



      for (int k = 0; k < len; k++)
 q->L[0]->m[i][j][k] = (char *) prealloc (q->L[0]->m[i][j][k], (newlen)*sizeof(char),__func__,__FILE__,__LINE__);



  for (int l = 1; l <= depth(); l++) {
    Layer * L = q->L[l];
    size_t len = L->len;
    Mempool * oldpool = L->pool;
    L->pool = mempool_new (poolsize (l, newlen), (1 << 3)*newlen);
    for (int i = 0; i < len; i += 2)
      if (L->m[i]) {
#line 1065 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"
 for (int j = 0; j < len; j += 2)
   if (L->m[i][j]) {
#line 1076 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"
     for (int k = 0; k < len; k += 2)
       if (L->m[i][j][k]) {
  char * new = (char *) mempool_alloc (L->pool);
  for (int l = 0; l < 2; l++)
    for (int m = 0; m < 2; m++)
      for (int n = 0; n < 2; n++) {
        memcpy (new, L->m[i+l][j+m][k+n], oldlen);
        L->m[i+l][j+m][k+n] = new;
        new += newlen;
      }
       }

   }

      }
    mempool_destroy (oldpool);
  }
}



#define VN v.x
#define VT v.y
#define VR v.z




#if 1
# define disable_fpe_for_mpi() disable_fpe (FE_DIVBYZERO|FE_INVALID)
# define enable_fpe_for_mpi() enable_fpe (FE_DIVBYZERO|FE_INVALID)
#else
# define disable_fpe_for_mpi()
# define enable_fpe_for_mpi()
#endif

static inline void no_restriction (Point point, scalar s);

static bool normal_neighbor (Point point, scalar * scalars, vector * vectors)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1115 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"

  for (int k = 1; k <= 2; k++)
    {
#line 1117

      for (int i = -k; i <= k; i += 2*k)
 if ((allocated(i,0,0) && !(neighbor(i,0,0).pid < 0))) {
   Point neighbor = neighborp(i,0,0);
   int id = (- cell.pid - 1);
   if (scalars) for (scalar s = *scalars, *_i12 = scalars; ((scalar *)&s)->i >= 0; s = *++_i12)
     val(s,0,0,0) = _attribute[s.i].boundary[id](neighbor, point, s);
   if (vectors) for (vector v = *vectors, *_i13 = vectors; ((scalar *)&v)->i >= 0; v = *++_i13) {
     scalar vn = VN;
     val(v.x,0,0,0) = _attribute[vn.i].boundary[id](neighbor, point, v.x);

     scalar vt = VT;
     val(v.y,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.y);


     scalar vr = VR;
     val(v.z,0,0,0) = _attribute[vr.i].boundary[id](neighbor, point, v.z);

   }
   return true;
 }
#line 1117

      for (int i = -k; i <= k; i += 2*k)
 if ((allocated(0,i,0) && !(neighbor(0,i,0).pid < 0))) {
   Point neighbor = neighborp(0,i,0);
   int id = (- cell.pid - 1);
   if (scalars) for (scalar s = *scalars, *_i12 = scalars; ((scalar *)&s)->i >= 0; s = *++_i12)
     val(s,0,0,0) = _attribute[s.i].boundary[id](neighbor, point, s);
   if (vectors) for (vector v = *vectors, *_i13 = vectors; ((scalar *)&v)->i >= 0; v = *++_i13) {
     scalar vn = VN;
     val(v.y,0,0,0) = _attribute[vn.i].boundary[id](neighbor, point, v.y);

     scalar vt = VT;
     val(v.z,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.z);


     scalar vr = VR;
     val(v.x,0,0,0) = _attribute[vr.i].boundary[id](neighbor, point, v.x);

   }
   return true;
 }
#line 1117

      for (int i = -k; i <= k; i += 2*k)
 if ((allocated(0,0,i) && !(neighbor(0,0,i).pid < 0))) {
   Point neighbor = neighborp(0,0,i);
   int id = (- cell.pid - 1);
   if (scalars) for (scalar s = *scalars, *_i12 = scalars; ((scalar *)&s)->i >= 0; s = *++_i12)
     val(s,0,0,0) = _attribute[s.i].boundary[id](neighbor, point, s);
   if (vectors) for (vector v = *vectors, *_i13 = vectors; ((scalar *)&v)->i >= 0; v = *++_i13) {
     scalar vn = VN;
     val(v.z,0,0,0) = _attribute[vn.i].boundary[id](neighbor, point, v.z);

     scalar vt = VT;
     val(v.x,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.x);


     scalar vr = VR;
     val(v.y,0,0,0) = _attribute[vr.i].boundary[id](neighbor, point, v.y);

   }
   return true;
 }}
  return false;
}

static bool diagonal_neighbor_2D (Point point,
      scalar * scalars, vector * vectors)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1143 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"


  for (int k = 1; k <= 2; k++)

    {
#line 1147


      for (int i = -k; i <= k; i += 2*k)
 for (int j = -k; j <= k; j += 2*k)
   if (allocated(i,j,0) && (allocated(i,j,0) && !(neighbor(i,j,0).pid < 0)) &&
       allocated(i,0,0) && (neighbor(i,0,0).pid < 0) &&
       allocated(0,j,0) && (neighbor(0,j,0).pid < 0)) {
     Point n = neighborp(i,j,0),
       n1 = neighborp(i,0,0), n2 = neighborp(0,j,0);
     int id1 = (- neighbor(i,0,0).pid - 1), id2 = (- neighbor(0,j,0).pid - 1);
     if (scalars) for (scalar s = *scalars, *_i14 = scalars; ((scalar *)&s)->i >= 0; s = *++_i14)
       val(s,0,0,0) = (_attribute[s.i].boundary[id1](n,n1,s) + _attribute[s.i].boundary[id2](n,n2,s) -
       val(s,i,j,0));
     if (vectors) for (vector v = *vectors, *_i15 = vectors; ((scalar *)&v)->i >= 0; v = *++_i15) {
       scalar vt = VT, vn = VN;
       val(v.x,0,0,0) = (_attribute[vt.i].boundary[id1](n,n1,v.x) +
         _attribute[vn.i].boundary[id2](n,n2,v.x) -
         val(v.x,i,j,0));
       val(v.y,0,0,0) = (_attribute[vn.i].boundary[id1](n,n1,v.y) +
         _attribute[vt.i].boundary[id2](n,n2,v.y) -
         val(v.y,i,j,0));

       scalar vr = VR;
       val(v.z,0,0,0) = (_attribute[vr.i].boundary[id1](n,n1,v.z) +
         _attribute[vr.i].boundary[id2](n,n2,v.z) -
         val(v.z,i,j,0));

     }
     return true;
   }
#line 1147


      for (int i = -k; i <= k; i += 2*k)
 for (int j = -k; j <= k; j += 2*k)
   if (allocated(0,i,j) && (allocated(0,i,j) && !(neighbor(0,i,j).pid < 0)) &&
       allocated(0,i,0) && (neighbor(0,i,0).pid < 0) &&
       allocated(0,0,j) && (neighbor(0,0,j).pid < 0)) {
     Point n = neighborp(0,i,j),
       n1 = neighborp(0,i,0), n2 = neighborp(0,0,j);
     int id1 = (- neighbor(0,i,0).pid - 1), id2 = (- neighbor(0,0,j).pid - 1);
     if (scalars) for (scalar s = *scalars, *_i14 = scalars; ((scalar *)&s)->i >= 0; s = *++_i14)
       val(s,0,0,0) = (_attribute[s.i].boundary[id1](n,n1,s) + _attribute[s.i].boundary[id2](n,n2,s) -
       val(s,0,i,j));
     if (vectors) for (vector v = *vectors, *_i15 = vectors; ((scalar *)&v)->i >= 0; v = *++_i15) {
       scalar vt = VT, vn = VN;
       val(v.y,0,0,0) = (_attribute[vt.i].boundary[id1](n,n1,v.y) +
         _attribute[vn.i].boundary[id2](n,n2,v.y) -
         val(v.y,0,i,j));
       val(v.z,0,0,0) = (_attribute[vn.i].boundary[id1](n,n1,v.z) +
         _attribute[vt.i].boundary[id2](n,n2,v.z) -
         val(v.z,0,i,j));

       scalar vr = VR;
       val(v.x,0,0,0) = (_attribute[vr.i].boundary[id1](n,n1,v.x) +
         _attribute[vr.i].boundary[id2](n,n2,v.x) -
         val(v.x,0,i,j));

     }
     return true;
   }
#line 1147


      for (int i = -k; i <= k; i += 2*k)
 for (int j = -k; j <= k; j += 2*k)
   if (allocated(j,0,i) && (allocated(j,0,i) && !(neighbor(j,0,i).pid < 0)) &&
       allocated(0,0,i) && (neighbor(0,0,i).pid < 0) &&
       allocated(j,0,0) && (neighbor(j,0,0).pid < 0)) {
     Point n = neighborp(j,0,i),
       n1 = neighborp(0,0,i), n2 = neighborp(j,0,0);
     int id1 = (- neighbor(0,0,i).pid - 1), id2 = (- neighbor(j,0,0).pid - 1);
     if (scalars) for (scalar s = *scalars, *_i14 = scalars; ((scalar *)&s)->i >= 0; s = *++_i14)
       val(s,0,0,0) = (_attribute[s.i].boundary[id1](n,n1,s) + _attribute[s.i].boundary[id2](n,n2,s) -
       val(s,j,0,i));
     if (vectors) for (vector v = *vectors, *_i15 = vectors; ((scalar *)&v)->i >= 0; v = *++_i15) {
       scalar vt = VT, vn = VN;
       val(v.z,0,0,0) = (_attribute[vt.i].boundary[id1](n,n1,v.z) +
         _attribute[vn.i].boundary[id2](n,n2,v.z) -
         val(v.z,j,0,i));
       val(v.x,0,0,0) = (_attribute[vn.i].boundary[id1](n,n1,v.x) +
         _attribute[vt.i].boundary[id2](n,n2,v.x) -
         val(v.x,j,0,i));

       scalar vr = VR;
       val(v.y,0,0,0) = (_attribute[vr.i].boundary[id1](n,n1,v.y) +
         _attribute[vr.i].boundary[id2](n,n2,v.y) -
         val(v.y,j,0,i));

     }
     return true;
   }}

  return false;
}

static bool diagonal_neighbor_3D (Point point,
      scalar * scalars, vector * vectors)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1183 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"


  for (int n = 1; n <= 2; n++)
    for (int i = -n; i <= n; i += 2*n)
      for (int j = -n; j <= n; j += 2*n)
 for (int k = -n; k <= n; k += 2*n)
   if ((allocated(i,j,k) && !(neighbor(i,j,k).pid < 0)) &&
       (neighbor(i,j,0).pid < 0) &&
       (neighbor(i,0,k).pid < 0) &&
       (neighbor(0,j,k).pid < 0)) {
     Point
       n0 = neighborp(i,j,k),
       n1 = neighborp(i,j,0),
       n2 = neighborp(i,0,k),
       n3 = neighborp(0,j,k);
     int
       id1 = (- neighbor(i,j,0).pid - 1),
       id2 = (- neighbor(i,0,k).pid - 1),
       id3 = (- neighbor(0,j,k).pid - 1);
     if (scalars) for (scalar s = *scalars, *_i16 = scalars; ((scalar *)&s)->i >= 0; s = *++_i16)
       val(s,0,0,0) = (_attribute[s.i].boundary[id1](n0,n1,s) +
       _attribute[s.i].boundary[id2](n0,n2,s) +
       _attribute[s.i].boundary[id3](n0,n3,s) -
       2.*val(s,i,j,k));
     if (vectors) for (vector v = *vectors, *_i17 = vectors; ((scalar *)&v)->i >= 0; v = *++_i17) {
       scalar vt = VT, vn = VN, vr = VR;
       val(v.x,0,0,0) = (_attribute[vt.i].boundary[id1](n0,n1,v.x) +
         _attribute[vt.i].boundary[id2](n0,n2,v.x) +
         _attribute[vn.i].boundary[id3](n0,n3,v.x) -
         2.*val(v.x,i,j,k));
       val(v.y,0,0,0) = (_attribute[vt.i].boundary[id1](n0,n1,v.y) +
         _attribute[vn.i].boundary[id2](n0,n2,v.y) +
         _attribute[vt.i].boundary[id3](n0,n3,v.y) -
         2.*val(v.y,i,j,k));
       val(v.z,0,0,0) = (_attribute[vn.i].boundary[id1](n0,n1,v.z) +
         _attribute[vr.i].boundary[id2](n0,n2,v.z) +
         _attribute[vr.i].boundary[id3](n0,n3,v.z) -
         2.*val(v.z,i,j,k));
     }
     return true;
   }

  return false;
}



#line 1229

static Point tangential_neighbor_x (Point point, bool * zn)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1231 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"

  for (int k = 1; k <= 2; k++)
    for (int j = -k; j <= k; j += 2*k) {
      if ((allocated(0,j,0) && !(neighbor(0,j,0).pid < 0)) || (allocated(-1,j,0) && !(neighbor(-1,j,0).pid < 0))) {
 *zn = false;
 return neighborp(0,j,0);
      }


      if ((allocated(0,0,j) && !(neighbor(0,0,j).pid < 0)) || (allocated(-1,0,j) && !(neighbor(-1,0,j).pid < 0))) {
 *zn = true;
 return neighborp(0,0,j);
      }

    }
  return (Point){.level = -1};
}
#line 1229

static Point tangential_neighbor_y (Point point, bool * zn)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1231 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"

  for (int k = 1; k <= 2; k++)
    for (int j = -k; j <= k; j += 2*k) {
      if ((allocated(0,0,j) && !(neighbor(0,0,j).pid < 0)) || (allocated(0,-1,j) && !(neighbor(0,-1,j).pid < 0))) {
 *zn = false;
 return neighborp(0,0,j);
      }


      if ((allocated(j,0,0) && !(neighbor(j,0,0).pid < 0)) || (allocated(j,-1,0) && !(neighbor(j,-1,0).pid < 0))) {
 *zn = true;
 return neighborp(j,0,0);
      }

    }
  return (Point){.level = -1};
}
#line 1229

static Point tangential_neighbor_z (Point point, bool * zn)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1231 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"

  for (int k = 1; k <= 2; k++)
    for (int j = -k; j <= k; j += 2*k) {
      if ((allocated(j,0,0) && !(neighbor(j,0,0).pid < 0)) || (allocated(j,0,-1) && !(neighbor(j,0,-1).pid < 0))) {
 *zn = false;
 return neighborp(j,0,0);
      }


      if ((allocated(0,j,0) && !(neighbor(0,j,0).pid < 0)) || (allocated(0,j,-1) && !(neighbor(0,j,-1).pid < 0))) {
 *zn = true;
 return neighborp(0,j,0);
      }

    }
  return (Point){.level = -1};
}


static inline bool is_boundary_point (Point point) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1250 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"

  return (cell.pid < 0);
}

static void box_boundary_level (const Boundary * b, scalar * list, int l)
{
  disable_fpe_for_mpi();
  scalar * scalars = NULL;
  vector * vectors = NULL, * faces = NULL;
  if (list) for (scalar s = *list, *_i18 = list; ((scalar *)&s)->i >= 0; s = *++_i18)
    if (!is_constant(s) && _attribute[s.i].refine != no_restriction) {
      if (_attribute[s.i].v.x.i == s.i) {
 if (_attribute[s.i].face)
   faces = vectors_add (faces, _attribute[s.i].v);
 else
   vectors = vectors_add (vectors, _attribute[s.i].v);
      }
      else if (_attribute[s.i].v.x.i < 0)
 scalars = list_add (scalars, s);
    }

   { foreach_boundary_level (l){

#line 1271 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"
 {
    if (!normal_neighbor (point, scalars, vectors) &&
 !diagonal_neighbor_2D (point, scalars, vectors) &&
 !diagonal_neighbor_3D (point, scalars, vectors)) {

      if (scalars) for (scalar s = *scalars, *_i19 = scalars; ((scalar *)&s)->i >= 0; s = *++_i19)
 val(s,0,0,0) = undefined;
      if (vectors) for (vector v = *vectors, *_i20 = vectors; ((scalar *)&v)->i >= 0; v = *++_i20)
 {
#line 1279

   val(v.x,0,0,0) = undefined;
#line 1279

   val(v.y,0,0,0) = undefined;
#line 1279

   val(v.z,0,0,0) = undefined;}
    }
    if (faces) {
      int id = (- cell.pid - 1);
      {
#line 1284

 for (int i = -1; i <= 1; i += 2) {

   if ((allocated(i,0,0) && !(neighbor(i,0,0).pid < 0))) {
     Point neighbor = neighborp(i,0,0);
     if (faces) for (vector v = *faces, *_i21 = faces; ((scalar *)&v)->i >= 0; v = *++_i21) {
       scalar vn = VN;
       if (_attribute[vn.i].boundary[id])
  val(v.x,(i + 1)/2,0,0) = _attribute[vn.i].boundary[id](neighbor, point, v.x);
     }
   }

   else if (i == -1) {

     bool zn;
     Point neighbor = tangential_neighbor_x (point, &zn);
     if (neighbor.level >= 0) {
       int id = is_boundary_point (neighbor) ?
  (- neighbor(-1,0,0).pid - 1) : (- cell.pid - 1);
       if (faces) for (vector v = *faces, *_i22 = faces; ((scalar *)&v)->i >= 0; v = *++_i22) {



  scalar vt = zn ? VT : VR;

  val(v.x,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.x);
       }
     }
     else

       if (faces) for (vector v = *faces, *_i23 = faces; ((scalar *)&v)->i >= 0; v = *++_i23)
  val(v.x,0,0,0) = 0.;
   }

 }
#line 1284

 for (int i = -1; i <= 1; i += 2) {

   if ((allocated(0,i,0) && !(neighbor(0,i,0).pid < 0))) {
     Point neighbor = neighborp(0,i,0);
     if (faces) for (vector v = *faces, *_i21 = faces; ((scalar *)&v)->i >= 0; v = *++_i21) {
       scalar vn = VN;
       if (_attribute[vn.i].boundary[id])
  val(v.y,0,(i + 1)/2,0) = _attribute[vn.i].boundary[id](neighbor, point, v.y);
     }
   }

   else if (i == -1) {

     bool zn;
     Point neighbor = tangential_neighbor_y (point, &zn);
     if (neighbor.level >= 0) {
       int id = is_boundary_point (neighbor) ?
  (- neighbor(0,-1,0).pid - 1) : (- cell.pid - 1);
       if (faces) for (vector v = *faces, *_i22 = faces; ((scalar *)&v)->i >= 0; v = *++_i22) {



  scalar vt = zn ? VT : VR;

  val(v.y,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.y);
       }
     }
     else

       if (faces) for (vector v = *faces, *_i23 = faces; ((scalar *)&v)->i >= 0; v = *++_i23)
  val(v.y,0,0,0) = 0.;
   }

 }
#line 1284

 for (int i = -1; i <= 1; i += 2) {

   if ((allocated(0,0,i) && !(neighbor(0,0,i).pid < 0))) {
     Point neighbor = neighborp(0,0,i);
     if (faces) for (vector v = *faces, *_i21 = faces; ((scalar *)&v)->i >= 0; v = *++_i21) {
       scalar vn = VN;
       if (_attribute[vn.i].boundary[id])
  val(v.z,0,0,(i + 1)/2) = _attribute[vn.i].boundary[id](neighbor, point, v.z);
     }
   }

   else if (i == -1) {

     bool zn;
     Point neighbor = tangential_neighbor_z (point, &zn);
     if (neighbor.level >= 0) {
       int id = is_boundary_point (neighbor) ?
  (- neighbor(0,0,-1).pid - 1) : (- cell.pid - 1);
       if (faces) for (vector v = *faces, *_i22 = faces; ((scalar *)&v)->i >= 0; v = *++_i22) {



  scalar vt = zn ? VT : VR;

  val(v.z,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.z);
       }
     }
     else

       if (faces) for (vector v = *faces, *_i23 = faces; ((scalar *)&v)->i >= 0; v = *++_i23)
  val(v.z,0,0,0) = 0.;
   }

 }}
    }
  } } end_foreach_boundary_level(); }

  pfree (scalars,__func__,__FILE__,__LINE__);
  pfree (vectors,__func__,__FILE__,__LINE__);
  pfree (faces,__func__,__FILE__,__LINE__);
  enable_fpe_for_mpi();
}



#undef VN
#undef VT
#define VN _attribute[s.i].v.x
#define VT _attribute[s.i].v.y

static double masked_average (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1336 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"

  double sum = 0., n = 0.;
   { foreach_child()
    if (!(cell.pid < 0) && val(s,0,0,0) != nodata)
      sum += val(s,0,0,0), n++; end_foreach_child(); }
  return n ? sum/n : nodata;
}


#line 1344

static double masked_average_x (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1346 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"

  double sum = 0., n = 0.;
   { foreach_child()
    if (child.x < 0 && (!(cell.pid < 0) || !(neighbor(1,0,0).pid < 0)) &&
 val(s,1,0,0) != nodata)
      sum += val(s,1,0,0), n++; end_foreach_child(); }
  return n ? sum/n : nodata;
}
#line 1344

static double masked_average_y (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1346 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"

  double sum = 0., n = 0.;
   { foreach_child()
    if (child.y < 0 && (!(cell.pid < 0) || !(neighbor(0,1,0).pid < 0)) &&
 val(s,0,1,0) != nodata)
      sum += val(s,0,1,0), n++; end_foreach_child(); }
  return n ? sum/n : nodata;
}
#line 1344

static double masked_average_z (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1346 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"

  double sum = 0., n = 0.;
   { foreach_child()
    if (child.z < 0 && (!(cell.pid < 0) || !(neighbor(0,0,1).pid < 0)) &&
 val(s,0,0,1) != nodata)
      sum += val(s,0,0,1), n++; end_foreach_child(); }
  return n ? sum/n : nodata;
}

static void masked_boundary_restriction (const Boundary * b,
      scalar * list, int l)
{
  scalar * scalars = NULL;
  vector * faces = NULL;
  if (list) for (scalar s = *list, *_i24 = list; ((scalar *)&s)->i >= 0; s = *++_i24)
    if (!is_constant(s) && _attribute[s.i].refine != no_restriction) {
      if (_attribute[s.i].v.x.i == s.i && _attribute[s.i].face)
 faces = vectors_add (faces, _attribute[s.i].v);
      else
 scalars = list_add (scalars, s);
    }

   { foreach_halo (restriction, l){

#line 1368 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"
 {
    if (scalars) for (scalar s = *scalars, *_i25 = scalars; ((scalar *)&s)->i >= 0; s = *++_i25)
      val(s,0,0,0) = masked_average (parent, s);
    if (faces) for (vector v = *faces, *_i26 = faces; ((scalar *)&v)->i >= 0; v = *++_i26)
      {
#line 1372
 {
 double average = masked_average_x (parent, v.x);
 if ((neighbor(-1,0,0).pid < 0))
   val(v.x,0,0,0) = average;
 if ((neighbor(1,0,0).pid < 0))
   val(v.x,1,0,0) = average;
      }
#line 1372
 {
 double average = masked_average_y (parent, v.y);
 if ((neighbor(0,-1,0).pid < 0))
   val(v.y,0,0,0) = average;
 if ((neighbor(0,1,0).pid < 0))
   val(v.y,0,1,0) = average;
      }
#line 1372
 {
 double average = masked_average_z (parent, v.z);
 if ((neighbor(0,0,-1).pid < 0))
   val(v.z,0,0,0) = average;
 if ((neighbor(0,0,1).pid < 0))
   val(v.z,0,0,1) = average;
      }}
  } } end_foreach_halo(); }

  pfree (scalars,__func__,__FILE__,__LINE__);
  pfree (faces,__func__,__FILE__,__LINE__);
}
#line 1410 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"
static double periodic_bc (Point, Point, scalar);


#line 1412

static void periodic_boundary_level_x (const Boundary * b, scalar * list, int l)
{
  scalar * list1 = NULL;
  if (list) for (scalar s = *list, *_i27 = list; ((scalar *)&s)->i >= 0; s = *++_i27)
    if (!is_constant(s)) {
      if (_attribute[s.i].face) {

 scalar vt = VT;
 if (_attribute[vt.i].boundary[right] == periodic_bc)
   list1 = list_add (list1, s);

      }
      else if (_attribute[s.i].boundary[right] == periodic_bc) {
 if (_attribute[s.i].v.x.i >= 0)
   {
#line 1427

     list1 = list_add (list1, _attribute[s.i].v.x);
#line 1427

     list1 = list_add (list1, _attribute[s.i].v.y);
#line 1427

     list1 = list_add (list1, _attribute[s.i].v.z);}
 else
   list1 = list_add (list1, s);
      }
    }
  if (!list1)
    return;

  if (l == 0) {
    Point point = {0};  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1437 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"

    int n = 1 << point.level;

    for (int j = 0; j < n + 2*2; j++)


      for (int k = 0; k < n + 2*2; k++)

 {
   for (int i = 0; i < 2; i++)
     if (list1) for (scalar s = *list1, *_i28 = list1; ((scalar *)&s)->i >= 0; s = *++_i28)
       val(s,i,j,k) = val(s,2,2,2);
   for (int i = n + 2; i < n + 2*2; i++)
     if (list1) for (scalar s = *list1, *_i29 = list1; ((scalar *)&s)->i >= 0; s = *++_i29)
       val(s,i,j,k) = val(s,2,2,2);
 }
  }
  else {
    OMP_PARALLEL() {



      Point point = {0,0,0};  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1459 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"


      point.level = l < 0 ? depth() : l;
      int n = 1 << point.level;

      int j;
      OMP(omp for schedule(static))
 for (j = 0; j < n + 2*2; j++)


   for (int k = 0; k < n + 2*2; k++)

     {
       for (int i = 0; i < 2; i++)
  if (allocated(i + n,j,k))
    if (list1) for (scalar s = *list1, *_i30 = list1; ((scalar *)&s)->i >= 0; s = *++_i30)
      val(s,i,j,k) = val(s,i + n,j,k);
       for (int i = n + 2; i < n + 2*2; i++)
  if (allocated(i - n,j,k))
    if (list1) for (scalar s = *list1, *_i31 = list1; ((scalar *)&s)->i >= 0; s = *++_i31)
      val(s,i,j,k) = val(s,i - n,j,k);
     }
    }
  }

  pfree (list1,__func__,__FILE__,__LINE__);
}
#line 1412

static void periodic_boundary_level_y (const Boundary * b, scalar * list, int l)
{
  scalar * list1 = NULL;
  if (list) for (scalar s = *list, *_i27 = list; ((scalar *)&s)->i >= 0; s = *++_i27)
    if (!is_constant(s)) {
      if (_attribute[s.i].face) {

 scalar vt = VT;
 if (_attribute[vt.i].boundary[top] == periodic_bc)
   list1 = list_add (list1, s);

      }
      else if (_attribute[s.i].boundary[top] == periodic_bc) {
 if (_attribute[s.i].v.y.i >= 0)
   {
#line 1427

     list1 = list_add (list1, _attribute[s.i].v.y);
#line 1427

     list1 = list_add (list1, _attribute[s.i].v.z);
#line 1427

     list1 = list_add (list1, _attribute[s.i].v.x);}
 else
   list1 = list_add (list1, s);
      }
    }
  if (!list1)
    return;

  if (l == 0) {
    Point point = {0};  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1437 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"

    int n = 1 << point.level;

    for (int j = 0; j < n + 2*2; j++)


      for (int k = 0; k < n + 2*2; k++)

 {
   for (int i = 0; i < 2; i++)
     if (list1) for (scalar s = *list1, *_i28 = list1; ((scalar *)&s)->i >= 0; s = *++_i28)
       val(s,k,i,j) = val(s,2,2,2);
   for (int i = n + 2; i < n + 2*2; i++)
     if (list1) for (scalar s = *list1, *_i29 = list1; ((scalar *)&s)->i >= 0; s = *++_i29)
       val(s,k,i,j) = val(s,2,2,2);
 }
  }
  else {
    OMP_PARALLEL() {



      Point point = {0,0,0};  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1459 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"


      point.level = l < 0 ? depth() : l;
      int n = 1 << point.level;

      int j;
      OMP(omp for schedule(static))
 for (j = 0; j < n + 2*2; j++)


   for (int k = 0; k < n + 2*2; k++)

     {
       for (int i = 0; i < 2; i++)
  if (allocated(k,i + n,j))
    if (list1) for (scalar s = *list1, *_i30 = list1; ((scalar *)&s)->i >= 0; s = *++_i30)
      val(s,k,i,j) = val(s,k,i + n,j);
       for (int i = n + 2; i < n + 2*2; i++)
  if (allocated(k,i - n,j))
    if (list1) for (scalar s = *list1, *_i31 = list1; ((scalar *)&s)->i >= 0; s = *++_i31)
      val(s,k,i,j) = val(s,k,i - n,j);
     }
    }
  }

  pfree (list1,__func__,__FILE__,__LINE__);
}
#line 1412

static void periodic_boundary_level_z (const Boundary * b, scalar * list, int l)
{
  scalar * list1 = NULL;
  if (list) for (scalar s = *list, *_i27 = list; ((scalar *)&s)->i >= 0; s = *++_i27)
    if (!is_constant(s)) {
      if (_attribute[s.i].face) {

 scalar vt = VT;
 if (_attribute[vt.i].boundary[front] == periodic_bc)
   list1 = list_add (list1, s);

      }
      else if (_attribute[s.i].boundary[front] == periodic_bc) {
 if (_attribute[s.i].v.z.i >= 0)
   {
#line 1427

     list1 = list_add (list1, _attribute[s.i].v.z);
#line 1427

     list1 = list_add (list1, _attribute[s.i].v.x);
#line 1427

     list1 = list_add (list1, _attribute[s.i].v.y);}
 else
   list1 = list_add (list1, s);
      }
    }
  if (!list1)
    return;

  if (l == 0) {
    Point point = {0};  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1437 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"

    int n = 1 << point.level;

    for (int j = 0; j < n + 2*2; j++)


      for (int k = 0; k < n + 2*2; k++)

 {
   for (int i = 0; i < 2; i++)
     if (list1) for (scalar s = *list1, *_i28 = list1; ((scalar *)&s)->i >= 0; s = *++_i28)
       val(s,j,k,i) = val(s,2,2,2);
   for (int i = n + 2; i < n + 2*2; i++)
     if (list1) for (scalar s = *list1, *_i29 = list1; ((scalar *)&s)->i >= 0; s = *++_i29)
       val(s,j,k,i) = val(s,2,2,2);
 }
  }
  else {
    OMP_PARALLEL() {



      Point point = {0,0,0};  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1459 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"


      point.level = l < 0 ? depth() : l;
      int n = 1 << point.level;

      int j;
      OMP(omp for schedule(static))
 for (j = 0; j < n + 2*2; j++)


   for (int k = 0; k < n + 2*2; k++)

     {
       for (int i = 0; i < 2; i++)
  if (allocated(j,k,i + n))
    if (list1) for (scalar s = *list1, *_i30 = list1; ((scalar *)&s)->i >= 0; s = *++_i30)
      val(s,j,k,i) = val(s,j,k,i + n);
       for (int i = n + 2; i < n + 2*2; i++)
  if (allocated(j,k,i - n))
    if (list1) for (scalar s = *list1, *_i31 = list1; ((scalar *)&s)->i >= 0; s = *++_i31)
      val(s,j,k,i) = val(s,j,k,i - n);
     }
    }
  }

  pfree (list1,__func__,__FILE__,__LINE__);
}

#undef VN
#undef VT

static void free_cache (CacheLevel * c)
{
  for (int l = 0; l <= depth(); l++)
    pfree (c[l].p,__func__,__FILE__,__LINE__);
  pfree (c,__func__,__FILE__,__LINE__);
}

void free_grid (void)
{
  if (!grid)
    return;
  free_boundaries();
  Tree * q = ((Tree *)grid);
  pfree (q->leaves.p,__func__,__FILE__,__LINE__);
  pfree (q->faces.p,__func__,__FILE__,__LINE__);
  pfree (q->vertices.p,__func__,__FILE__,__LINE__);
  pfree (q->refined.p,__func__,__FILE__,__LINE__);


  Layer * L = q->L[0];
#line 1527 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"
  for (int i = 0; i < L->len; i++) {
    for (int j = 0; j < L->len; j++) {
      for (int k = 0; k < L->len; k++)
 pfree (L->m[i][j][k],__func__,__FILE__,__LINE__);
      pfree (L->m[i][j],__func__,__FILE__,__LINE__);
    }
    pfree (L->m[i],__func__,__FILE__,__LINE__);
  }

  for (int l = 1; l <= depth(); l++) {
    Layer * L = q->L[l];
    for (int i = 0; i < L->len; i++)
      if (L->m[i]) {
 for (int j = 0; j < L->len; j++)
   if (L->m[i][j])
     pfree (L->m[i][j],__func__,__FILE__,__LINE__);
 pfree (L->m[i],__func__,__FILE__,__LINE__);
      }
  }

  for (int l = 0; l <= depth(); l++)
    destroy_layer (q->L[l]);
  q->L = &(q->L[-1]);
  pfree (q->L,__func__,__FILE__,__LINE__);
  free_cache (q->active);
  free_cache (q->prolongation);
  free_cache (q->boundary);
  free_cache (q->restriction);
  pfree (q,__func__,__FILE__,__LINE__);
  grid = NULL;
}

static void refine_level (int depth);


void init_grid (int n)
{ trace ("init_grid", "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h", 1563);

  assert (sizeof(Cell) % 8 == 0);

  free_grid();
  int depth = 0;
  while (n > 1) {
    if (n % 2) {
      fprintf (qstderr(), "tree: N must be a power-of-two\n");
      exit (1);
    }
    n /= 2;
    depth++;
  }
  Tree * q = ((Tree *) pcalloc (1, sizeof(Tree),__func__,__FILE__,__LINE__));
  grid = (Grid *) q;
  grid->depth = 0;


  q->L = ((Layer * *) pmalloc ((2)*sizeof(Layer *),__func__,__FILE__,__LINE__));

  q->L[0] = NULL; q->L = &(q->L[1]);

  Layer * L = new_layer (0);
  q->L[0] = L;
#line 1619 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"
  for (int i = 0; i < L->len; i++)
    for (int j = 0; j < L->len; j++) {
      layer_add_row (L, i, j);
      for (int k = 0; k < L->len; k++)
 L->m[i][j][k] = (char *) pcalloc (1, sizeof(Cell) + datasize,__func__,__FILE__,__LINE__);
    }
  CELL(L->m[2][2][2]).flags |= leaf;
  if (pid() == 0)
    CELL(L->m[2][2][2]).flags |= active;
  for (int k = -2; k <= 2; k++)
    for (int l = -2; l <= 2; l++)
      for (int n = -2; n <= 2; n++)
 CELL(L->m[2 +k][2 +l][2 +n]).pid =
   (k > 0 ? -1 - right :
    k < 0 ? -1 - left :
    l > 0 ? -1 - top :
    l < 0 ? -1 - bottom :
    n > 0 ? -1 - front :
    n < 0 ? -1 - back :
    0);

  q->active = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,__LINE__));
  q->prolongation = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,__LINE__));
  q->boundary = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,__LINE__));
  q->restriction = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,__LINE__));
  q->dirty = true;
  N = 1 << depth;
#if 1
  void mpi_boundary_new();
  mpi_boundary_new();
#endif
  { if (((Tree *)grid)->dirty) update_cache_f(); };

  Boundary * b = ((Boundary *) pcalloc (1, sizeof(Boundary),__func__,__FILE__,__LINE__));
  b->level = box_boundary_level;
  b->restriction = masked_boundary_restriction;
  add_boundary (b);

  {
#line 1657
 {
    Boundary * b = ((Boundary *) pcalloc (1, sizeof(Boundary),__func__,__FILE__,__LINE__));
    b->level = periodic_boundary_level_x;
    add_boundary (b);
  }
#line 1657
 {
    Boundary * b = ((Boundary *) pcalloc (1, sizeof(Boundary),__func__,__FILE__,__LINE__));
    b->level = periodic_boundary_level_y;
    add_boundary (b);
  }
#line 1657
 {
    Boundary * b = ((Boundary *) pcalloc (1, sizeof(Boundary),__func__,__FILE__,__LINE__));
    b->level = periodic_boundary_level_z;
    add_boundary (b);
  }}
  refine_level (depth);
  reset (all, 0.);
 end_trace("init_grid", "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h", 1664); }
#line 1700 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"
struct _locate { double x, y, z; };

Point locate (struct _locate p)
{
  for (int l = depth(); l >= 0; l--) {
    Point point = { .level = l };
    int n = 1 << point.level;
    point.i = (p.x - X0)/L0*n + 2;

    point.j = (p.y - Y0)/L0*n + 2;


    point.k = (p.z - Z0)/L0*n + 2;

    if (point.i >= 0 && point.i < n + 2*2

 && point.j >= 0 && point.j < n + 2*2


 && point.k >= 0 && point.k < n + 2*2

 ) {
      if (allocated(0,0,0) && is_local(cell) && is_leaf(cell))
 return point;
    }
    else
      break;
  }
  Point point = { .level = -1 };
  return point;
}



bool tree_is_full()
{
  { if (((Tree *)grid)->dirty) update_cache_f(); };
  return (grid->tn == 1L << grid->maxdepth*3);
}

#line 1 "grid/tree-common.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h"



#line 1 "grid/multigrid-common.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/grid/multigrid-common.h"


#line 1 "grid/cartesian-common.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/grid/cartesian-common.h"
#line 1 "grid/events.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/grid/events.h"







static void event_error (Event * ev, const char * s)
{
  fprintf (qstderr(), "%s:%d: error: %s\n", ev->file, ev->line, s);
  exit (1);
}

static void init_event (Event * ev)
{
  if (ev->arrayi || ev->arrayt) {
    ev->i = ev->t = -1;
    if (ev->arrayi)
      ev->i = ev->arrayi[0];
    else
      ev->t = ev->arrayt[0];
    ev->a = 1;
    ev->expr[1] = NULL;
  }
  else {
    if (ev->nexpr > 0) {
      Expr init = NULL, cond = NULL, inc = NULL;
      for (int j = 0; j < ev->nexpr; j++) {
 int i = -123456; double t = i;
 (* ev->expr[j]) (&i, &t, ev);
 if (i == -123456 && t == -123456) {

   if (cond)
     event_error (ev, "events can only use a single condition");
   cond = ev->expr[j];
 }
 else {

   int i1 = i; double t1 = t;
   (* ev->expr[j]) (&i1, &t1, ev);
   if (i1 == i && t1 == t) {


     if (init)
       event_error (ev, "events can only use a single initialisation");
     init = ev->expr[j];
   }
   else {

     if (inc)
       event_error (ev, "events can only use a single increment");
     inc = ev->expr[j];
   }
 }
      }
      ev->expr[0] = init;
      ev->expr[1] = cond;
      ev->expr[2] = inc;
      ev->nexpr = 0;
    }
    ev->i = ev->t = -1;
    if (ev->expr[0]) {
      (* ev->expr[0]) (&ev->i, &ev->t, ev);
      if (ev->i == 1234567890 || ev->t == 1234567890) {
 ev->i = 1234567890; ev->t = -1;
      }
    }
    else if (ev->expr[2]) {
      (* ev->expr[2]) (&ev->i, &ev->t, ev);
      if (ev->i != -1)
 ev->i = 0;
      if (ev->t != -1)
 ev->t = 0;
    }
  }
}

enum { event_done, event_alive, event_stop };

static int event_finished (Event * ev)
{
  ev->t = ev->i = -1;
  return event_done;
}

void event_register (Event event) {
  assert (Events);
  assert (!event.last);
  int n = 0, parent = -1;
  for (Event * ev = Events; !ev->last; ev++) {
    if (!strcmp (event.name, ev->name)) {
      assert (parent < 0);
      parent = n;
    }
    n++;
  }
  if (parent < 0) {
    Events = (Event *) prealloc (Events, (n + 2)*sizeof(Event),__func__,__FILE__,__LINE__);
    Events[n] = event;
    Events[n].next = NULL;
    Events[n + 1].last = true;
    init_event (&Events[n]);
  }
  else {
    Event * ev = ((Event *) pcalloc (1, sizeof(Event),__func__,__FILE__,__LINE__));
    *ev = Events[parent];
    Events[parent] = event;
    Events[parent].next = ev;
    init_event (&Events[parent]);
  }
}

static int event_cond (Event * ev, int i, double t)
{
  if (!ev->expr[1])
    return true;
  return (* ev->expr[1]) (&i, &t, ev);
}
#line 131 "/home/miro3/Documents/Programming/basilisk/src/grid/events.h"
static int event_do (Event * ev, bool action)
{
  if ((iter > ev->i && t > ev->t) || !event_cond (ev, iter, t))
    return event_finished (ev);
  if (iter == ev->i || fabs (t - ev->t) <= 1e-9) {
    if (action) {
      bool finished = false;
      for (Event * e = ev; e; e = e->next) {



 if ((* e->action) (iter, t, e))
   finished = true;
      }
      if (finished) {
 event_finished (ev);
 return event_stop;
      }
    }
    if (ev->arrayi) {
      ev->i = ev->arrayi[ev->a++];
      if (ev->i < 0)
 return event_finished (ev);
    }
    if (ev->arrayt) {
      ev->t = ev->arrayt[ev->a++];
      if (ev->t < 0)
 return event_finished (ev);
    }
    else if (ev->expr[2]) {
      int i0 = ev->i;
      (* ev->expr[2]) (&ev->i, &ev->t, ev);
      if (i0 == -1 && ev->i != i0)
 ev->i += iter + 1;
      if (!event_cond (ev, iter + 1, ev->t))
 return event_finished (ev);
    }
    else if (ev->expr[0] && !ev->expr[1])
      return event_finished (ev);
  }
  return event_alive;
}

static void end_event_do (bool action)
{




  for (Event * ev = Events; !ev->last; ev++)
    if (ev->i == 1234567890 && action)
      for (Event * e = ev; e; e = e->next) {



 e->action (iter, t, e);
      }
}

int events (bool action)
{





  if (iter == 0)
    for (Event * ev = Events; !ev->last; ev++)
      init_event (ev);

  int cond = 0, cond1 = 0;
  inext = 1234567890; tnext = HUGE;
  for (Event * ev = Events; !ev->last && !cond; ev++)
    if (ev->i != 1234567890 &&
 (ev->expr[1] || (ev->expr[0] && !ev->expr[1] && !ev->expr[2]) || ev->arrayi || ev->arrayt))
      cond = 1;
  for (Event * ev = Events; !ev->last; ev++) {
    int status = event_do (ev, action);
    if (status == event_stop) {
      end_event_do (action);
      return 0;
    }
    if (status == event_alive && ev->i != 1234567890 &&
 (ev->expr[1] || (ev->expr[0] && !ev->expr[1] && !ev->expr[2]) || ev->arrayi || ev->arrayt))
      cond1 = 1;
    if (ev->t > t && ev->t < tnext)
      tnext = ev->t;
    if (ev->i > iter && ev->i < inext)
      inext = ev->i;
  }
  if ((!cond || cond1) && (tnext != HUGE || inext != 1234567890)) {
    inext = iter + 1;
    return 1;
  }
  end_event_do (action);
  return 0;
}

void event (const char * name)
{
  for (Event * ev = Events; !ev->last; ev++)
    if (!strcmp (ev->name, name))
      for (Event * e = ev; e; e = e->next) {



 (* e->action) (0, 0, e);
      }
}

double dtnext (double dt)
{
  if (tnext != HUGE && tnext > t) {
    unsigned int n = (tnext - t)/dt;
    assert (n < INT_MAX);
    if (n == 0)
      dt = tnext - t;
    else {
      double dt1 = (tnext - t)/n;
      if (dt1 > dt + 1e-9)
 dt = (tnext - t)/(n + 1);
      else if (dt1 < dt)
 dt = dt1;
      tnext = t + dt;
    }
  }
  else
    tnext = t + dt;
  return dt;
}
#line 2 "/home/miro3/Documents/Programming/basilisk/src/grid/cartesian-common.h"

void (* debug) (Point);

#define _val_constant(a,k,l,m) ((const double) _constant[a.i -_NVARMAX])

#undef VARIABLES
#define VARIABLES\
  double Delta = L0*(1./(1 << point.level));\
  \
    double Delta_x = Delta;\
    double Delta_y = Delta;\
    double Delta_z = Delta;\
\
  double x = (ig/2. + (point.i - 2) + 0.5)*Delta + X0; NOT_UNUSED(x);\
\
  double y = (jg/2. + (point.j - 2) + 0.5)*Delta + Y0;\
\
\
\
 NOT_UNUSED(y);\
\
  double z = (kg/2. + (point.k - 2) + 0.5)*Delta + Z0;\
\
\
\
  NOT_UNUSED(z);\
\
  NOT_UNUSED(Delta);\
  \
    NOT_UNUSED(Delta_x);\
    NOT_UNUSED(Delta_y);\
    NOT_UNUSED(Delta_z);\
\
  ;\

#line 32


#line 1 "grid/fpe.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/grid/fpe.h"


#include <signal.h>
#include <unistd.h>

static void gdb()
{
  if (last_point.level >= 0) {
    debug (last_point);
    fputc ('\n', qstderr());
    fflush (qstderr());
  }
  char command[80];
  sprintf (command, "exec xterm -e 'gdb -p %d' & xterm -e 'gnuplot plot -'",
    getpid());
  system (command);
}

static void caught_abort (int sig)
{
  fprintf (qstderr(), "Caught signal %d (Aborted)\n", sig);
  gdb();
}

static void caught_fpe (int sig)
{
  fprintf (qstderr(), "Caught signal %d (Floating Point Exception)\n", sig);
  gdb();
  exit (1);
}

static void caught_segfault (int sig)
{
  fprintf (qstderr(), "Caught signal %d (Segmentation Fault)\n", sig);
  gdb();
  exit (2);
}

void catch_fpe (void)
{
  struct sigaction act;
  act.sa_handler = caught_fpe;
  sigemptyset (&act.sa_mask);
  act.sa_flags = 0;
  last_point.level = -1;
  sigaction (8, &act, NULL);
  act.sa_handler = caught_segfault;
  sigaction (11, &act, NULL);
  act.sa_handler = caught_abort;
  act.sa_flags = SA_RESETHAND;
  sigaction (6, &act, NULL);
}
#line 35 "/home/miro3/Documents/Programming/basilisk/src/grid/cartesian-common.h"

#define end_foreach_face()

scalar new_scalar (const char * name)
{
  int nvar = datasize/sizeof(double);
  scalar s;
  for (s.i = 0; s.i < nvar; s.i++)
    if (!list_lookup (all, s)) {
      init_scalar (s, name);
      trash (((scalar []){s, {-1}}));
      all = list_append (all, s);
      return s;
    }


  assert (nvar < _NVARMAX);
  datasize += sizeof(double); nvar++;
  _attribute = (_Attributes *) prealloc (_attribute, (nvar)*sizeof(_Attributes),__func__,__FILE__,__LINE__);
  memset (&_attribute[nvar-1], 0, sizeof (_Attributes));
  s = (scalar){nvar - 1};
  realloc_scalar();
  init_scalar (s, name);
  trash (((scalar []){s, {-1}}));
  all = list_append (all, s);
  return s;
}

scalar new_vertex_scalar (const char * name)
{
  scalar s = new_scalar (name);
  {
#line 66

    _attribute[s.i].d.x = -1;
#line 66

    _attribute[s.i].d.y = -1;
#line 66

    _attribute[s.i].d.z = -1;}
  return s;
}

static vector alloc_vector (const char * name)
{
  vector v;
  char cname[strlen(name) + 3];
  struct { char * x, * y, * z; } ext = {"%s.x", "%s.y", "%s.z"};
  {
#line 76
 {
    sprintf (cname, ext.x, name);
    v.x = new_scalar (cname);
  }
#line 76
 {
    sprintf (cname, ext.y, name);
    v.y = new_scalar (cname);
  }
#line 76
 {
    sprintf (cname, ext.z, name);
    v.z = new_scalar (cname);
  }}
  return v;
}

vector new_vector (const char * name)
{
  vector v = alloc_vector (name);
  init_vector (v, NULL);
  return v;
}

vector new_face_vector (const char * name)
{
  vector v = alloc_vector (name);
  init_face_vector (v, NULL);
  return v;
}

tensor new_tensor (const char * name)
{
  char cname[strlen(name) + 3];
  struct { char * x, * y, * z; } ext = {"%s.x", "%s.y", "%s.z"};
  tensor t;
  {
#line 102
 {
    sprintf (cname, ext.x, name);
    t.x = new_vector (cname);
  }
#line 102
 {
    sprintf (cname, ext.y, name);
    t.y = new_vector (cname);
  }
#line 102
 {
    sprintf (cname, ext.z, name);
    t.z = new_vector (cname);
  }}
  init_tensor (t, NULL);
  return t;
}

tensor new_symmetric_tensor (const char * name)
{
  char cname[strlen(name) + 5];
  struct { char * x, * y, * z; } ext = {"%s.x.x", "%s.y.y", "%s.z.z"};
  tensor t;
  {
#line 115
 {
    sprintf (cname, ext.x, name);
    t.x.x = new_scalar(cname);
  }
#line 115
 {
    sprintf (cname, ext.y, name);
    t.y.y = new_scalar(cname);
  }
#line 115
 {
    sprintf (cname, ext.z, name);
    t.z.z = new_scalar(cname);
  }}

    sprintf (cname, "%s.x.y", name);
    t.x.y = new_scalar(cname);
    t.y.x = t.x.y;


    sprintf (cname, "%s.x.z", name);
    t.x.z = new_scalar(cname);
    t.z.x = t.x.z;
    sprintf (cname, "%s.y.z", name);
    t.y.z = new_scalar(cname);
    t.z.y = t.y.z;




  init_tensor (t, NULL);
  return t;
}

static int nconst = 0;

void init_const_scalar (scalar s, const char * name, double val)
{
  if (s.i - _NVARMAX >= nconst) {
    nconst = s.i - _NVARMAX + 1;
    _constant = (double *) prealloc (_constant, (nconst)*sizeof(double),__func__,__FILE__,__LINE__);
  }
  _constant[s.i - _NVARMAX] = val;
}

scalar new_const_scalar (const char * name, int i, double val)
{
  scalar s = (scalar){i + _NVARMAX};
  init_const_scalar (s, name, val);
  return s;
}

void init_const_vector (vector v, const char * name, double * val)
{
  {
#line 159

    init_const_scalar (v.x, name, *val++);
#line 159

    init_const_scalar (v.y, name, *val++);
#line 159

    init_const_scalar (v.z, name, *val++);}
}

vector new_const_vector (const char * name, int i, double * val)
{
  vector v;
  {
#line 166

    v.x.i = _NVARMAX + i++;
#line 166

    v.y.i = _NVARMAX + i++;
#line 166

    v.z.i = _NVARMAX + i++;}
  init_const_vector (v, name, val);
  return v;
}

void scalar_clone (scalar a, scalar b)
{
  char * name = _attribute[a.i].name;
  double (** boundary) (Point, Point, scalar) = _attribute[a.i].boundary;
  double (** boundary_homogeneous) (Point, Point, scalar) =
    _attribute[a.i].boundary_homogeneous;
  _attribute[a.i] = _attribute[b.i];
  _attribute[a.i].name = name;
  _attribute[a.i].boundary = boundary;
  _attribute[a.i].boundary_homogeneous = boundary_homogeneous;
  for (int i = 0; i < nboundary; i++) {
    _attribute[a.i].boundary[i] = _attribute[b.i].boundary[i];
    _attribute[a.i].boundary_homogeneous[i] = _attribute[b.i].boundary_homogeneous[i];
  }
}

scalar * list_clone (scalar * l)
{
  scalar * list = NULL;
  int nvar = datasize/sizeof(double), map[nvar];
  for (int i = 0; i < nvar; i++)
    map[i] = -1;
  if (l) for (scalar s = *l, *_i32 = l; ((scalar *)&s)->i >= 0; s = *++_i32) {
    scalar c = new_scalar("c");
    scalar_clone (c, s);
    map[s.i] = c.i;
    list = list_append (list, c);
  }
  if (list) for (scalar s = *list, *_i33 = list; ((scalar *)&s)->i >= 0; s = *++_i33)
    {
#line 201

      if (_attribute[s.i].v.x.i >= 0 && map[_attribute[s.i].v.x.i] >= 0)
 _attribute[s.i].v.x.i = map[_attribute[s.i].v.x.i];
#line 201

      if (_attribute[s.i].v.y.i >= 0 && map[_attribute[s.i].v.y.i] >= 0)
 _attribute[s.i].v.y.i = map[_attribute[s.i].v.y.i];
#line 201

      if (_attribute[s.i].v.z.i >= 0 && map[_attribute[s.i].v.z.i] >= 0)
 _attribute[s.i].v.z.i = map[_attribute[s.i].v.z.i];}
  return list;
}

void delete (scalar * list)
{
  if (all == NULL)
    return;

  if (list) for (scalar f = *list, *_i34 = list; ((scalar *)&f)->i >= 0; f = *++_i34) {
    if (_attribute[f.i].delete)
      _attribute[f.i].delete (f);
    pfree (_attribute[f.i].name,__func__,__FILE__,__LINE__); _attribute[f.i].name = NULL;
    pfree (_attribute[f.i].boundary,__func__,__FILE__,__LINE__); _attribute[f.i].boundary = NULL;
    pfree (_attribute[f.i].boundary_homogeneous,__func__,__FILE__,__LINE__); _attribute[f.i].boundary_homogeneous = NULL;
  }

  if (list == all) {
    all[0].i = -1;
    return;
  }

  trash (list);
  if (list) for (scalar f = *list, *_i35 = list; ((scalar *)&f)->i >= 0; f = *++_i35) {
    scalar * s = all;
    for (; s->i >= 0 && s->i != f.i; s++);
    if (s->i == f.i)
      for (; s->i >= 0; s++)
 s[0] = s[1];
  }
}

typedef void (* free_solver_func) (void);

static Array * free_solver_funcs = NULL;

void free_solver_func_add (free_solver_func func)
{
  if (!free_solver_funcs)
    free_solver_funcs = array_new();
  array_append (free_solver_funcs, &func, sizeof(free_solver_func));
}

void free_solver()
{
  if (free_solver_funcs) {
    free_solver_func * a = (free_solver_func *) free_solver_funcs->p;
    for (int i = 0; i < free_solver_funcs->len/sizeof(free_solver_func); i++)
      a[i] ();
    array_free (free_solver_funcs);
  }

  delete (all);
  pfree (all,__func__,__FILE__,__LINE__); all = NULL;
  for (Event * ev = Events; !ev->last; ev++) {
    Event * e = ev->next;
    while (e) {
      Event * next = e->next;
      pfree (e,__func__,__FILE__,__LINE__);
      e = next;
    }
  }

  pfree (Events,__func__,__FILE__,__LINE__); Events = NULL;
  pfree (_attribute,__func__,__FILE__,__LINE__); _attribute = NULL;
  pfree (_constant,__func__,__FILE__,__LINE__); _constant = NULL;
  free_grid();
  qpclose_all();
#if 2
  trace_off();
#endif
#if MTRACE
  pmuntrace();
#endif
#if _CADNA
  cadna_end();
#endif
}



void (* boundary_level) (scalar *, int l);
void (* boundary_flux) (vector *);


void boundary (scalar * list)
{ trace ("boundary", "/home/miro3/Documents/Programming/basilisk/src/grid/cartesian-common.h", 289);
  if (list == NULL)
    { ; end_trace("boundary", "/home/miro3/Documents/Programming/basilisk/src/grid/cartesian-common.h", 291);  return; }
  vector * listf = NULL;
  if (list) for (scalar s = *list, *_i36 = list; ((scalar *)&s)->i >= 0; s = *++_i36)
    if (!is_constant(s) && _attribute[s.i].face)
      listf = vectors_add (listf, _attribute[s.i].v);
  if (listf) {
    boundary_flux (listf);
    pfree (listf,__func__,__FILE__,__LINE__);
  }
  boundary_level (list, -1);
 end_trace("boundary", "/home/miro3/Documents/Programming/basilisk/src/grid/cartesian-common.h", 301); }

void cartesian_boundary_level (scalar * list, int l)
{
  { Boundary ** _i = boundaries, * _b; while ((_b = *_i++)) if (_b->level) _b->level (_b, list, l); };
}

void cartesian_boundary_flux (vector * list)
{

}

static double symmetry (Point point, Point neighbor, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 314 "/home/miro3/Documents/Programming/basilisk/src/grid/cartesian-common.h"

  return val(s,0,0,0);
}

static double antisymmetry (Point point, Point neighbor, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 319 "/home/miro3/Documents/Programming/basilisk/src/grid/cartesian-common.h"

  return -val(s,0,0,0);
}

double (* default_scalar_bc[]) (Point, Point, scalar) = {
  symmetry, symmetry, symmetry, symmetry, symmetry, symmetry
};

scalar cartesian_init_scalar (scalar s, const char * name)
{

  char * pname;
  if (name) {
    pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
    pname = pstrdup (name,__func__,__FILE__,__LINE__);
  }
  else
    pname = _attribute[s.i].name;
  pfree (_attribute[s.i].boundary,__func__,__FILE__,__LINE__);
  pfree (_attribute[s.i].boundary_homogeneous,__func__,__FILE__,__LINE__);

  _attribute[s.i] = (const _Attributes){0};
  _attribute[s.i].name = pname;

  _attribute[s.i].boundary = (double (**)(Point, Point, scalar))
    pmalloc (nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  _attribute[s.i].boundary_homogeneous = (double (**)(Point, Point, scalar))
    pmalloc (nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  for (int b = 0; b < nboundary; b++)
    _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b] =
      b < 2*3 ? default_scalar_bc[b] : symmetry;
  _attribute[s.i].gradient = NULL;
  {
#line 351
 {
    _attribute[s.i].d.x = 0;
    _attribute[s.i].v.x.i = -1;
  }
#line 351
 {
    _attribute[s.i].d.y = 0;
    _attribute[s.i].v.y.i = -1;
  }
#line 351
 {
    _attribute[s.i].d.z = 0;
    _attribute[s.i].v.z.i = -1;
  }}
  _attribute[s.i].face = false;
  return s;
}

double (* default_vector_bc[]) (Point, Point, scalar) = {
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry
};

vector cartesian_init_vector (vector v, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  {
#line 368
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.x);
      init_scalar (v.x, cname);
    }
    else
      init_scalar (v.x, NULL);
    _attribute[v.x.i].v = v;
  }
#line 368
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.y);
      init_scalar (v.y, cname);
    }
    else
      init_scalar (v.y, NULL);
    _attribute[v.y.i].v = v;
  }
#line 368
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.z);
      init_scalar (v.z, cname);
    }
    else
      init_scalar (v.z, NULL);
    _attribute[v.z.i].v = v;
  }}

  for (int d = 0; d < nboundary; d++)
    _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] =
      d < 2*3 ? default_vector_bc[d] : antisymmetry;
  return v;
}

vector cartesian_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_vector (v, name);
  {
#line 388
 {
    _attribute[v.x.i].d.x = -1;
    _attribute[v.x.i].face = true;
  }
#line 388
 {
    _attribute[v.y.i].d.y = -1;
    _attribute[v.y.i].face = true;
  }
#line 388
 {
    _attribute[v.z.i].d.z = -1;
    _attribute[v.z.i].face = true;
  }}
  for (int d = 0; d < nboundary; d++)
    _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] = NULL;
  return v;
}

tensor cartesian_init_tensor (tensor t, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  {
#line 400
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.x);
      init_vector (t.x, cname);
    }
    else
      init_vector (t.x, NULL);
  }
#line 400
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.y);
      init_vector (t.y, cname);
    }
    else
      init_vector (t.y, NULL);
  }
#line 400
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.z);
      init_vector (t.z, cname);
    }
    else
      init_vector (t.z, NULL);
  }}
#line 424 "/home/miro3/Documents/Programming/basilisk/src/grid/cartesian-common.h"
    assert (false);

  return t;
}

void output_cells (FILE * fp)
{
   { foreach(){

#line 431 "/home/miro3/Documents/Programming/basilisk/src/grid/cartesian-common.h"
 {
    Delta /= 2.;
#line 443 "/home/miro3/Documents/Programming/basilisk/src/grid/cartesian-common.h"
      for (int i = -1; i <= 1; i += 2) {
 fprintf (fp, "%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n\n",
   x - Delta, y - Delta, z + i*Delta,
   x - Delta, y + Delta, z + i*Delta,
   x + Delta, y + Delta, z + i*Delta,
   x + Delta, y - Delta, z + i*Delta,
   x - Delta, y - Delta, z + i*Delta);
 for (int j = -1; j <= 1; j += 2)
   fprintf (fp, "%g %g %g\n%g %g %g\n\n",
     x + i*Delta, y + j*Delta, z - Delta,
     x + i*Delta, y + j*Delta, z + Delta);
      }

  } } end_foreach(); }
  fflush (fp);
}

static char * replace_ (const char * vname)
{
  char * name = pstrdup (vname,__func__,__FILE__,__LINE__), * c = name;
  while (*c != '\0') {
    if (*c == '.')
      *c = '_';
    c++;
  }
  return name;
}

static void debug_plot (FILE * fp, const char * name, const char * cells,
   const char * stencil)
{
  char * vname = replace_ (name);
  fprintf (fp,
    "  load 'debug.plot'\n"
    "  v=%s\n"







    "  splot '%s' w l lc 0, "
    "'%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 1"
           " title columnhead(4+4*v)",

    vname, cells, stencil);
  pfree (vname,__func__,__FILE__,__LINE__);
}

void cartesian_debug (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 494 "/home/miro3/Documents/Programming/basilisk/src/grid/cartesian-common.h"

  char name[80] = "cells";
  if (pid() > 0)
    sprintf (name, "cells-%d", pid());
  FILE * fp = fopen (name, "w");
  output_cells (fp);
  fclose (fp);

  char stencil[80] = "stencil";
  if (pid() > 0)
    sprintf (stencil, "stencil-%d", pid());
  fp = fopen (stencil, "w");
  if (all) for (scalar v = *all, *_i37 = all; ((scalar *)&v)->i >= 0; v = *++_i37)





    fprintf (fp, "x y z %s ", _attribute[v.i].name);

  fputc ('\n', fp);
#line 541 "/home/miro3/Documents/Programming/basilisk/src/grid/cartesian-common.h"
    for (int k = -2; k <= 2; k++)
      for (int l = -2; l <= 2; l++)
 for (int m = -2; m <= 2; m++) {
   if (all) for (scalar v = *all, *_i38 = all; ((scalar *)&v)->i >= 0; v = *++_i38) {
     fprintf (fp, "%g %g %g ",
       x + k*Delta + _attribute[v.i].d.x*Delta/2.,
       y + l*Delta + _attribute[v.i].d.y*Delta/2.,
       z + m*Delta + _attribute[v.i].d.z*Delta/2.);
     if (allocated(k,l,m))
       fprintf (fp, "%g ", val(v,k,l,m));
     else
       fputs ("n/a ", fp);
   }
   fputc ('\n', fp);
 }

  fclose (fp);

  fp = fopen ("debug.plot", "w");
  fprintf (fp,
    "set term x11\n"
    "set size ratio -1\n"
    "set key outside\n");
  if (all) for (scalar s = *all, *_i39 = all; ((scalar *)&s)->i >= 0; s = *++_i39) {
    char * name = replace_ (_attribute[s.i].name);
    fprintf (fp, "%s = %d\n", name, s.i);
    pfree (name,__func__,__FILE__,__LINE__);
  }
  fclose (fp);

  fprintf (qstderr(), "Last point stencils can be displayed using (in gnuplot)\n");
  debug_plot (qstderr(), _attribute[0].name, name, stencil);
  fflush (qstderr());

  fp = fopen ("plot", "w");
  debug_plot (fp, _attribute[0].name, name, stencil);
  fclose (fp);
}

void cartesian_methods()
{
  init_scalar = cartesian_init_scalar;
  init_vector = cartesian_init_vector;
  init_tensor = cartesian_init_tensor;
  init_face_vector = cartesian_init_face_vector;
  boundary_level = cartesian_boundary_level;
  boundary_flux = cartesian_boundary_flux;
  debug = cartesian_debug;
}

struct _interpolate {
  scalar v;
  double x, y, z;
};

static double interpolate_linear (Point point, struct _interpolate p)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 597 "/home/miro3/Documents/Programming/basilisk/src/grid/cartesian-common.h"

  scalar v = p.v;
#line 614 "/home/miro3/Documents/Programming/basilisk/src/grid/cartesian-common.h"
  x = (p.x - x)/Delta - _attribute[v.i].d.x/2.;
  y = (p.y - y)/Delta - _attribute[v.i].d.y/2.;
  z = (p.z - z)/Delta - _attribute[v.i].d.z/2.;
  int i = sign(x), j = sign(y), k = sign(z);
  x = fabs(x); y = fabs(y); z = fabs(z);

  return (((val(v,0,0,0)*(1. - x) + val(v,i,0,0)*x)*(1. - y) +
    (val(v,0,j,0)*(1. - x) + val(v,i,j,0)*x)*y)*(1. - z) +
   ((val(v,0,0,k)*(1. - x) + val(v,i,0,k)*x)*(1. - y) +
    (val(v,0,j,k)*(1. - x) + val(v,i,j,k)*x)*y)*z);

}


double interpolate (struct _interpolate p)
{ trace ("interpolate", "/home/miro3/Documents/Programming/basilisk/src/grid/cartesian-common.h", 629);
  Point point = locate ((struct _locate){p.x, p.y, p.z});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 630 "/home/miro3/Documents/Programming/basilisk/src/grid/cartesian-common.h"

  if (point.level < 0)
    { double _ret =  nodata; end_trace("interpolate", "/home/miro3/Documents/Programming/basilisk/src/grid/cartesian-common.h", 632);  return _ret; }
  { double _ret =  interpolate_linear (point, p); end_trace("interpolate", "/home/miro3/Documents/Programming/basilisk/src/grid/cartesian-common.h", 633);  return _ret; }
 end_trace("interpolate", "/home/miro3/Documents/Programming/basilisk/src/grid/cartesian-common.h", 634); }


void interpolate_array (scalar * list, coord * a, int n, double * v, bool linear)
{ trace ("interpolate_array", "/home/miro3/Documents/Programming/basilisk/src/grid/cartesian-common.h", 638);
  int j = 0;
  for (int i = 0; i < n; i++) {
    Point point = locate ((struct _locate){a[i].x, a[i].y, a[i].z});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 641 "/home/miro3/Documents/Programming/basilisk/src/grid/cartesian-common.h"

    if (point.level >= 0) {
      if (list) for (scalar s = *list, *_i40 = list; ((scalar *)&s)->i >= 0; s = *++_i40)
 v[j++] = !linear ? val(s,0,0,0) :
   interpolate_linear (point,
         (struct _interpolate){s, a[i].x, a[i].y, a[i].z});
    }
    else
      if (list) for (scalar s = *list, *_i41 = list; ((scalar *)&s)->i >= 0; s = *++_i41)
 v[j++] = nodata;
  }
#if 1
  if (pid() == 0)
    MPI_Reduce (MPI_IN_PLACE, v, n*list_len(list), MPI_DOUBLE,
  MPI_MIN, 0, MPI_COMM_WORLD);
  else
    MPI_Reduce (v, v, n*list_len(list), MPI_DOUBLE,
  MPI_MIN, 0, MPI_COMM_WORLD);
#endif
 end_trace("interpolate_array", "/home/miro3/Documents/Programming/basilisk/src/grid/cartesian-common.h", 660); }



typedef int bid;

bid new_bid()
{
  int b = nboundary++;
  if (all) for (scalar s = *all, *_i42 = all; ((scalar *)&s)->i >= 0; s = *++_i42) {
    _attribute[s.i].boundary = (double (**)(Point, Point, scalar))
      prealloc (_attribute[s.i].boundary, nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
    _attribute[s.i].boundary_homogeneous = (double (**)(Point, Point, scalar))
      prealloc (_attribute[s.i].boundary_homogeneous, nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  }
  if (all) for (scalar s = *all, *_i43 = all; ((scalar *)&s)->i >= 0; s = *++_i43) {
    if (_attribute[s.i].v.x.i < 0)
      _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b] = symmetry;
    else if (_attribute[s.i].v.x.i == s.i) {
      vector v = _attribute[s.i].v;
      {
#line 680

 _attribute[v.y.i].boundary[b] = _attribute[v.y.i].boundary_homogeneous[b] = symmetry;
#line 680

 _attribute[v.z.i].boundary[b] = _attribute[v.z.i].boundary_homogeneous[b] = symmetry;
#line 680

 _attribute[v.x.i].boundary[b] = _attribute[v.x.i].boundary_homogeneous[b] = symmetry;}
      _attribute[v.x.i].boundary[b] = _attribute[v.x.i].boundary_homogeneous[b] =
 _attribute[v.x.i].face ? NULL : antisymmetry;
    }
  }
  return b;
}



static double periodic_bc (Point point, Point neighbor, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 692 "/home/miro3/Documents/Programming/basilisk/src/grid/cartesian-common.h"

  return nodata;
}

static void periodic_boundary (int d)
{

  if (all) for (scalar s = *all, *_i44 = all; ((scalar *)&s)->i >= 0; s = *++_i44)
    _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = periodic_bc;

  if (all) for (scalar s = *all, *_i45 = all; ((scalar *)&s)->i >= 0; s = *++_i45)
    if (_attribute[s.i].face) {
      vector v = _attribute[s.i].v;
      _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] = NULL;
    }

  default_scalar_bc[d] = periodic_bc;
  default_vector_bc[d] = periodic_bc;
}

int Periods[3];

void periodic (int dir)
{





    assert (dir <= back);


  int c = dir/2;
  periodic_boundary (2*c);
  periodic_boundary (2*c + 1);
  Periods[c] = true;
}
#line 4 "/home/miro3/Documents/Programming/basilisk/src/grid/multigrid-common.h"

#ifndef foreach_level_or_leaf
# define foreach_level_or_leaf foreach_level
# define end_foreach_level_or_leaf end_foreach_level
#endif

#ifndef foreach_coarse_level
# define foreach_coarse_level foreach_level
# define end_foreach_coarse_level end_foreach_level
#endif










void (* restriction) (scalar *);

static inline void restriction_average (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 27 "/home/miro3/Documents/Programming/basilisk/src/grid/multigrid-common.h"

  double sum = 0.;
   { foreach_child()
    sum += val(s,0,0,0); end_foreach_child(); }
  val(s,0,0,0) = sum/(1 << 3);
}

static inline void restriction_volume_average (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 35 "/home/miro3/Documents/Programming/basilisk/src/grid/multigrid-common.h"

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 35

  double sum = 0.;
   { foreach_child()
    sum += val_cm(cm,0,0,0)*val(s,0,0,0); end_foreach_child(); }
  val(s,0,0,0) = sum/(1 << 3)/val_cm(cm,0,0,0);
 }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 35

  double sum = 0.;
   { foreach_child()
    sum += val_cm(cm,0,0,0)*val(s,0,0,0); end_foreach_child(); }
  val(s,0,0,0) = sum/(1 << 3)/val_cm(cm,0,0,0);
 }}

static inline void face_average (Point point, vector v)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 43 "/home/miro3/Documents/Programming/basilisk/src/grid/multigrid-common.h"

  {
#line 44
 {







      val(v.x,0,0,0) = (fine(v.x,0,0,0) + fine(v.x,0,1,0) +
        fine(v.x,0,0,1) + fine(v.x,0,1,1))/4.;
      val(v.x,1,0,0) = (fine(v.x,2,0,0) + fine(v.x,2,1,0) +
  fine(v.x,2,0,1) + fine(v.x,2,1,1))/4.;

  }
#line 44
 {







      val(v.y,0,0,0) = (fine(v.y,0,0,0) + fine(v.y,0,0,1) +
        fine(v.y,1,0,0) + fine(v.y,1,0,1))/4.;
      val(v.y,0,1,0) = (fine(v.y,0,2,0) + fine(v.y,0,2,1) +
  fine(v.y,1,2,0) + fine(v.y,1,2,1))/4.;

  }
#line 44
 {







      val(v.z,0,0,0) = (fine(v.z,0,0,0) + fine(v.z,1,0,0) +
        fine(v.z,0,1,0) + fine(v.z,1,1,0))/4.;
      val(v.z,0,0,1) = (fine(v.z,0,0,2) + fine(v.z,1,0,2) +
  fine(v.z,0,1,2) + fine(v.z,1,1,2))/4.;

  }}
}

static inline void restriction_face (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 61 "/home/miro3/Documents/Programming/basilisk/src/grid/multigrid-common.h"

  face_average (point, _attribute[s.i].v);
}

static inline void no_restriction (Point point, scalar s) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 64 "/home/miro3/Documents/Programming/basilisk/src/grid/multigrid-common.h"
}

static inline void no_data (Point point, scalar s) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 67 "/home/miro3/Documents/Programming/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    val(s,0,0,0) = nodata; end_foreach_child(); }
}

void wavelet (scalar s, scalar w)
{
  restriction (((scalar []){s,{-1}}));
  for (int l = depth() - 1; l >= 0; l--)
     { foreach_coarse_level (l){

#line 76 "/home/miro3/Documents/Programming/basilisk/src/grid/multigrid-common.h"
 {
      double sc[1 << 3];
      int c = 0;
       { foreach_child()
 sc[c++] = val(s,0,0,0); end_foreach_child(); }
      _attribute[s.i].prolongation (point, s);
      c = 0;
       { foreach_child() {

 val(w,0,0,0) = sc[c] - val(s,0,0,0);
 val(s,0,0,0) = sc[c++];
      } end_foreach_child(); }
    } } end_foreach_coarse_level(); }

   { foreach_level(0){

#line 90 "/home/miro3/Documents/Programming/basilisk/src/grid/multigrid-common.h"
 val(w,0,0,0) = 0.; } end_foreach_level(); }
}

static inline double bilinear (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 94 "/home/miro3/Documents/Programming/basilisk/src/grid/multigrid-common.h"








    return (27.*coarse(s,0,0,0) +
     9.*(coarse(s,child.x,0,0) + coarse(s,0,child.y,0) +
  coarse(s,0,0,child.z)) +
     3.*(coarse(s,child.x,child.y,0) + coarse(s,child.x,0,child.z) +
  coarse(s,0,child.y,child.z)) +
     coarse(s,child.x,child.y,child.z))/64.;

}

static inline void refine_bilinear (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 112 "/home/miro3/Documents/Programming/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    val(s,0,0,0) = bilinear (point, s); end_foreach_child(); }
}

static inline double quadratic (double a, double b, double c)
{
  return (30.*a + 5.*b - 3.*c)/32.;
}

static inline double biquadratic (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 123 "/home/miro3/Documents/Programming/basilisk/src/grid/multigrid-common.h"

#line 138 "/home/miro3/Documents/Programming/basilisk/src/grid/multigrid-common.h"
  assert (false);
  return 0.;

}

static inline double biquadratic_vertex (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 144 "/home/miro3/Documents/Programming/basilisk/src/grid/multigrid-common.h"







  assert (false);
  return 0.;

}

static inline void refine_biquadratic (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 157 "/home/miro3/Documents/Programming/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    val(s,0,0,0) = biquadratic (point, s); end_foreach_child(); }
}

static inline void refine_linear (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 163 "/home/miro3/Documents/Programming/basilisk/src/grid/multigrid-common.h"

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 163

  coord g;
  if (_attribute[s.i].gradient)
    {
#line 166

      g.x = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0));
#line 166

      g.y = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0));
#line 166

      g.z = _attribute[s.i].gradient (val(s,0,0,-1), val(s,0,0,0), val(s,0,0,1));}
  else
    {
#line 169

      g.x = (val(s,1,0,0) - val(s,-1,0,0))/2.;
#line 169

      g.y = (val(s,0,1,0) - val(s,0,-1,0))/2.;
#line 169

      g.z = (val(s,0,0,1) - val(s,0,0,-1))/2.;}

  double sc = val(s,0,0,0), cmc = 4.*val_cm(cm,0,0,0), sum = val_cm(cm,0,0,0)*(1 << 3);
   { foreach_child() {
    val(s,0,0,0) = sc;
    {
#line 175

      val(s,0,0,0) += child.x*g.x*val_cm(cm,-child.x,0,0)/cmc;
#line 175

      val(s,0,0,0) += child.y*g.y*val_cm(cm,0,-child.y,0)/cmc;
#line 175

      val(s,0,0,0) += child.z*g.z*val_cm(cm,0,0,-child.z)/cmc;}
    sum -= val_cm(cm,0,0,0);
  } end_foreach_child(); }
  assert (fabs(sum) < 1e-10);
 }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 163

  coord g;
  if (_attribute[s.i].gradient)
    {
#line 166

      g.x = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0));
#line 166

      g.y = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0));
#line 166

      g.z = _attribute[s.i].gradient (val(s,0,0,-1), val(s,0,0,0), val(s,0,0,1));}
  else
    {
#line 169

      g.x = (val(s,1,0,0) - val(s,-1,0,0))/2.;
#line 169

      g.y = (val(s,0,1,0) - val(s,0,-1,0))/2.;
#line 169

      g.z = (val(s,0,0,1) - val(s,0,0,-1))/2.;}

  double sc = val(s,0,0,0), cmc = 4.*val_cm(cm,0,0,0), sum = val_cm(cm,0,0,0)*(1 << 3);
   { foreach_child() {
    val(s,0,0,0) = sc;
    {
#line 175

      val(s,0,0,0) += child.x*g.x*val_cm(cm,-child.x,0,0)/cmc;
#line 175

      val(s,0,0,0) += child.y*g.y*val_cm(cm,0,-child.y,0)/cmc;
#line 175

      val(s,0,0,0) += child.z*g.z*val_cm(cm,0,0,-child.z)/cmc;}
    sum -= val_cm(cm,0,0,0);
  } end_foreach_child(); }
  assert (fabs(sum) < 1e-10);
 }}

static inline void refine_reset (Point point, scalar v)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 183 "/home/miro3/Documents/Programming/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    val(v,0,0,0) = 0.; end_foreach_child(); }
}

static inline void refine_injection (Point point, scalar v)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 189 "/home/miro3/Documents/Programming/basilisk/src/grid/multigrid-common.h"

  double val = val(v,0,0,0);
   { foreach_child()
    val(v,0,0,0) = val; end_foreach_child(); }
}

static scalar multigrid_init_scalar (scalar s, const char * name)
{
  s = cartesian_init_scalar (s, name);
  _attribute[s.i].prolongation = refine_bilinear;
  _attribute[s.i].restriction = restriction_average;
  return s;
}

static vector multigrid_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_face_vector (v, name);
  {
#line 206

    _attribute[v.y.i].restriction = no_restriction;
#line 206

    _attribute[v.z.i].restriction = no_restriction;
#line 206

    _attribute[v.x.i].restriction = no_restriction;}
  _attribute[v.x.i].restriction = restriction_face;
  return v;
}

void multigrid_debug (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 213 "/home/miro3/Documents/Programming/basilisk/src/grid/multigrid-common.h"

  cartesian_debug (point);

  FILE * plot = fopen ("plot", "a");
  if (point.level > 0) {
    char name[80] = "coarse";
    if (pid() > 0)
      sprintf (name, "coarse-%d", pid());
    FILE * fp = fopen (name, "w");
#line 246 "/home/miro3/Documents/Programming/basilisk/src/grid/multigrid-common.h"
      double xc = x - child.x*Delta/2., yc = y - child.y*Delta/2.;
      double zc = z - child.z*Delta/2.;
      for (int k = 0; k <= 1; k++)
 for (int l = 0; l <= 1; l++)
   for (int m = 0; m <= 1; m++) {
     if (all) for (scalar v = *all, *_i46 = all; ((scalar *)&v)->i >= 0; v = *++_i46)
       fprintf (fp, "%g %g %g %g ",
         xc + k*child.x*Delta*2. + _attribute[v.i].d.x*Delta,
         yc + l*child.y*Delta*2. + _attribute[v.i].d.y*Delta,
         zc + m*child.z*Delta*2. + _attribute[v.i].d.z*Delta,
         coarse(v,k*child.x,l*child.y,m*child.z));
     fputc ('\n', fp);
   }
      fprintf (qstderr(), ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 3 t ''",
        name);
      fprintf (plot, ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 3 t ''",
        name);

    fclose (fp);
  }

  if (is_coarse()) {
    char name[80] = "fine";
    if (pid() > 0)
      sprintf (name, "fine-%d", pid());
    FILE * fp = fopen (name, "w");
#line 303 "/home/miro3/Documents/Programming/basilisk/src/grid/multigrid-common.h"
      double xf = x - Delta/4., yf = y - Delta/4., zf = z - Delta/4.;
      for (int k = -2; k <= 3; k++)
 for (int l = -2; l <= 3; l++)
   for (int m = -2; m <= 3; m++) {
     if (all) for (scalar v = *all, *_i47 = all; ((scalar *)&v)->i >= 0; v = *++_i47) {
       fprintf (fp, "%g %g %g ",
         xf + k*Delta/2. + _attribute[v.i].d.x*Delta/4.,
         yf + l*Delta/2. + _attribute[v.i].d.y*Delta/4.,
         zf + m*Delta/2. + _attribute[v.i].d.z*Delta/4.);
       if (allocated_child(k,l,m))
  fprintf (fp, "%g ", fine(v,k,l,m));
       else
  fputs ("n/a ", fp);
     }
     fputc ('\n', fp);
   }
      fprintf (qstderr(), ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 2 t ''",
        name);
      fprintf (plot, ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 2 t ''",
        name);

    fclose (fp);
  }
  fflush (qstderr());
  fclose (plot);
}

static void multigrid_restriction (scalar * list)
{
  scalar * listdef = NULL, * listc = NULL, * list2 = NULL;
  if (list) for (scalar s = *list, *_i48 = list; ((scalar *)&s)->i >= 0; s = *++_i48)
    if (!is_constant (s)) {
      if (_attribute[s.i].restriction == restriction_average) {
 listdef = list_add (listdef, s);
 list2 = list_add (list2, s);
      }
      else if (_attribute[s.i].restriction != no_restriction) {
 listc = list_add (listc, s);
 if (_attribute[s.i].face)
   {
#line 342

     list2 = list_add (list2, _attribute[s.i].v.x);
#line 342

     list2 = list_add (list2, _attribute[s.i].v.y);
#line 342

     list2 = list_add (list2, _attribute[s.i].v.z);}
 else
   list2 = list_add (list2, s);
      }
    }

  if (listdef || listc) {
    for (int l = depth() - 1; l >= 0; l--) {
       { foreach_coarse_level(l){

#line 351 "/home/miro3/Documents/Programming/basilisk/src/grid/multigrid-common.h"
 {
 if (listdef) for (scalar s = *listdef, *_i49 = listdef; ((scalar *)&s)->i >= 0; s = *++_i49)
   restriction_average (point, s);
 if (listc) for (scalar s = *listc, *_i50 = listc; ((scalar *)&s)->i >= 0; s = *++_i50)
   _attribute[s.i].restriction (point, s);
      } } end_foreach_coarse_level(); }
      { Boundary ** _i = boundaries, * _b; while ((_b = *_i++)) if (_b->level) _b->level (_b, list2, l); };
    }
    pfree (listdef,__func__,__FILE__,__LINE__);
    pfree (listc,__func__,__FILE__,__LINE__);
    pfree (list2,__func__,__FILE__,__LINE__);
  }
}

void multigrid_methods()
{
  cartesian_methods();
  debug = multigrid_debug;
  init_scalar = multigrid_init_scalar;
  init_face_vector = multigrid_init_face_vector;
  restriction = multigrid_restriction;
}







void subtree_size (scalar size, bool leaves)
{




   { foreach(){

#line 386 "/home/miro3/Documents/Programming/basilisk/src/grid/multigrid-common.h"

    val(size,0,0,0) = 1; } end_foreach(); }





  { Boundary ** _i = boundaries, * _b; while ((_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar []){size,{-1}}), depth()); };
  for (int l = depth() - 1; l >= 0; l--) {
     { foreach_coarse_level(l){

#line 395 "/home/miro3/Documents/Programming/basilisk/src/grid/multigrid-common.h"
 {
      double sum = !leaves;
       { foreach_child()
 sum += val(size,0,0,0); end_foreach_child(); }
      val(size,0,0,0) = sum;
    } } end_foreach_coarse_level(); }
    { Boundary ** _i = boundaries, * _b; while ((_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar []){size,{-1}}), l); };
  }
}
#line 5 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h"












int refine_cell (Point point, scalar * list, int flag, Cache * refined)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 18 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h"

  int nr = 0;


  if (level > 0)
    for (int k = 0; k != 2*child.x; k += child.x)

      for (int l = 0; l != 2*child.y; l += child.y)


 for (int m = 0; m != 2*child.z; m += child.z)

   if (aparent(k,l,m).pid >= 0 && is_leaf(aparent(k,l,m))) {
     Point p = point;


     p.level = point.level - 1;
     p.i = (point.i + 2)/2 + k;

       p.j = (point.j + 2)/2 + l;


       p.k = (point.k + 2)/2 + m;

       nr += refine_cell (p, list, flag, refined);
     aparent(k,l,m).flags |= flag;
   }



  cell.flags &= ~leaf;


  increment_neighbors (point);

  int cflag = is_active(cell) ? (active|leaf) : leaf;
   { foreach_child()
    cell.flags |= cflag; end_foreach_child(); }


  if (list) for (scalar s = *list, *_i51 = list; ((scalar *)&s)->i >= 0; s = *++_i51)
    if (is_local(cell) || _attribute[s.i].face)
      _attribute[s.i].refine (point, s);

#if 1
  if (is_border(cell)) {
     { foreach_child() {
      bool bord = false;
       { foreach_neighbor() {
 if (!is_local(cell) || (level > 0 && !is_local(aparent(0,0,0))))
   bord = true, foreach_neighbor_break();
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
    { foreach_child()
     if (!is_local(cell))
       bord = true, foreach_child_break(); end_foreach_child(); }
 if (bord)
   foreach_neighbor_break();
      } end_foreach_neighbor(); }
      if (bord)
 cell.flags |= border;
    } end_foreach_child(); }
    if (refined)
      cache_append (refined, point, cell.flags);
    nr++;
  }
#endif
  return nr;
}





bool coarsen_cell (Point point, scalar * list)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 92 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h"




  int pid = cell.pid;
   { foreach_child()
    if (cell.neighbors || (cell.pid < 0 && cell.pid != pid))
      return false; end_foreach_child(); }



  if (list) for (scalar s = *list, *_i52 = list; ((scalar *)&s)->i >= 0; s = *++_i52) {
    _attribute[s.i].restriction (point, s);
    if (_attribute[s.i].coarsen)
      _attribute[s.i].coarsen (point, s);
  }


  cell.flags |= leaf;


  decrement_neighbors (point);

#if 1
  if (!is_local(cell)) {
    cell.flags &= ~(active|border);
     { foreach_neighbor(1)
      if (cell.neighbors)
  { foreach_child()
   if (allocated(0,0,0) && is_local(cell) && !is_border(cell))
     cell.flags |= border; end_foreach_child(); } end_foreach_neighbor(); }
  }
#endif

  return true;
}

void coarsen_cell_recursive (Point point, scalar * list)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 130 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h"



   { foreach_child()
    if (cell.neighbors)
       { foreach_neighbor(1)
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
   coarsen_cell_recursive (point, list); end_foreach_neighbor(); } end_foreach_child(); }

  assert (coarsen_cell (point, list));
}

void mpi_boundary_refine (scalar *);
void mpi_boundary_coarsen (int, int);
void mpi_boundary_update (scalar *);

typedef struct {
  int nc, nf;
} astats;

struct Adapt {
  scalar * slist;
  double * max;
  int maxlevel;
  int minlevel;
  scalar * list;
};


astats adapt_wavelet (struct Adapt p)
{ trace ("adapt_wavelet", "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h", 160);
  if (p.list == NULL)
    p.list = all;
  if (is_constant(cm))
    restriction (p.slist);
  else {
    scalar * listr = list_concat (((scalar []){cm,{-1}}), p.slist);
    restriction (listr);
    pfree (listr,__func__,__FILE__,__LINE__);
  }

  astats st = {0, 0};
  scalar * listc = NULL;
  if (p.list) for (scalar s = *p.list, *_i53 = p.list; ((scalar *)&s)->i >= 0; s = *++_i53)
    if (!is_constant(s) && _attribute[s.i].restriction != no_restriction)
      listc = list_add (listc, s);


  if (p.minlevel < 1)
    p.minlevel = 1;
  ((Tree *)grid)->refined.n = 0;
  static const int refined = 1 << user, too_fine = 1 << (user + 1);
   { foreach_cell(){

#line 182 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h"
 {
    if (is_active(cell)) {
      static const int too_coarse = 1 << (user + 2);
      if (is_leaf (cell)) {
 if (cell.flags & too_coarse) {
   cell.flags &= ~too_coarse;
   refine_cell (point, listc, refined, &((Tree *)grid)->refined);
   st.nf++;
 }
 continue;
      }
      else {
 if (cell.flags & refined) {

   cell.flags &= ~too_coarse;
   continue;
 }

 bool local = is_local(cell);
 if (!local)
    { foreach_child()
     if (is_local(cell))
       local = true, foreach_child_break(); end_foreach_child(); }
 if (local) {
   int i = 0;
   static const int just_fine = 1 << (user + 3);
   if (p.slist) for (scalar s = *p.slist, *_i54 = p.slist; ((scalar *)&s)->i >= 0; s = *++_i54) {
     double max = p.max[i++], sc[1 << 3];
     int c = 0;
      { foreach_child()
       sc[c++] = val(s,0,0,0); end_foreach_child(); }
     _attribute[s.i].prolongation (point, s);
     c = 0;
      { foreach_child() {
       double e = fabs(sc[c] - val(s,0,0,0));
       if (e > max && level < p.maxlevel) {
  cell.flags &= ~too_fine;
  cell.flags |= too_coarse;
       }
       else if ((e <= max/1.5 || level > p.maxlevel) &&
         !(cell.flags & (too_coarse|just_fine))) {
  if (level >= p.minlevel)
    cell.flags |= too_fine;
       }
       else if (!(cell.flags & too_coarse)) {
  cell.flags &= ~too_fine;
  cell.flags |= just_fine;
       }
       val(s,0,0,0) = sc[c++];
     } end_foreach_child(); }
   }
    { foreach_child() {
     cell.flags &= ~just_fine;
     if (!is_leaf(cell)) {
       cell.flags &= ~too_coarse;
       if (level >= p.maxlevel)
  cell.flags |= too_fine;
     }
     else if (!is_active(cell))
       cell.flags &= ~too_coarse;
   } end_foreach_child(); }
 }
      }
    }
    else
      continue;
  } } end_foreach_cell(); }
  mpi_boundary_refine (listc);



  for (int l = depth(); l >= 0; l--) {
     { foreach_cell(){

#line 254 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h"

      if (!(cell.pid < 0)) {
 if (level == l) {
   if (!is_leaf(cell)) {
     if (cell.flags & refined)

       cell.flags &= ~(refined|too_fine);
     else if (cell.flags & too_fine) {
       if (is_local(cell) && coarsen_cell (point, listc))
  st.nc++;
       cell.flags &= ~too_fine;
     }
   }
   if (cell.flags & too_fine)
     cell.flags &= ~too_fine;
   else if (level > 0 && (aparent(0,0,0).flags & too_fine))
     aparent(0,0,0).flags &= ~too_fine;
   continue;
 }
 else if (is_leaf(cell))
   continue;
      } } end_foreach_cell(); }
    mpi_boundary_coarsen (l, too_fine);
  }
  pfree (listc,__func__,__FILE__,__LINE__);

  mpi_all_reduce (st.nf, MPI_INT, MPI_SUM);
  mpi_all_reduce (st.nc, MPI_INT, MPI_SUM);
  if (st.nc || st.nf)
    mpi_boundary_update (p.list);

  { astats _ret =  st; end_trace("adapt_wavelet", "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h", 285);  return _ret; }
 end_trace("adapt_wavelet", "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h", 286); }
#line 307 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h"
static void refine_level (int depth)
{
  int refined;
  do {
    refined = 0;
    ((Tree *)grid)->refined.n = 0;
     { foreach_leaf(){

#line 313 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h"

      if (level < depth) {
 refine_cell (point, NULL, 0, &((Tree *)grid)->refined);
 refined++;
 continue;
      } } end_foreach_leaf(); }
    mpi_all_reduce (refined, MPI_INT, MPI_SUM);
    if (refined) {
      mpi_boundary_refine (NULL);
      mpi_boundary_update (NULL);
    }
  } while (refined);
}
#line 352 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h"
static void halo_flux (vector * list)
{
  vector * listv = NULL;
  if (list) for (vector v = *list, *_i55 = list; ((scalar *)&v)->i >= 0; v = *++_i55)
    if (!is_constant(v.x))
      listv = vectors_add (listv, v);

  if (listv) {
    for (int l = depth() - 1; l >= 0; l--)
       { foreach_halo (prolongation, l){

#line 361 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h"

 {
#line 362
 {
#line 378 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h"
   if ((!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0))
     if (listv) for (vector f = *listv, *_i56 = listv; ((scalar *)&f)->i >= 0; f = *++_i56)
       val(f.x,0,0,0) = (fine(f.x,0,0,0) + fine(f.x,0,1,0) +
         fine(f.x,0,0,1) + fine(f.x,0,1,1))/4.;
   if ((!is_leaf (neighbor(1,0,0)) && neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0))
     if (listv) for (vector f = *listv, *_i57 = listv; ((scalar *)&f)->i >= 0; f = *++_i57)
       val(f.x,1,0,0) = (fine(f.x,2,0,0) + fine(f.x,2,1,0) +
   fine(f.x,2,0,1) + fine(f.x,2,1,1))/4.;

      }
#line 362
 {
#line 378 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h"
   if ((!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0))
     if (listv) for (vector f = *listv, *_i56 = listv; ((scalar *)&f)->i >= 0; f = *++_i56)
       val(f.y,0,0,0) = (fine(f.y,0,0,0) + fine(f.y,0,0,1) +
         fine(f.y,1,0,0) + fine(f.y,1,0,1))/4.;
   if ((!is_leaf (neighbor(0,1,0)) && neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0))
     if (listv) for (vector f = *listv, *_i57 = listv; ((scalar *)&f)->i >= 0; f = *++_i57)
       val(f.y,0,1,0) = (fine(f.y,0,2,0) + fine(f.y,0,2,1) +
   fine(f.y,1,2,0) + fine(f.y,1,2,1))/4.;

      }
#line 362
 {
#line 378 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h"
   if ((!is_leaf (neighbor(0,0,-1)) && neighbor(0,0,-1).neighbors && neighbor(0,0,-1).pid >= 0))
     if (listv) for (vector f = *listv, *_i56 = listv; ((scalar *)&f)->i >= 0; f = *++_i56)
       val(f.z,0,0,0) = (fine(f.z,0,0,0) + fine(f.z,1,0,0) +
         fine(f.z,0,1,0) + fine(f.z,1,1,0))/4.;
   if ((!is_leaf (neighbor(0,0,1)) && neighbor(0,0,1).neighbors && neighbor(0,0,1).pid >= 0))
     if (listv) for (vector f = *listv, *_i57 = listv; ((scalar *)&f)->i >= 0; f = *++_i57)
       val(f.z,0,0,1) = (fine(f.z,0,0,2) + fine(f.z,1,0,2) +
   fine(f.z,0,1,2) + fine(f.z,1,1,2))/4.;

      }} } end_foreach_halo(); }
    pfree (listv,__func__,__FILE__,__LINE__);
  }
}



static scalar tree_init_scalar (scalar s, const char * name)
{
  s = multigrid_init_scalar (s, name);
  _attribute[s.i].refine = _attribute[s.i].prolongation;
  return s;
}


#line 401

static void refine_face_x (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 403 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h"

  assert (3 > 1);
  vector v = _attribute[s.i].v;
  if (!(!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) &&
      (is_local(cell) || is_local(neighbor(-1,0,0)))) {
    double g1 = (val(v.x,0,+1,0) - val(v.x,0,-1,0))/8.;
    double g2 = (val(v.x,0,0,+1) - val(v.x,0,0,-1))/8.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.x,0,j,k) = val(v.x,0,0,0) + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
  if (!(!is_leaf (neighbor(1,0,0)) && neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0) && neighbor(1,0,0).neighbors &&
      (is_local(cell) || is_local(neighbor(1,0,0)))) {
    double g1 = (val(v.x,1,+1,0) - val(v.x,1,-1,0))/8.;
    double g2 = (val(v.x,1,0,+1) - val(v.x,1,0,-1))/8.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.x,2,j,k) = val(v.x,1,0,0) + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
  if (is_local(cell)) {
    double g1 = (val(v.x,0,+1,0) + val(v.x,1,+1,0) - val(v.x,0,-1,0) - val(v.x,1,-1,0))/16.;
    double g2 = (val(v.x,0,0,+1) + val(v.x,1,0,+1) - val(v.x,0,0,-1) - val(v.x,1,0,-1))/16.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.x,1,j,k) = (val(v.x,0,0,0) + val(v.x,1,0,0))/2. + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
}
#line 401

static void refine_face_y (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 403 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h"

  assert (3 > 1);
  vector v = _attribute[s.i].v;
  if (!(!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) &&
      (is_local(cell) || is_local(neighbor(0,-1,0)))) {
    double g1 = (val(v.y,0,0,+1) - val(v.y,0,0,-1))/8.;
    double g2 = (val(v.y,+1,0,0) - val(v.y,-1,0,0))/8.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.y,k,0,j) = val(v.y,0,0,0) + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
  if (!(!is_leaf (neighbor(0,1,0)) && neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0) && neighbor(0,1,0).neighbors &&
      (is_local(cell) || is_local(neighbor(0,1,0)))) {
    double g1 = (val(v.y,0,1,+1) - val(v.y,0,1,-1))/8.;
    double g2 = (val(v.y,+1,1,0) - val(v.y,-1,1,0))/8.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.y,k,2,j) = val(v.y,0,1,0) + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
  if (is_local(cell)) {
    double g1 = (val(v.y,0,0,+1) + val(v.y,0,1,+1) - val(v.y,0,0,-1) - val(v.y,0,1,-1))/16.;
    double g2 = (val(v.y,+1,0,0) + val(v.y,+1,1,0) - val(v.y,-1,0,0) - val(v.y,-1,1,0))/16.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.y,k,1,j) = (val(v.y,0,0,0) + val(v.y,0,1,0))/2. + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
}
#line 401

static void refine_face_z (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 403 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h"

  assert (3 > 1);
  vector v = _attribute[s.i].v;
  if (!(!is_leaf (neighbor(0,0,-1)) && neighbor(0,0,-1).neighbors && neighbor(0,0,-1).pid >= 0) &&
      (is_local(cell) || is_local(neighbor(0,0,-1)))) {
    double g1 = (val(v.z,+1,0,0) - val(v.z,-1,0,0))/8.;
    double g2 = (val(v.z,0,+1,0) - val(v.z,0,-1,0))/8.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.z,j,k,0) = val(v.z,0,0,0) + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
  if (!(!is_leaf (neighbor(0,0,1)) && neighbor(0,0,1).neighbors && neighbor(0,0,1).pid >= 0) && neighbor(0,0,1).neighbors &&
      (is_local(cell) || is_local(neighbor(0,0,1)))) {
    double g1 = (val(v.z,+1,0,1) - val(v.z,-1,0,1))/8.;
    double g2 = (val(v.z,0,+1,1) - val(v.z,0,-1,1))/8.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.z,j,k,2) = val(v.z,0,0,1) + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
  if (is_local(cell)) {
    double g1 = (val(v.z,+1,0,0) + val(v.z,+1,0,1) - val(v.z,-1,0,0) - val(v.z,-1,0,1))/16.;
    double g2 = (val(v.z,0,+1,0) + val(v.z,0,+1,1) - val(v.z,0,-1,0) - val(v.z,0,-1,1))/16.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.z,j,k,1) = (val(v.z,0,0,0) + val(v.z,0,0,1))/2. + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
}

void refine_face (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 432 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h"

  vector v = _attribute[s.i].v;
  {
#line 434

    _attribute[v.x.i].prolongation (point, v.x);
#line 434

    _attribute[v.y.i].prolongation (point, v.y);
#line 434

    _attribute[v.z.i].prolongation (point, v.z);}
}

void refine_face_solenoidal (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 439 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h"

  refine_face (point, s);

  if (is_local(cell)) {

    vector v = _attribute[s.i].v;
    double d[1 << 3], p[1 << 3];
    int i = 0;
     { foreach_child() {
      d[i] = 0.;
      {
#line 449

 d[i] += val(v.x,1,0,0) - val(v.x,0,0,0);
#line 449

 d[i] += val(v.y,0,1,0) - val(v.y,0,0,0);
#line 449

 d[i] += val(v.z,0,0,1) - val(v.z,0,0,0);}
      i++;
    } end_foreach_child(); }
#line 463 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h"
    static double m[7][7] = {
      {7./12,5./24,3./8,5./24,3./8,1./4,1./3},
      {5./24,7./12,3./8,5./24,1./4,3./8,1./3},
      {3./8,3./8,3./4,1./4,3./8,3./8,1./2},
      {5./24,5./24,1./4,7./12,3./8,3./8,1./3},
      {3./8,1./4,3./8,3./8,3./4,3./8,1./2},
      {1./4,3./8,3./8,3./8,3./8,3./4,1./2},
      {1./3,1./3,1./2,1./3,1./2,1./2,5./6}
    };
    p[0] = 0.;
    for (int i = 0; i < 7; i++) {
      p[i + 1] = 0.;
      for (int j = 0; j < 7; j++)
 p[i + 1] += m[i][j]*d[j+1];
    }
    for (int k = 0; k <= 1; k++) {
      fine(v.x,1,0,k) += p[4+k] - p[0+k];
      fine(v.x,1,1,k) += p[6+k] - p[2+k];
      fine(v.y,0,1,k) += p[2+k] - p[0+k];
      fine(v.y,1,1,k) += p[6+k] - p[4+k];
    }
    fine(v.z,0,0,1) += p[1] - p[0];
    fine(v.z,0,1,1) += p[3] - p[2];
    fine(v.z,1,0,1) += p[5] - p[4];
    fine(v.z,1,1,1) += p[7] - p[6];

  }

}

vector tree_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_face_vector (v, name);
  {
#line 496

    _attribute[v.x.i].restriction = _attribute[v.x.i].refine = no_restriction;
#line 496

    _attribute[v.y.i].restriction = _attribute[v.y.i].refine = no_restriction;
#line 496

    _attribute[v.z.i].restriction = _attribute[v.z.i].refine = no_restriction;}
  _attribute[v.x.i].restriction = restriction_face;
  _attribute[v.x.i].refine = refine_face;
  {
#line 500

    _attribute[v.x.i].prolongation = refine_face_x;
#line 500

    _attribute[v.y.i].prolongation = refine_face_y;
#line 500

    _attribute[v.z.i].prolongation = refine_face_z;}
  return v;
}

static void tree_boundary_level (scalar * list, int l)
{
  int depth = l < 0 ? depth() : l;

  if (tree_is_full()) {
    { Boundary ** _i = boundaries, * _b; while ((_b = *_i++)) if (_b->level) _b->level (_b, list, depth); };
    return;
  }

  scalar * listdef = NULL, * listc = NULL, * list2 = NULL;
  if (list) for (scalar s = *list, *_i58 = list; ((scalar *)&s)->i >= 0; s = *++_i58)
    if (!is_constant (s)) {
      if (_attribute[s.i].restriction == restriction_average) {
 listdef = list_add (listdef, s);
 list2 = list_add (list2, s);
      }
      else if (_attribute[s.i].restriction != no_restriction) {
 listc = list_add (listc, s);
 if (_attribute[s.i].face)
   {
#line 524

     list2 = list_add (list2, _attribute[s.i].v.x);
#line 524

     list2 = list_add (list2, _attribute[s.i].v.y);
#line 524

     list2 = list_add (list2, _attribute[s.i].v.z);}
 else
   list2 = list_add (list2, s);
      }
    }

  if (listdef || listc) {
    { Boundary ** _i = boundaries, * _b; while ((_b = *_i++)) if (_b->restriction) _b->restriction (_b, list2, depth); };
    for (int l = depth - 1; l >= 0; l--) {
       { foreach_coarse_level(l){

#line 534 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h"
 {
 if (listdef) for (scalar s = *listdef, *_i59 = listdef; ((scalar *)&s)->i >= 0; s = *++_i59)
   restriction_average (point, s);
 if (listc) for (scalar s = *listc, *_i60 = listc; ((scalar *)&s)->i >= 0; s = *++_i60)
   _attribute[s.i].restriction (point, s);
      } } end_foreach_coarse_level(); }
      { Boundary ** _i = boundaries, * _b; while ((_b = *_i++)) if (_b->restriction) _b->restriction (_b, list2, l); };
    }
    pfree (listdef,__func__,__FILE__,__LINE__);
    pfree (listc,__func__,__FILE__,__LINE__);
    pfree (list2,__func__,__FILE__,__LINE__);
  }

  scalar * listr = NULL;
  vector * listf = NULL;
  if (list) for (scalar s = *list, *_i61 = list; ((scalar *)&s)->i >= 0; s = *++_i61)
    if (!is_constant (s) && _attribute[s.i].refine != no_restriction) {
      if (_attribute[s.i].face)
 listf = vectors_add (listf, _attribute[s.i].v);
      else
 listr = list_add (listr, s);
    }

  if (listr || listf) {
    { Boundary ** _i = boundaries, * _b; while ((_b = *_i++)) if (_b->level) _b->level (_b, list, 0); };
    for (int i = 0; i < depth; i++) {
       { foreach_halo (prolongation, i){

#line 560 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h"
 {
 if (listr) for (scalar s = *listr, *_i62 = listr; ((scalar *)&s)->i >= 0; s = *++_i62)
          _attribute[s.i].prolongation (point, s);
 if (listf) for (vector v = *listf, *_i63 = listf; ((scalar *)&v)->i >= 0; v = *++_i63)
   {
#line 564

     _attribute[v.x.i].prolongation (point, v.x);
#line 564

     _attribute[v.y.i].prolongation (point, v.y);
#line 564

     _attribute[v.z.i].prolongation (point, v.z);}
      } } end_foreach_halo(); }
      { Boundary ** _i = boundaries, * _b; while ((_b = *_i++)) if (_b->level) _b->level (_b, list, i + 1); };
    }
    pfree (listr,__func__,__FILE__,__LINE__);
    pfree (listf,__func__,__FILE__,__LINE__);
  }
}

double treex (Point point) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 574 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h"

  if (level == 0)
    return 0;




  assert (false);
  double i = 0;

  return treex(parent) + i/(1 << 2*(level - 1));
}

double treey (Point point) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 587 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h"

  if (level == 0)
    return 0;
  return treey(parent) + 4./(1 << 2*(level - 1));
}

void output_tree (FILE * fp)
{
   { foreach_cell(){

#line 595 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h"

    if (cell.neighbors)
       { foreach_child()
 if (is_local(cell))
   fprintf (fp, "%g %g\n%g %g\n\n",
     treex(parent), treey(parent), treex(point), treey(point)); end_foreach_child(); }; } end_foreach_cell(); }
}

void tree_check()
{


  long nleaves = 0, nactive = 0;
   { foreach_cell_all(){

#line 608 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h"
 {
    if (is_leaf(cell)) {
      assert (cell.pid >= 0);
      nleaves++;
    }
    if (is_local(cell))
      assert (is_active(cell) || (!is_leaf(cell) && !cell.neighbors && cell.pid >= 0));
    if (is_active(cell))
      nactive++;

    int neighbors = 0;
     { foreach_neighbor(1)
      if (allocated(0,0,0) && (!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
 neighbors++; end_foreach_neighbor(); }
    assert (cell.neighbors == neighbors);


    if (!cell.neighbors)
      assert (!allocated_child(0,0,0));
  } } end_foreach_cell_all(); }


  long reachable = 0;
   { foreach_cell(){

#line 631 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h"
 {
    if (is_active(cell))
      reachable++;
    else
      continue;
  } } end_foreach_cell(); }
  assert (nactive == reachable);


  reachable = 0;
   { foreach_cell(){

#line 641 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-common.h"

    if (is_leaf(cell)) {
      reachable++;
      continue;
    } } end_foreach_cell(); }
  assert (nleaves == reachable);
}

static void tree_restriction (scalar * list) {
  if (tree_is_full())
    multigrid_restriction (list);

}

void tree_methods()
{
  multigrid_methods();
  init_scalar = tree_init_scalar;
  init_face_vector = tree_init_face_vector;
  boundary_level = tree_boundary_level;
  boundary_flux = halo_flux;
  restriction = tree_restriction;
}
#line 1741 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"

#if 1
#line 1 "grid/tree-mpi.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"

int debug_iteration = -1;

void debug_mpi (FILE * fp1);

typedef struct {
  CacheLevel * halo;
  void * buf;
  MPI_Request r;
  int depth;
  int pid;
  int maxdepth;
} Rcv;

typedef struct {
  Rcv * rcv;
  char * name;
  int npid;
} RcvPid;

typedef struct {
  RcvPid * rcv, * snd;
} SndRcv;

typedef struct {
  Boundary parent;

  SndRcv mpi_level, mpi_level_root, restriction;
  Array * send, * receive;
} MpiBoundary;

static void cache_level_init (CacheLevel * c)
{
  c->p = NULL;
  c->n = c->nm = 0;
}

static void rcv_append (Point point, Rcv * rcv)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 39 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"

  if (level > rcv->depth) {
    rcv->halo = (CacheLevel *) prealloc (rcv->halo, (level + 1)*sizeof(CacheLevel),__func__,__FILE__,__LINE__);
    for (int j = rcv->depth + 1; j <= level; j++)
      cache_level_init (&rcv->halo[j]);
    rcv->depth = level;
  }
  cache_level_append (&rcv->halo[level], point);
  if (level > rcv->maxdepth)
    rcv->maxdepth = level;
}

void rcv_print (Rcv * rcv, FILE * fp, const char * prefix)
{
  for (int l = 0; l <= rcv->depth; l++)
    if (rcv->halo[l].n > 0)
       { foreach_cache_level(rcv->halo[l], l){

#line 55 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"

 fprintf (fp, "%s%g %g %g %d %d\n", prefix, x, y, z, rcv->pid, level); } end_foreach_cache_level(); }
}

static void rcv_free_buf (Rcv * rcv)
{
  if (rcv->buf) {
    prof_start ("rcv_pid_receive");
    MPI_Wait (&rcv->r, MPI_STATUS_IGNORE);
    pfree (rcv->buf,__func__,__FILE__,__LINE__);
    rcv->buf = NULL;
    prof_stop();
  }
}

static void rcv_destroy (Rcv * rcv)
{
  rcv_free_buf (rcv);
  for (int i = 0; i <= rcv->depth; i++)
    if (rcv->halo[i].n > 0)
      pfree (rcv->halo[i].p,__func__,__FILE__,__LINE__);
  pfree (rcv->halo,__func__,__FILE__,__LINE__);
}

static RcvPid * rcv_pid_new (const char * name)
{
  RcvPid * r = ((RcvPid *) pcalloc (1, sizeof(RcvPid),__func__,__FILE__,__LINE__));
  r->name = pstrdup (name,__func__,__FILE__,__LINE__);
  return r;
}

static Rcv * rcv_pid_pointer (RcvPid * p, int pid)
{
  assert (pid >= 0 && pid < npe());

  int i;
  for (i = 0; i < p->npid; i++)
    if (pid == p->rcv[i].pid)
      break;

  if (i == p->npid) {
    p->rcv = (Rcv *) prealloc (p->rcv, (++p->npid)*sizeof(Rcv),__func__,__FILE__,__LINE__);
    Rcv * rcv = &p->rcv[p->npid-1];
    rcv->pid = pid;
    rcv->depth = rcv->maxdepth = 0;
    rcv->halo = ((CacheLevel *) pmalloc ((1)*sizeof(CacheLevel),__func__,__FILE__,__LINE__));
    rcv->buf = NULL;
    cache_level_init (&rcv->halo[0]);
  }
  return &p->rcv[i];
}

static void rcv_pid_append (RcvPid * p, int pid, Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 108 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"

  rcv_append (point, rcv_pid_pointer (p, pid));
}

static void rcv_pid_append_pids (RcvPid * p, Array * pids)
{

  for (int i = 0; i < p->npid; i++) {
    int pid = p->rcv[i].pid, j, * a;
    for (j = 0, a = pids->p; j < pids->len/sizeof(int); j++,a++)
      if (*a == pid)
 break;
    if (j == pids->len/sizeof(int))
      array_append (pids, &pid, sizeof(int));
  }
}

void rcv_pid_write (RcvPid * p, const char * name)
{
  for (int i = 0; i < p->npid; i++) {
    Rcv * rcv = &p->rcv[i];
    char fname[80];
    sprintf (fname, "%s-%d-%d", name, pid(), rcv->pid);
    FILE * fp = fopen (fname, "w");
    rcv_print (rcv, fp, "");
    fclose (fp);
  }
}

static void rcv_pid_print (RcvPid * p, FILE * fp, const char * prefix)
{
  for (int i = 0; i < p->npid; i++)
    rcv_print (&p->rcv[i], fp, prefix);
}

static void rcv_pid_destroy (RcvPid * p)
{
  for (int i = 0; i < p->npid; i++)
    rcv_destroy (&p->rcv[i]);
  pfree (p->rcv,__func__,__FILE__,__LINE__);
  pfree (p->name,__func__,__FILE__,__LINE__);
  pfree (p,__func__,__FILE__,__LINE__);
}

static Boundary * mpi_boundary = NULL;






void debug_mpi (FILE * fp1);

static void apply_bc (Rcv * rcv, scalar * list, vector * listf, int l,
        MPI_Status s)
{
  double * b = rcv->buf;
   { foreach_cache_level(rcv->halo[l], l){

#line 165 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"
 {
    if (list) for (scalar s = *list, *_i64 = list; ((scalar *)&s)->i >= 0; s = *++_i64)
      val(s,0,0,0) = *b++;
    if (listf) for (vector v = *listf, *_i65 = listf; ((scalar *)&v)->i >= 0; v = *++_i65) {
      {
#line 169
 {
 val(v.x,0,0,0) = *b++;
 if (allocated(1,0,0))
   val(v.x,1,0,0) = *b;
 b++;
      }
#line 169
 {
 val(v.y,0,0,0) = *b++;
 if (allocated(0,1,0))
   val(v.y,0,1,0) = *b;
 b++;
      }
#line 169
 {
 val(v.z,0,0,0) = *b++;
 if (allocated(0,0,1))
   val(v.z,0,0,1) = *b;
 b++;
      }}
    }
  } } end_foreach_cache_level(); }
  size_t size = b - (double *) rcv->buf;
  pfree (rcv->buf,__func__,__FILE__,__LINE__);
  rcv->buf = NULL;

  int rlen;
  MPI_Get_count (&s, MPI_DOUBLE, &rlen);
  if (rlen != size) {
    fprintf (qstderr(),
      "rlen (%d) != size (%ld), %d receiving from %d at level %d\n"
      "Calling debug_mpi(NULL)...\n"
      "Aborting...\n",
      rlen, size, pid(), rcv->pid, l);
    fflush (qstderr());
    debug_mpi (NULL);
    MPI_Abort (MPI_COMM_WORLD, -2);
  }
}
#line 215 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"
static void mpi_recv_check (void * buf, int count, MPI_Datatype datatype,
       int source, int tag,
       MPI_Comm comm, MPI_Status * status,
       const char * name)
{
#line 250 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"
  int errorcode = MPI_Recv (buf, count, datatype, source, tag, comm, status);
  if (errorcode != MPI_SUCCESS) {
    char string[MPI_MAX_ERROR_STRING];
    int resultlen;
    MPI_Error_string (errorcode, string, &resultlen);
    fprintf (qstderr(),
      "ERROR MPI_Recv \"%s\" (count = %d, source = %d, tag = %d):\n%s\n"
      "Calling debug_mpi(NULL)...\n"
      "Aborting...\n",
      name, count, source, tag, string);
    fflush (qstderr());
    debug_mpi (NULL);
    MPI_Abort (MPI_COMM_WORLD, -1);
  }





}

static void rcv_pid_receive (RcvPid * m, scalar * list, vector * listf, int l)
{
  if (m->npid == 0)
    return;

  prof_start ("rcv_pid_receive");

  int len = list_len (list) + 2*3*vectors_len (listf);

  MPI_Request r[m->npid];
  Rcv * rrcv[m->npid];
  int nr = 0;
  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0) {
      assert (!rcv->buf);
      rcv->buf = pmalloc (sizeof (double)*rcv->halo[l].n*len,__func__,__FILE__,__LINE__);






      MPI_Irecv (rcv->buf, rcv->halo[l].n*len, MPI_DOUBLE, rcv->pid,
   (l), MPI_COMM_WORLD, &r[nr]);
      rrcv[nr++] = rcv;






    }
  }


  if (nr > 0) {
    int i;
    MPI_Status s;
    MPI_Waitany (nr, r, &i, &s);
    while (i != MPI_UNDEFINED) {
      Rcv * rcv = rrcv[i];
      assert (l <= rcv->depth && rcv->halo[l].n > 0);
      assert (rcv->buf);
      apply_bc (rcv, list, listf, l, s);
      MPI_Waitany (nr, r, &i, &s);
    }
  }

  prof_stop();
}

static void rcv_pid_wait (RcvPid * m)
{

  for (int i = 0; i < m->npid; i++)
    rcv_free_buf (&m->rcv[i]);
}

static void rcv_pid_send (RcvPid * m, scalar * list, vector * listf, int l)
{
  if (m->npid == 0)
    return;

  prof_start ("rcv_pid_send");

  int len = list_len (list) + 2*3*vectors_len (listf);


  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0) {
      assert (!rcv->buf);
      rcv->buf = pmalloc (sizeof (double)*rcv->halo[l].n*len,__func__,__FILE__,__LINE__);
      double * b = rcv->buf;
       { foreach_cache_level(rcv->halo[l], l){

#line 346 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"
 {
 if (list) for (scalar s = *list, *_i66 = list; ((scalar *)&s)->i >= 0; s = *++_i66)
   *b++ = val(s,0,0,0);
 if (listf) for (vector v = *listf, *_i67 = listf; ((scalar *)&v)->i >= 0; v = *++_i67)
   {
#line 350
 {
     *b++ = val(v.x,0,0,0);
     *b++ = allocated(1,0,0) ? val(v.x,1,0,0) : undefined;
   }
#line 350
 {
     *b++ = val(v.y,0,0,0);
     *b++ = allocated(0,1,0) ? val(v.y,0,1,0) : undefined;
   }
#line 350
 {
     *b++ = val(v.z,0,0,0);
     *b++ = allocated(0,0,1) ? val(v.z,0,0,1) : undefined;
   }}
      } } end_foreach_cache_level(); }





      MPI_Isend (rcv->buf, (b - (double *) rcv->buf),
   MPI_DOUBLE, rcv->pid, (l), MPI_COMM_WORLD,
   &rcv->r);
    }
  }

  prof_stop();
}

static void rcv_pid_sync (SndRcv * m, scalar * list, int l)
{
  scalar * listr = NULL;
  vector * listf = NULL;
  if (list) for (scalar s = *list, *_i68 = list; ((scalar *)&s)->i >= 0; s = *++_i68)
    if (!is_constant(s)) {
      if (_attribute[s.i].face)
 listf = vectors_add (listf, _attribute[s.i].v);
      else
 listr = list_add (listr, s);
    }
  rcv_pid_send (m->snd, listr, listf, l);
  rcv_pid_receive (m->rcv, listr, listf, l);
  rcv_pid_wait (m->snd);
  pfree (listr,__func__,__FILE__,__LINE__);
  pfree (listf,__func__,__FILE__,__LINE__);
}

static void snd_rcv_destroy (SndRcv * m)
{
  rcv_pid_destroy (m->rcv);
  rcv_pid_destroy (m->snd);
}

static void snd_rcv_init (SndRcv * m, const char * name)
{
  char s[strlen(name) + 5];
  strcpy (s, name);
  strcat (s, ".rcv");
  m->rcv = rcv_pid_new (s);
  strcpy (s, name);
  strcat (s, ".snd");
  m->snd = rcv_pid_new (s);
}

static void mpi_boundary_destroy (Boundary * b)
{
  MpiBoundary * m = (MpiBoundary *) b;
  snd_rcv_destroy (&m->mpi_level);
  snd_rcv_destroy (&m->mpi_level_root);
  snd_rcv_destroy (&m->restriction);
  array_free (m->send);
  array_free (m->receive);
  pfree (m,__func__,__FILE__,__LINE__);
}


static void mpi_boundary_level (const Boundary * b, scalar * list, int l)
{ trace ("mpi_boundary_level", "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h", 417);
  MpiBoundary * m = (MpiBoundary *) b;
  rcv_pid_sync (&m->mpi_level, list, l);
  rcv_pid_sync (&m->mpi_level_root, list, l);
 end_trace("mpi_boundary_level", "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h", 421); }


static void mpi_boundary_restriction (const Boundary * b, scalar * list, int l)
{ trace ("mpi_boundary_restriction", "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h", 425);
  MpiBoundary * m = (MpiBoundary *) b;
  rcv_pid_sync (&m->restriction, list, l);
 end_trace("mpi_boundary_restriction", "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h", 428); }

void mpi_boundary_new()
{
  mpi_boundary = (Boundary *) ((MpiBoundary *) pcalloc (1, sizeof(MpiBoundary),__func__,__FILE__,__LINE__));
  mpi_boundary->destroy = mpi_boundary_destroy;
  mpi_boundary->level = mpi_boundary_level;
  mpi_boundary->restriction = mpi_boundary_restriction;
  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;
  snd_rcv_init (&mpi->mpi_level, "mpi_level");
  snd_rcv_init (&mpi->mpi_level_root, "mpi_level_root");
  snd_rcv_init (&mpi->restriction, "restriction");
  mpi->send = array_new();
  mpi->receive = array_new();
  add_boundary (mpi_boundary);
}

static FILE * fopen_prefix (FILE * fp, const char * name, char * prefix)
{
  if (fp) {
    sprintf (prefix, "%s-%d ", name, pid());
    return fp;
  }
  else {
    strcpy (prefix, "");
    char fname[80];
    if (debug_iteration >= 0)
      sprintf (fname, "%s-%d-%d", name, debug_iteration, pid());
    else
      sprintf (fname, "%s-%d", name, pid());
    return fopen (fname, "w");
  }
}

void debug_mpi (FILE * fp1)
{
  void output_cells (FILE * fp);

  char prefix[80];
  FILE * fp;


  if (fp1 == NULL) {
    char name[80];
    sprintf (name, "halo-%d", pid()); remove (name);
    sprintf (name, "cells-%d", pid()); remove (name);
    sprintf (name, "faces-%d", pid()); remove (name);
    sprintf (name, "neighbors-%d", pid()); remove (name);
    sprintf (name, "mpi-level-rcv-%d", pid()); remove (name);
    sprintf (name, "mpi-level-snd-%d", pid()); remove (name);
    sprintf (name, "mpi-level-root-rcv-%d", pid()); remove (name);
    sprintf (name, "mpi-level-root-snd-%d", pid()); remove (name);
    sprintf (name, "mpi-restriction-rcv-%d", pid()); remove (name);
    sprintf (name, "mpi-restriction-snd-%d", pid()); remove (name);
    sprintf (name, "mpi-border-%d", pid()); remove (name);
    sprintf (name, "exterior-%d", pid()); remove (name);
    sprintf (name, "depth-%d", pid()); remove (name);
    sprintf (name, "refined-%d", pid()); remove (name);
  }


  fp = fopen_prefix (fp1, "halo", prefix);
  for (int l = 0; l < depth(); l++)
     { foreach_halo (prolongation, l){

#line 491 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"

       { foreach_child()
      fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level); end_foreach_child(); }; } end_foreach_halo(); }
  if (!fp1)
    fclose (fp);

  if (!fp1) {
    fp = fopen_prefix (fp1, "cells", prefix);
    output_cells (fp);
    fclose (fp);
  }

  fp = fopen_prefix (fp1, "faces", prefix);
   { foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 504
{

#line 504 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"

    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 504
{

#line 504 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"

    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 504
{

#line 504 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"

    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level); }  }}  end_foreach_face_generic()
#line 505
 end_foreach_face(); }
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "neighbors", prefix);
   { foreach(){

#line 510 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"
 {
    int n = 0;
     { foreach_neighbor(1)
      if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
 n++; end_foreach_neighbor(); }
    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, cell.neighbors);
    assert (cell.neighbors == n);
  } } end_foreach(); }
  if (!fp1)
    fclose (fp);

  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;

  fp = fopen_prefix (fp1, "mpi-level-rcv", prefix);
  rcv_pid_print (mpi->mpi_level.rcv, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-level-root-rcv", prefix);
  rcv_pid_print (mpi->mpi_level_root.rcv, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-restriction-rcv", prefix);
  rcv_pid_print (mpi->restriction.rcv, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-level-snd", prefix);
  rcv_pid_print (mpi->mpi_level.snd, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-level-root-snd", prefix);
  rcv_pid_print (mpi->mpi_level_root.snd, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-restriction-snd", prefix);
  rcv_pid_print (mpi->restriction.snd, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-border", prefix);
   { foreach_cell(){

#line 554 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"
 {
    if (is_border(cell))
      fprintf (fp, "%s%g %g %g %d %d %d\n",
        prefix, x, y, z, level, cell.neighbors, cell.pid);
    else
      continue;
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "exterior", prefix);
   { foreach_cell(){

#line 567 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"
 {
    if (!is_local(cell))
      fprintf (fp, "%s%g %g %g %d %d %d %d\n",
        prefix, x, y, z, level, cell.neighbors,
        cell.pid, cell.flags & leaf);






  } } end_foreach_cell(); }
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "depth", prefix);
  fprintf (fp, "depth: %d %d\n", pid(), depth());
  fprintf (fp, "======= mpi_level.snd ======\n");
  RcvPid * snd = mpi->mpi_level.snd;
  for (int i = 0; i < snd->npid; i++)
    fprintf (fp, "%d %d %d\n", pid(), snd->rcv[i].pid, snd->rcv[i].maxdepth);
  fprintf (fp, "======= mpi_level.rcv ======\n");
  snd = mpi->mpi_level.rcv;
  for (int i = 0; i < snd->npid; i++)
    fprintf (fp, "%d %d %d\n", pid(), snd->rcv[i].pid, snd->rcv[i].maxdepth);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "refined", prefix);
   { foreach_cache (((Tree *)grid)->refined){

#line 596 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"

    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level); } end_foreach_cache(); }
  if (!fp1)
    fclose (fp);
}

static void snd_rcv_free (SndRcv * p)
{
  char name[strlen(p->rcv->name) + 1];
  strcpy (name, p->rcv->name);
  rcv_pid_destroy (p->rcv);
  p->rcv = rcv_pid_new (name);
  strcpy (name, p->snd->name);
  rcv_pid_destroy (p->snd);
  p->snd = rcv_pid_new (name);
}

static bool is_root (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 614 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"

  if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
     { foreach_child()
      if (is_local(cell))
 return true; end_foreach_child(); }
  return false;
}


static bool is_local_prolongation (Point point, Point p)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 624 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"




  struct { int x, y, z; } dp = {p.i - point.i, p.j - point.j, p.k - point.k};

  {
#line 630
 {
    if (dp.x == 0 && ((!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) || (!is_leaf (neighbor(1,0,0)) && neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0)))
      return true;
    if ((!is_leaf (neighbor(dp.x,0,0)) && neighbor(dp.x,0,0).neighbors && neighbor(dp.x,0,0).pid >= 0))
      return true;
  }
#line 630
 {
    if (dp.y == 0 && ((!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) || (!is_leaf (neighbor(0,1,0)) && neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0)))
      return true;
    if ((!is_leaf (neighbor(0,dp.y,0)) && neighbor(0,dp.y,0).neighbors && neighbor(0,dp.y,0).pid >= 0))
      return true;
  }
#line 630
 {
    if (dp.z == 0 && ((!is_leaf (neighbor(0,0,-1)) && neighbor(0,0,-1).neighbors && neighbor(0,0,-1).pid >= 0) || (!is_leaf (neighbor(0,0,1)) && neighbor(0,0,1).neighbors && neighbor(0,0,1).pid >= 0)))
      return true;
    if ((!is_leaf (neighbor(0,0,dp.z)) && neighbor(0,0,dp.z).neighbors && neighbor(0,0,dp.z).pid >= 0))
      return true;
  }}
  return false;
}



static void append_pid (Array * pids, int pid)
{
  for (int i = 0, * p = (int *) pids->p; i < pids->len/sizeof(int); i++, p++)
    if (*p == pid)
      return;
  array_append (pids, &pid, sizeof(int));
}

static int locals_pids (Point point, Array * pids)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 650 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"

  if (is_leaf(cell)) {
    if (is_local(cell)) {
      Point p = point;
       { foreach_neighbor(1) {
 if ((cell.pid >= 0 && cell.pid != pid()) &&
     ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) || is_local_prolongation (point, p)))
   append_pid (pids, cell.pid);
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
    { foreach_child()
     if ((cell.pid >= 0 && cell.pid != pid()))
       append_pid (pids, cell.pid); end_foreach_child(); }
      } end_foreach_neighbor(); }
    }
  }
  else
     { foreach_neighbor(1) {
      if ((cell.pid >= 0 && cell.pid != pid()))
 append_pid (pids, cell.pid);
      if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
  { foreach_child()
   if ((cell.pid >= 0 && cell.pid != pid()))
     append_pid (pids, cell.pid); end_foreach_child(); }
    } end_foreach_neighbor(); }
  return pids->len/sizeof(int);
}

static int root_pids (Point point, Array * pids)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 678 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"

   { foreach_child()
    if ((cell.pid >= 0 && cell.pid != pid()))
      append_pid (pids, cell.pid); end_foreach_child(); }
  return pids->len/sizeof(int);
}







static void rcv_pid_row (RcvPid * m, int l, int * row)
{
  for (int i = 0; i < npe(); i++)
    row[i] = 0;
  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0)
      row[rcv->pid] = rcv->halo[l].n;
  }
}

void check_snd_rcv_matrix (SndRcv * sndrcv, const char * name)
{
  int maxlevel = depth();
  mpi_all_reduce (maxlevel, MPI_INT, MPI_MAX);
  int * row = ((int *) pmalloc ((npe())*sizeof(int),__func__,__FILE__,__LINE__));
  for (int l = 0; l <= maxlevel; l++) {
    int status = 0;
    if (pid() == 0) {


      int ** send = matrix_new (npe(), npe(), sizeof(int));
      int ** receive = matrix_new (npe(), npe(), sizeof(int));
      rcv_pid_row (sndrcv->snd, l, row);
      MPI_Gather (row, npe(), MPI_INT, &send[0][0], npe(), MPI_INT, 0,
    MPI_COMM_WORLD);
      rcv_pid_row (sndrcv->rcv, l, row);
      MPI_Gather (row, npe(), MPI_INT, &receive[0][0], npe(), MPI_INT, 0,
    MPI_COMM_WORLD);

      int * astatus = ((int *) pmalloc ((npe())*sizeof(int),__func__,__FILE__,__LINE__));
      for (int i = 0; i < npe(); i++)
 astatus[i] = 0;
      for (int i = 0; i < npe(); i++)
 for (int j = 0; j < npe(); j++)
   if (send[i][j] != receive[j][i]) {
     fprintf (qstderr(), "%s: %d sends    %d to   %d at level %d\n",
       name, i, send[i][j], j, l);
     fprintf (qstderr(), "%s: %d receives %d from %d at level %d\n",
       name, j, receive[j][i], i, l);
     fflush (qstderr());
     for (int k = i - 2; k <= i + 2; k++)
       if (k >= 0 && k < npe())
  astatus[k] = 1;
     for (int k = j - 2; k <= j + 2; k++)
       if (k >= 0 && k < npe())
  astatus[k] = 1;
   }
      MPI_Scatter (astatus, 1, MPI_INT, &status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      pfree (astatus,__func__,__FILE__,__LINE__);

      matrix_free (send);
      matrix_free (receive);
    }
    else {
      rcv_pid_row (sndrcv->snd, l, row);
      MPI_Gather (row, npe(), MPI_INT, NULL, npe(), MPI_INT, 0, MPI_COMM_WORLD);
      rcv_pid_row (sndrcv->rcv, l, row);
      MPI_Gather (row, npe(), MPI_INT, NULL, npe(), MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Scatter (NULL, 1, MPI_INT, &status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    if (status) {
      fprintf (qstderr(),
        "check_snd_rcv_matrix \"%s\" failed\n"
        "Calling debug_mpi(NULL)...\n"
        "Aborting...\n",
        name);
      fflush (qstderr());
      debug_mpi (NULL);
      MPI_Abort (MPI_COMM_WORLD, -3);
    }
  }
  pfree (row,__func__,__FILE__,__LINE__);
}

static bool has_local_child (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 767 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"

   { foreach_child()
    if (is_local(cell))
      return true; end_foreach_child(); }
  return false;
}


void mpi_boundary_update_buffers()
{ trace ("mpi_boundary_update_buffers", "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h", 776);
  if (npe() == 1)
    { ; end_trace("mpi_boundary_update_buffers", "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h", 778);  return; }

  prof_start ("mpi_boundary_update_buffers");

  MpiBoundary * m = (MpiBoundary *) mpi_boundary;
  SndRcv * mpi_level = &m->mpi_level;
  SndRcv * mpi_level_root = &m->mpi_level_root;
  SndRcv * restriction = &m->restriction;

  snd_rcv_free (mpi_level);
  snd_rcv_free (mpi_level_root);
  snd_rcv_free (restriction);

  static const unsigned short used = 1 << user;
   { foreach_cell(){

#line 792 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"
 {
    if (is_active(cell) && !is_border(cell))



      continue;

    if (cell.neighbors) {

      Array pids = {NULL, 0, 0};
      int n = locals_pids (point, &pids);
      if (n) {
  { foreach_child()
   if (is_local(cell))
     for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
       rcv_pid_append (mpi_level->snd, *p, point); end_foreach_child(); }
 pfree (pids.p,__func__,__FILE__,__LINE__);
      }

      bool locals = false;
      if (is_leaf(cell)) {
 if ((cell.pid >= 0 && cell.pid != pid())) {
   Point p = point;
    { foreach_neighbor(1)
     if ((is_local(cell) &&
   ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) || is_local_prolongation (point, p))) ||
  is_root(point))
       locals = true, foreach_neighbor_break(); end_foreach_neighbor(); }
 }
      }
      else
  { foreach_neighbor(1)
   if (is_local(cell) || is_root(point))
     locals = true, foreach_neighbor_break(); end_foreach_neighbor(); }
      if (locals)
  { foreach_child()
   if ((cell.pid >= 0 && cell.pid != pid()))
            rcv_pid_append (mpi_level->rcv, cell.pid, point),
       cell.flags |= used; end_foreach_child(); }


      if (!is_leaf(cell)) {

 if (is_local(cell)) {
   Array pids = {NULL, 0, 0};

   int n = root_pids (point, &pids);
   if (n) {
      { foreach_neighbor()
       for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
  if (cell.pid >= 0 && cell.pid != *p)
    rcv_pid_append (mpi_level_root->snd, *p, point); end_foreach_neighbor(); }

     for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
       rcv_pid_append (restriction->snd, *p, point);
     pfree (pids.p,__func__,__FILE__,__LINE__);
   }
 }

 else if ((cell.pid >= 0 && cell.pid != pid())) {
   bool root = false;
    { foreach_child()
     if (is_local(cell))
       root = true, foreach_child_break(); end_foreach_child(); }
   if (root) {
     int pid = cell.pid;
      { foreach_neighbor()
       if ((cell.pid >= 0 && cell.pid != pid()))
  rcv_pid_append (mpi_level_root->rcv, pid, point),
    cell.flags |= used; end_foreach_neighbor(); }

     rcv_pid_append (restriction->rcv, pid, point);
   }
 }
      }
    }


    if (level > 0) {
      if (is_local(cell)) {

 Array pids = {NULL, 0, 0};
 if ((aparent(0,0,0).pid >= 0 && aparent(0,0,0).pid != pid()))
   append_pid (&pids, aparent(0,0,0).pid);
 int n = root_pids (parent, &pids);
 if (n) {
   for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
     rcv_pid_append (restriction->snd, *p, point);
   pfree (pids.p,__func__,__FILE__,__LINE__);
 }
      }
      else if ((cell.pid >= 0 && cell.pid != pid())) {

 if (is_local(aparent(0,0,0)) || has_local_child (parent))
   rcv_pid_append (restriction->rcv, cell.pid, point);
      }
    }
  } } end_foreach_cell(); }





  static const unsigned short keep = 1 << (user + 1);
  for (int l = depth(); l >= 0; l--)
     { foreach_cell(){

#line 897 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"

      if (level == l) {
 if (level > 0 && (cell.pid < 0 || is_local(cell) || (cell.flags & used)))
   aparent(0,0,0).flags |= keep;
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) && !(cell.flags & keep))
   coarsen_cell (point, NULL);
 cell.flags &= ~(used|keep);
 continue;
      } } end_foreach_cell(); }


  m->send->len = m->receive->len = 0;
  rcv_pid_append_pids (mpi_level->snd, m->send);
  rcv_pid_append_pids (mpi_level_root->snd, m->send);
  rcv_pid_append_pids (mpi_level->rcv, m->receive);
  rcv_pid_append_pids (mpi_level_root->rcv, m->receive);

  prof_stop();
#line 929 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"
 end_trace("mpi_boundary_update_buffers", "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h", 929); }


void mpi_boundary_refine (scalar * list)
{ trace ("mpi_boundary_refine", "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h", 933);
  prof_start ("mpi_boundary_refine");

  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;


  Array * snd = mpi->send;
  MPI_Request r[2*snd->len/sizeof(int)];
  int nr = 0;
  for (int i = 0, * dest = snd->p; i < snd->len/sizeof(int); i++,dest++) {
    int len = ((Tree *)grid)->refined.n;
    MPI_Isend (&((Tree *)grid)->refined.n, 1, MPI_INT, *dest,
        (128), MPI_COMM_WORLD, &r[nr++]);
    if (len > 0)
      MPI_Isend (((Tree *)grid)->refined.p, sizeof(Index)/sizeof(int)*len,
   MPI_INT, *dest, (128), MPI_COMM_WORLD, &r[nr++]);
  }



  Array * rcv = mpi->receive;
  Cache rerefined = {NULL, 0, 0};
  for (int i = 0, * source = rcv->p; i < rcv->len/sizeof(int); i++,source++) {
    int len;
    mpi_recv_check (&len, 1, MPI_INT, *source, (128),
      MPI_COMM_WORLD, MPI_STATUS_IGNORE,
      "mpi_boundary_refine (len)");
    if (len > 0) {
      Index p[len];
      mpi_recv_check (p, sizeof(Index)/sizeof(int)*len,
        MPI_INT, *source, (128),
        MPI_COMM_WORLD, MPI_STATUS_IGNORE,
        "mpi_boundary_refine (p)");
      Cache refined = {p, len, len};
       { foreach_cache (refined){

#line 967 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"

 if (level <= depth() && allocated(0,0,0)) {
   if (is_leaf(cell)) {
     bool neighbors = false;
      { foreach_neighbor()
       if (allocated(0,0,0) && (is_active(cell) || is_local(aparent(0,0,0))))
  neighbors = true, foreach_neighbor_break(); end_foreach_neighbor(); }

     if (neighbors)
       refine_cell (point, list, 0, &rerefined);
   }
 } } end_foreach_cache(); }
    }
  }


  if (nr)
    MPI_Waitall (nr, r, MPI_STATUSES_IGNORE);


  pfree (((Tree *)grid)->refined.p,__func__,__FILE__,__LINE__);
  ((Tree *)grid)->refined = rerefined;

  prof_stop();



  mpi_all_reduce (rerefined.n, MPI_INT, MPI_SUM);
  if (rerefined.n)
    mpi_boundary_refine (list);
 end_trace("mpi_boundary_refine", "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h", 997); }

static void check_depth()
{
#line 1032 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"
}

typedef struct {
  int refined, leaf;
} Remote;




void mpi_boundary_coarsen (int l, int too_fine)
{ trace ("mpi_boundary_coarsen", "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h", 1042);
  if (npe() == 1)
    { ; end_trace("mpi_boundary_coarsen", "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h", 1044);  return; }

  check_depth();

  assert (sizeof(Remote) == sizeof(double));

  scalar remote= new_scalar("remote");
   { foreach_cell(){

#line 1051 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"
 {
    if (level == l) {
      if (is_local(cell)) {
 ((Remote *)&val(remote,0,0,0))->refined = (!is_leaf (cell) && cell.neighbors && cell.pid >= 0);
 ((Remote *)&val(remote,0,0,0))->leaf = is_leaf(cell);
      }
      else {
 ((Remote *)&val(remote,0,0,0))->refined = true;
 ((Remote *)&val(remote,0,0,0))->leaf = false;
      }
      continue;
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }
  mpi_boundary_level (mpi_boundary, ((scalar []){remote,{-1}}), l);

   { foreach_cell(){

#line 1068 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"
 {
    if (level == l) {
      if (!is_local(cell)) {
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) && !((Remote *)&val(remote,0,0,0))->refined)
   coarsen_cell_recursive (point, NULL);
 else if (is_leaf(cell) && cell.neighbors && ((Remote *)&val(remote,0,0,0))->leaf) {
   int pid = cell.pid;
    { foreach_child()
     cell.pid = pid; end_foreach_child(); }
 }
      }
      continue;
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

  check_depth();

  if (l > 0) {
     { foreach_cell(){

#line 1088 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"
 {
      if (level == l) {
 val(remote,0,0,0) = is_local(cell) ? cell.neighbors : 0;
 continue;
      }
      if (is_leaf(cell))
 continue;
    } } end_foreach_cell(); }
    mpi_boundary_level (mpi_boundary, ((scalar []){remote,{-1}}), l);
     { foreach_cell(){

#line 1097 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"
 {
      if (level == l)
 if (!is_local(cell) && is_local(aparent(0,0,0)) && val(remote,0,0,0)) {
   aparent(0,0,0).flags &= ~too_fine;
   continue;
 }
      if (is_leaf(cell))
 continue;
    } } end_foreach_cell(); }
  }
 delete (((scalar []){remote,{-1}}));  end_trace("mpi_boundary_coarsen", "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h", 1107); }

static void flag_border_cells()
{
   { foreach_cell(){

#line 1111 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"
 {
    if (is_active(cell)) {
      short flags = cell.flags & ~border;
       { foreach_neighbor() {
 if (!is_local(cell) || (level > 0 && !is_local(aparent(0,0,0))))
   flags |= border, foreach_neighbor_break();

 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
    { foreach_child()
     if (!is_local(cell))
       flags |= border, foreach_child_break(); end_foreach_child(); }
 if (flags & border)
   foreach_neighbor_break();
      } end_foreach_neighbor(); }
      cell.flags = flags;
    }
    else {
      cell.flags &= ~border;

    }
    if (is_leaf(cell)) {
      if (cell.neighbors) {
  { foreach_child()
   cell.flags &= ~border; end_foreach_child(); }
 if (is_border(cell)) {
   bool remote = false;
    { foreach_neighbor (2/2)
     if (!is_local(cell))
       remote = true, foreach_neighbor_break(); end_foreach_neighbor(); }
   if (remote)
      { foreach_child()
       cell.flags |= border; end_foreach_child(); }
 }
      }
      continue;
    }
  } } end_foreach_cell(); }
}

static int balanced_pid (long index, long nt, int nproc)
{
  long ne = max(1, nt/nproc), nr = nt % nproc;
  int pid = index < nr*(ne + 1) ?
    index/(ne + 1) :
    nr + (index - nr*(ne + 1))/ne;
  return min(nproc - 1, pid);
}



void mpi_partitioning()
{ trace ("mpi_partitioning", "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h", 1162);
  prof_start ("mpi_partitioning");

  long nt = 0;
   { foreach(){

#line 1166 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"

    nt++; } end_foreach(); }


  long i = 0;
  ((Tree *)grid)->dirty = true;
   { foreach_cell_post (is_active (cell)){

#line 1172 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"

    if (is_active (cell)) {
      if (is_leaf (cell)) {
 cell.pid = balanced_pid (i++, nt, npe());
 if (cell.neighbors > 0) {
   int pid = cell.pid;
    { foreach_child()
     cell.pid = pid; end_foreach_child(); }
 }
 if (!is_local(cell))
   cell.flags &= ~active;
      }
      else {
 cell.pid = child(0,0,0).pid;
 bool inactive = true;
  { foreach_child()
   if (is_active(cell))
     inactive = false, foreach_child_break(); end_foreach_child(); }
 if (inactive)
   cell.flags &= ~active;
      }
    } } end_foreach_cell_post(); }

  flag_border_cells();

  prof_stop();

  mpi_boundary_update_buffers();
 end_trace("mpi_partitioning", "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h", 1200); }

void restore_mpi (FILE * fp, scalar * list1)
{
  long index = 0, nt = 0, start = ftell (fp);
  scalar size= new_scalar("size"), * list = list_concat (((scalar []){size,{-1}}), list1);;
  long offset = sizeof(double)*list_len(list);


  static const unsigned short set = 1 << user;
  scalar * listm = is_constant(cm) ? NULL : (scalar *)((vector []){{fm.x,fm.y,fm.z},{{-1},{-1},{-1}}});
   { foreach_cell(){

#line 1211 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"

    if (balanced_pid (index, nt, npe()) <= pid()) {
      unsigned flags;
      if (fread (&flags, sizeof(unsigned), 1, fp) != 1) {
 fprintf (qstderr(), "restore(): error: expecting 'flags'\n");
 exit (1);
      }
      if (list) for (scalar s = *list, *_i69 = list; ((scalar *)&s)->i >= 0; s = *++_i69) {
 double val;
 if (fread (&val, sizeof(double), 1, fp) != 1) {
   fprintf (qstderr(), "restore(): error: expecting scalar\n");
   exit (1);
 }
 if (s.i != INT_MAX)
   val(s,0,0,0) = val;
      }
      if (level == 0)
 nt = val(size,0,0,0);
      cell.pid = balanced_pid (index, nt, npe());
      cell.flags |= set;
      if (!(flags & leaf) && is_leaf(cell)) {
 if (balanced_pid (index + val(size,0,0,0) - 1, nt, npe()) < pid()) {
   fseek (fp, (sizeof(unsigned) + offset)*(val(size,0,0,0) - 1), SEEK_CUR);
   index += val(size,0,0,0);
   continue;
 }
 refine_cell (point, listm, 0, NULL);
      }
      index++;
      if (is_leaf(cell))
 continue;
    } } end_foreach_cell(); }


  fseek (fp, start, SEEK_SET);
  index = 0;
   { foreach_cell(){

#line 1247 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"
 {
    unsigned flags;
    if (fread (&flags, sizeof(unsigned), 1, fp) != 1) {
      fprintf (qstderr(), "restore(): error: expecting 'flags'\n");
      exit (1);
    }
    if (cell.flags & set)
      fseek (fp, offset, SEEK_CUR);
    else {
      if (list) for (scalar s = *list, *_i70 = list; ((scalar *)&s)->i >= 0; s = *++_i70) {
 double val;
 if (fread (&val, sizeof(double), 1, fp) != 1) {
   fprintf (qstderr(), "restore(): error: expecting a scalar\n");
   exit (1);
 }
 if (s.i != INT_MAX)
   val(s,0,0,0) = val;
      }
      cell.pid = balanced_pid (index, nt, npe());
      if (is_leaf(cell) && cell.neighbors) {
 int pid = cell.pid;
  { foreach_child()
   cell.pid = pid; end_foreach_child(); }
      }
    }
    if (!(flags & leaf) && is_leaf(cell)) {
      bool locals = false;
       { foreach_neighbor(1)
 if ((cell.flags & set) && (is_local(cell) || is_root(point)))
   locals = true, foreach_neighbor_break(); end_foreach_neighbor(); }
      if (locals)
 refine_cell (point, listm, 0, NULL);
      else {
 fseek (fp, (sizeof(unsigned) + offset)*(val(size,0,0,0) - 1), SEEK_CUR);
 index += val(size,0,0,0);
 continue;
      }
    }
    index++;
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }


   { foreach_cell_post (is_active (cell)){

#line 1291 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"
 {
    cell.flags &= ~set;
    if (is_active (cell)) {
      if (is_leaf (cell)) {
 if (cell.neighbors > 0) {
   int pid = cell.pid;
    { foreach_child()
     cell.pid = pid; end_foreach_child(); }
 }
 if (!is_local(cell))
   cell.flags &= ~active;
      }
      else if (!is_local(cell)) {
 bool inactive = true;
  { foreach_child()
   if (is_active(cell))
     inactive = false, foreach_child_break(); end_foreach_child(); }
 if (inactive)
   cell.flags &= ~active;
      }
    }
  } } end_foreach_cell_post(); }

  flag_border_cells();

  mpi_boundary_update (list);
  pfree (list,__func__,__FILE__,__LINE__);
 delete (((scalar []){size,{-1}})); }
#line 1340 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"

double z_indexing (scalar index, bool leaves)
{ trace ("z_indexing", "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h", 1342);



  scalar size= new_scalar("size");
  subtree_size (size, leaves);






  double maxi = -1.;
  if (pid() == 0)
     { foreach_level(0){

#line 1356 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"

      maxi = val(size,0,0,0) - 1.; } end_foreach_level(); }




   { foreach_level(0){

#line 1362 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"

    val(index,0,0,0) = 0; } end_foreach_level(); }
  for (int l = 0; l < depth(); l++) {
    { Boundary ** _i = boundaries, * _b; while ((_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar []){index,{-1}}), l); };
     { foreach_cell(){

#line 1366 "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h"
 {
      if (level == l) {
 if (is_leaf(cell)) {
   if (is_local(cell) && cell.neighbors) {
     int i = val(index,0,0,0);
      { foreach_child()
       val(index,0,0,0) = i; end_foreach_child(); }
   }
 }
 else {
   bool loc = is_local(cell);
   if (!loc)
      { foreach_child()
       if (is_local(cell))
  loc = true, foreach_child_break(); end_foreach_child(); }
   if (loc) {
     int i = val(index,0,0,0) + !leaves;
      { foreach_child() {
       val(index,0,0,0) = i;
       i += val(size,0,0,0);
     } end_foreach_child(); }
   }
 }
 continue;
      }
      if (is_leaf(cell))
 continue;
    } } end_foreach_cell(); }
  }
  { Boundary ** _i = boundaries, * _b; while ((_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar []){index,{-1}}), depth()); };

  { double _ret =  maxi; delete (((scalar []){size,{-1}}));  end_trace("z_indexing", "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h", 1397);  return _ret; }
 delete (((scalar []){size,{-1}}));  end_trace("z_indexing", "/home/miro3/Documents/Programming/basilisk/src/grid/tree-mpi.h", 1398); }
#line 1744 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"
#line 1 "grid/balance.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/grid/balance.h"


typedef struct {
  short leaf, prolongation;
  int pid;
} NewPid;



#if TRASH
# define is_newpid() (!isnan(val(newpid,0,0,0)) && ((NewPid *)&val(newpid,0,0,0))->pid > 0)
#else
# define is_newpid() (((NewPid *)&val(newpid,0,0,0))->pid > 0)
#endif

Array * linear_tree (size_t size, scalar newpid)
{
  const unsigned short sent = 1 << user, next = 1 << (user + 1);
  Array * a = array_new();

   { foreach_cell_post_all (true){

#line 21 "/home/miro3/Documents/Programming/basilisk/src/grid/balance.h"

    if (level > 0 && (cell.flags & (sent|next)))
      aparent(0,0,0).flags |= next; } end_foreach_cell_post_all(); }

  bool empty = true;
   { foreach_cell_all(){

#line 26 "/home/miro3/Documents/Programming/basilisk/src/grid/balance.h"
 {
    if (cell.flags & sent) {
      array_append (a, &cell, size);
      cell.flags &= ~sent;
      empty = false;
    }
    else {
      if (cell.pid >= 0 && ((NewPid *)&val(newpid,0,0,0))->leaf)
 assert (is_leaf(cell));
      if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0)) {


 bool prolo = false;
  { foreach_child()
   if (((NewPid *)&val(newpid,0,0,0))->prolongation)
     prolo = true; end_foreach_child(); }
 if (prolo) {

   cell.flags |= leaf;
   array_append (a, &cell, sizeof(Cell));
   cell.flags &= ~leaf;
 }
 else
   array_append (a, &cell, sizeof(Cell));
      }
      else
 array_append (a, &cell, sizeof(Cell));
    }
    if (cell.flags & next)
      cell.flags &= ~next;
    else
      continue;
  } } end_foreach_cell_all(); }

  if (empty)
    a->len = 0;
  return a;
}

#define foreach_tree(t, size, list)\
{\
  const unsigned short _sent = 1 << user, _next = 1 << (user + 1);\
  scalar * _list = list;\
  char * _i = (char *) (t)->p;\
  foreach_cell_all() {\
    Cell * c = (Cell *) _i;\
    if (c->flags & _sent) {\
      _i += size;\

#line 74


#define end_foreach_tree()\
    }\
    else\
      _i += sizeof(Cell);\
    if (c->flags & _next) {\
      assert (c->neighbors);\
      if (!(c->flags & leaf) && is_leaf(cell) &&\
   (!is_newpid() || !((NewPid *)&val(newpid,0,0,0))->leaf))\
\
 refine_cell (point, _list, 0, NULL);\
      else if (!cell.neighbors)\
\
 alloc_children (point);\
    }\
    else\
      continue;\
  } end_foreach_cell_all();\
}\

#line 94


Array * neighborhood (scalar newpid, int nextpid, FILE * fp)
{
  const unsigned short sent = 1 << user;
   { foreach_cell(){

#line 99 "/home/miro3/Documents/Programming/basilisk/src/grid/balance.h"
 {

    bool root = false;
    if ((!is_local(cell) || ((NewPid *)&val(newpid,0,0,0))->pid - 1 != nextpid) && (!is_leaf (cell) && cell.neighbors && cell.pid >= 0)) {
       { foreach_child()
 if (is_local(cell) && ((NewPid *)&val(newpid,0,0,0))->pid - 1 == nextpid)
   root = true, foreach_child_break(); end_foreach_child(); }
      if (root && cell.pid != nextpid) {
  { foreach_neighbor()
   if (cell.pid != nextpid && is_newpid()) {
     if (fp)
       fprintf (fp, "%g %g %g %d %d root\n",
         x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
     cell.flags |= sent;
   } end_foreach_neighbor(); }
      }
    }

    if ((is_local(cell) && ((NewPid *)&val(newpid,0,0,0))->pid - 1 == nextpid) || root) {
       { foreach_neighbor(1)
 if (cell.neighbors && cell.pid != nextpid)
    { foreach_child()
     if (cell.pid != nextpid && is_newpid()) {
       if (fp)
  fprintf (fp, "%g %g %g %d %d nextpid\n",
    x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
       cell.flags |= sent;
     } end_foreach_child(); } end_foreach_neighbor(); }
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

  return linear_tree (sizeof(Cell) + datasize, newpid);
}

static void send_tree (Array * a, int to, MPI_Request * r)
{
  MPI_Isend (&a->len, 1, MPI_LONG, to, (256), MPI_COMM_WORLD, &r[0]);
  if (a->len > 0) {
    MPI_Isend (a->p, a->len, MPI_BYTE, to, (256), MPI_COMM_WORLD, &r[1]);
    ((Tree *)grid)->dirty = true;
  }
}

static void receive_tree (int from, scalar newpid, FILE * fp)
{
  Array a;
  mpi_recv_check (&a.len, 1, MPI_LONG, from, (256),
    MPI_COMM_WORLD, MPI_STATUS_IGNORE, "receive_tree (len)");
  if (a.len > 0) {
    a.p = pmalloc (a.len,__func__,__FILE__,__LINE__);
    if (fp)
      fprintf (fp, "receiving %ld from %d\n", a.len, from);
    mpi_recv_check (a.p, a.len, MPI_BYTE, from, (256),
      MPI_COMM_WORLD, MPI_STATUS_IGNORE, "receive_tree (p)");

     { foreach_tree (&a, sizeof(Cell) + datasize, NULL){

#line 156 "/home/miro3/Documents/Programming/basilisk/src/grid/balance.h"
 {
      memcpy (((char *)&cell) + sizeof(Cell), ((char *)c) + sizeof(Cell),
       datasize);
      assert (((NewPid *)&val(newpid,0,0,0))->pid > 0);
      if (fp)
 fprintf (fp, "%g %g %g %d %d %d %d %d %d recv\n",
   x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid,
   c->flags & leaf,
   cell.flags & leaf, from, ((NewPid *)&val(newpid,0,0,0))->leaf);
    } } end_foreach_tree(); }
    pfree (a.p,__func__,__FILE__,__LINE__);
    ((Tree *)grid)->dirty = true;
  }
}

static void wait_tree (Array * a, MPI_Request * r)
{
  MPI_Wait (&r[0], MPI_STATUS_IGNORE);
  if (a->len > 0)
    MPI_Wait (&r[1], MPI_STATUS_IGNORE);
}

static void check_flags()
{







}

struct {
  int min;
  bool leaves;

  int npe;
} mpi = {
  1,
  true
};


bool balance()
{ trace ("balance", "/home/miro3/Documents/Programming/basilisk/src/grid/balance.h", 201);
  if (npe() == 1)
    { bool _ret =  false; end_trace("balance", "/home/miro3/Documents/Programming/basilisk/src/grid/balance.h", 203);  return _ret; }

  assert (sizeof(NewPid) == sizeof(double));

  check_flags();

  long nl = 0, nt = 0;
   { foreach_cell(){

#line 210 "/home/miro3/Documents/Programming/basilisk/src/grid/balance.h"
 {
    if (is_local(cell)) {
      nt++;
      if (is_leaf(cell))
 nl++;
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

  grid->n = grid->tn = nl;
  grid->maxdepth = depth();
  long nmin = nl, nmax = nl;

  mpi_all_reduce (nmax, MPI_LONG, MPI_MAX);
  mpi_all_reduce (nmin, MPI_LONG, MPI_MIN);
  mpi_all_reduce (grid->tn, MPI_LONG, MPI_SUM);
  mpi_all_reduce (grid->maxdepth, MPI_INT, MPI_MAX);
  if (mpi.leaves)
    nt = grid->tn;
  else
    mpi_all_reduce (nt, MPI_LONG, MPI_SUM);

  long ne = max(1, nt/npe());

  if (ne < mpi.min) {
    mpi.npe = max(1, nt/mpi.min);
    ne = max(1, nt/mpi.npe);
  }
  else
    mpi.npe = npe();

  if (nmax - nmin <= 1)
    { bool _ret =  false; end_trace("balance", "/home/miro3/Documents/Programming/basilisk/src/grid/balance.h", 243);  return _ret; }

  scalar newpid= new_scalar("newpid");
  double zn = z_indexing (newpid, mpi.leaves);
  if (pid() == 0)
    assert (zn + 1 == nt);

  FILE * fp = NULL;
#line 260 "/home/miro3/Documents/Programming/basilisk/src/grid/balance.h"
  bool next = false, prev = false;
   { foreach_cell_all(){

#line 261 "/home/miro3/Documents/Programming/basilisk/src/grid/balance.h"
 {
    if (is_local(cell)) {
      int pid = balanced_pid (val(newpid,0,0,0), nt, mpi.npe);
      pid = clamp (pid, cell.pid - 1, cell.pid + 1);
      if (pid == pid() + 1)
 next = true;
      else if (pid == pid() - 1)
 prev = true;
      ((NewPid *)&val(newpid,0,0,0))->pid = pid + 1;
      ((NewPid *)&val(newpid,0,0,0))->leaf = is_leaf(cell);
      ((NewPid *)&val(newpid,0,0,0))->prolongation = (!is_leaf(cell) && !cell.neighbors && cell.pid >= 0);
      if (fp)
 fprintf (fp, "%g %g %d %d newpid\n", x, y, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
    }
    else
      val(newpid,0,0,0) = 0;
  } } end_foreach_cell_all(); }
  for (int l = 0; l <= depth(); l++)
    { Boundary ** _i = boundaries, * _b; while ((_b = *_i++)) if (_b->level) _b->level (_b, ((scalar []){newpid,{-1}}), l); };
#line 304 "/home/miro3/Documents/Programming/basilisk/src/grid/balance.h"
  Array * anext = next ? neighborhood (newpid, pid() + 1, fp) : array_new();
  Array * aprev = prev ? neighborhood (newpid, pid() - 1, fp) : array_new();

  if (fp)
    fflush (fp);

  check_flags();


  MPI_Request rprev[2], rnext[2];
  if (pid() > 0)
    send_tree (aprev, pid() - 1, rprev);
  if (pid() < npe() - 1)
    send_tree (anext, pid() + 1, rnext);


  if (pid() < npe() - 1)
    receive_tree (pid() + 1, newpid, fp);
  if (pid() > 0)
    receive_tree (pid() - 1, newpid, fp);


  if (pid() > 0)
    wait_tree (aprev, rprev);
  array_free (aprev);
  if (pid() < npe() - 1)
    wait_tree (anext, rnext);
  array_free (anext);

  if (fp)
    fflush (fp);


  int pid_changed = false;
   { foreach_cell_all(){

#line 338 "/home/miro3/Documents/Programming/basilisk/src/grid/balance.h"
 {
    if (cell.pid >= 0) {
      if (is_newpid()) {
 if (fp)
   fprintf (fp, "%g %g %g %d %d %d %d %d new\n",
     x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid,
     is_leaf(cell), cell.neighbors, ((NewPid *)&val(newpid,0,0,0))->leaf);
 if (cell.pid != ((NewPid *)&val(newpid,0,0,0))->pid - 1) {
   cell.pid = ((NewPid *)&val(newpid,0,0,0))->pid - 1;
   cell.flags &= ~(active|border);
   if (is_local(cell))
     cell.flags |= active;
   pid_changed = true;
 }
 if (((NewPid *)&val(newpid,0,0,0))->leaf && !is_leaf(cell) && cell.neighbors)
   coarsen_cell_recursive (point, NULL);
      }
      else if (level > 0 && ((NewPid *)&coarse(newpid,0,0,0))->leaf)
 cell.pid = aparent(0,0,0).pid;
    }

    if (!cell.neighbors && allocated_child(0,0,0)) {
      if (fp)
 fprintf (fp, "%g %g %g %d %d freechildren\n",
   x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
      free_children (point);
    }
  } } end_foreach_cell_all(); }

  if (((Tree *)grid)->dirty || pid_changed) {


     { foreach_cell_post (!is_leaf (cell)){

#line 370 "/home/miro3/Documents/Programming/basilisk/src/grid/balance.h"

      if (!is_leaf(cell) && !is_local(cell)) {
 unsigned short flags = cell.flags & ~active;
  { foreach_child()
   if (is_active(cell))
     flags |= active, foreach_child_break(); end_foreach_child(); }
 cell.flags = flags;
      } } end_foreach_cell_post(); }

    flag_border_cells();
    pid_changed = true;
  }

  if (fp)
    fclose (fp);

  mpi_all_reduce (pid_changed, MPI_INT, MPI_MAX);
  if (pid_changed)
    mpi_boundary_update_buffers();

  { bool _ret =  pid_changed; delete (((scalar []){newpid,{-1}}));  end_trace("balance", "/home/miro3/Documents/Programming/basilisk/src/grid/balance.h", 390);  return _ret; }
 delete (((scalar []){newpid,{-1}}));  end_trace("balance", "/home/miro3/Documents/Programming/basilisk/src/grid/balance.h", 391); }

void mpi_boundary_update (scalar * list)
{
  mpi_boundary_update_buffers();
  grid->tn = 0;
  boundary (list);
  while (balance());
}
#line 1745 "/home/miro3/Documents/Programming/basilisk/src/grid/tree.h"
#else
void mpi_boundary_refine (scalar * list){}
void mpi_boundary_coarsen (int a, int b){}
void mpi_boundary_update (scalar * list) {
  boundary (list);
}
#endif
#line 4 "/home/miro3/Documents/Programming/basilisk/src/grid/octree.h"

void octree_methods() {
  tree_methods();
}
#line 2 "milk_crown.c"
#line 1 "fractions.h"
#line 1 "./fractions.h"
#line 11 "./fractions.h"
#line 1 "geometry.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/geometry.h"
#line 28 "/home/miro3/Documents/Programming/basilisk/src/geometry.h"
double line_alpha (double c, coord n)
{
  double alpha, n1, n2;

  n1 = fabs (n.x); n2 = fabs (n.y);
  if (n1 > n2)
    swap (double, n1, n2);

  c = clamp (c, 0., 1.);
  double v1 = n1/2.;
  if (c <= v1/n2)
    alpha = sqrt (2.*c*n1*n2);
  else if (c <= 1. - v1/n2)
    alpha = c*n2 + v1;
  else
    alpha = n1 + n2 - sqrt (2.*n1*n2*(1. - c));

  if (n.x < 0.)
    alpha += n.x;
  if (n.y < 0.)
    alpha += n.y;

  return alpha - (n.x + n.y)/2.;
}



double plane_alpha (double c, coord n)
{
  double alpha;
  coord n1;

  n1.x = fabs (n.x); n1.y = fabs (n.y); n1.z = fabs (n.z);

  double m1, m2, m3;
  m1 = min(n1.x, n1.y);
  m3 = max(n1.x, n1.y);
  m2 = n1.z;
  if (m2 < m1) {
    double tmp = m1;
    m1 = m2;
    m2 = tmp;
  }
  else if (m2 > m3) {
    double tmp = m3;
    m3 = m2;
    m2 = tmp;
  }
  double m12 = m1 + m2;
  double pr = max(6.*m1*m2*m3, 1e-50);
  double V1 = m1*m1*m1/pr;
  double V2 = V1 + (m2 - m1)/(2.*m3), V3;
  double mm;
  if (m3 < m12) {
    mm = m3;
    V3 = (m3*m3*(3.*m12 - m3) + m1*m1*(m1 - 3.*m3) + m2*m2*(m2 - 3.*m3))/pr;
  }
  else {
    mm = m12;
    V3 = mm/(2.*m3);
  }

  c = clamp (c, 0., 1.);
  double ch = min(c, 1. - c);
  if (ch < V1)
    alpha = pow (pr*ch, 1./3.);
  else if (ch < V2)
    alpha = (m1 + sqrt(m1*m1 + 8.*m2*m3*(ch - V1)))/2.;
  else if (ch < V3) {
    double p12 = sqrt (2.*m1*m2);
    double q = 3.*(m12 - 2.*m3*ch)/(4.*p12);
    double teta = acos(clamp(q,-1.,1.))/3.;
    double cs = cos(teta);
    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + m12;
  }
  else if (m12 < m3)
    alpha = m3*ch + mm/2.;
  else {
    double p = m1*(m2 + m3) + m2*m3 - 1./4., p12 = sqrt(p);
    double q = 3.*m1*m2*m3*(1./2. - ch)/(2.*p*p12);
    double teta = acos(clamp(q,-1.,1.))/3.;
    double cs = cos(teta);
    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + 1./2.;
  }
  if (c > 1./2.) alpha = 1. - alpha;

  if (n.x < 0.)
    alpha += n.x;
  if (n.y < 0.)
    alpha += n.y;
  if (n.z < 0.)
    alpha += n.z;

  return alpha - (n.x + n.y + n.z)/2.;;
}
#line 133 "/home/miro3/Documents/Programming/basilisk/src/geometry.h"
double line_area (double nx, double ny, double alpha)
{
  double a, v, area;

  alpha += (nx + ny)/2.;
  if (nx < 0.) {
    alpha -= nx;
    nx = - nx;
  }
  if (ny < 0.) {
    alpha -= ny;
    ny = - ny;
  }

  if (alpha <= 0.)
    return 0.;

  if (alpha >= nx + ny)
    return 1.;

  if (nx < 1e-10)
    area = alpha/ny;
  else if (ny < 1e-10)
    area = alpha/nx;
  else {
    v = sq(alpha);

    a = alpha - nx;
    if (a > 0.)
      v -= a*a;

    a = alpha - ny;
    if (a > 0.)
      v -= a*a;

    area = v/(2.*nx*ny);
  }

  return clamp (area, 0., 1.);
}



double plane_volume (coord n, double alpha)
{
  double al = alpha + (n.x + n.y + n.z)/2. +
    max(0., -n.x) + max(0., -n.y) + max(0., -n.z);
  if (al <= 0.)
    return 0.;
  double tmp = fabs(n.x) + fabs(n.y) + fabs(n.z);
  if (al >= tmp)
    return 1.;
  if (tmp < 1e-10)
    return 0.;
  double n1 = fabs(n.x)/tmp;
  double n2 = fabs(n.y)/tmp;
  double n3 = fabs(n.z)/tmp;
  al = max(0., min(1., al/tmp));
  double al0 = min(al, 1. - al);
  double b1 = min(n1, n2);
  double b3 = max(n1, n2);
  double b2 = n3;
  if (b2 < b1) {
    tmp = b1;
    b1 = b2;
    b2 = tmp;
  }
  else if (b2 > b3) {
    tmp = b3;
    b3 = b2;
    b2 = tmp;
  }
  double b12 = b1 + b2;
  double bm = min(b12, b3);
  double pr = max(6.*b1*b2*b3, 1e-50);
  if (al0 < b1)
    tmp = al0*al0*al0/pr;
  else if (al0 < b2)
    tmp = 0.5*al0*(al0 - b1)/(b2*b3) + b1*b1*b1/pr;
  else if (al0 < bm)
    tmp = (al0*al0*(3.*b12 - al0) + b1*b1*(b1 - 3.*al0) +
    b2*b2*(b2 - 3.*al0))/pr;
  else if (b12 < b3)
    tmp = (al0 - 0.5*bm)/b3;
  else
    tmp = (al0*al0*(3. - 2.*al0) + b1*b1*(b1 - 3.*al0) +
    b2*b2*(b2 - 3.*al0) + b3*b3*(b3 - 3.*al0))/pr;

  double volume = al <= 0.5 ? tmp : 1. - tmp;
  return clamp (volume, 0., 1.);
}
#line 237 "/home/miro3/Documents/Programming/basilisk/src/geometry.h"
double rectangle_fraction (coord n, double alpha, coord a, coord b)
{
  coord n1;
  {
#line 240
 {
    alpha -= n.x*(b.x + a.x)/2.;
    n1.x = n.x*(b.x - a.x);
  }
#line 240
 {
    alpha -= n.y*(b.y + a.y)/2.;
    n1.y = n.y*(b.y - a.y);
  }
#line 240
 {
    alpha -= n.z*(b.z + a.z)/2.;
    n1.z = n.z*(b.z - a.z);
  }}
  return plane_volume (n1, alpha);
}
#line 277 "/home/miro3/Documents/Programming/basilisk/src/geometry.h"
static coord cube_edge[12][2] = {
  {{0.,0.,0.},{1.,0.,0.}},{{0.,0.,1.},{1.,0.,1.}},
  {{0.,1.,1.},{1.,1.,1.}},{{0.,1.,0.},{1.,1.,0.}},
  {{0.,0.,0.},{0.,1.,0.}},{{0.,0.,1.},{0.,1.,1.}},
  {{1.,0.,1.},{1.,1.,1.}},{{1.,0.,0.},{1.,1.,0.}},
  {{0.,0.,0.},{0.,0.,1.}},{{1.,0.,0.},{1.,0.,1.}},
  {{1.,1.,0.},{1.,1.,1.}},{{0.,1.,0.},{0.,1.,1.}}
};




static int cube_connect[12][2][4] = {
  {{9, 1, 8}, {4, 3, 7}},
  {{6, 2, 5}, {8, 0, 9}},
  {{10, 3, 11}, {5, 1, 6}},
  {{7, 0, 4}, {11, 2, 10}},
  {{3, 7, 0}, {8, 5, 11}},
  {{11, 4, 8}, {1, 6, 2}},
  {{2, 5, 1}, {9, 7, 10}},
  {{10, 6, 9}, {0, 4, 3}},
  {{5, 11, 4}, {0, 9, 1}},
  {{1, 8, 0}, {7, 10, 6}},
  {{6, 9, 7}, {3, 11, 2}},
  {{2, 10, 3}, {4, 8, 5}}
};

int facets (coord n, double alpha, coord v[12], double h)
{
  coord a[12];
  int orient[12];

  for (int i = 0; i < 12; i++) {
    coord e, d;
    double den = 0., t = alpha;
    {
#line 312
 {
      d.x = h*(cube_edge[i][0].x - 0.5);
      e.x = h*(cube_edge[i][1].x - 0.5);
      den += n.x*(e.x - d.x);
      t -= n.x*d.x;
    }
#line 312
 {
      d.y = h*(cube_edge[i][0].y - 0.5);
      e.y = h*(cube_edge[i][1].y - 0.5);
      den += n.y*(e.y - d.y);
      t -= n.y*d.y;
    }
#line 312
 {
      d.z = h*(cube_edge[i][0].z - 0.5);
      e.z = h*(cube_edge[i][1].z - 0.5);
      den += n.z*(e.z - d.z);
      t -= n.z*d.z;
    }}
    orient[i] = -1;
    if (fabs (den) > 1e-10) {
      t /= den;
      if (t >= 0. && t < 1.) {
 double s = - alpha;
 {
#line 323
 {
   a[i].x = d.x + t*(e.x - d.x);
   s += n.x*e.x;
 }
#line 323
 {
   a[i].y = d.y + t*(e.y - d.y);
   s += n.y*e.y;
 }
#line 323
 {
   a[i].z = d.z + t*(e.z - d.z);
   s += n.z*e.z;
 }}
 orient[i] = (s > 0.);
      }
    }
  }

  for (int i = 0; i < 12; i++) {
    int nv = 0, e = i;
    while (orient[e] >= 0) {
      int m = 0, * ne = cube_connect[e][orient[e]];
      v[nv++] = a[e];
      orient[e] = -1;
      while (m < 3 && orient[e] < 0)
 e = ne[m++];
    }
    if (nv > 2)
      return nv;
  }
  return 0;
}






double line_length_center (coord m, double alpha, coord * p)
{
  alpha += (m.x + m.y)/2.;

  coord n = m;
  if (n.x < 0.) {
    alpha -= n.x;
    n.x = - n.x;
  }
  if (n.y < 0.) {
    alpha -= n.y;
    n.y = - n.y;
  }

  p->x = p->y = 0.;

  if (alpha <= 0. || alpha >= n.x + n.y)
    return 0.;

  {
#line 371

    if (n.x < 1e-4) {
      p->x = 0.;
      p->y = (m.y < 0. ? 1. - alpha : alpha) - 0.5;
      return 1.;
    }
#line 371

    if (n.y < 1e-4) {
      p->y = 0.;
      p->x = (m.x < 0. ? 1. - alpha : alpha) - 0.5;
      return 1.;
    }}

  if (alpha >= n.x) {
    p->x += 1.;
    p->y += (alpha - n.x)/n.y;
  }
  else
    p->x += alpha/n.x;

  double ax = p->x, ay = p->y;
  if (alpha >= n.y) {
    p->y += 1.;
    ay -= 1.;
    p->x += (alpha - n.y)/n.x;
    ax -= (alpha - n.y)/n.x;
  }
  else {
    p->y += alpha/n.y;
    ay -= alpha/n.y;
  }

  {
#line 397
 {
    p->x /= 2.;
    p->x = clamp (p->x, 0., 1.);
    if (m.x < 0.)
      p->x = 1. - p->x;
    p->x -= 0.5;
  }
#line 397
 {
    p->y /= 2.;
    p->y = clamp (p->y, 0., 1.);
    if (m.y < 0.)
      p->y = 1. - p->y;
    p->y -= 0.5;
  }}

  return sqrt (ax*ax + ay*ay);
}




double plane_area_center (coord m, double alpha, coord * p)
{
  if (fabs (m.x) < 1e-4) {
    coord n, q;
    n.x = m.y;
    n.y = m.z;
    double length = line_length_center (n, alpha, &q);
    p->x = 0.;
    p->y = q.x;
    p->z = q.y;
    return sq(length);
  }
  if (fabs (m.y) < 1e-4) {
    coord n, q;
    n.x = m.z;
    n.y = m.x;
    double length = line_length_center (n, alpha, &q);
    p->x = q.y;
    p->y = 0.;
    p->z = q.x;
    return sq(length);
  }
  if (fabs (m.z) < 1e-4) {
    double length = line_length_center (m, alpha, p);
    p->z = 0.;
    return sq(length);
  }

  alpha += (m.x + m.y + m.z)/2.;
  coord n = m;
  {
#line 441

    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }
#line 441

    if (n.y < 0.) {
      alpha -= n.y;
      n.y = - n.y;
    }
#line 441

    if (n.z < 0.) {
      alpha -= n.z;
      n.z = - n.z;
    }}

  double amax = n.x + n.y + n.z;
  if (alpha < 0. || alpha > amax) {
    p->x = p->y = p->z = 0.;
    return 0.;
  }

  double area = sq(alpha);
  p->x = p->y = p->z = area*alpha;

  {
#line 456
 {
    double b = alpha - n.x;
    if (b > 0.) {
      area -= b*b;
      p->x -= b*b*(2.*n.x + alpha);
      p->y -= b*b*b;
      p->z -= b*b*b;
    }
  }
#line 456
 {
    double b = alpha - n.y;
    if (b > 0.) {
      area -= b*b;
      p->y -= b*b*(2.*n.y + alpha);
      p->z -= b*b*b;
      p->x -= b*b*b;
    }
  }
#line 456
 {
    double b = alpha - n.z;
    if (b > 0.) {
      area -= b*b;
      p->z -= b*b*(2.*n.z + alpha);
      p->x -= b*b*b;
      p->y -= b*b*b;
    }
  }}

  amax = alpha - amax;
  {
#line 467
 {
    double b = amax + n.x;
    if (b > 0.) {
      area += b*b;
      p->y += b*b*(2.*n.y + alpha - n.z);
      p->z += b*b*(2.*n.z + alpha - n.y);
      p->x += b*b*b;
    }
  }
#line 467
 {
    double b = amax + n.y;
    if (b > 0.) {
      area += b*b;
      p->z += b*b*(2.*n.z + alpha - n.x);
      p->x += b*b*(2.*n.x + alpha - n.z);
      p->y += b*b*b;
    }
  }
#line 467
 {
    double b = amax + n.z;
    if (b > 0.) {
      area += b*b;
      p->x += b*b*(2.*n.x + alpha - n.y);
      p->y += b*b*(2.*n.y + alpha - n.x);
      p->z += b*b*b;
    }
  }}

  area *= 3.;
  {
#line 478
 {
    if (area) {
      p->x /= area*n.x;
      p->x = clamp (p->x, 0., 1.);
    }
    else
      p->x = 0.;
    if (m.x < 0.) p->x = 1. - p->x;
    p->x -= 0.5;
  }
#line 478
 {
    if (area) {
      p->y /= area*n.y;
      p->y = clamp (p->y, 0., 1.);
    }
    else
      p->y = 0.;
    if (m.y < 0.) p->y = 1. - p->y;
    p->y -= 0.5;
  }
#line 478
 {
    if (area) {
      p->z /= area*n.z;
      p->z = clamp (p->z, 0., 1.);
    }
    else
      p->z = 0.;
    if (m.z < 0.) p->z = 1. - p->z;
    p->z -= 0.5;
  }}

  return area*sqrt (1./(sq(n.x)*sq(n.y)) +
      1./(sq(n.x)*sq(n.z)) +
      1./(sq(n.z)*sq(n.y)))/6.;
}
#line 12 "./fractions.h"
#line 20 "./fractions.h"
#line 1 "myc.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/myc.h"
#line 16 "/home/miro3/Documents/Programming/basilisk/src/myc.h"
coord mycs (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 17 "/home/miro3/Documents/Programming/basilisk/src/myc.h"

  double m1,m2,m[4][3],t0,t1,t2;
  int cn;



  m1 = val(c,-1,0,-1) + val(c,-1,0,1) + val(c,-1,-1,0) + val(c,-1,1,0) +
       val(c,-1,0,0);
  m2 = val(c,1,0,-1) + val(c,1,0,1) + val(c,1,-1,0) + val(c,1,1,0) +
       val(c,1,0,0);
  m[0][0] = m1 > m2 ? 1. : -1.;

  m1 = val(c,-1,-1,0)+ val(c,1,-1,0)+ val(c,0,-1,0);
  m2 = val(c,-1,1,0)+ val(c,1,1,0)+ val(c,0,1,0);
  m[0][1] = 0.5*(m1-m2);

  m1 = val(c,-1,0,-1)+ val(c,1,0,-1)+ val(c,0,0,-1);
  m2 = val(c,-1,0,1)+ val(c,1,0,1)+ val(c,0,0,1);
  m[0][2] = 0.5*(m1-m2);



  m1 = val(c,-1,-1,0) + val(c,-1,1,0) + val(c,-1,0,0);
  m2 = val(c,1,-1,0) + val(c,1,1,0) + val(c,1,0,0);
  m[1][0] = 0.5*(m1-m2);

  m1 = val(c,0,-1,-1) + val(c,0,-1,1) + val(c,1,-1,0) + val(c,-1,-1,0) +
       val(c,0,-1,0);
  m2 = val(c,0,1,-1) + val(c,0,1,1) + val(c,1,1,0) + val(c,-1,1,0) +
       val(c,0,1,0);
  m[1][1] = m1 > m2 ? 1. : -1.;

  m1 = val(c,0,-1,-1)+ val(c,0,0,-1)+ val(c,0,1,-1);
  m2 = val(c,0,-1,1)+ val(c,0,0,1)+ val(c,0,1,1);
  m[1][2] = 0.5*(m1-m2);




  m1 = val(c,-1,0,-1)+ val(c,-1,0,1)+ val(c,-1,0,0);
  m2 = val(c,1,0,-1)+ val(c,1,0,1)+ val(c,1,0,0);
  m[2][0] = 0.5*(m1-m2);

  m1 = val(c,0,-1,-1)+ val(c,0,-1,1)+ val(c,0,-1,0);
  m2 = val(c,0,1,-1)+ val(c,0,1,1)+ val(c,0,1,0);
  m[2][1] = 0.5*(m1-m2);

  m1 = val(c,-1,0,-1) + val(c,1,0,-1) + val(c,0,-1,-1) + val(c,0,1,-1) +
       val(c,0,0,-1);
  m2 = val(c,-1,0,1) + val(c,1,0,1) + val(c,0,-1,1) + val(c,0,1,1) +
       val(c,0,0,1);
  m[2][2] = m1 > m2 ? 1. : -1.;


  t0 = fabs(m[0][0]) + fabs(m[0][1]) + fabs(m[0][2]);
  m[0][0] /= t0;
  m[0][1] /= t0;
  m[0][2] /= t0;

  t0 = fabs(m[1][0]) + fabs(m[1][1]) + fabs(m[1][2]);
  m[1][0] /= t0;
  m[1][1] /= t0;
  m[1][2] /= t0;

  t0 = fabs(m[2][0]) + fabs(m[2][1]) + fabs(m[2][2]);
  m[2][0] /= t0;
  m[2][1] /= t0;
  m[2][2] /= t0;


  t0 = fabs(m[0][0]);
  t1 = fabs(m[1][1]);
  t2 = fabs(m[2][2]);

  cn = 0;
  if (t1 > t0) {
    t0 = t1;
    cn = 1;
  }
  if (t2 > t0)
    cn = 2;


  m1 = val(c,-1,-1,-1) + val(c,-1,1,-1) + val(c,-1,-1,1) + val(c,-1,1,1) +
       2.*(val(c,-1,-1,0) + val(c,-1,1,0) + val(c,-1,0,-1) + val(c,-1,0,1)) +
       4.*val(c,-1,0,0);
  m2 = val(c,1,-1,-1) + val(c,1,1,-1) + val(c,1,-1,1) + val(c,1,1,1) +
       2.*(val(c,1,-1,0) + val(c,1,1,0) + val(c,1,0,-1) + val(c,1,0,1)) +
       4.*val(c,1,0,0);
  m[3][0] = m1 - m2;

  m1 = val(c,-1,-1,-1) + val(c,-1,-1,1) + val(c,1,-1,-1) + val(c,1,-1,1) +
       2.*( val(c,-1,-1,0) + val(c,1,-1,0) + val(c,0,-1,-1) + val(c,0,-1,1)) +
       4.*val(c,0,-1,0);
  m2 = val(c,-1,1,-1) + val(c,-1,1,1) + val(c,1,1,-1) + val(c,1,1,1) +
       2.*(val(c,-1,1,0) + val(c,1,1,0) + val(c,0,1,-1) + val(c,0,1,1)) +
       4.*val(c,0,1,0);
  m[3][1] = m1 - m2;

  m1 = val(c,-1,-1,-1) + val(c,-1,1,-1) + val(c,1,-1,-1) + val(c,1,1,-1) +
       2.*(val(c,-1,0,-1) + val(c,1,0,-1) + val(c,0,-1,-1) + val(c,0,1,-1)) +
       4.*val(c,0,0,-1);
  m2 = val(c,-1,-1,1) + val(c,-1,1,1) + val(c,1,-1,1) + val(c,1,1,1) +
       2.*(val(c,-1,0,1) + val(c,1,0,1) + val(c,0,-1,1) + val(c,0,1,1)) +
       4.*val(c,0,0,1);
  m[3][2] = m1 - m2;


  t0 = fabs(m[3][0]) + fabs(m[3][1]) + fabs(m[3][2]);
  if (t0 < 1e-30) {
    coord mxyz = {1., 0., 0.};
    return mxyz;
  }

  m[3][0] /= t0;
  m[3][1] /= t0;
  m[3][2] /= t0;


  t0 = fabs (m[3][0]);
  t1 = fabs (m[3][1]);
  t2 = fabs (m[3][2]);
  if (t1 > t0)
    t0 = t1;
  if (t2 > t0)
    t0 = t2;

  if (fabs(m[cn][cn]) > t0)
    cn = 3;


  coord mxyz = {m[cn][0], m[cn][1], m[cn][2]};
  return mxyz;
}
#line 21 "./fractions.h"
#line 32 "./fractions.h"
void fraction_refine (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 33 "./fractions.h"






  double cc = val(c,0,0,0);
  if (cc <= 0. || cc >= 1.)
     { foreach_child()
      val(c,0,0,0) = cc; end_foreach_child(); }
  else {




    coord n = mycs (point, c);
    double alpha = plane_alpha (cc, n);






    coord a, b;
    {
#line 57
 {
      a.x = 0.; b.x = 0.5;
    }
#line 57
 {
      a.y = 0.; b.y = 0.5;
    }
#line 57
 {
      a.z = 0.; b.z = 0.5;
    }}

     { foreach_child() {
      coord nc;
      {
#line 63

 nc.x = child.x*n.x;
#line 63

 nc.y = child.y*n.y;
#line 63

 nc.z = child.z*n.z;}
      val(c,0,0,0) = rectangle_fraction (nc, alpha, a, b);
    } end_foreach_child(); }
  }
}











static void alpha_refine (Point point, scalar alpha)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 81 "./fractions.h"

  vector n = _attribute[alpha.i].n;
  double alphac = 2.*val(alpha,0,0,0);
  coord m;
  {
#line 85

    m.x = val(n.x,0,0,0);
#line 85

    m.y = val(n.y,0,0,0);
#line 85

    m.z = val(n.z,0,0,0);}
   { foreach_child() {
    val(alpha,0,0,0) = alphac;
    {
#line 89

      val(alpha,0,0,0) -= child.x*m.x/2.;
#line 89

      val(alpha,0,0,0) -= child.y*m.y/2.;
#line 89

      val(alpha,0,0,0) -= child.z*m.z/2.;}
  } end_foreach_child(); }
}
#line 116 "./fractions.h"
struct Fractions {
  scalar Phi;
  scalar c;
  vector s;
};


void fractions (struct Fractions a)
{ trace ("fractions", "./fractions.h", 124);
  scalar Phi = a.Phi;
  scalar c = a.c;
  vector s = (a.s).x.i ? (a.s) : new_face_vector("s");
#line 136 "./fractions.h"
  vector p= new_vector("p");
#line 148 "./fractions.h"
   { foreach_vertex(){

#line 148 "./fractions.h"
 {
#line 148
 if (neighbor(1,0,0).flags & vertex) {





    if (val(Phi,0,0,0)*val(Phi,1,0,0) < 0.) {






      val(p.x,0,0,0) = val(Phi,0,0,0)/(val(Phi,0,0,0) - val(Phi,1,0,0));
      if (val(Phi,0,0,0) < 0.)
 val(p.x,0,0,0) = 1. - val(p.x,0,0,0);
    }
#line 173 "./fractions.h"
    else
      val(p.x,0,0,0) = (val(Phi,0,0,0) > 0. || val(Phi,1,0,0) > 0.);
  }
#line 148
 if (neighbor(0,1,0).flags & vertex) {





    if (val(Phi,0,0,0)*val(Phi,0,1,0) < 0.) {






      val(p.y,0,0,0) = val(Phi,0,0,0)/(val(Phi,0,0,0) - val(Phi,0,1,0));
      if (val(Phi,0,0,0) < 0.)
 val(p.y,0,0,0) = 1. - val(p.y,0,0,0);
    }
#line 173 "./fractions.h"
    else
      val(p.y,0,0,0) = (val(Phi,0,0,0) > 0. || val(Phi,0,1,0) > 0.);
  }
#line 148
 if (neighbor(0,0,1).flags & vertex) {





    if (val(Phi,0,0,0)*val(Phi,0,0,1) < 0.) {






      val(p.z,0,0,0) = val(Phi,0,0,0)/(val(Phi,0,0,0) - val(Phi,0,0,1));
      if (val(Phi,0,0,0) < 0.)
 val(p.z,0,0,0) = 1. - val(p.z,0,0,0);
    }
#line 173 "./fractions.h"
    else
      val(p.z,0,0,0) = (val(Phi,0,0,0) > 0. || val(Phi,0,0,1) > 0.);
  }} } end_foreach_vertex(); }
#line 187 "./fractions.h"
  scalar s_x = s.x, s_y = s.y, s_z = s.z;
   { foreach_face_generic() { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 188
{

#line 188 "./fractions.h"






  {
#line 226 "./fractions.h"
    coord n;
    double nn = 0.;
    {
#line 228
 {
      n.x = val(p.y,0,0,0) - val(p.y,1,0,0);
      nn += fabs(n.x);
    }
#line 228
 {
      n.y = val(p.x,0,0,0) - val(p.x,0,1,0);
      nn += fabs(n.y);
    }}





    if (nn == 0.)
      val(s_z,0,0,0) = val(p.x,0,0,0);
    else {





      {
#line 245

 n.x /= nn;
#line 245

 n.y /= nn;}






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 {
#line 255

   if (val(p.x,0,i,0) > 0. && val(p.x,0,i,0) < 1.) {
     double a = sign(val(Phi,0,i,0))*(val(p.x,0,i,0) - 0.5);
     alpha += n.x*a + n.y*(i - 0.5);
     ni++;
   }
#line 255

   if (val(p.y,i,0,0) > 0. && val(p.y,i,0,0) < 1.) {
     double a = sign(val(Phi,i,0,0))*(val(p.y,i,0,0) - 0.5);
     alpha += n.y*a + n.x*(i - 0.5);
     ni++;
   }}







      val(s_z,0,0,0) = ni ? line_area (n.x, n.y, alpha/ni) : val(p.x,0,0,0);
    }
  } }  }}  { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 188
{

#line 188 "./fractions.h"






  {
#line 226 "./fractions.h"
    coord n;
    double nn = 0.;
    {
#line 228
 {
      n.y = val(p.z,0,0,0) - val(p.z,0,1,0);
      nn += fabs(n.y);
    }
#line 228
 {
      n.z = val(p.y,0,0,0) - val(p.y,0,0,1);
      nn += fabs(n.z);
    }}





    if (nn == 0.)
      val(s_x,0,0,0) = val(p.y,0,0,0);
    else {





      {
#line 245

 n.y /= nn;
#line 245

 n.z /= nn;}






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 {
#line 255

   if (val(p.y,0,0,i) > 0. && val(p.y,0,0,i) < 1.) {
     double a = sign(val(Phi,0,0,i))*(val(p.y,0,0,i) - 0.5);
     alpha += n.y*a + n.z*(i - 0.5);
     ni++;
   }
#line 255

   if (val(p.z,0,i,0) > 0. && val(p.z,0,i,0) < 1.) {
     double a = sign(val(Phi,0,i,0))*(val(p.z,0,i,0) - 0.5);
     alpha += n.z*a + n.y*(i - 0.5);
     ni++;
   }}







      val(s_x,0,0,0) = ni ? line_area (n.y, n.z, alpha/ni) : val(p.y,0,0,0);
    }
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 188
{

#line 188 "./fractions.h"






  {
#line 226 "./fractions.h"
    coord n;
    double nn = 0.;
    {
#line 228
 {
      n.z = val(p.x,0,0,0) - val(p.x,0,0,1);
      nn += fabs(n.z);
    }
#line 228
 {
      n.x = val(p.z,0,0,0) - val(p.z,1,0,0);
      nn += fabs(n.x);
    }}





    if (nn == 0.)
      val(s_y,0,0,0) = val(p.z,0,0,0);
    else {





      {
#line 245

 n.z /= nn;
#line 245

 n.x /= nn;}






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 {
#line 255

   if (val(p.z,i,0,0) > 0. && val(p.z,i,0,0) < 1.) {
     double a = sign(val(Phi,i,0,0))*(val(p.z,i,0,0) - 0.5);
     alpha += n.z*a + n.x*(i - 0.5);
     ni++;
   }
#line 255

   if (val(p.x,0,0,i) > 0. && val(p.x,0,0,i) < 1.) {
     double a = sign(val(Phi,0,0,i))*(val(p.x,0,0,i) - 0.5);
     alpha += n.x*a + n.z*(i - 0.5);
     ni++;
   }}







      val(s_y,0,0,0) = ni ? line_area (n.z, n.x, alpha/ni) : val(p.z,0,0,0);
    }
  } }  }}  end_foreach_face_generic()
#line 270
 end_foreach_face(); }







  boundary_flux (((vector []){{s.x,s.y,s.z},{{-1},{-1},{-1}}}));
   { foreach(){

#line 279 "./fractions.h"
 {




    coord n;
    double nn = 0.;
    {
#line 286
 {
      n.x = val(s.x,0,0,0) - val(s.x,1,0,0);
      nn += fabs(n.x);
    }
#line 286
 {
      n.y = val(s.y,0,0,0) - val(s.y,0,1,0);
      nn += fabs(n.y);
    }
#line 286
 {
      n.z = val(s.z,0,0,0) - val(s.z,0,0,1);
      nn += fabs(n.z);
    }}
    if (nn == 0.)
      val(c,0,0,0) = val(s.x,0,0,0);
    else {
      {
#line 293

 n.x /= nn;
#line 293

 n.y /= nn;
#line 293

 n.z /= nn;}






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 for (int j = 0; j <= 1; j++)
   {
#line 304

     if (val(p.x,0,i,j) > 0. && val(p.x,0,i,j) < 1.) {
       double a = sign(val(Phi,0,i,j))*(val(p.x,0,i,j) - 0.5);
       alpha += n.x*a + n.y*(i - 0.5) + n.z*(j - 0.5);
       ni++;
     }
#line 304

     if (val(p.y,j,0,i) > 0. && val(p.y,j,0,i) < 1.) {
       double a = sign(val(Phi,j,0,i))*(val(p.y,j,0,i) - 0.5);
       alpha += n.y*a + n.z*(i - 0.5) + n.x*(j - 0.5);
       ni++;
     }
#line 304

     if (val(p.z,i,j,0) > 0. && val(p.z,i,j,0) < 1.) {
       double a = sign(val(Phi,i,j,0))*(val(p.z,i,j,0) - 0.5);
       alpha += n.z*a + n.x*(i - 0.5) + n.y*(j - 0.5);
       ni++;
     }}




      val(c,0,0,0) = ni ? plane_volume (n, alpha/ni) : val(s.x,0,0,0);
    }
  } } end_foreach(); }





  boundary (((scalar []){c,{-1}}));
 delete (((scalar []){p.x,p.y,p.z,{-1}}));  { if (!(a.s).x.i) delete (((scalar []){s.x,s.y,s.z,{-1}})); }  end_trace("fractions", "./fractions.h", 323); }
#line 348 "./fractions.h"
coord youngs_normal (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 349 "./fractions.h"

  coord n;
  double nn = 0.;
  assert (3 == 2);
  {
#line 353
 {
    n.x = (val(c,-1,1,0) + 2.*val(c,-1,0,0) + val(c,-1,-1,0) -
    val(c,+1,1,0) - 2.*val(c,+1,0,0) - val(c,+1,-1,0));
    nn += fabs(n.x);
  }
#line 353
 {
    n.y = (val(c,0,-1,1) + 2.*val(c,0,-1,0) + val(c,0,-1,-1) -
    val(c,0,+1,1) - 2.*val(c,0,+1,0) - val(c,0,+1,-1));
    nn += fabs(n.y);
  }
#line 353
 {
    n.z = (val(c,1,0,-1) + 2.*val(c,0,0,-1) + val(c,-1,0,-1) -
    val(c,1,0,+1) - 2.*val(c,0,0,+1) - val(c,-1,0,+1));
    nn += fabs(n.z);
  }}

  if (nn > 0.)
    {
#line 360

      n.x /= nn;
#line 360

      n.y /= nn;
#line 360

      n.z /= nn;}
  else
    n.x = 1.;
  return n;
}
#line 374 "./fractions.h"

void reconstruction (const scalar c, vector n, scalar alpha)
{ trace ("reconstruction", "./fractions.h", 376);
   { foreach(){

#line 377 "./fractions.h"
 {





    if (val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) {
      val(alpha,0,0,0) = 0.;
      {
#line 385

 val(n.x,0,0,0) = 0.;
#line 385

 val(n.y,0,0,0) = 0.;
#line 385

 val(n.z,0,0,0) = 0.;}
    }
    else {






      coord m = mycs (point, c);

      {
#line 397

 val(n.x,0,0,0) = m.x;
#line 397

 val(n.y,0,0,0) = m.y;
#line 397

 val(n.z,0,0,0) = m.z;}
      val(alpha,0,0,0) = plane_alpha (val(c,0,0,0), m);
    }
  } } end_foreach(); }
#line 410 "./fractions.h"
  {
#line 410

    _attribute[n.x.i].refine = _attribute[n.x.i].prolongation = refine_injection;
#line 410

    _attribute[n.y.i].refine = _attribute[n.y.i].prolongation = refine_injection;
#line 410

    _attribute[n.z.i].refine = _attribute[n.z.i].prolongation = refine_injection;}




  _attribute[alpha.i].n = n;
  _attribute[alpha.i].refine = _attribute[alpha.i].prolongation = alpha_refine;







  boundary (((scalar []){n.x,n.y,n.z,alpha,{-1}}));
 end_trace("reconstruction", "./fractions.h", 426); }
#line 450 "./fractions.h"
struct OutputFacets {
  scalar c;
  FILE * fp;
  vector s;
};


void output_facets (struct OutputFacets p)
{ trace ("output_facets", "./fractions.h", 458);
  scalar c = p.c;
  vector s = p.s;
  if (!p.fp) p.fp = qstdout();

   { foreach(){

#line 463 "./fractions.h"

    if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n;
      if (!s.x.i)

 n = mycs (point, c);
      else {

 double nn = 0.;
 {
#line 472
 {
   n.x = val(s.x,0,0,0) - val(s.x,1,0,0);
   nn += fabs(n.x);
 }
#line 472
 {
   n.y = val(s.y,0,0,0) - val(s.y,0,1,0);
   nn += fabs(n.y);
 }
#line 472
 {
   n.z = val(s.z,0,0,0) - val(s.z,0,0,1);
   nn += fabs(n.z);
 }}
 assert (nn > 0.);
 {
#line 477

   n.x /= nn;
#line 477

   n.y /= nn;
#line 477

   n.z /= nn;}
      }
      double alpha = plane_alpha (val(c,0,0,0), n);







      coord v[12];
      int m = facets (n, alpha, v, 1.);
      for (int i = 0; i < m; i++)
 fprintf (p.fp, "%g %g %g\n",
   x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
      if (m > 0)
 fputc ('\n', p.fp);

    } } end_foreach(); }

  fflush (p.fp);
 end_trace("output_facets", "./fractions.h", 499); }



void output_ply (struct OutputFacets p)
{
  scalar c = p.c;
  vector s = p.s;
  if (!p.fp) p.fp = qstdout();


  fputs ("ply\n", p.fp);
  fputs ("format ascii 1.0\n", p.fp);

  int nverts = 0;
  int nfacets = 0;

   { foreach(){

#line 516 "./fractions.h"

    if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n;
      if (!s.x.i)

 n = mycs (point, c);
      else {

 double nn = 0.;
 {
#line 525
 {
   n.x = val(s.x,0,0,0) - val(s.x,1,0,0);
   nn += fabs(n.x);
 }
#line 525
 {
   n.y = val(s.y,0,0,0) - val(s.y,0,1,0);
   nn += fabs(n.y);
 }
#line 525
 {
   n.z = val(s.z,0,0,0) - val(s.z,0,0,1);
   nn += fabs(n.z);
 }}
 assert (nn > 0.);
 {
#line 530

   n.x /= nn;
#line 530

   n.y /= nn;
#line 530

   n.z /= nn;}
      }
      double alpha = plane_alpha (val(c,0,0,0), n);

      coord v[12];
      int m = facets (n, alpha, v, 1.);
      for (int i = 0; i < m; i++) {
        nverts ++;
        }
      if (m > 0) {
 nfacets ++;
        }
    } } end_foreach(); }

  fprintf (p.fp, "element vertex %i\n", nverts);
  fputs ("property float x\n", p.fp);
  fputs ("property float y\n", p.fp);
  fputs ("property float z\n", p.fp);
  fprintf (p.fp, "element face %i\n", nfacets);
  fputs ("property list uchar int vertex_index\n", p.fp);
  fputs ("end_header\n", p.fp);

  int facet_num[nfacets];

  int ifacet = 0;
  int ivert = 0;

   { foreach(){

#line 558 "./fractions.h"

    if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n;
      if (!s.x.i)

 n = mycs (point, c);
      else {

 double nn = 0.;
 {
#line 567
 {
   n.x = val(s.x,0,0,0) - val(s.x,1,0,0);
   nn += fabs(n.x);
 }
#line 567
 {
   n.y = val(s.y,0,0,0) - val(s.y,0,1,0);
   nn += fabs(n.y);
 }
#line 567
 {
   n.z = val(s.z,0,0,0) - val(s.z,0,0,1);
   nn += fabs(n.z);
 }}
 assert (nn > 0.);
 {
#line 572

   n.x /= nn;
#line 572

   n.y /= nn;
#line 572

   n.z /= nn;}
      }
      double alpha = plane_alpha (val(c,0,0,0), n);

      coord v[12];
      int m = facets (n, alpha, v, 1.);
      for (int i = 0; i < m; i++) {
 fprintf (p.fp, "%g %g %g\n",
   x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
        }
      if (m > 0) {
 facet_num[ifacet] = m;
 ifacet ++;
        }
    } } end_foreach(); }


  for (ifacet = 0; ifacet < nfacets; ifacet++) {
    fprintf (p.fp, "%i ", facet_num[ifacet]);
    for (int iv = 0; iv < facet_num[ifacet]; iv ++) {
      fprintf (p.fp, "%i ", ivert);
      ivert ++;
      }
    fputc ('\n', p.fp);
    }

  fflush (p.fp);
}


void output_vtu (struct OutputFacets p)
{
  scalar c = p.c;
  vector s = p.s;
  if (!p.fp) p.fp = qstdout();


  fputs ("<?xml version=\"1.0\"?>\n", p.fp);
  fputs ("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n", p.fp);
  fputs ("  <UnstructuredGrid>\n", p.fp);

  int nverts = 0;
  int nfacets = 0;

   { foreach(){

#line 617 "./fractions.h"

    if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n;
      if (!s.x.i)

 n = mycs (point, c);
      else {

 double nn = 0.;
 {
#line 626
 {
   n.x = val(s.x,0,0,0) - val(s.x,1,0,0);
   nn += fabs(n.x);
 }
#line 626
 {
   n.y = val(s.y,0,0,0) - val(s.y,0,1,0);
   nn += fabs(n.y);
 }
#line 626
 {
   n.z = val(s.z,0,0,0) - val(s.z,0,0,1);
   nn += fabs(n.z);
 }}
 assert (nn > 0.);
 {
#line 631

   n.x /= nn;
#line 631

   n.y /= nn;
#line 631

   n.z /= nn;}
      }
      double alpha = plane_alpha (val(c,0,0,0), n);

      coord v[12];
      int m = facets (n, alpha, v, 1.);
      for (int i = 0; i < m; i++) {
        nverts ++;
        }
      if (m > 0) {
 nfacets ++;
        }
    } } end_foreach(); }

  fprintf (p.fp, "    <Piece NumberOfPoints=\"%i\" NumberOfCells=\"%i\">\n", nverts, nfacets);
  fputs ("      <Points>\n", p.fp);
  fputs ("        <DataArray type=\"Float32\" Name=\"vertices\" NumberOfComponents=\"3\" format=\"ascii\">\n", p.fp);

  int offsets[nfacets];

  int ifacet = 0;
  int offset = 0;

   { foreach(){

#line 655 "./fractions.h"

    if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n;
      if (!s.x.i)

 n = mycs (point, c);
      else {

 double nn = 0.;
 {
#line 664
 {
   n.x = val(s.x,0,0,0) - val(s.x,1,0,0);
   nn += fabs(n.x);
 }
#line 664
 {
   n.y = val(s.y,0,0,0) - val(s.y,0,1,0);
   nn += fabs(n.y);
 }
#line 664
 {
   n.z = val(s.z,0,0,0) - val(s.z,0,0,1);
   nn += fabs(n.z);
 }}
 assert (nn > 0.);
 {
#line 669

   n.x /= nn;
#line 669

   n.y /= nn;
#line 669

   n.z /= nn;}
      }
      double alpha = plane_alpha (val(c,0,0,0), n);

      coord v[12];
      int m = facets (n, alpha, v, 1.);
      for (int i = 0; i < m; i++) {
 fprintf (p.fp, "%g %g %g ",
   x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
        }
      if (m > 0) {
        offset += m;
        offsets[ifacet] = offset;
 ifacet ++;
        }
    } } end_foreach(); }


  fputs ("        </DataArray>\n", p.fp);
  fputs ("      </Points>\n", p.fp);
  fputs ("      <Cells>\n", p.fp);

  fputs ("        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n", p.fp);


  for (int ivert = 0; ivert < nverts; ivert++)
    fprintf (p.fp, "%i ", ivert);

  fputs ("        </DataArray>\n", p.fp);
  fputs ("        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n", p.fp);


  for (ifacet = 0; ifacet < nfacets; ifacet++)
    fprintf (p.fp, "%i ", offsets[ifacet]);

  fputs ("        </DataArray>\n", p.fp);
  fputs ("        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n", p.fp);


  for (ifacet = 0; ifacet < nfacets; ifacet++)
    fprintf (p.fp, "7 ");

  fputs ("        </DataArray>\n", p.fp);
  fputs ("      </Cells>\n", p.fp);
  fputs ("      <PointData>\n", p.fp);
  fputs ("      </PointData>\n", p.fp);
  fputs ("      <CellData>\n", p.fp);
  fputs ("      </CellData>\n", p.fp);
  fputs ("    </Piece>\n", p.fp);
  fputs ("  </UnstructuredGrid>\n", p.fp);
  fputs ("</VTKFile>\n", p.fp);

  fflush (p.fp);
}








double interface_area (scalar c)
{ trace ("interface_area", "./fractions.h", 733);
  double area = 0.;
   { foreach(){

#line 735 "./fractions.h"

    if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n = mycs (point, c), p;
      double alpha = plane_alpha (val(c,0,0,0), n);
      area += pow(Delta, 3 - 1)*plane_area_center (n, alpha, &p);
    } } end_foreach(); }
  { double _ret =  area; end_trace("interface_area", "./fractions.h", 741);  return _ret; }
 end_trace("interface_area", "./fractions.h", 742); }
#line 3 "milk_crown.c"

#line 1 "navier-stokes/centered.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
#line 24 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
#line 1 "./run.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/run.h"
#line 9 "/home/miro3/Documents/Programming/basilisk/src/run.h"
double dt = 1.;

#line 1 "./utils.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/utils.h"







double DT = 1e10, CFL = 0.5;




struct {

  long nc;

  long tnc;

  double t;

  double speed;

  timer gt;
} perf;





void update_perf() {
  perf.nc += grid->n;
  perf.tnc += grid->tn;
  perf.t = timer_elapsed (perf.gt);
  perf.speed = perf.tnc/perf.t;
}






typedef struct {
  double cpu;
  double real;
  double speed;
  double min;
  double avg;
  double max;
  size_t tnc;
  long mem;
} timing;






timing timer_timing (timer t, int i, size_t tnc, double * mpi)
{
  timing s;
#if 1
  s.avg = mpi_time - t.tm;
#endif
  clock_t end = clock();
  s.cpu = ((double) (end - t.c))/CLOCKS_PER_SEC;
  s.real = timer_elapsed (t);
  if (tnc == 0) {
    double n = 0;
     { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _n = n; 
#line 69
foreach(){

#line 69 "/home/miro3/Documents/Programming/basilisk/src/utils.h"
 _n++; } end_foreach();OMP(omp critical) n += _n;
mpi_all_reduce_double (n, MPI_SUM);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 69
 }
    s.tnc = n;
    tnc = n*i;
  }
  else
    s.tnc = tnc;
#if _GNU_SOURCE
  struct rusage usage;
  getrusage (RUSAGE_SELF, &usage);
  s.mem = usage.ru_maxrss;
#else
  s.mem = 0;
#endif
#if 1
  if (mpi)
    MPI_Allgather (&s.avg, 1, MPI_DOUBLE, mpi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
  s.max = s.min = s.avg;
  mpi_all_reduce (s.max, MPI_DOUBLE, MPI_MAX);
  mpi_all_reduce (s.min, MPI_DOUBLE, MPI_MIN);
  mpi_all_reduce (s.avg, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (s.real, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (s.mem, MPI_LONG, MPI_SUM);
  s.real /= npe();
  s.avg /= npe();
  s.mem /= npe();
#else
  s.min = s.max = s.avg = 0.;
#endif
  s.speed = s.real > 0. ? tnc/s.real : -1.;
  return s;
}




void timer_print (timer t, int i, size_t tnc)
{
  timing s = timer_timing (t, i, tnc, NULL);
  fprintf (fout,
    "\n# " "Octree"
    ", %d steps, %g CPU, %.4g real, %.3g points.step/s, %d var\n",
    i, s.cpu, s.real, s.speed, (int) (datasize/sizeof(double)));
#if 1
  fprintf (fout,
    "# %d procs, MPI: min %.2g (%.2g%%) "
    "avg %.2g (%.2g%%) max %.2g (%.2g%%)\n",
    npe(),
    s.min, 100.*s.min/s.real,
    s.avg, 100.*s.avg/s.real,
    s.max, 100.*s.max/s.real);
#endif
}







typedef struct {
  double avg, rms, max, volume;
} norm;

norm normf (scalar f)
{
  double avg = 0., rms = 0., max = 0., volume = 0.;
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _max = max; double _avg = avg; double _rms = rms; double _volume = volume; 
#line 135

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 135
foreach(){

#line 136 "/home/miro3/Documents/Programming/basilisk/src/utils.h"

    if (val(f,0,0,0) != nodata) {
      double v = fabs(val(f,0,0,0));
      if (v > _max) _max = v;
      _volume += (cube(Delta)*val_cm(cm,0,0,0));
      _avg += (cube(Delta)*val_cm(cm,0,0,0))*v;
      _rms += (cube(Delta)*val_cm(cm,0,0,0))*sq(v);
    } } end_foreach(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 135
foreach(){

#line 136 "/home/miro3/Documents/Programming/basilisk/src/utils.h"

    if (val(f,0,0,0) != nodata) {
      double v = fabs(val(f,0,0,0));
      if (v > _max) _max = v;
      _volume += (cube(Delta)*val_cm(cm,0,0,0));
      _avg += (cube(Delta)*val_cm(cm,0,0,0))*v;
      _rms += (cube(Delta)*val_cm(cm,0,0,0))*sq(v);
    } } end_foreach(); }OMP(omp critical) if (_max > max) max = _max;
mpi_all_reduce_double (max, MPI_MAX);
OMP(omp critical) avg += _avg;
mpi_all_reduce_double (avg, MPI_SUM);
OMP(omp critical) rms += _rms;
mpi_all_reduce_double (rms, MPI_SUM);
OMP(omp critical) volume += _volume;
mpi_all_reduce_double (volume, MPI_SUM);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 143
 }
  norm n;
  n.avg = avg/volume;
  n.rms = sqrt(rms/volume);
  n.max = max;
  n.volume = volume;
  return n;
}





typedef struct {
  double min, max, sum, stddev, volume;
} stats;

stats statsf (scalar f)
{
  double min = 1e100, max = -1e100, sum = 0., sum2 = 0., volume = 0.;
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _sum = sum; double _sum2 = sum2; double _volume = volume; double _max = max; double _min = min; 
#line 163

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 163
foreach(){

#line 164 "/home/miro3/Documents/Programming/basilisk/src/utils.h"

    if (val(f,0,0,0) != nodata) {
      _volume += (cube(Delta)*val_cm(cm,0,0,0));
      _sum += (cube(Delta)*val_cm(cm,0,0,0))*val(f,0,0,0);
      _sum2 += (cube(Delta)*val_cm(cm,0,0,0))*sq(val(f,0,0,0));
      if (val(f,0,0,0) > _max) _max = val(f,0,0,0);
      if (val(f,0,0,0) < _min) _min = val(f,0,0,0);
    } } end_foreach(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 163
foreach(){

#line 164 "/home/miro3/Documents/Programming/basilisk/src/utils.h"

    if (val(f,0,0,0) != nodata) {
      _volume += (cube(Delta)*val_cm(cm,0,0,0));
      _sum += (cube(Delta)*val_cm(cm,0,0,0))*val(f,0,0,0);
      _sum2 += (cube(Delta)*val_cm(cm,0,0,0))*sq(val(f,0,0,0));
      if (val(f,0,0,0) > _max) _max = val(f,0,0,0);
      if (val(f,0,0,0) < _min) _min = val(f,0,0,0);
    } } end_foreach(); }OMP(omp critical) sum += _sum;
mpi_all_reduce_double (sum, MPI_SUM);
OMP(omp critical) sum2 += _sum2;
mpi_all_reduce_double (sum2, MPI_SUM);
OMP(omp critical) volume += _volume;
mpi_all_reduce_double (volume, MPI_SUM);
OMP(omp critical) if (_max > max) max = _max;
mpi_all_reduce_double (max, MPI_MAX);
OMP(omp critical) if (_min < min) min = _min;
mpi_all_reduce_double (min, MPI_MIN);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 171
 }
  stats s;
  s.min = min, s.max = max, s.sum = sum, s.volume = volume;
  sum2 -= sum*sum/volume;
  s.stddev = sum2 > 0. ? sqrt(sum2/volume) : 0.;
  return s;
}
#line 186 "/home/miro3/Documents/Programming/basilisk/src/utils.h"
static double generic_limiter (double r, double beta)
{
  double v1 = min (r, beta), v2 = min (beta*r, 1.);
  v1 = max (0., v1);
  return max (v1, v2);
}

double minmod (double s0, double s1, double s2) {
  return generic_limiter ((s2 - s1)/(s1 - s0), 1.)*(s1 - s0);
}

double superbee (double s0, double s1, double s2) {
  return generic_limiter ((s2 - s1)/(s1 - s0), 2.)*(s1 - s0);
}

double sweby (double s0, double s1, double s2) {
  return generic_limiter ((s2 - s1)/(s1 - s0), 1.5)*(s1 - s0);
}
#line 212 "/home/miro3/Documents/Programming/basilisk/src/utils.h"
double theta = 1.3;

double minmod2 (double s0, double s1, double s2)
{
  if (s0 < s1 && s1 < s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 < d1) d1 = d2;
    return min(d1, d3);
  }
  if (s0 > s1 && s1 > s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 > d1) d1 = d2;
    return max(d1, d3);
  }
  return 0.;
}
#line 236 "/home/miro3/Documents/Programming/basilisk/src/utils.h"
void gradients (scalar * f, vector * g)
{
  assert (list_len(f) == vectors_len(g));
   { foreach(){

#line 239 "/home/miro3/Documents/Programming/basilisk/src/utils.h"
 {
    scalar s; vector v;
    scalar * _i0 = f; vector * _i1 = g; if (f) for (s = *f, v = *g; ((scalar *)&s)->i >= 0; s = *++_i0, v = *++_i1) {
      if (_attribute[s.i].gradient)
 {
#line 243

   val(v.x,0,0,0) = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0))/Delta;
#line 243

   val(v.y,0,0,0) = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0))/Delta;
#line 243

   val(v.z,0,0,0) = _attribute[s.i].gradient (val(s,0,0,-1), val(s,0,0,0), val(s,0,0,1))/Delta;}
      else
 {
#line 246

   val(v.x,0,0,0) = (val(s,1,0,0) - val(s,-1,0,0))/(2.*Delta);
#line 246

   val(v.y,0,0,0) = (val(s,0,1,0) - val(s,0,-1,0))/(2.*Delta);
#line 246

   val(v.z,0,0,0) = (val(s,0,0,1) - val(s,0,0,-1))/(2.*Delta);}
    }
  } } end_foreach(); }
  boundary ((scalar *) g);
}
#line 262 "/home/miro3/Documents/Programming/basilisk/src/utils.h"
void vorticity (const vector u, scalar omega)
{

     { foreach(){

#line 265 "/home/miro3/Documents/Programming/basilisk/src/utils.h"

      val(omega,0,0,0) = (val(u.y,1,0,0) - val(u.y,-1,0,0) + val(u.x,0,-1,0) - val(u.x,0,1,0))/(2.*Delta); } end_foreach(); }
    boundary (((scalar []){omega,{-1}}));

}





double change (scalar v, scalar vn)
{
  double max = 0.;
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _max = max; 
#line 278
foreach(){

#line 278 "/home/miro3/Documents/Programming/basilisk/src/utils.h"
 {
    double dv = fabs (val(v,0,0,0) - val(vn,0,0,0));
    if (dv > _max)
      _max = dv;
    val(vn,0,0,0) = val(v,0,0,0);
  } } end_foreach();OMP(omp critical) if (_max > max) max = _max;
mpi_all_reduce_double (max, MPI_MAX);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 283
 }
  return max;
}

#line 1 "./output.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/output.h"
#line 33 "/home/miro3/Documents/Programming/basilisk/src/output.h"
struct OutputField {
  scalar * list;
  FILE * fp;
  int n;
  bool linear;
};


void output_field (struct OutputField p)
{ trace ("output_field", "/home/miro3/Documents/Programming/basilisk/src/output.h", 42);
  if (!p.list) p.list = all;
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = qstdout();

  int len = list_len(p.list);
  double ** field = (double **) matrix_new (p.n, p.n, len*sizeof(double));

  double Delta = L0/p.n;
  for (int i = 0; i < p.n; i++) {
    double xp = Delta*i + X0 + Delta/2.;
    for (int j = 0; j < p.n; j++) {
      double yp = Delta*j + Y0 + Delta/2.;
      if (p.linear) {
 int k = 0;
 if (p.list) for (scalar s = *p.list, *_i71 = p.list; ((scalar *)&s)->i >= 0; s = *++_i71)
   field[i][len*j + k++] = interpolate ((struct _interpolate){s, xp, yp});
      }
      else {
 Point point = locate ((struct _locate){xp, yp});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 61 "/home/miro3/Documents/Programming/basilisk/src/output.h"

 int k = 0;
 if (p.list) for (scalar s = *p.list, *_i72 = p.list; ((scalar *)&s)->i >= 0; s = *++_i72)
   field[i][len*j + k++] = point.level >= 0 ? val(s,0,0,0) : nodata;
      }
    }
  }

  if (pid() == 0) {
#if 1
    MPI_Reduce (MPI_IN_PLACE, field[0], len*p.n*p.n, MPI_DOUBLE, MPI_MIN, 0,
  MPI_COMM_WORLD);
#endif
    fprintf (p.fp, "# 1:x 2:y");
    int i = 3;
    if (p.list) for (scalar s = *p.list, *_i73 = p.list; ((scalar *)&s)->i >= 0; s = *++_i73)
      fprintf (p.fp, " %d:%s", i++, _attribute[s.i].name);
    fputc('\n', p.fp);
    for (int i = 0; i < p.n; i++) {
      double xp = Delta*i + X0 + Delta/2.;
      for (int j = 0; j < p.n; j++) {
 double yp = Delta*j + Y0 + Delta/2.;
 fprintf (p.fp, "%g %g", xp, yp);
 int k = 0;
 if (p.list) for (scalar s = *p.list, *_i74 = p.list; ((scalar *)&s)->i >= 0; s = *++_i74)
   fprintf (p.fp, " %g", field[i][len*j + k++]);
 fputc ('\n', p.fp);
      }
      fputc ('\n', p.fp);
    }
    fflush (p.fp);
  }
#if 1
  else
    MPI_Reduce (field[0], NULL, len*p.n*p.n, MPI_DOUBLE, MPI_MIN, 0,
  MPI_COMM_WORLD);
#endif

  matrix_free (field);
 end_trace("output_field", "/home/miro3/Documents/Programming/basilisk/src/output.h", 100); }
#line 128 "/home/miro3/Documents/Programming/basilisk/src/output.h"
struct OutputMatrix {
  scalar f;
  FILE * fp;
  int n;
  bool linear;
};


void output_matrix (struct OutputMatrix p)
{ trace ("output_matrix", "/home/miro3/Documents/Programming/basilisk/src/output.h", 137);
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = qstdout();
  float fn = p.n;
  float Delta = (float) L0/fn;
  fwrite (&fn, sizeof(float), 1, p.fp);
  for (int j = 0; j < p.n; j++) {
    float yp = (float) (Delta*j + X0 + Delta/2.);
    fwrite (&yp, sizeof(float), 1, p.fp);
  }
  for (int i = 0; i < p.n; i++) {
    float xp = (float) (Delta*i + X0 + Delta/2.);
    fwrite (&xp, sizeof(float), 1, p.fp);
    for (int j = 0; j < p.n; j++) {
      float yp = (float)(Delta*j + Y0 + Delta/2.), v;
      if (p.linear)
 v = interpolate ((struct _interpolate){p.f, xp, yp});
      else {
 Point point = locate ((struct _locate){xp, yp});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 155 "/home/miro3/Documents/Programming/basilisk/src/output.h"

 assert (point.level >= 0);
 v = val(p.f,0,0,0);
      }
      fwrite (&v, sizeof(float), 1, p.fp);
    }
  }
  fflush (p.fp);
 end_trace("output_matrix", "/home/miro3/Documents/Programming/basilisk/src/output.h", 163); }
#line 172 "/home/miro3/Documents/Programming/basilisk/src/output.h"
typedef void (* colormap) (double cmap[127][3]);

void jet (double cmap[127][3])
{
  for (int i = 0; i < 127; i++) {
    cmap[i][0] =
      i <= 46 ? 0. :
      i >= 111 ? -0.03125*(i - 111) + 1. :
      i >= 78 ? 1. :
      0.03125*(i - 46);
    cmap[i][1] =
      i <= 14 || i >= 111 ? 0. :
      i >= 79 ? -0.03125*(i - 111) :
      i <= 46 ? 0.03125*(i - 14) :
      1.;
    cmap[i][2] =
      i >= 79 ? 0. :
      i >= 47 ? -0.03125*(i - 79) :
      i <= 14 ? 0.03125*(i - 14) + 1.:
      1.;
  }
}

void cool_warm (double cmap[127][3])
{






  static double basemap[33][3] = {
    {0.2298057, 0.298717966, 0.753683153},
    {0.26623388, 0.353094838, 0.801466763},
    {0.30386891, 0.406535296, 0.84495867},
    {0.342804478, 0.458757618, 0.883725899},
    {0.38301334, 0.50941904, 0.917387822},
    {0.424369608, 0.558148092, 0.945619588},
    {0.46666708, 0.604562568, 0.968154911},
    {0.509635204, 0.648280772, 0.98478814},
    {0.552953156, 0.688929332, 0.995375608},
    {0.596262162, 0.726149107, 0.999836203},
    {0.639176211, 0.759599947, 0.998151185},
    {0.681291281, 0.788964712, 0.990363227},
    {0.722193294, 0.813952739, 0.976574709},
    {0.761464949, 0.834302879, 0.956945269},
    {0.798691636, 0.849786142, 0.931688648},
    {0.833466556, 0.860207984, 0.901068838},
    {0.865395197, 0.86541021, 0.865395561},
    {0.897787179, 0.848937047, 0.820880546},
    {0.924127593, 0.827384882, 0.774508472},
    {0.944468518, 0.800927443, 0.726736146},
    {0.958852946, 0.769767752, 0.678007945},
    {0.96732803, 0.734132809, 0.628751763},
    {0.969954137, 0.694266682, 0.579375448},
    {0.966811177, 0.650421156, 0.530263762},
    {0.958003065, 0.602842431, 0.481775914},
    {0.943660866, 0.551750968, 0.434243684},
    {0.923944917, 0.49730856, 0.387970225},
    {0.89904617, 0.439559467, 0.343229596},
    {0.869186849, 0.378313092, 0.300267182},
    {0.834620542, 0.312874446, 0.259301199},
    {0.795631745, 0.24128379, 0.220525627},
    {0.752534934, 0.157246067, 0.184115123},
    {0.705673158, 0.01555616, 0.150232812}
  };

  for (int i = 0; i < 127; i++) {
    double x = i*(32 - 1e-10)/(127 - 1);
    int j = x; x -= j;
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (1. - x)*basemap[j][k] + x*basemap[j+1][k];
  }
}

void gray (double cmap[127][3])
{
  for (int i = 0; i < 127; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = i/(127 - 1.);
}

void randomap (double cmap[127][3])
{
  srand(0);
  for (int i = 0; i < 127; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (noise() + 1.)/2.;
}





typedef struct {
  unsigned char r, g, b;
} color;

color colormap_color (double cmap[127][3],
        double val, double min, double max)
{
  color c;
  if (val == nodata) {
    c.r = c.g = c.b = 0;
    return c;
  }
  val = val <= min ? 0. : val >= max ? 0.9999 : (val - min)/(max - min);
  int i = val*(127 - 1);
  double coef = val*(127 - 1) - i;
  assert (i < 127 - 1);
  unsigned char * c1 = (unsigned char *) &c;
  for (int j = 0; j < 3; j++)
    c1[j] = 255*(cmap[i][j]*(1. - coef) + cmap[i + 1][j]*coef);
  return c;
}
#line 347 "/home/miro3/Documents/Programming/basilisk/src/output.h"
struct OutputPPM {
  scalar f;
  FILE * fp;
  int n;
  char * file;
  double min, max, spread, z;
  bool linear;
  double box[2][2];
  scalar mask;
  colormap map;
};


void output_ppm (struct OutputPPM p)
{ trace ("output_ppm", "/home/miro3/Documents/Programming/basilisk/src/output.h", 361);

  if (p.n == 0) p.n = N;
  if (p.min == 0 && p.max == 0) {
    stats s = statsf (p.f);
    double avg = s.sum/s.volume, spread = (p.spread ? p.spread : 5.)*s.stddev;
    p.min = avg - spread; p.max = avg + spread;
  }
  if (p.box[0][0] == 0. && p.box[0][1] == 0. &&
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0; p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
  }
  if (!p.map)
    p.map = jet;

  double fn = p.n;
  double Delta = (p.box[1][0] - p.box[0][0])/fn;
  int ny = (p.box[1][1] - p.box[0][1])/Delta;

  color ** ppm = (color **) matrix_new (ny, p.n, sizeof(color));
  double cmap[127][3];
  p.map (cmap);
  OMP_PARALLEL() {
    OMP(omp for schedule(static))
      for (int j = 0; j < ny; j++) {
 double yp = Delta*j + p.box[0][1] + Delta/2.;
 for (int i = 0; i < p.n; i++) {
   double xp = Delta*i + p.box[0][0] + Delta/2., v;
   if (p.mask.i) {
     if (p.linear) {
       double m = interpolate ((struct _interpolate){p.mask, xp, yp, p.z});
       if (m < 0.)
  v = nodata;
       else
  v = interpolate ((struct _interpolate){p.f, xp, yp, p.z});
     }
     else {
       Point point = locate ((struct _locate){xp, yp, p.z});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 399 "/home/miro3/Documents/Programming/basilisk/src/output.h"

       if (point.level < 0 || val(p.mask,0,0,0) < 0.)
  v = nodata;
       else
  v = val(p.f,0,0,0);
     }
   }
   else if (p.linear)
     v = interpolate ((struct _interpolate){p.f, xp, yp, p.z});
   else {
     Point point = locate ((struct _locate){xp, yp, p.z});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 409 "/home/miro3/Documents/Programming/basilisk/src/output.h"

     v = point.level >= 0 ? val(p.f,0,0,0) : nodata;
   }
   ppm[ny - 1 - j][i] = colormap_color (cmap, v, p.min, p.max);
 }
      }
  }

  if (pid() == 0) {
#if 1
    MPI_Reduce (MPI_IN_PLACE, ppm[0], 3*ny*p.n, MPI_UNSIGNED_CHAR, MPI_MAX, 0,
  MPI_COMM_WORLD);
#endif
    if (!p.fp) p.fp = qstdout();
    if (p.file) {
      char command[strlen ("convert ppm:- -transparent black ") +
     strlen (p.file) + 1];
      strcpy (command, "convert ppm:- -transparent black ");
      strcat (command, p.file);
      p.fp = qpopen (command, "w");
    }

    fprintf (p.fp, "P6\n%u %u 255\n", p.n, ny);
    fwrite (((void **) ppm)[0], sizeof(color), ny*p.n, p.fp);

    if (p.file)
      qpclose (p.fp);
    else
      fflush (p.fp);
  }
#if 1
  else
    MPI_Reduce (ppm[0], NULL, 3*ny*p.n, MPI_UNSIGNED_CHAR, MPI_MAX, 0,
  MPI_COMM_WORLD);
#endif

  matrix_free (ppm);
 end_trace("output_ppm", "/home/miro3/Documents/Programming/basilisk/src/output.h", 446); }
#line 487 "/home/miro3/Documents/Programming/basilisk/src/output.h"
struct OutputGRD {
  scalar f;
  FILE * fp;
  double Delta;
  bool linear;
  double box[2][2];
  scalar mask;
};


void output_grd (struct OutputGRD p)
{ trace ("output_grd", "/home/miro3/Documents/Programming/basilisk/src/output.h", 498);

  if (!p.fp) p.fp = qstdout();
  if (p.box[0][0] == 0. && p.box[0][1] == 0. &&
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0; p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
    if (p.Delta == 0) p.Delta = L0/N;
  }

  double Delta = p.Delta;
  int nx = (p.box[1][0] - p.box[0][0])/Delta;
  int ny = (p.box[1][1] - p.box[0][1])/Delta;


  fprintf (p.fp, "ncols          %d\n", nx);
  fprintf (p.fp, "nrows          %d\n", ny);
  fprintf (p.fp, "xllcorner      %g\n", p.box[0][0]);
  fprintf (p.fp, "yllcorner      %g\n", p.box[0][1]);
  fprintf (p.fp, "cellsize       %g\n", Delta);
  fprintf (p.fp, "nodata_value   -9999\n");


  for (int j = ny-1; j >= 0; j--) {
    double yp = Delta*j + p.box[0][1] + Delta/2.;
    for (int i = 0; i < nx; i++) {
      double xp = Delta*i + p.box[0][0] + Delta/2., v;
      if (p.mask.i) {
 if (p.linear) {
   double m = interpolate ((struct _interpolate){p.mask, xp, yp});
   if (m < 0.)
     v = nodata;
   else
     v = interpolate ((struct _interpolate){p.f, xp, yp});
 }
 else {
   Point point = locate ((struct _locate){xp, yp});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 534 "/home/miro3/Documents/Programming/basilisk/src/output.h"

   if (point.level < 0 || val(p.mask,0,0,0) < 0.)
     v = nodata;
   else
     v = val(p.f,0,0,0);
 }
      }
      else if (p.linear)
 v = interpolate ((struct _interpolate){p.f, xp, yp});
      else {
 Point point = locate ((struct _locate){xp, yp});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 544 "/home/miro3/Documents/Programming/basilisk/src/output.h"

 v = point.level >= 0 ? val(p.f,0,0,0) : nodata;
      }
      if (v == nodata)
 fprintf (p.fp, "-9999 ");
      else
 fprintf (p.fp, "%f ", v);
    }
    fprintf (p.fp, "\n");
  }

  fflush (p.fp);
 end_trace("output_grd", "/home/miro3/Documents/Programming/basilisk/src/output.h", 556); }
#line 585 "/home/miro3/Documents/Programming/basilisk/src/output.h"
struct OutputGfs {
  FILE * fp;
  scalar * list;
  double t;
  char * file;
  bool translate;
};

static char * replace (const char * input, int target, int with,
         bool translate)
{
  if (translate) {
    if (!strcmp (input, "u.x"))
      return pstrdup ("U",__func__,__FILE__,__LINE__);
    if (!strcmp (input, "u.y"))
      return pstrdup ("V",__func__,__FILE__,__LINE__);
    if (!strcmp (input, "u.z"))
      return pstrdup ("W",__func__,__FILE__,__LINE__);
  }
  char * name = pstrdup (input,__func__,__FILE__,__LINE__), * i = name;
  while (*i != '\0') {
    if (*i == target)
      *i = with;
    i++;
  }
  return name;
}


void output_gfs (struct OutputGfs p)
{ trace ("output_gfs", "/home/miro3/Documents/Programming/basilisk/src/output.h", 615);
  char * fname = p.file;

#if 1
  FILE * fp = p.fp;
  if (p.file == NULL) {
    long pid = getpid();
    MPI_Bcast (&pid, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    fname = ((char *) pmalloc ((80)*sizeof(char),__func__,__FILE__,__LINE__));
    snprintf (fname, 80, ".output-%ld", pid);
    p.fp = NULL;
  }
#endif

  bool opened = false;
  if (p.fp == NULL) {
    if (fname == NULL)
      p.fp = qstdout();
    else if (!(p.fp = fopen (fname, "w"))) {
      perror (fname);
      exit (1);
    }
    else
      opened = true;
  }

  scalar * list = p.list ? p.list : list_copy (all);

  fprintf (p.fp,
    "1 0 GfsSimulation GfsBox GfsGEdge { binary = 1"
    " x = %g y = %g ",
    0.5 + X0/L0, 0.5 + Y0/L0);

  fprintf (p.fp, "z = %g ", 0.5 + Z0/L0);


  if (list != NULL && list[0].i != -1) {
    scalar s = list[0];
    char * name = replace (_attribute[s.i].name, '.', '_', p.translate);
    fprintf (p.fp, "variables = %s", name);
    pfree (name,__func__,__FILE__,__LINE__);
    for (int i = 1; i < list_len(list); i++) {
      scalar s = list[i];
      if (_attribute[s.i].name) {
 char * name = replace (_attribute[s.i].name, '.', '_', p.translate);
 fprintf (p.fp, ",%s", name);
 pfree (name,__func__,__FILE__,__LINE__);
      }
    }
    fprintf (p.fp, " ");
  }
  fprintf (p.fp, "} {\n");
  fprintf (p.fp, "  Time { t = %g }\n", t);
  if (L0 != 1.)
    fprintf (p.fp, "  PhysicalParams { L = %g }\n", L0);
  fprintf (p.fp, "  VariableTracerVOF f\n");
  fprintf (p.fp, "}\nGfsBox { x = 0 y = 0 z = 0 } {\n");

#if 1
  long header;
  if ((header = ftell (p.fp)) < 0) {
    perror ("output_gfs(): error in header");
    exit (1);
  }
  int cell_size = sizeof(unsigned) + sizeof(double);
  if (list) for (scalar s = *list, *_i75 = list; ((scalar *)&s)->i >= 0; s = *++_i75)
    if (_attribute[s.i].name)
      cell_size += sizeof(double);
  scalar index = new_scalar("index");
  size_t total_size = header + (z_indexing (index, false) + 1)*cell_size;
#endif



   { foreach_cell(){

#line 689 "/home/miro3/Documents/Programming/basilisk/src/output.h"
 {
#if 1
    if (is_local(cell))
#endif
    {
#if 1
      if (fseek (p.fp, header + val(index,0,0,0)*cell_size, SEEK_SET) < 0) {
 perror ("output_gfs(): error while seeking");
 exit (1);
      }
#endif
      unsigned flags =
 level == 0 ? 0 :
#line 710 "/home/miro3/Documents/Programming/basilisk/src/output.h"
      child.x == -1 && child.y == -1 && child.z == -1 ? 0 :
 child.x == -1 && child.y == -1 && child.z == 1 ? 1 :
 child.x == -1 && child.y == 1 && child.z == -1 ? 2 :
 child.x == -1 && child.y == 1 && child.z == 1 ? 3 :
 child.x == 1 && child.y == -1 && child.z == -1 ? 4 :
 child.x == 1 && child.y == -1 && child.z == 1 ? 5 :
 child.x == 1 && child.y == 1 && child.z == -1 ? 6 :
 7;

      if (is_leaf(cell))
 flags |= (1 << 4);
      fwrite (&flags, sizeof (unsigned), 1, p.fp);
      double a = -1;
      fwrite (&a, sizeof (double), 1, p.fp);
      if (list) for (scalar s = *list, *_i76 = list; ((scalar *)&s)->i >= 0; s = *++_i76)
 if (_attribute[s.i].name) {
   if (_attribute[s.i].v.x.i >= 0) {




     if (_attribute[s.i].v.x.i == s.i) {
       s = _attribute[s.i].v.y;
       a = is_local(cell) && !isnan(val(s,0,0,0)) && val(s,0,0,0) != nodata ? val(s,0,0,0) :
  (double) DBL_MAX;
     }
     else if (_attribute[s.i].v.y.i == s.i) {
       s = _attribute[s.i].v.x;
       a = is_local(cell) && !isnan(val(s,0,0,0)) && val(s,0,0,0) != nodata ?
  - val(s,0,0,0) : (double) DBL_MAX;
     }


     else
       a = is_local(cell) && !isnan(val(s,0,0,0)) && val(s,0,0,0) != nodata ? val(s,0,0,0) :
  (double) DBL_MAX;

   }
   else
     a = is_local(cell) && !isnan(val(s,0,0,0)) && val(s,0,0,0) != nodata ? val(s,0,0,0) :
       (double) DBL_MAX;
   fwrite (&a, sizeof (double), 1, p.fp);
 }
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

#if 1
  delete (((scalar []){index,{-1}}));
  if (!pid() && fseek (p.fp, total_size, SEEK_SET) < 0) {
    perror ("output_gfs(): error while finishing");
    exit (1);
  }
  if (!pid())
#endif
    fputs ("}\n", p.fp);
  fflush (p.fp);

  if (!p.list)
    pfree (list,__func__,__FILE__,__LINE__);
  if (opened)
    fclose (p.fp);

#if 1
  if (p.file == NULL) {
    MPI_Barrier (MPI_COMM_WORLD);
    if (pid() == 0) {
      if (fp == NULL)
 fp = qstdout();
      p.fp = fopen (fname, "r");
      size_t l;
      unsigned char buffer[8192];
      while ((l = fread (buffer, 1, 8192, p.fp)) > 0)
 fwrite (buffer, 1, l, fp);
      fflush (fp);
      remove (fname);
    }
    pfree (fname,__func__,__FILE__,__LINE__);
  }
#endif
 end_trace("output_gfs", "/home/miro3/Documents/Programming/basilisk/src/output.h", 791); }
#line 813 "/home/miro3/Documents/Programming/basilisk/src/output.h"
struct Dump {
  char * file;
  scalar * list;
  FILE * fp;
};

struct DumpHeader {
  double t;
  long len;
  int i, depth, npe, version;
  coord n;
};

static const int dump_version =

  170901;

static scalar * dump_list (scalar * lista)
{
  scalar * list = is_constant(cm) ? NULL : list_concat (((scalar []){cm,{-1}}), NULL);
  if (lista) for (scalar s = *lista, *_i77 = lista; ((scalar *)&s)->i >= 0; s = *++_i77)
    if (!_attribute[s.i].face && !_attribute[s.i].nodump && s.i != cm.i)
      list = list_add (list, s);
  return list;
}

static void dump_header (FILE * fp, struct DumpHeader * header, scalar * list)
{
  if (fwrite (header, sizeof(struct DumpHeader), 1, fp) < 1) {
    perror ("dump(): error while writing header");
    exit (1);
  }
  if (list) for (scalar s = *list, *_i78 = list; ((scalar *)&s)->i >= 0; s = *++_i78) {
    unsigned len = strlen(_attribute[s.i].name);
    if (fwrite (&len, sizeof(unsigned), 1, fp) < 1) {
      perror ("dump(): error while writing len");
      exit (1);
    }
    if (fwrite (_attribute[s.i].name, sizeof(char), len, fp) < len) {
      perror ("dump(): error while writing s.name");
      exit (1);
    }
  }
  double o[4] = {X0,Y0,Z0,L0};
  if (fwrite (o, sizeof(double), 4, fp) < 4) {
    perror ("dump(): error while writing coordinates");
    exit (1);
  }
}

#if !1

void dump (struct Dump p)
{ trace ("dump", "/home/miro3/Documents/Programming/basilisk/src/output.h", 866);
  FILE * fp = p.fp;
  char * file = p.file;

  if (file && (fp = fopen (file, "w")) == NULL) {
    perror (file);
    exit (1);
  }
  assert (fp);

  scalar * dlist = dump_list (p.list ? p.list : all);
  scalar size= new_scalar("size");
  scalar * list = list_concat (((scalar []){size,{-1}}), dlist); pfree (dlist,__func__,__FILE__,__LINE__);
  struct DumpHeader header = { t, list_len(list), iter, depth(), npe(),
          dump_version };
  dump_header (fp, &header, list);

  subtree_size (size, false);

   { foreach_cell(){

#line 885 "/home/miro3/Documents/Programming/basilisk/src/output.h"
 {
    unsigned flags = is_leaf(cell) ? leaf : 0;
    if (fwrite (&flags, sizeof(unsigned), 1, fp) < 1) {
      perror ("dump(): error while writing flags");
      exit (1);
    }
    if (list) for (scalar s = *list, *_i79 = list; ((scalar *)&s)->i >= 0; s = *++_i79)
      if (fwrite (&val(s,0,0,0), sizeof(double), 1, fp) < 1) {
 perror ("dump(): error while writing scalars");
 exit (1);
      }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

  pfree (list,__func__,__FILE__,__LINE__);
  if (file)
    fclose (fp);
 delete (((scalar []){size,{-1}}));  end_trace("dump", "/home/miro3/Documents/Programming/basilisk/src/output.h", 903); }
#else

void dump (struct Dump p)
{ trace ("dump", "/home/miro3/Documents/Programming/basilisk/src/output.h", 907);
  FILE * fp = p.fp;
  char * file = p.file;

  if (fp != NULL || file == NULL) {
    fprintf (qstderr(), "dump(): must specify a file name when using MPI\n");
    exit(1);
  }

  FILE * fh = fopen (file, "w");

  scalar * dlist = dump_list (p.list ? p.list : all);
  scalar size= new_scalar("size");
  scalar * list = list_concat (((scalar []){size,{-1}}), dlist); pfree (dlist,__func__,__FILE__,__LINE__);
  struct DumpHeader header = { t, list_len(list), iter, depth(), npe(),
          dump_version };







  if (pid() == 0)
    dump_header (fh, &header, list);

  scalar index = {-1};

  index = new_scalar("index");
  z_indexing (index, false);
  int cell_size = sizeof(unsigned) + header.len*sizeof(double);
  int sizeofheader = sizeof(header) + 4*sizeof(double);
  if (list) for (scalar s = *list, *_i80 = list; ((scalar *)&s)->i >= 0; s = *++_i80)
    sizeofheader += sizeof(unsigned) + sizeof(char)*strlen(_attribute[s.i].name);

  subtree_size (size, false);

   { foreach_cell(){

#line 944 "/home/miro3/Documents/Programming/basilisk/src/output.h"
 {

    if (is_local(cell)) {
      long offset = sizeofheader + val(index,0,0,0)*cell_size;
      fseek (fh, offset, SEEK_SET);
      unsigned flags = is_leaf(cell) ? leaf : 0;
      fwrite (&flags, 1, sizeof(unsigned), fh);
      if (list) for (scalar s = *list, *_i81 = list; ((scalar *)&s)->i >= 0; s = *++_i81)
 fwrite (&val(s,0,0,0), 1, sizeof(double), fh);
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

  delete (((scalar []){index,{-1}}));

  pfree (list,__func__,__FILE__,__LINE__);
  fclose (fh);
 delete (((scalar []){size,{-1}}));  end_trace("dump", "/home/miro3/Documents/Programming/basilisk/src/output.h", 962); }
#endif


bool restore (struct Dump p)
{ trace ("restore", "/home/miro3/Documents/Programming/basilisk/src/output.h", 967);
  FILE * fp = p.fp;
  char * file = p.file;
  if (file && (fp = fopen (file, "r")) == NULL)
    { bool _ret =  false; end_trace("restore", "/home/miro3/Documents/Programming/basilisk/src/output.h", 971);  return _ret; }
  assert (fp);

  struct DumpHeader header;
  if (fread (&header, sizeof(header), 1, fp) < 1) {
    fprintf (qstderr(), "restore(): error: expecting header\n");
    exit (1);
  }


  init_grid (1);
   { foreach_cell(){

#line 982 "/home/miro3/Documents/Programming/basilisk/src/output.h"
 {
    cell.pid = pid();
    cell.flags |= active;
  } } end_foreach_cell(); }
  ((Tree *)grid)->dirty = true;
#line 1006 "/home/miro3/Documents/Programming/basilisk/src/output.h"
  bool restore_all = (p.list == all);
  scalar * list = dump_list (p.list ? p.list : all);
  if (header.version == 161020) {
    if (header.len - 1 != list_len (list)) {
      fprintf (qstderr(),
        "restore(): error: the list lengths don't match: "
        "%ld (file) != %d (code)\n",
        header.len - 1, list_len (list));
      exit (1);
    }
  }
  else {
    if (header.version != dump_version) {
      fprintf (qstderr(),
        "restore(): error: file version mismatch: "
        "%d (file) != %d (code)\n",
        header.version, dump_version);
      exit (1);
    }

    scalar * input = NULL;
    for (int i = 0; i < header.len; i++) {
      unsigned len;
      if (fread (&len, sizeof(unsigned), 1, fp) < 1) {
 fprintf (qstderr(), "restore(): error: expecting len\n");
 exit (1);
      }
      char name[len + 1];
      if (fread (name, sizeof(char), len, fp) < 1) {
 fprintf (qstderr(), "restore(): error: expecting s.name\n");
 exit (1);
      }
      name[len] = '\0';

      if (i > 0) {
 bool found = false;
 if (list) for (scalar s = *list, *_i82 = list; ((scalar *)&s)->i >= 0; s = *++_i82)
   if (!strcmp (_attribute[s.i].name, name)) {
     input = list_append (input, s);
     found = true; break;
   }
 if (!found) {
   if (restore_all) {
     scalar s = new_scalar("s");
     pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
     _attribute[s.i].name = pstrdup (name,__func__,__FILE__,__LINE__);
     input = list_append (input, s);
   }
   else
     input = list_append (input, (scalar){INT_MAX});
 }
      }
    }
    pfree (list,__func__,__FILE__,__LINE__);
    list = input;

    double o[4];
    if (fread (o, sizeof(double), 4, fp) < 4) {
      fprintf (qstderr(), "restore(): error: expecting coordinates\n");
      exit (1);
    }
    origin ((struct _origin){o[0], o[1], o[2]});
    size (o[3]);
  }
#line 1081 "/home/miro3/Documents/Programming/basilisk/src/output.h"
  scalar * listm = is_constant(cm) ? NULL : (scalar *)((vector []){{fm.x,fm.y,fm.z},{{-1},{-1},{-1}}});

  restore_mpi (fp, list);
#line 1109 "/home/miro3/Documents/Programming/basilisk/src/output.h"
  boundary (listm);

  scalar * other = NULL;
  if (all) for (scalar s = *all, *_i83 = all; ((scalar *)&s)->i >= 0; s = *++_i83)
    if (!list_lookup (list, s) && !list_lookup (listm, s))
      other = list_append (other, s);
  reset (other, 0.);
  pfree (other,__func__,__FILE__,__LINE__);

  pfree (list,__func__,__FILE__,__LINE__);
  if (file)
    fclose (fp);


  while (iter < header.i && events (false))
    iter = inext;
  events (false);
  while (t < header.t && events (false))
    t = tnext;
  t = header.t;
  events (false);

  { bool _ret =  true; end_trace("restore", "/home/miro3/Documents/Programming/basilisk/src/output.h", 1131);  return _ret; }
 end_trace("restore", "/home/miro3/Documents/Programming/basilisk/src/output.h", 1132); }
#line 288 "/home/miro3/Documents/Programming/basilisk/src/utils.h"
#line 12 "/home/miro3/Documents/Programming/basilisk/src/run.h"


void run (void)
{ trace ("run", "/home/miro3/Documents/Programming/basilisk/src/run.h", 15);
  iter = 0, t = 0., dt = 1.;
  init_grid (N);

  perf.nc = perf.tnc = 0;
  perf.gt = timer_start();
  while (events (true)) {





    update_perf();
    iter = inext, t = tnext;
  }




  timer_print (perf.gt, iter, perf.tnc);

  free_grid();
 end_trace("run", "/home/miro3/Documents/Programming/basilisk/src/run.h", 37); }
#line 25 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
#line 1 "./timestep.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/timestep.h"

double timestep (const vector u, double dtmax)
{
  static double previous = 0.;
  dtmax /= CFL;
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _dtmax = dtmax; 
#line 6

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 6
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 6
{

#line 6 "/home/miro3/Documents/Programming/basilisk/src/timestep.h"

    if (val(u.x,0,0,0) != 0.) {
      double dt = Delta*val_cm(cm,0,0,0)/fabs(val(u.x,0,0,0));
      if (dt < _dtmax) _dtmax = dt;
    } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 6
{

#line 6 "/home/miro3/Documents/Programming/basilisk/src/timestep.h"

    if (val(u.y,0,0,0) != 0.) {
      double dt = Delta*val_cm(cm,0,0,0)/fabs(val(u.y,0,0,0));
      if (dt < _dtmax) _dtmax = dt;
    } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 6
{

#line 6 "/home/miro3/Documents/Programming/basilisk/src/timestep.h"

    if (val(u.z,0,0,0) != 0.) {
      double dt = Delta*val_cm(cm,0,0,0)/fabs(val(u.z,0,0,0));
      if (dt < _dtmax) _dtmax = dt;
    } }  }}  end_foreach_face_generic()
#line 10
 end_foreach_face(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 6
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 6
{

#line 6 "/home/miro3/Documents/Programming/basilisk/src/timestep.h"

    if (val(u.x,0,0,0) != 0.) {
      double dt = Delta*val_cm(cm,0,0,0)/fabs(val(u.x,0,0,0));
      if (dt < _dtmax) _dtmax = dt;
    } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 6
{

#line 6 "/home/miro3/Documents/Programming/basilisk/src/timestep.h"

    if (val(u.y,0,0,0) != 0.) {
      double dt = Delta*val_cm(cm,0,0,0)/fabs(val(u.y,0,0,0));
      if (dt < _dtmax) _dtmax = dt;
    } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 6
{

#line 6 "/home/miro3/Documents/Programming/basilisk/src/timestep.h"

    if (val(u.z,0,0,0) != 0.) {
      double dt = Delta*val_cm(cm,0,0,0)/fabs(val(u.z,0,0,0));
      if (dt < _dtmax) _dtmax = dt;
    } }  }}  end_foreach_face_generic()
#line 10
 end_foreach_face(); }OMP(omp critical) if (_dtmax < dtmax) dtmax = _dtmax;
mpi_all_reduce_double (dtmax, MPI_MIN);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 10
 }
  dtmax *= CFL;
  if (dtmax > previous)
    dtmax = (previous + 0.1*dtmax)/1.1;
  previous = dtmax;
  return dtmax;
}
#line 26 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
#line 1 "./bcg.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/bcg.h"
#line 11 "/home/miro3/Documents/Programming/basilisk/src/bcg.h"
void tracer_fluxes (scalar f,
      vector uf,
      vector flux,
      double dt,
       scalar src)
{





  vector g= new_vector("g");
  gradients (((scalar []){f,{-1}}), ((vector []){{g.x,g.y,g.z},{{-1},{-1},{-1}}}));




   { 
if (!is_constant(fm.x) && !is_constant(src)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_src
#define val_src(a,i,j,k) val(a,i,j,k)
#undef fine_src
#define fine_src(a,i,j,k) fine(a,i,j,k)
#undef coarse_src
#define coarse_src(a,i,j,k) coarse(a,i,j,k)
#line 28
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 28
{

#line 28 "/home/miro3/Documents/Programming/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.x,0,0,0)/(val_fm_x(fm.x,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val_src(src,0,0,0) + val_src(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





      double vn = val(uf.y,i,0,0)/val_fm_y(fm.y,i,0,0) + val(uf.y,i,1,0)/val_fm_y(fm.y,i,1,0);
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.z,i,0,0)/val_fm_z(fm.z,i,0,0) + val(uf.z,i,0,1)/val_fm_z(fm.z,i,0,1);
      double fzz = wn < 0. ? val(f,i,0,1) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,0,-1);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 28
{

#line 28 "/home/miro3/Documents/Programming/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.y,0,0,0)/(val_fm_y(fm.y,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val_src(src,0,0,0) + val_src(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





      double vn = val(uf.z,0,i,0)/val_fm_z(fm.z,0,i,0) + val(uf.z,0,i,1)/val_fm_z(fm.z,0,i,1);
      double fyy = vn < 0. ? val(f,0,i,1) - val(f,0,i,0) : val(f,0,i,0) - val(f,0,i,-1);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.x,0,i,0)/val_fm_x(fm.x,0,i,0) + val(uf.x,1,i,0)/val_fm_x(fm.x,1,i,0);
      double fzz = wn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 28
{

#line 28 "/home/miro3/Documents/Programming/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.z,0,0,0)/(val_fm_z(fm.z,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,0,i) + (val_src(src,0,0,0) + val_src(src,0,0,-1))*dt/4. + s*(1. - s*un)*val(g.z,0,0,i)*Delta/2.;





      double vn = val(uf.x,0,0,i)/val_fm_x(fm.x,0,0,i) + val(uf.x,1,0,i)/val_fm_x(fm.x,1,0,i);
      double fyy = vn < 0. ? val(f,1,0,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,-1,0,i);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.y,0,0,i)/val_fm_y(fm.y,0,0,i) + val(uf.y,0,1,i)/val_fm_y(fm.y,0,1,i);
      double fzz = wn < 0. ? val(f,0,1,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,0,-1,i);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.z,0,0,0) = f2*val(uf.z,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 55
 end_foreach_face(); }
if (is_constant(fm.x) && !is_constant(src)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_src
#define val_src(a,i,j,k) val(a,i,j,k)
#undef fine_src
#define fine_src(a,i,j,k) fine(a,i,j,k)
#undef coarse_src
#define coarse_src(a,i,j,k) coarse(a,i,j,k)
#line 28
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 28
{

#line 28 "/home/miro3/Documents/Programming/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.x,0,0,0)/(val_fm_x(fm.x,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val_src(src,0,0,0) + val_src(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





      double vn = val(uf.y,i,0,0)/val_fm_y(fm.y,i,0,0) + val(uf.y,i,1,0)/val_fm_y(fm.y,i,1,0);
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.z,i,0,0)/val_fm_z(fm.z,i,0,0) + val(uf.z,i,0,1)/val_fm_z(fm.z,i,0,1);
      double fzz = wn < 0. ? val(f,i,0,1) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,0,-1);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 28
{

#line 28 "/home/miro3/Documents/Programming/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.y,0,0,0)/(val_fm_y(fm.y,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val_src(src,0,0,0) + val_src(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





      double vn = val(uf.z,0,i,0)/val_fm_z(fm.z,0,i,0) + val(uf.z,0,i,1)/val_fm_z(fm.z,0,i,1);
      double fyy = vn < 0. ? val(f,0,i,1) - val(f,0,i,0) : val(f,0,i,0) - val(f,0,i,-1);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.x,0,i,0)/val_fm_x(fm.x,0,i,0) + val(uf.x,1,i,0)/val_fm_x(fm.x,1,i,0);
      double fzz = wn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 28
{

#line 28 "/home/miro3/Documents/Programming/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.z,0,0,0)/(val_fm_z(fm.z,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,0,i) + (val_src(src,0,0,0) + val_src(src,0,0,-1))*dt/4. + s*(1. - s*un)*val(g.z,0,0,i)*Delta/2.;





      double vn = val(uf.x,0,0,i)/val_fm_x(fm.x,0,0,i) + val(uf.x,1,0,i)/val_fm_x(fm.x,1,0,i);
      double fyy = vn < 0. ? val(f,1,0,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,-1,0,i);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.y,0,0,i)/val_fm_y(fm.y,0,0,i) + val(uf.y,0,1,i)/val_fm_y(fm.y,0,1,i);
      double fzz = wn < 0. ? val(f,0,1,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,0,-1,i);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.z,0,0,0) = f2*val(uf.z,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 55
 end_foreach_face(); }
if (!is_constant(fm.x) && is_constant(src)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const double _const_src = _constant[src.i -_NVARMAX];
NOT_UNUSED(_const_src);
#undef val_src
#define val_src(a,i,j,k) _const_src
#undef fine_src
#define fine_src(a,i,j,k) _const_src
#undef coarse_src
#define coarse_src(a,i,j,k) _const_src
#line 28
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 28
{

#line 28 "/home/miro3/Documents/Programming/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.x,0,0,0)/(val_fm_x(fm.x,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val_src(src,0,0,0) + val_src(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





      double vn = val(uf.y,i,0,0)/val_fm_y(fm.y,i,0,0) + val(uf.y,i,1,0)/val_fm_y(fm.y,i,1,0);
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.z,i,0,0)/val_fm_z(fm.z,i,0,0) + val(uf.z,i,0,1)/val_fm_z(fm.z,i,0,1);
      double fzz = wn < 0. ? val(f,i,0,1) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,0,-1);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 28
{

#line 28 "/home/miro3/Documents/Programming/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.y,0,0,0)/(val_fm_y(fm.y,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val_src(src,0,0,0) + val_src(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





      double vn = val(uf.z,0,i,0)/val_fm_z(fm.z,0,i,0) + val(uf.z,0,i,1)/val_fm_z(fm.z,0,i,1);
      double fyy = vn < 0. ? val(f,0,i,1) - val(f,0,i,0) : val(f,0,i,0) - val(f,0,i,-1);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.x,0,i,0)/val_fm_x(fm.x,0,i,0) + val(uf.x,1,i,0)/val_fm_x(fm.x,1,i,0);
      double fzz = wn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 28
{

#line 28 "/home/miro3/Documents/Programming/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.z,0,0,0)/(val_fm_z(fm.z,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,0,i) + (val_src(src,0,0,0) + val_src(src,0,0,-1))*dt/4. + s*(1. - s*un)*val(g.z,0,0,i)*Delta/2.;





      double vn = val(uf.x,0,0,i)/val_fm_x(fm.x,0,0,i) + val(uf.x,1,0,i)/val_fm_x(fm.x,1,0,i);
      double fyy = vn < 0. ? val(f,1,0,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,-1,0,i);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.y,0,0,i)/val_fm_y(fm.y,0,0,i) + val(uf.y,0,1,i)/val_fm_y(fm.y,0,1,i);
      double fzz = wn < 0. ? val(f,0,1,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,0,-1,i);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.z,0,0,0) = f2*val(uf.z,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 55
 end_foreach_face(); }
if (is_constant(fm.x) && is_constant(src)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const double _const_src = _constant[src.i -_NVARMAX];
NOT_UNUSED(_const_src);
#undef val_src
#define val_src(a,i,j,k) _const_src
#undef fine_src
#define fine_src(a,i,j,k) _const_src
#undef coarse_src
#define coarse_src(a,i,j,k) _const_src
#line 28
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 28
{

#line 28 "/home/miro3/Documents/Programming/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.x,0,0,0)/(val_fm_x(fm.x,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val_src(src,0,0,0) + val_src(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





      double vn = val(uf.y,i,0,0)/val_fm_y(fm.y,i,0,0) + val(uf.y,i,1,0)/val_fm_y(fm.y,i,1,0);
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.z,i,0,0)/val_fm_z(fm.z,i,0,0) + val(uf.z,i,0,1)/val_fm_z(fm.z,i,0,1);
      double fzz = wn < 0. ? val(f,i,0,1) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,0,-1);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 28
{

#line 28 "/home/miro3/Documents/Programming/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.y,0,0,0)/(val_fm_y(fm.y,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val_src(src,0,0,0) + val_src(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





      double vn = val(uf.z,0,i,0)/val_fm_z(fm.z,0,i,0) + val(uf.z,0,i,1)/val_fm_z(fm.z,0,i,1);
      double fyy = vn < 0. ? val(f,0,i,1) - val(f,0,i,0) : val(f,0,i,0) - val(f,0,i,-1);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.x,0,i,0)/val_fm_x(fm.x,0,i,0) + val(uf.x,1,i,0)/val_fm_x(fm.x,1,i,0);
      double fzz = wn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 28
{

#line 28 "/home/miro3/Documents/Programming/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.z,0,0,0)/(val_fm_z(fm.z,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,0,i) + (val_src(src,0,0,0) + val_src(src,0,0,-1))*dt/4. + s*(1. - s*un)*val(g.z,0,0,i)*Delta/2.;





      double vn = val(uf.x,0,0,i)/val_fm_x(fm.x,0,0,i) + val(uf.x,1,0,i)/val_fm_x(fm.x,1,0,i);
      double fyy = vn < 0. ? val(f,1,0,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,-1,0,i);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.y,0,0,i)/val_fm_y(fm.y,0,0,i) + val(uf.y,0,1,i)/val_fm_y(fm.y,0,1,i);
      double fzz = wn < 0. ? val(f,0,1,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,0,-1,i);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.z,0,0,0) = f2*val(uf.z,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 55
 end_foreach_face(); } }





  boundary_flux (((vector []){{flux.x,flux.y,flux.z},{{-1},{-1},{-1}}}));
 delete (((scalar []){g.x,g.y,g.z,{-1}})); }






struct Advection {
  scalar * tracers;
  vector u;
  double dt;
  scalar * src;
};

void advection (struct Advection p)
{




  scalar * lsrc = p.src;
  if (!lsrc) {
    scalar zero= new_const_scalar("zero", 8,  0.);
    if (p.tracers) for (scalar s = *p.tracers, *_i84 = p.tracers; ((scalar *)&s)->i >= 0; s = *++_i84)
      lsrc = list_append (lsrc, zero);
  }

  assert (list_len(p.tracers) == list_len(lsrc));
  scalar f, src;
  scalar * _i2 = p.tracers; scalar * _i3 = lsrc; if (p.tracers) for (f = *p.tracers, src = *lsrc; ((scalar *)&f)->i >= 0; f = *++_i2, src = *++_i3) {
    vector flux= new_face_vector("flux");
    tracer_fluxes (f, p.u, flux, p.dt, src);
     { 
if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 94
foreach(){

#line 94 "/home/miro3/Documents/Programming/basilisk/src/bcg.h"

      {
#line 95

        val(f,0,0,0) += p.dt*(val(flux.x,0,0,0) - val(flux.x,1,0,0))/(Delta*val_cm(cm,0,0,0));
#line 95

        val(f,0,0,0) += p.dt*(val(flux.y,0,0,0) - val(flux.y,0,1,0))/(Delta*val_cm(cm,0,0,0));
#line 95

        val(f,0,0,0) += p.dt*(val(flux.z,0,0,0) - val(flux.z,0,0,1))/(Delta*val_cm(cm,0,0,0));}; } end_foreach(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 94
foreach(){

#line 94 "/home/miro3/Documents/Programming/basilisk/src/bcg.h"

      {
#line 95

        val(f,0,0,0) += p.dt*(val(flux.x,0,0,0) - val(flux.x,1,0,0))/(Delta*val_cm(cm,0,0,0));
#line 95

        val(f,0,0,0) += p.dt*(val(flux.y,0,0,0) - val(flux.y,0,1,0))/(Delta*val_cm(cm,0,0,0));
#line 95

        val(f,0,0,0) += p.dt*(val(flux.z,0,0,0) - val(flux.z,0,0,1))/(Delta*val_cm(cm,0,0,0));}; } end_foreach(); } }
   delete (((scalar []){flux.x,flux.y,flux.z,{-1}})); }
  boundary (p.tracers);

  if (!p.src)
    pfree (lsrc,__func__,__FILE__,__LINE__);
}
#line 27 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
#line 1 "./viscosity.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"
#line 1 "./poisson.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"
#line 32 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"
void mg_cycle (scalar * a, scalar * res, scalar * da,
        void (* relax) (scalar * da, scalar * res,
          int depth, void * data),
        void * data,
        int nrelax, int minlevel, int maxlevel)
{




  restriction (res);





  for (int l = minlevel; l <= maxlevel; l++) {




    if (l == minlevel)
       { foreach_level_or_leaf (l){

#line 54 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"

 if (da) for (scalar s = *da, *_i85 = da; ((scalar *)&s)->i >= 0; s = *++_i85)
   val(s,0,0,0) = 0.; } end_foreach_level_or_leaf(); }





    else
       { foreach_level (l){

#line 63 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"

 if (da) for (scalar s = *da, *_i86 = da; ((scalar *)&s)->i >= 0; s = *++_i86)
   val(s,0,0,0) = bilinear (point, s); } end_foreach_level(); }





    boundary_level (da, l);
    for (int i = 0; i < nrelax; i++) {
      relax (da, res, l, data);
      boundary_level (da, l);
    }
  }




   { foreach(){

#line 81 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"
 {
    scalar s, ds;
    scalar * _i4 = a; scalar * _i5 = da; if (a) for (s = *a, ds = *da; ((scalar *)&s)->i >= 0; s = *++_i4, ds = *++_i5)
      val(s,0,0,0) += val(ds,0,0,0);
  } } end_foreach(); }
  boundary (a);
}
#line 99 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"
int NITERMAX = 100, NITERMIN = 1;
double TOLERANCE = 1e-3;




typedef struct {
  int i;
  double resb, resa;
  double sum;
  int nrelax;
} mgstats;
#line 120 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"
struct MGSolve {
  scalar * a, * b;
  double (* residual) (scalar * a, scalar * b, scalar * res,
         void * data);
  void (* relax) (scalar * da, scalar * res, int depth,
    void * data);
  void * data;

  int nrelax;
  scalar * res;
};

mgstats mg_solve (struct MGSolve p)
{





  scalar * da = list_clone (p.a), * res = p.res;
  if (!res)
    if (p.a) for (scalar s = *p.a, *_i87 = p.a; ((scalar *)&s)->i >= 0; s = *++_i87) {
      scalar r = new_scalar("r");
      res = list_append (res, r);
    }






  for (int b = 0; b < nboundary; b++)
    if (da) for (scalar s = *da, *_i88 = da; ((scalar *)&s)->i >= 0; s = *++_i88)
      _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b];




  mgstats s = {0};
  double sum = 0.;
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _sum = sum; 
#line 160
foreach (){

#line 160 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"

    if (p.b) for (scalar s = *p.b, *_i89 = p.b; ((scalar *)&s)->i >= 0; s = *++_i89)
      _sum += val(s,0,0,0); } end_foreach();OMP(omp critical) sum += _sum;
mpi_all_reduce_double (sum, MPI_SUM);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 162
 }
  s.sum = sum;
  s.nrelax = p.nrelax > 0 ? p.nrelax : 4;




  double resb;
  resb = s.resb = s.resa = p.residual (p.a, p.b, res, p.data);






  for (s.i = 0;
       s.i < NITERMAX && (s.i < NITERMIN || s.resa > TOLERANCE);
       s.i++) {
    mg_cycle (p.a, res, da, p.relax, p.data, s.nrelax, 0, grid->maxdepth);
    s.resa = p.residual (p.a, p.b, res, p.data);







    if (s.resa > TOLERANCE) {
      if (resb/s.resa < 2 && s.nrelax < 100)
 s.nrelax++;
      else if (resb/s.resa > 20 && s.nrelax > 1)
 s.nrelax--;
    }
    resb = s.resa;
  }




  if (s.resa > TOLERANCE)
    fprintf (ferr,
      "WARNING: convergence not reached after %d iterations\n"
      "  res: %g sum: %g nrelax: %d\n",
      s.i, s.resa, s.sum, s.nrelax), fflush (ferr);




  if (!p.res)
    delete (res), pfree (res,__func__,__FILE__,__LINE__);
  delete (da), pfree (da,__func__,__FILE__,__LINE__);

  return s;
}
#line 237 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"
struct Poisson {
  scalar a, b;
   vector alpha;
   scalar lambda;
  double tolerance;
  int nrelax;
  scalar * res;
};





static void relax (scalar * al, scalar * bl, int l, void * data)
{
  scalar a = al[0], b = bl[0];
  struct Poisson * p = (struct Poisson *) data;
   vector alpha = p->alpha;
   scalar lambda = p->lambda;
#line 272 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"
  scalar c = a;






   { 
if (!is_constant(lambda) && !is_constant(alpha.x)) {
#undef val_lambda
#define val_lambda(a,i,j,k) val(a,i,j,k)
#undef fine_lambda
#define fine_lambda(a,i,j,k) fine(a,i,j,k)
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 279
foreach_level_or_leaf (l){

#line 279 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"
 {
    double n = - sq(Delta)*val(b,0,0,0), d = - val_lambda(lambda,0,0,0)*sq(Delta);
    {
#line 281
 {
      n += val_alpha_x(alpha.x,1,0,0)*val(a,1,0,0) + val_alpha_x(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val_alpha_x(alpha.x,1,0,0) + val_alpha_x(alpha.x,0,0,0);
    }
#line 281
 {
      n += val_alpha_y(alpha.y,0,1,0)*val(a,0,1,0) + val_alpha_y(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val_alpha_y(alpha.y,0,1,0) + val_alpha_y(alpha.y,0,0,0);
    }
#line 281
 {
      n += val_alpha_z(alpha.z,0,0,1)*val(a,0,0,1) + val_alpha_z(alpha.z,0,0,0)*val(a,0,0,-1);
      d += val_alpha_z(alpha.z,0,0,1) + val_alpha_z(alpha.z,0,0,0);
    }}
    val(c,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); }
if (is_constant(lambda) && !is_constant(alpha.x)) {
const double _const_lambda = _constant[lambda.i -_NVARMAX];
NOT_UNUSED(_const_lambda);
#undef val_lambda
#define val_lambda(a,i,j,k) _const_lambda
#undef fine_lambda
#define fine_lambda(a,i,j,k) _const_lambda
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) _const_lambda
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 279
foreach_level_or_leaf (l){

#line 279 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"
 {
    double n = - sq(Delta)*val(b,0,0,0), d = - val_lambda(lambda,0,0,0)*sq(Delta);
    {
#line 281
 {
      n += val_alpha_x(alpha.x,1,0,0)*val(a,1,0,0) + val_alpha_x(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val_alpha_x(alpha.x,1,0,0) + val_alpha_x(alpha.x,0,0,0);
    }
#line 281
 {
      n += val_alpha_y(alpha.y,0,1,0)*val(a,0,1,0) + val_alpha_y(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val_alpha_y(alpha.y,0,1,0) + val_alpha_y(alpha.y,0,0,0);
    }
#line 281
 {
      n += val_alpha_z(alpha.z,0,0,1)*val(a,0,0,1) + val_alpha_z(alpha.z,0,0,0)*val(a,0,0,-1);
      d += val_alpha_z(alpha.z,0,0,1) + val_alpha_z(alpha.z,0,0,0);
    }}
    val(c,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); }
if (!is_constant(lambda) && is_constant(alpha.x)) {
#undef val_lambda
#define val_lambda(a,i,j,k) val(a,i,j,k)
#undef fine_lambda
#define fine_lambda(a,i,j,k) fine(a,i,j,k)
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 279
foreach_level_or_leaf (l){

#line 279 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"
 {
    double n = - sq(Delta)*val(b,0,0,0), d = - val_lambda(lambda,0,0,0)*sq(Delta);
    {
#line 281
 {
      n += val_alpha_x(alpha.x,1,0,0)*val(a,1,0,0) + val_alpha_x(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val_alpha_x(alpha.x,1,0,0) + val_alpha_x(alpha.x,0,0,0);
    }
#line 281
 {
      n += val_alpha_y(alpha.y,0,1,0)*val(a,0,1,0) + val_alpha_y(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val_alpha_y(alpha.y,0,1,0) + val_alpha_y(alpha.y,0,0,0);
    }
#line 281
 {
      n += val_alpha_z(alpha.z,0,0,1)*val(a,0,0,1) + val_alpha_z(alpha.z,0,0,0)*val(a,0,0,-1);
      d += val_alpha_z(alpha.z,0,0,1) + val_alpha_z(alpha.z,0,0,0);
    }}
    val(c,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); }
if (is_constant(lambda) && is_constant(alpha.x)) {
const double _const_lambda = _constant[lambda.i -_NVARMAX];
NOT_UNUSED(_const_lambda);
#undef val_lambda
#define val_lambda(a,i,j,k) _const_lambda
#undef fine_lambda
#define fine_lambda(a,i,j,k) _const_lambda
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) _const_lambda
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 279
foreach_level_or_leaf (l){

#line 279 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"
 {
    double n = - sq(Delta)*val(b,0,0,0), d = - val_lambda(lambda,0,0,0)*sq(Delta);
    {
#line 281
 {
      n += val_alpha_x(alpha.x,1,0,0)*val(a,1,0,0) + val_alpha_x(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val_alpha_x(alpha.x,1,0,0) + val_alpha_x(alpha.x,0,0,0);
    }
#line 281
 {
      n += val_alpha_y(alpha.y,0,1,0)*val(a,0,1,0) + val_alpha_y(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val_alpha_y(alpha.y,0,1,0) + val_alpha_y(alpha.y,0,0,0);
    }
#line 281
 {
      n += val_alpha_z(alpha.z,0,0,1)*val(a,0,0,1) + val_alpha_z(alpha.z,0,0,0)*val(a,0,0,-1);
      d += val_alpha_z(alpha.z,0,0,1) + val_alpha_z(alpha.z,0,0,0);
    }}
    val(c,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); } }
#line 304 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"
}






static double residual (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  struct Poisson * p = (struct Poisson *) data;
   vector alpha = p->alpha;
   scalar lambda = p->lambda;
  double maxres = 0.;


  vector g= new_face_vector("g");
   { 
if (!is_constant(alpha.x)) {
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 321
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 321
{

#line 321 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"

    val(g.x,0,0,0) = val_alpha_x(alpha.x,0,0,0)*(val(a,0,0,0) - val(a,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 321
{

#line 321 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"

    val(g.y,0,0,0) = val_alpha_y(alpha.y,0,0,0)*(val(a,0,0,0) - val(a,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 321
{

#line 321 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"

    val(g.z,0,0,0) = val_alpha_z(alpha.z,0,0,0)*(val(a,0,0,0) - val(a,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 322
 end_foreach_face(); }
if (is_constant(alpha.x)) {
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 321
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 321
{

#line 321 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"

    val(g.x,0,0,0) = val_alpha_x(alpha.x,0,0,0)*(val(a,0,0,0) - val(a,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 321
{

#line 321 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"

    val(g.y,0,0,0) = val_alpha_y(alpha.y,0,0,0)*(val(a,0,0,0) - val(a,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 321
{

#line 321 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"

    val(g.z,0,0,0) = val_alpha_z(alpha.z,0,0,0)*(val(a,0,0,0) - val(a,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 322
 end_foreach_face(); } }
  boundary_flux (((vector []){{g.x,g.y,g.z},{{-1},{-1},{-1}}}));
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _maxres = maxres; 
#line 324

if (!is_constant(lambda)) {
#undef val_lambda
#define val_lambda(a,i,j,k) val(a,i,j,k)
#undef fine_lambda
#define fine_lambda(a,i,j,k) fine(a,i,j,k)
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) coarse(a,i,j,k)
#line 324
foreach (){

#line 324 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"
 {
    val(res,0,0,0) = val(b,0,0,0) - val_lambda(lambda,0,0,0)*val(a,0,0,0);
    {
#line 326

      val(res,0,0,0) += (val(g.x,0,0,0) - val(g.x,1,0,0))/Delta;
#line 326

      val(res,0,0,0) += (val(g.y,0,0,0) - val(g.y,0,1,0))/Delta;
#line 326

      val(res,0,0,0) += (val(g.z,0,0,0) - val(g.z,0,0,1))/Delta;}
    if (fabs (val(res,0,0,0)) > _maxres)
      _maxres = fabs (val(res,0,0,0));
  } } end_foreach(); }
if (is_constant(lambda)) {
const double _const_lambda = _constant[lambda.i -_NVARMAX];
NOT_UNUSED(_const_lambda);
#undef val_lambda
#define val_lambda(a,i,j,k) _const_lambda
#undef fine_lambda
#define fine_lambda(a,i,j,k) _const_lambda
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) _const_lambda
#line 324
foreach (){

#line 324 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"
 {
    val(res,0,0,0) = val(b,0,0,0) - val_lambda(lambda,0,0,0)*val(a,0,0,0);
    {
#line 326

      val(res,0,0,0) += (val(g.x,0,0,0) - val(g.x,1,0,0))/Delta;
#line 326

      val(res,0,0,0) += (val(g.y,0,0,0) - val(g.y,0,1,0))/Delta;
#line 326

      val(res,0,0,0) += (val(g.z,0,0,0) - val(g.z,0,0,1))/Delta;}
    if (fabs (val(res,0,0,0)) > _maxres)
      _maxres = fabs (val(res,0,0,0));
  } } end_foreach(); }OMP(omp critical) if (_maxres > maxres) maxres = _maxres;
mpi_all_reduce_double (maxres, MPI_MAX);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 330
 }
#line 342 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"
  boundary (resl);
  { double _ret =  maxres; delete (((scalar []){g.x,g.y,g.z,{-1}}));  return _ret; }
 delete (((scalar []){g.x,g.y,g.z,{-1}})); }
#line 355 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"
mgstats poisson (struct Poisson p)
{






  if (!p.alpha.x.i) {
    vector alpha= new_const_vector("alpha", 9, (double []) {1.,1.,1.});
    p.alpha = alpha;
  }
  if (!p.lambda.i) {
    scalar lambda= new_const_scalar("lambda", 12,  0.);
    p.lambda = lambda;
  }




  vector alpha = p.alpha;
  scalar lambda = p.lambda;
  restriction (((scalar []){alpha.x,alpha.y,alpha.z,lambda,{-1}}));





  double defaultol = TOLERANCE;
  if (p.tolerance)
    TOLERANCE = p.tolerance;

  scalar a = p.a, b = p.b;
  mgstats s = mg_solve ((struct MGSolve){((scalar []){a,{-1}}), ((scalar []){b,{-1}}), residual, relax, &p, p.nrelax, p.res});




  if (p.tolerance)
    TOLERANCE = defaultol;

  return s;
}
#line 416 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"
struct Project {
  vector u;
  scalar p;
  vector alpha;
  double dt;
  int nrelax;
};


mgstats project (struct Project q)
{ trace ("project", "/home/miro3/Documents/Programming/basilisk/src/poisson.h", 426);
  vector u = q.u;
  scalar p = q.p;
   vector alpha = q.alpha.x.i ? q.alpha : unityf;
  double dt = q.dt ? q.dt : 1.;
  int nrelax = q.nrelax ? q.nrelax : 4;






  scalar div= new_scalar("div");
   { foreach(){

#line 439 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"
 {
    val(div,0,0,0) = 0.;
    {
#line 441

      val(div,0,0,0) += val(u.x,1,0,0) - val(u.x,0,0,0);
#line 441

      val(div,0,0,0) += val(u.y,0,1,0) - val(u.y,0,0,0);
#line 441

      val(div,0,0,0) += val(u.z,0,0,1) - val(u.z,0,0,0);}
    val(div,0,0,0) /= dt*Delta;
  } } end_foreach(); }
#line 455 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"
  mgstats mgp = poisson ((struct Poisson){p, div, alpha,
    .tolerance = TOLERANCE/sq(dt), .nrelax = nrelax});




   { 
if (!is_constant(alpha.x)) {
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 461
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 461
{

#line 461 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"

    val(u.x,0,0,0) -= dt*val_alpha_x(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 461
{

#line 461 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"

    val(u.y,0,0,0) -= dt*val_alpha_y(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 461
{

#line 461 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"

    val(u.z,0,0,0) -= dt*val_alpha_z(alpha.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 462
 end_foreach_face(); }
if (is_constant(alpha.x)) {
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 461
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 461
{

#line 461 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"

    val(u.x,0,0,0) -= dt*val_alpha_x(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 461
{

#line 461 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"

    val(u.y,0,0,0) -= dt*val_alpha_y(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 461
{

#line 461 "/home/miro3/Documents/Programming/basilisk/src/poisson.h"

    val(u.z,0,0,0) -= dt*val_alpha_z(alpha.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 462
 end_foreach_face(); } }
  boundary ((scalar *)((vector []){{u.x,u.y,u.z},{{-1},{-1},{-1}}}));

  { mgstats _ret =  mgp; delete (((scalar []){div,{-1}}));  end_trace("project", "/home/miro3/Documents/Programming/basilisk/src/poisson.h", 465);  return _ret; }
 delete (((scalar []){div,{-1}}));  end_trace("project", "/home/miro3/Documents/Programming/basilisk/src/poisson.h", 466); }
#line 2 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"

struct Viscosity {
  vector u;
  vector mu;
  scalar rho;
  double dt;
  int nrelax;
  scalar * res;
};
#line 25 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"
static void relax_viscosity (scalar * a, scalar * b, int l, void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
   vector mu = p->mu;
   scalar rho = p->rho;
  double dt = p->dt;
  vector u = (*((vector *)&(a[0]))), r = (*((vector *)&(b[0])));




  vector w = u;


   { 
if (!is_constant(rho) && !is_constant(mu.x)) {
#undef val_rho
#define val_rho(a,i,j,k) val(a,i,j,k)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,i,j,k)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_x
#define val_mu_x(a,i,j,k) val(a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_y
#define val_mu_y(a,i,j,k) val(a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_z
#define val_mu_z(a,i,j,k) val(a,i,j,k)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,i,j,k)
#line 39
foreach_level_or_leaf (l){

#line 39 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"
 {
    {
#line 40

      val(w.x,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)*val(u.x,-1,0,0)

      + val_mu_y(mu.y,0,1,0)*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


      + val_mu_z(mu.z,0,0,1)*(val(u.x,0,0,1) +
       (val(u.z,1,0,0) + val(u.z,1,0,1))/4. -
       (val(u.z,-1,0,0) + val(u.z,-1,0,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.x,0,0,-1) +
         (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
         (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).x + dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)

          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)


          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)

        ));
#line 40

      val(w.y,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)*val(u.y,0,-1,0)

      + val_mu_z(mu.z,0,0,1)*(val(u.y,0,0,1) +
     (val(u.z,0,1,0) + val(u.z,0,1,1))/4. -
     (val(u.z,0,-1,0) + val(u.z,0,-1,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.y,0,0,-1) +
         (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
         (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


      + val_mu_x(mu.x,1,0,0)*(val(u.y,1,0,0) +
       (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
       (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).y + dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)

          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)


          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)

        ));
#line 40

      val(w.z,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1)*val(u.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)*val(u.z,0,0,-1)

      + val_mu_x(mu.x,1,0,0)*(val(u.z,1,0,0) +
     (val(u.x,0,0,1) + val(u.x,1,0,1))/4. -
     (val(u.x,0,0,-1) + val(u.x,1,0,-1))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.z,-1,0,0) +
         (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
         (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


      + val_mu_y(mu.y,0,1,0)*(val(u.z,0,1,0) +
       (val(u.y,0,0,1) + val(u.y,0,1,1))/4. -
       (val(u.y,0,0,-1) + val(u.y,0,1,-1))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.z,0,-1,0) +
         (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
         (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

      ) + val(r.z,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).z + dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)

          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)


          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)

        ));}
  } } end_foreach_level_or_leaf(); }
if (is_constant(rho) && !is_constant(mu.x)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,i,j,k) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
#undef val_mu_x
#define val_mu_x(a,i,j,k) val(a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_y
#define val_mu_y(a,i,j,k) val(a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_z
#define val_mu_z(a,i,j,k) val(a,i,j,k)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,i,j,k)
#line 39
foreach_level_or_leaf (l){

#line 39 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"
 {
    {
#line 40

      val(w.x,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)*val(u.x,-1,0,0)

      + val_mu_y(mu.y,0,1,0)*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


      + val_mu_z(mu.z,0,0,1)*(val(u.x,0,0,1) +
       (val(u.z,1,0,0) + val(u.z,1,0,1))/4. -
       (val(u.z,-1,0,0) + val(u.z,-1,0,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.x,0,0,-1) +
         (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
         (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).x + dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)

          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)


          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)

        ));
#line 40

      val(w.y,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)*val(u.y,0,-1,0)

      + val_mu_z(mu.z,0,0,1)*(val(u.y,0,0,1) +
     (val(u.z,0,1,0) + val(u.z,0,1,1))/4. -
     (val(u.z,0,-1,0) + val(u.z,0,-1,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.y,0,0,-1) +
         (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
         (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


      + val_mu_x(mu.x,1,0,0)*(val(u.y,1,0,0) +
       (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
       (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).y + dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)

          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)


          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)

        ));
#line 40

      val(w.z,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1)*val(u.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)*val(u.z,0,0,-1)

      + val_mu_x(mu.x,1,0,0)*(val(u.z,1,0,0) +
     (val(u.x,0,0,1) + val(u.x,1,0,1))/4. -
     (val(u.x,0,0,-1) + val(u.x,1,0,-1))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.z,-1,0,0) +
         (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
         (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


      + val_mu_y(mu.y,0,1,0)*(val(u.z,0,1,0) +
       (val(u.y,0,0,1) + val(u.y,0,1,1))/4. -
       (val(u.y,0,0,-1) + val(u.y,0,1,-1))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.z,0,-1,0) +
         (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
         (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

      ) + val(r.z,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).z + dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)

          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)


          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)

        ));}
  } } end_foreach_level_or_leaf(); }
if (!is_constant(rho) && is_constant(mu.x)) {
#undef val_rho
#define val_rho(a,i,j,k) val(a,i,j,k)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,i,j,k)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_mu = {_constant[mu.x.i -_NVARMAX], _constant[mu.y.i - _NVARMAX], _constant[mu.z.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_x
#define val_mu_x(a,i,j,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,i,j,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_z
#define val_mu_z(a,i,j,k) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#line 39
foreach_level_or_leaf (l){

#line 39 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"
 {
    {
#line 40

      val(w.x,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)*val(u.x,-1,0,0)

      + val_mu_y(mu.y,0,1,0)*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


      + val_mu_z(mu.z,0,0,1)*(val(u.x,0,0,1) +
       (val(u.z,1,0,0) + val(u.z,1,0,1))/4. -
       (val(u.z,-1,0,0) + val(u.z,-1,0,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.x,0,0,-1) +
         (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
         (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).x + dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)

          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)


          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)

        ));
#line 40

      val(w.y,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)*val(u.y,0,-1,0)

      + val_mu_z(mu.z,0,0,1)*(val(u.y,0,0,1) +
     (val(u.z,0,1,0) + val(u.z,0,1,1))/4. -
     (val(u.z,0,-1,0) + val(u.z,0,-1,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.y,0,0,-1) +
         (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
         (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


      + val_mu_x(mu.x,1,0,0)*(val(u.y,1,0,0) +
       (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
       (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).y + dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)

          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)


          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)

        ));
#line 40

      val(w.z,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1)*val(u.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)*val(u.z,0,0,-1)

      + val_mu_x(mu.x,1,0,0)*(val(u.z,1,0,0) +
     (val(u.x,0,0,1) + val(u.x,1,0,1))/4. -
     (val(u.x,0,0,-1) + val(u.x,1,0,-1))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.z,-1,0,0) +
         (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
         (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


      + val_mu_y(mu.y,0,1,0)*(val(u.z,0,1,0) +
       (val(u.y,0,0,1) + val(u.y,0,1,1))/4. -
       (val(u.y,0,0,-1) + val(u.y,0,1,-1))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.z,0,-1,0) +
         (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
         (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

      ) + val(r.z,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).z + dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)

          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)


          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)

        ));}
  } } end_foreach_level_or_leaf(); }
if (is_constant(rho) && is_constant(mu.x)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,i,j,k) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
const struct { double x, y, z; } _const_mu = {_constant[mu.x.i -_NVARMAX], _constant[mu.y.i - _NVARMAX], _constant[mu.z.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_x
#define val_mu_x(a,i,j,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,i,j,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_z
#define val_mu_z(a,i,j,k) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#line 39
foreach_level_or_leaf (l){

#line 39 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"
 {
    {
#line 40

      val(w.x,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)*val(u.x,-1,0,0)

      + val_mu_y(mu.y,0,1,0)*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


      + val_mu_z(mu.z,0,0,1)*(val(u.x,0,0,1) +
       (val(u.z,1,0,0) + val(u.z,1,0,1))/4. -
       (val(u.z,-1,0,0) + val(u.z,-1,0,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.x,0,0,-1) +
         (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
         (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).x + dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)

          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)


          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)

        ));
#line 40

      val(w.y,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)*val(u.y,0,-1,0)

      + val_mu_z(mu.z,0,0,1)*(val(u.y,0,0,1) +
     (val(u.z,0,1,0) + val(u.z,0,1,1))/4. -
     (val(u.z,0,-1,0) + val(u.z,0,-1,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.y,0,0,-1) +
         (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
         (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


      + val_mu_x(mu.x,1,0,0)*(val(u.y,1,0,0) +
       (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
       (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).y + dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)

          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)


          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)

        ));
#line 40

      val(w.z,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1)*val(u.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)*val(u.z,0,0,-1)

      + val_mu_x(mu.x,1,0,0)*(val(u.z,1,0,0) +
     (val(u.x,0,0,1) + val(u.x,1,0,1))/4. -
     (val(u.x,0,0,-1) + val(u.x,1,0,-1))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.z,-1,0,0) +
         (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
         (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


      + val_mu_y(mu.y,0,1,0)*(val(u.z,0,1,0) +
       (val(u.y,0,0,1) + val(u.y,0,1,1))/4. -
       (val(u.y,0,0,-1) + val(u.y,0,1,-1))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.z,0,-1,0) +
         (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
         (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

      ) + val(r.z,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).z + dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)

          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)


          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)

        ));}
  } } end_foreach_level_or_leaf(); } }
#line 85 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"
}

static double residual_viscosity (scalar * a, scalar * b, scalar * resl,
      void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
   vector mu = p->mu;
   scalar rho = p->rho;
  double dt = p->dt;
  vector u = (*((vector *)&(a[0]))), r = (*((vector *)&(b[0]))), res = (*((vector *)&(resl[0])));
  double maxres = 0.;


  {
#line 98
 {
    vector taux= new_face_vector("taux");
     { 
if (!is_constant(mu.x)) {
#undef val_mu_x
#define val_mu_x(a,i,j,k) val(a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_y
#define val_mu_y(a,i,j,k) val(a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_z
#define val_mu_z(a,i,j,k) val(a,i,j,k)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,i,j,k)
#line 100
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 100
{

#line 100 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"

      val(taux.x,0,0,0) = 2.*val_mu_x(mu.x,0,0,0)*(val(u.x,0,0,0) - val(u.x,-1,0,0))/Delta; }  }}  end_foreach_face_generic()
#line 101
 end_foreach_face(); }
if (is_constant(mu.x)) {
const struct { double x, y, z; } _const_mu = {_constant[mu.x.i -_NVARMAX], _constant[mu.y.i - _NVARMAX], _constant[mu.z.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_x
#define val_mu_x(a,i,j,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,i,j,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_z
#define val_mu_z(a,i,j,k) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#line 100
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 100
{

#line 100 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"

      val(taux.x,0,0,0) = 2.*val_mu_x(mu.x,0,0,0)*(val(u.x,0,0,0) - val(u.x,-1,0,0))/Delta; }  }}  end_foreach_face_generic()
#line 101
 end_foreach_face(); } }

       { 
if (!is_constant(mu.x)) {
#undef val_mu_x
#define val_mu_x(a,i,j,k) val(a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_y
#define val_mu_y(a,i,j,k) val(a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_z
#define val_mu_z(a,i,j,k) val(a,i,j,k)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,i,j,k)
#line 103
foreach_face_generic() { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 103
{

#line 103 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"

 val(taux.y,0,0,0) = val_mu_y(mu.y,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,-1,0) +
      (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
      (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 106
 end_foreach_face(); }
if (is_constant(mu.x)) {
const struct { double x, y, z; } _const_mu = {_constant[mu.x.i -_NVARMAX], _constant[mu.y.i - _NVARMAX], _constant[mu.z.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_x
#define val_mu_x(a,i,j,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,i,j,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_z
#define val_mu_z(a,i,j,k) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#line 103
foreach_face_generic() { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 103
{

#line 103 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"

 val(taux.y,0,0,0) = val_mu_y(mu.y,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,-1,0) +
      (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
      (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 106
 end_foreach_face(); } }


       { 
if (!is_constant(mu.x)) {
#undef val_mu_x
#define val_mu_x(a,i,j,k) val(a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_y
#define val_mu_y(a,i,j,k) val(a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_z
#define val_mu_z(a,i,j,k) val(a,i,j,k)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,i,j,k)
#line 109
foreach_face_generic() { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 109
{

#line 109 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"

 val(taux.z,0,0,0) = val_mu_z(mu.z,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,0,-1) +
      (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
      (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 112
 end_foreach_face(); }
if (is_constant(mu.x)) {
const struct { double x, y, z; } _const_mu = {_constant[mu.x.i -_NVARMAX], _constant[mu.y.i - _NVARMAX], _constant[mu.z.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_x
#define val_mu_x(a,i,j,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,i,j,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_z
#define val_mu_z(a,i,j,k) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#line 109
foreach_face_generic() { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 109
{

#line 109 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"

 val(taux.z,0,0,0) = val_mu_z(mu.z,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,0,-1) +
      (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
      (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 112
 end_foreach_face(); } }

    boundary_flux (((vector []){{taux.x,taux.y,taux.z},{{-1},{-1},{-1}}}));
     { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _maxres = maxres; 
#line 115

if (!is_constant(rho)) {
#undef val_rho
#define val_rho(a,i,j,k) val(a,i,j,k)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,i,j,k)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,i,j,k)
#line 115
foreach (){

#line 115 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"
 {
      double d = 0.;
      {
#line 117

 d += val(taux.x,1,0,0) - val(taux.x,0,0,0);
#line 117

 d += val(taux.y,0,1,0) - val(taux.y,0,0,0);
#line 117

 d += val(taux.z,0,0,1) - val(taux.z,0,0,0);}
      val(res.x,0,0,0) = val(r.x,0,0,0) - ((coord){1.,1.,1.}).x*val(u.x,0,0,0) + dt/val_rho(rho,0,0,0)*d/Delta;
      if (fabs (val(res.x,0,0,0)) > _maxres)
 _maxres = fabs (val(res.x,0,0,0));
    } } end_foreach(); }
if (is_constant(rho)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,i,j,k) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
#line 115
foreach (){

#line 115 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"
 {
      double d = 0.;
      {
#line 117

 d += val(taux.x,1,0,0) - val(taux.x,0,0,0);
#line 117

 d += val(taux.y,0,1,0) - val(taux.y,0,0,0);
#line 117

 d += val(taux.z,0,0,1) - val(taux.z,0,0,0);}
      val(res.x,0,0,0) = val(r.x,0,0,0) - ((coord){1.,1.,1.}).x*val(u.x,0,0,0) + dt/val_rho(rho,0,0,0)*d/Delta;
      if (fabs (val(res.x,0,0,0)) > _maxres)
 _maxres = fabs (val(res.x,0,0,0));
    } } end_foreach(); }OMP(omp critical) if (_maxres > maxres) maxres = _maxres;
mpi_all_reduce_double (maxres, MPI_MAX);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 122
 }
   delete (((scalar []){taux.x,taux.y,taux.z,{-1}})); }
#line 98
 {
    vector taux= new_face_vector("taux");
     { 
if (!is_constant(mu.y)) {
#undef val_mu_y
#define val_mu_y(a,k,i,j) val(a,k,i,j)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,k,i,j)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,k,i,j)
#undef val_mu_z
#define val_mu_z(a,k,i,j) val(a,k,i,j)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,k,i,j)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,k,i,j)
#undef val_mu_x
#define val_mu_x(a,k,i,j) val(a,k,i,j)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,k,i,j)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,k,i,j)
#line 100
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_y()) {
#line 100
{

#line 100 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"

      val(taux.y,0,0,0) = 2.*val_mu_y(mu.y,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,-1,0))/Delta; }  }}  end_foreach_face_generic()
#line 101
 end_foreach_face(); }
if (is_constant(mu.y)) {
const struct { double x, y, z; } _const_mu = {_constant[mu.y.i -_NVARMAX], _constant[mu.z.i - _NVARMAX], _constant[mu.x.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_y
#define val_mu_y(a,k,i,j) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_z
#define val_mu_z(a,k,i,j) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#undef val_mu_x
#define val_mu_x(a,k,i,j) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#line 100
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_y()) {
#line 100
{

#line 100 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"

      val(taux.y,0,0,0) = 2.*val_mu_y(mu.y,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,-1,0))/Delta; }  }}  end_foreach_face_generic()
#line 101
 end_foreach_face(); } }

       { 
if (!is_constant(mu.y)) {
#undef val_mu_y
#define val_mu_y(a,k,i,j) val(a,k,i,j)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,k,i,j)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,k,i,j)
#undef val_mu_z
#define val_mu_z(a,k,i,j) val(a,k,i,j)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,k,i,j)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,k,i,j)
#undef val_mu_x
#define val_mu_x(a,k,i,j) val(a,k,i,j)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,k,i,j)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,k,i,j)
#line 103
foreach_face_generic() { int jg = -1; VARIABLES;  if (is_face_z()) {
#line 103
{

#line 103 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"

 val(taux.z,0,0,0) = val_mu_z(mu.z,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,0,-1) +
      (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
      (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 106
 end_foreach_face(); }
if (is_constant(mu.y)) {
const struct { double x, y, z; } _const_mu = {_constant[mu.y.i -_NVARMAX], _constant[mu.z.i - _NVARMAX], _constant[mu.x.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_y
#define val_mu_y(a,k,i,j) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_z
#define val_mu_z(a,k,i,j) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#undef val_mu_x
#define val_mu_x(a,k,i,j) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#line 103
foreach_face_generic() { int jg = -1; VARIABLES;  if (is_face_z()) {
#line 103
{

#line 103 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"

 val(taux.z,0,0,0) = val_mu_z(mu.z,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,0,-1) +
      (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
      (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 106
 end_foreach_face(); } }


       { 
if (!is_constant(mu.y)) {
#undef val_mu_y
#define val_mu_y(a,k,i,j) val(a,k,i,j)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,k,i,j)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,k,i,j)
#undef val_mu_z
#define val_mu_z(a,k,i,j) val(a,k,i,j)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,k,i,j)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,k,i,j)
#undef val_mu_x
#define val_mu_x(a,k,i,j) val(a,k,i,j)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,k,i,j)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,k,i,j)
#line 109
foreach_face_generic() { int kg = -1; VARIABLES;  if (is_face_x()) {
#line 109
{

#line 109 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"

 val(taux.x,0,0,0) = val_mu_x(mu.x,0,0,0)*(val(u.y,0,0,0) - val(u.y,-1,0,0) +
      (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
      (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 112
 end_foreach_face(); }
if (is_constant(mu.y)) {
const struct { double x, y, z; } _const_mu = {_constant[mu.y.i -_NVARMAX], _constant[mu.z.i - _NVARMAX], _constant[mu.x.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_y
#define val_mu_y(a,k,i,j) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_z
#define val_mu_z(a,k,i,j) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#undef val_mu_x
#define val_mu_x(a,k,i,j) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#line 109
foreach_face_generic() { int kg = -1; VARIABLES;  if (is_face_x()) {
#line 109
{

#line 109 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"

 val(taux.x,0,0,0) = val_mu_x(mu.x,0,0,0)*(val(u.y,0,0,0) - val(u.y,-1,0,0) +
      (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
      (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 112
 end_foreach_face(); } }

    boundary_flux (((vector []){{taux.x,taux.y,taux.z},{{-1},{-1},{-1}}}));
     { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _maxres = maxres; 
#line 115

if (!is_constant(rho)) {
#undef val_rho
#define val_rho(a,k,i,j) val(a,k,i,j)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,k,i,j)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,k,i,j)
#line 115
foreach (){

#line 115 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"
 {
      double d = 0.;
      {
#line 117

 d += val(taux.y,0,1,0) - val(taux.y,0,0,0);
#line 117

 d += val(taux.z,0,0,1) - val(taux.z,0,0,0);
#line 117

 d += val(taux.x,1,0,0) - val(taux.x,0,0,0);}
      val(res.y,0,0,0) = val(r.y,0,0,0) - ((coord){1.,1.,1.}).y*val(u.y,0,0,0) + dt/val_rho(rho,0,0,0)*d/Delta;
      if (fabs (val(res.y,0,0,0)) > _maxres)
 _maxres = fabs (val(res.y,0,0,0));
    } } end_foreach(); }
if (is_constant(rho)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,k,i,j) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
#line 115
foreach (){

#line 115 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"
 {
      double d = 0.;
      {
#line 117

 d += val(taux.y,0,1,0) - val(taux.y,0,0,0);
#line 117

 d += val(taux.z,0,0,1) - val(taux.z,0,0,0);
#line 117

 d += val(taux.x,1,0,0) - val(taux.x,0,0,0);}
      val(res.y,0,0,0) = val(r.y,0,0,0) - ((coord){1.,1.,1.}).y*val(u.y,0,0,0) + dt/val_rho(rho,0,0,0)*d/Delta;
      if (fabs (val(res.y,0,0,0)) > _maxres)
 _maxres = fabs (val(res.y,0,0,0));
    } } end_foreach(); }OMP(omp critical) if (_maxres > maxres) maxres = _maxres;
mpi_all_reduce_double (maxres, MPI_MAX);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 122
 }
   delete (((scalar []){taux.x,taux.y,taux.z,{-1}})); }
#line 98
 {
    vector taux= new_face_vector("taux");
     { 
if (!is_constant(mu.z)) {
#undef val_mu_z
#define val_mu_z(a,j,k,i) val(a,j,k,i)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,j,k,i)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,j,k,i)
#undef val_mu_x
#define val_mu_x(a,j,k,i) val(a,j,k,i)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,j,k,i)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,j,k,i)
#undef val_mu_y
#define val_mu_y(a,j,k,i) val(a,j,k,i)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,j,k,i)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,j,k,i)
#line 100
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_z()) {
#line 100
{

#line 100 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"

      val(taux.z,0,0,0) = 2.*val_mu_z(mu.z,0,0,0)*(val(u.z,0,0,0) - val(u.z,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 101
 end_foreach_face(); }
if (is_constant(mu.z)) {
const struct { double x, y, z; } _const_mu = {_constant[mu.z.i -_NVARMAX], _constant[mu.x.i - _NVARMAX], _constant[mu.y.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_z
#define val_mu_z(a,j,k,i) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#undef val_mu_x
#define val_mu_x(a,j,k,i) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,j,k,i) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#line 100
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_z()) {
#line 100
{

#line 100 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"

      val(taux.z,0,0,0) = 2.*val_mu_z(mu.z,0,0,0)*(val(u.z,0,0,0) - val(u.z,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 101
 end_foreach_face(); } }

       { 
if (!is_constant(mu.z)) {
#undef val_mu_z
#define val_mu_z(a,j,k,i) val(a,j,k,i)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,j,k,i)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,j,k,i)
#undef val_mu_x
#define val_mu_x(a,j,k,i) val(a,j,k,i)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,j,k,i)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,j,k,i)
#undef val_mu_y
#define val_mu_y(a,j,k,i) val(a,j,k,i)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,j,k,i)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,j,k,i)
#line 103
foreach_face_generic() { int jg = -1; VARIABLES;  if (is_face_x()) {
#line 103
{

#line 103 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"

 val(taux.x,0,0,0) = val_mu_x(mu.x,0,0,0)*(val(u.z,0,0,0) - val(u.z,-1,0,0) +
      (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
      (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 106
 end_foreach_face(); }
if (is_constant(mu.z)) {
const struct { double x, y, z; } _const_mu = {_constant[mu.z.i -_NVARMAX], _constant[mu.x.i - _NVARMAX], _constant[mu.y.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_z
#define val_mu_z(a,j,k,i) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#undef val_mu_x
#define val_mu_x(a,j,k,i) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,j,k,i) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#line 103
foreach_face_generic() { int jg = -1; VARIABLES;  if (is_face_x()) {
#line 103
{

#line 103 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"

 val(taux.x,0,0,0) = val_mu_x(mu.x,0,0,0)*(val(u.z,0,0,0) - val(u.z,-1,0,0) +
      (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
      (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 106
 end_foreach_face(); } }


       { 
if (!is_constant(mu.z)) {
#undef val_mu_z
#define val_mu_z(a,j,k,i) val(a,j,k,i)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,j,k,i)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,j,k,i)
#undef val_mu_x
#define val_mu_x(a,j,k,i) val(a,j,k,i)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,j,k,i)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,j,k,i)
#undef val_mu_y
#define val_mu_y(a,j,k,i) val(a,j,k,i)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,j,k,i)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,j,k,i)
#line 109
foreach_face_generic() { int kg = -1; VARIABLES;  if (is_face_y()) {
#line 109
{

#line 109 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"

 val(taux.y,0,0,0) = val_mu_y(mu.y,0,0,0)*(val(u.z,0,0,0) - val(u.z,0,-1,0) +
      (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
      (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 112
 end_foreach_face(); }
if (is_constant(mu.z)) {
const struct { double x, y, z; } _const_mu = {_constant[mu.z.i -_NVARMAX], _constant[mu.x.i - _NVARMAX], _constant[mu.y.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_z
#define val_mu_z(a,j,k,i) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#undef val_mu_x
#define val_mu_x(a,j,k,i) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,j,k,i) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#line 109
foreach_face_generic() { int kg = -1; VARIABLES;  if (is_face_y()) {
#line 109
{

#line 109 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"

 val(taux.y,0,0,0) = val_mu_y(mu.y,0,0,0)*(val(u.z,0,0,0) - val(u.z,0,-1,0) +
      (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
      (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 112
 end_foreach_face(); } }

    boundary_flux (((vector []){{taux.x,taux.y,taux.z},{{-1},{-1},{-1}}}));
     { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _maxres = maxres; 
#line 115

if (!is_constant(rho)) {
#undef val_rho
#define val_rho(a,j,k,i) val(a,j,k,i)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,j,k,i)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,j,k,i)
#line 115
foreach (){

#line 115 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"
 {
      double d = 0.;
      {
#line 117

 d += val(taux.z,0,0,1) - val(taux.z,0,0,0);
#line 117

 d += val(taux.x,1,0,0) - val(taux.x,0,0,0);
#line 117

 d += val(taux.y,0,1,0) - val(taux.y,0,0,0);}
      val(res.z,0,0,0) = val(r.z,0,0,0) - ((coord){1.,1.,1.}).z*val(u.z,0,0,0) + dt/val_rho(rho,0,0,0)*d/Delta;
      if (fabs (val(res.z,0,0,0)) > _maxres)
 _maxres = fabs (val(res.z,0,0,0));
    } } end_foreach(); }
if (is_constant(rho)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,j,k,i) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
#line 115
foreach (){

#line 115 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"
 {
      double d = 0.;
      {
#line 117

 d += val(taux.z,0,0,1) - val(taux.z,0,0,0);
#line 117

 d += val(taux.x,1,0,0) - val(taux.x,0,0,0);
#line 117

 d += val(taux.y,0,1,0) - val(taux.y,0,0,0);}
      val(res.z,0,0,0) = val(r.z,0,0,0) - ((coord){1.,1.,1.}).z*val(u.z,0,0,0) + dt/val_rho(rho,0,0,0)*d/Delta;
      if (fabs (val(res.z,0,0,0)) > _maxres)
 _maxres = fabs (val(res.z,0,0,0));
    } } end_foreach(); }OMP(omp critical) if (_maxres > maxres) maxres = _maxres;
mpi_all_reduce_double (maxres, MPI_MAX);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 122
 }
   delete (((scalar []){taux.x,taux.y,taux.z,{-1}})); }}
  boundary (resl);
#line 153 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"
  return maxres;
}



mgstats viscosity (struct Viscosity p)
{
  vector u = p.u, r= new_vector("r");
   { foreach(){

#line 161 "/home/miro3/Documents/Programming/basilisk/src/viscosity.h"

    {
#line 162

      val(r.x,0,0,0) = val(u.x,0,0,0);
#line 162

      val(r.y,0,0,0) = val(u.y,0,0,0);
#line 162

      val(r.z,0,0,0) = val(u.z,0,0,0);}; } end_foreach(); }

  vector mu = p.mu;
  scalar rho = p.rho;
  restriction (((scalar []){mu.x,mu.y,mu.z,rho,{-1}}));

  { mgstats _ret =  mg_solve ((struct MGSolve){(scalar *)((vector []){{u.x,u.y,u.z},{{-1},{-1},{-1}}}), (scalar *)((vector []){{r.x,r.y,r.z},{{-1},{-1},{-1}}}),
     residual_viscosity, relax_viscosity, &p, p.nrelax, p.res}); delete (((scalar []){r.x,r.y,r.z,{-1}}));  return _ret; }
 delete (((scalar []){r.x,r.y,r.z,{-1}})); }
#line 28 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
#line 37 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
scalar p= {0};
vector u= {{1},{2},{3}}, g= {{4},{5},{6}};
scalar pf= {7};
vector uf= {{8},{9},{10}};
#line 63 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
 vector mu = {{_NVARMAX + 0},{_NVARMAX + 1},{_NVARMAX + 2}}, a = {{_NVARMAX + 0},{_NVARMAX + 1},{_NVARMAX + 2}}, alpha = {{_NVARMAX + 3},{_NVARMAX + 4},{_NVARMAX + 5}};
 scalar rho = {(_NVARMAX + 6)};
mgstats mgp, mgpf, mgu;
bool stokes = false;
#line 77 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
static void _set_boundary0 (void) { _attribute[p.i].boundary[right] = _boundary0; _attribute[p.i].boundary_homogeneous[right] = _boundary0_homogeneous; } 
static void _set_boundary1 (void) { _attribute[p.i].boundary[left] = _boundary1; _attribute[p.i].boundary_homogeneous[left] = _boundary1_homogeneous; } 





static void _set_boundary2 (void) { _attribute[p.i].boundary[top] = _boundary2; _attribute[p.i].boundary_homogeneous[top] = _boundary2_homogeneous; } 
static void _set_boundary3 (void) { _attribute[p.i].boundary[bottom] = _boundary3; _attribute[p.i].boundary_homogeneous[bottom] = _boundary3_homogeneous; } 


static void _set_boundary4 (void) { _attribute[p.i].boundary[front] = _boundary4; _attribute[p.i].boundary_homogeneous[front] = _boundary4_homogeneous; } 
static void _set_boundary5 (void) { _attribute[p.i].boundary[back] = _boundary5; _attribute[p.i].boundary_homogeneous[back] = _boundary5_homogeneous; } 






static int defaults_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int defaults (const int i, const double t, Event * _ev) { trace ("defaults", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 96); 
{

  CFL = 0.8;




  _attribute[p.i].nodump = _attribute[pf.i].nodump = true;




  if (alpha.x.i == unityf.x.i) {
    alpha = fm;
    rho = cm;
  }
  else if (!is_constant(alpha.x)) {
    vector alphav = alpha;
     { 
if (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 115
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 115
{

#line 115 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

      val(alphav.x,0,0,0) = val_fm_x(fm.x,0,0,0); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 115
{

#line 115 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

      val(alphav.y,0,0,0) = val_fm_y(fm.y,0,0,0); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 115
{

#line 115 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

      val(alphav.z,0,0,0) = val_fm_z(fm.z,0,0,0); }  }}  end_foreach_face_generic()
#line 116
 end_foreach_face(); }
if (is_constant(fm.x)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 115
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 115
{

#line 115 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

      val(alphav.x,0,0,0) = val_fm_x(fm.x,0,0,0); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 115
{

#line 115 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

      val(alphav.y,0,0,0) = val_fm_y(fm.y,0,0,0); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 115
{

#line 115 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

      val(alphav.z,0,0,0) = val_fm_z(fm.z,0,0,0); }  }}  end_foreach_face_generic()
#line 116
 end_foreach_face(); } }
    boundary ((scalar *)((vector []){{alpha.x,alpha.y,alpha.z},{{-1},{-1},{-1}}}));
  }






  _attribute[uf.x.i].refine = refine_face_solenoidal;

 end_trace("defaults", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 127); } return 0; } 





double dtmax;

static int init_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int init (const int i, const double t, Event * _ev) { trace ("init", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 135); 
{
  boundary ((scalar *)((vector []){{u.x,u.y,u.z},{{-1},{-1},{-1}}}));
  trash (((vector []){{uf.x,uf.y,uf.z},{{-1},{-1},{-1}}}));
   { 
if (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 139
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 139
{

#line 139 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(val(u.x,0,0,0) + val(u.x,-1,0,0))/2.; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 139
{

#line 139 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(val(u.y,0,0,0) + val(u.y,0,-1,0))/2.; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 139
{

#line 139 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(uf.z,0,0,0) = val_fm_z(fm.z,0,0,0)*(val(u.z,0,0,0) + val(u.z,0,0,-1))/2.; }  }}  end_foreach_face_generic()
#line 140
 end_foreach_face(); }
if (is_constant(fm.x)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 139
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 139
{

#line 139 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(val(u.x,0,0,0) + val(u.x,-1,0,0))/2.; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 139
{

#line 139 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(val(u.y,0,0,0) + val(u.y,0,-1,0))/2.; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 139
{

#line 139 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(uf.z,0,0,0) = val_fm_z(fm.z,0,0,0)*(val(u.z,0,0,0) + val(u.z,0,0,-1))/2.; }  }}  end_foreach_face_generic()
#line 140
 end_foreach_face(); } }
  boundary ((scalar *)((vector []){{uf.x,uf.y,uf.z},{{-1},{-1},{-1}}}));




  event ("properties");





  dtmax = DT;
  event ("stability");
 end_trace("init", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 154); } return 0; } 
#line 163 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
static int set_dtmax_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int set_dtmax (const int i, const double t, Event * _ev) { trace ("set_dtmax", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 163);  dtmax = DT; end_trace("set_dtmax", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 163);  return 0; } 

static int stability_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int stability (const int i, const double t, Event * _ev) { trace ("stability", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 165);  {
  dt = dtnext (timestep (uf, dtmax));
 end_trace("stability", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 167); } return 0; } 







static int vof_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int vof (const int i, const double t, Event * _ev) { trace ("vof", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 175); ; end_trace("vof", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 175);  return 0; } 
static int tracer_advection_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int tracer_advection (const int i, const double t, Event * _ev) { trace ("tracer_advection", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 176); ; end_trace("tracer_advection", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 176);  return 0; } 
static int tracer_diffusion_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int tracer_diffusion (const int i, const double t, Event * _ev) { trace ("tracer_diffusion", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 177); ; end_trace("tracer_diffusion", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 177);  return 0; } 






static int properties_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int properties (const int i, const double t, Event * _ev) { trace ("properties", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 184);  {
  boundary (((scalar []){alpha.x,alpha.y,alpha.z,mu.x,mu.y,mu.z,rho,{-1}}));
 end_trace("properties", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 186); } return 0; } 
#line 198 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
void prediction()
{
  vector du;
  {
#line 201
 {
    scalar s = new_scalar("s");
    du.x = s;
  }
#line 201
 {
    scalar s = new_scalar("s");
    du.y = s;
  }
#line 201
 {
    scalar s = new_scalar("s");
    du.z = s;
  }}

  if (_attribute[u.x.i].gradient)
     { foreach(){

#line 207 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

      {
#line 208

        val(du.x,0,0,0) = _attribute[u.x.i].gradient (val(u.x,-1,0,0), val(u.x,0,0,0), val(u.x,1,0,0))/Delta;
#line 208

        val(du.y,0,0,0) = _attribute[u.y.i].gradient (val(u.y,0,-1,0), val(u.y,0,0,0), val(u.y,0,1,0))/Delta;
#line 208

        val(du.z,0,0,0) = _attribute[u.z.i].gradient (val(u.z,0,0,-1), val(u.z,0,0,0), val(u.z,0,0,1))/Delta;}; } end_foreach(); }
  else
     { foreach(){

#line 211 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

      {
#line 212

        val(du.x,0,0,0) = (val(u.x,1,0,0) - val(u.x,-1,0,0))/(2.*Delta);
#line 212

        val(du.y,0,0,0) = (val(u.y,0,1,0) - val(u.y,0,-1,0))/(2.*Delta);
#line 212

        val(du.z,0,0,0) = (val(u.z,0,0,1) - val(u.z,0,0,-1))/(2.*Delta);}; } end_foreach(); }
  boundary ((scalar *)((vector []){{du.x,du.y,du.z},{{-1},{-1},{-1}}}));

  trash (((vector []){{uf.x,uf.y,uf.z},{{-1},{-1},{-1}}}));
   { 
if (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 217
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 217
{

#line 217 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.x,0,0,0) + val(u.x,-1,0,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.x,0,0,0) = val(u.x,i,0,0) + (val(g.x,0,0,0) + val(g.x,-1,0,0))*dt/4. + s*(1. - s*un)*val(du.x,i,0,0)*Delta/2.;

      double fyy = val(u.y,i,0,0) < 0. ? val(u.x,i,1,0) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,-1,0);
      val(uf.x,0,0,0) -= dt*val(u.y,i,0,0)*fyy/(2.*Delta);


      double fzz = val(u.z,i,0,0) < 0. ? val(u.x,i,0,1) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,0,-1);
      val(uf.x,0,0,0) -= dt*val(u.z,i,0,0)*fzz/(2.*Delta);

    val(uf.x,0,0,0) *= val_fm_x(fm.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 217
{

#line 217 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.y,0,0,0) + val(u.y,0,-1,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.y,0,0,0) = val(u.y,0,i,0) + (val(g.y,0,0,0) + val(g.y,0,-1,0))*dt/4. + s*(1. - s*un)*val(du.y,0,i,0)*Delta/2.;

      double fyy = val(u.z,0,i,0) < 0. ? val(u.y,0,i,1) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,0,i,-1);
      val(uf.y,0,0,0) -= dt*val(u.z,0,i,0)*fyy/(2.*Delta);


      double fzz = val(u.x,0,i,0) < 0. ? val(u.y,1,i,0) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,-1,i,0);
      val(uf.y,0,0,0) -= dt*val(u.x,0,i,0)*fzz/(2.*Delta);

    val(uf.y,0,0,0) *= val_fm_y(fm.y,0,0,0);
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 217
{

#line 217 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.z,0,0,0) + val(u.z,0,0,-1))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.z,0,0,0) = val(u.z,0,0,i) + (val(g.z,0,0,0) + val(g.z,0,0,-1))*dt/4. + s*(1. - s*un)*val(du.z,0,0,i)*Delta/2.;

      double fyy = val(u.x,0,0,i) < 0. ? val(u.z,1,0,i) - val(u.z,0,0,i) : val(u.z,0,0,i) - val(u.z,-1,0,i);
      val(uf.z,0,0,0) -= dt*val(u.x,0,0,i)*fyy/(2.*Delta);


      double fzz = val(u.y,0,0,i) < 0. ? val(u.z,0,1,i) - val(u.z,0,0,i) : val(u.z,0,0,i) - val(u.z,0,-1,i);
      val(uf.z,0,0,0) -= dt*val(u.y,0,0,i)*fzz/(2.*Delta);

    val(uf.z,0,0,0) *= val_fm_z(fm.z,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 230
 end_foreach_face(); }
if (is_constant(fm.x)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 217
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 217
{

#line 217 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.x,0,0,0) + val(u.x,-1,0,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.x,0,0,0) = val(u.x,i,0,0) + (val(g.x,0,0,0) + val(g.x,-1,0,0))*dt/4. + s*(1. - s*un)*val(du.x,i,0,0)*Delta/2.;

      double fyy = val(u.y,i,0,0) < 0. ? val(u.x,i,1,0) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,-1,0);
      val(uf.x,0,0,0) -= dt*val(u.y,i,0,0)*fyy/(2.*Delta);


      double fzz = val(u.z,i,0,0) < 0. ? val(u.x,i,0,1) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,0,-1);
      val(uf.x,0,0,0) -= dt*val(u.z,i,0,0)*fzz/(2.*Delta);

    val(uf.x,0,0,0) *= val_fm_x(fm.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 217
{

#line 217 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.y,0,0,0) + val(u.y,0,-1,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.y,0,0,0) = val(u.y,0,i,0) + (val(g.y,0,0,0) + val(g.y,0,-1,0))*dt/4. + s*(1. - s*un)*val(du.y,0,i,0)*Delta/2.;

      double fyy = val(u.z,0,i,0) < 0. ? val(u.y,0,i,1) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,0,i,-1);
      val(uf.y,0,0,0) -= dt*val(u.z,0,i,0)*fyy/(2.*Delta);


      double fzz = val(u.x,0,i,0) < 0. ? val(u.y,1,i,0) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,-1,i,0);
      val(uf.y,0,0,0) -= dt*val(u.x,0,i,0)*fzz/(2.*Delta);

    val(uf.y,0,0,0) *= val_fm_y(fm.y,0,0,0);
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 217
{

#line 217 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.z,0,0,0) + val(u.z,0,0,-1))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.z,0,0,0) = val(u.z,0,0,i) + (val(g.z,0,0,0) + val(g.z,0,0,-1))*dt/4. + s*(1. - s*un)*val(du.z,0,0,i)*Delta/2.;

      double fyy = val(u.x,0,0,i) < 0. ? val(u.z,1,0,i) - val(u.z,0,0,i) : val(u.z,0,0,i) - val(u.z,-1,0,i);
      val(uf.z,0,0,0) -= dt*val(u.x,0,0,i)*fyy/(2.*Delta);


      double fzz = val(u.y,0,0,i) < 0. ? val(u.z,0,1,i) - val(u.z,0,0,i) : val(u.z,0,0,i) - val(u.z,0,-1,i);
      val(uf.z,0,0,0) -= dt*val(u.y,0,0,i)*fzz/(2.*Delta);

    val(uf.z,0,0,0) *= val_fm_z(fm.z,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 230
 end_foreach_face(); } }
  boundary ((scalar *)((vector []){{uf.x,uf.y,uf.z},{{-1},{-1},{-1}}}));

  delete ((scalar *)((vector []){{du.x,du.y,du.z},{{-1},{-1},{-1}}}));
}
#line 245 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
static int advection_term_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int advection_term (const int i, const double t, Event * _ev) { trace ("advection_term", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 245); 
{
  if (!stokes) {
    prediction();
    mgpf = project ((struct Project){uf, pf, alpha, dt/2., mgpf.nrelax});
    advection ((struct Advection){(scalar *)((vector []){{u.x,u.y,u.z},{{-1},{-1},{-1}}}), uf, dt, (scalar *)((vector []){{g.x,g.y,g.z},{{-1},{-1},{-1}}})});
  }
 end_trace("advection_term", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 252); } return 0; } 







static void correction (double dt)
{
   { foreach(){

#line 262 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    {
#line 263

      val(u.x,0,0,0) += dt*val(g.x,0,0,0);
#line 263

      val(u.y,0,0,0) += dt*val(g.y,0,0,0);
#line 263

      val(u.z,0,0,0) += dt*val(g.z,0,0,0);}; } end_foreach(); }
  boundary ((scalar *)((vector []){{u.x,u.y,u.z},{{-1},{-1},{-1}}}));
}
#line 275 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
static int viscous_term_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int viscous_term (const int i, const double t, Event * _ev) { trace ("viscous_term", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 275); 
{
  if (constant(mu.x) != 0.) {
    correction (dt);
    mgu = viscosity ((struct Viscosity){u, mu, rho, dt, mgu.nrelax});
    correction (-dt);
  }






  vector af = a;
  trash (((vector []){{uf.x,uf.y,uf.z},{af.x,af.y,af.z},{{-1},{-1},{-1}}}));
   { 
if (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 290
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 290
{

#line 290 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
 {
    val(uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(val(u.x,0,0,0) + val(u.x,-1,0,0))/2.;
    if (!is_constant(af.x))
      val(af.x,0,0,0) = 0.;
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 290
{

#line 290 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
 {
    val(uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(val(u.y,0,0,0) + val(u.y,0,-1,0))/2.;
    if (!is_constant(af.y))
      val(af.y,0,0,0) = 0.;
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 290
{

#line 290 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
 {
    val(uf.z,0,0,0) = val_fm_z(fm.z,0,0,0)*(val(u.z,0,0,0) + val(u.z,0,0,-1))/2.;
    if (!is_constant(af.z))
      val(af.z,0,0,0) = 0.;
  } }  }}  end_foreach_face_generic()
#line 294
 end_foreach_face(); }
if (is_constant(fm.x)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 290
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 290
{

#line 290 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
 {
    val(uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(val(u.x,0,0,0) + val(u.x,-1,0,0))/2.;
    if (!is_constant(af.x))
      val(af.x,0,0,0) = 0.;
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 290
{

#line 290 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
 {
    val(uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(val(u.y,0,0,0) + val(u.y,0,-1,0))/2.;
    if (!is_constant(af.y))
      val(af.y,0,0,0) = 0.;
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 290
{

#line 290 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
 {
    val(uf.z,0,0,0) = val_fm_z(fm.z,0,0,0)*(val(u.z,0,0,0) + val(u.z,0,0,-1))/2.;
    if (!is_constant(af.z))
      val(af.z,0,0,0) = 0.;
  } }  }}  end_foreach_face_generic()
#line 294
 end_foreach_face(); } }
 end_trace("viscous_term", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 295); } return 0; } 
#line 310 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
static int acceleration_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int acceleration (const int i, const double t, Event * _ev) { trace ("acceleration", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 310); 
{
  boundary ((scalar *)((vector []){{a.x,a.y,a.z},{{-1},{-1},{-1}}}));
   { 
if (!is_constant(fm.x) && !is_constant(a.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#line 313
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 313
{

#line 313 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) += dt*val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 313
{

#line 313 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) += dt*val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 313
{

#line 313 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(uf.z,0,0,0) += dt*val_fm_z(fm.z,0,0,0)*val_a_z(a.z,0,0,0); }  }}  end_foreach_face_generic()
#line 314
 end_foreach_face(); }
if (is_constant(fm.x) && !is_constant(a.x)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#line 313
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 313
{

#line 313 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) += dt*val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 313
{

#line 313 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) += dt*val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 313
{

#line 313 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(uf.z,0,0,0) += dt*val_fm_z(fm.z,0,0,0)*val_a_z(a.z,0,0,0); }  }}  end_foreach_face_generic()
#line 314
 end_foreach_face(); }
if (!is_constant(fm.x) && is_constant(a.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#line 313
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 313
{

#line 313 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) += dt*val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 313
{

#line 313 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) += dt*val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 313
{

#line 313 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(uf.z,0,0,0) += dt*val_fm_z(fm.z,0,0,0)*val_a_z(a.z,0,0,0); }  }}  end_foreach_face_generic()
#line 314
 end_foreach_face(); }
if (is_constant(fm.x) && is_constant(a.x)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#line 313
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 313
{

#line 313 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) += dt*val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 313
{

#line 313 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) += dt*val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 313
{

#line 313 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(uf.z,0,0,0) += dt*val_fm_z(fm.z,0,0,0)*val_a_z(a.z,0,0,0); }  }}  end_foreach_face_generic()
#line 314
 end_foreach_face(); } }
  boundary ((scalar *)((vector []){{uf.x,uf.y,uf.z},{{-1},{-1},{-1}}}));
 end_trace("acceleration", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 316); } return 0; } 
#line 325 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
static int projection_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int projection (const int i, const double t, Event * _ev) { trace ("projection", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 325); 
{
  mgp = project ((struct Project){uf, p, alpha, dt, mgp.nrelax});





  vector gf= new_face_vector("gf");
   { 
if (!is_constant(a.x) && !is_constant(alpha.x) && !is_constant(fm.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 334
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 334
{

#line 334 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 334
{

#line 334 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 334
{

#line 334 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 335
 end_foreach_face(); }
if (is_constant(a.x) && !is_constant(alpha.x) && !is_constant(fm.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 334
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 334
{

#line 334 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 334
{

#line 334 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 334
{

#line 334 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 335
 end_foreach_face(); }
if (!is_constant(a.x) && is_constant(alpha.x) && !is_constant(fm.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 334
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 334
{

#line 334 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 334
{

#line 334 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 334
{

#line 334 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 335
 end_foreach_face(); }
if (is_constant(a.x) && is_constant(alpha.x) && !is_constant(fm.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 334
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 334
{

#line 334 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 334
{

#line 334 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 334
{

#line 334 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 335
 end_foreach_face(); }
if (!is_constant(a.x) && !is_constant(alpha.x) && is_constant(fm.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 334
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 334
{

#line 334 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 334
{

#line 334 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 334
{

#line 334 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 335
 end_foreach_face(); }
if (is_constant(a.x) && !is_constant(alpha.x) && is_constant(fm.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 334
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 334
{

#line 334 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 334
{

#line 334 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 334
{

#line 334 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 335
 end_foreach_face(); }
if (!is_constant(a.x) && is_constant(alpha.x) && is_constant(fm.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 334
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 334
{

#line 334 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 334
{

#line 334 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 334
{

#line 334 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 335
 end_foreach_face(); }
if (is_constant(a.x) && is_constant(alpha.x) && is_constant(fm.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 334
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 334
{

#line 334 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 334
{

#line 334 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 334
{

#line 334 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 335
 end_foreach_face(); } }
  boundary_flux (((vector []){{gf.x,gf.y,gf.z},{{-1},{-1},{-1}}}));





  trash (((vector []){{g.x,g.y,g.z},{{-1},{-1},{-1}}}));
   { foreach(){

#line 343 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

    {
#line 344

      val(g.x,0,0,0) = (val(gf.x,0,0,0) + val(gf.x,1,0,0))/2.;
#line 344

      val(g.y,0,0,0) = (val(gf.y,0,0,0) + val(gf.y,0,1,0))/2.;
#line 344

      val(g.z,0,0,0) = (val(gf.z,0,0,0) + val(gf.z,0,0,1))/2.;}; } end_foreach(); }
  boundary ((scalar *)((vector []){{g.x,g.y,g.z},{{-1},{-1},{-1}}}));




  correction (dt);
 delete (((scalar []){gf.x,gf.y,gf.z,{-1}}));  end_trace("projection", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 352); } return 0; } 







static int adapt_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int adapt (const int i, const double t, Event * _ev) { trace ("adapt", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 360);  {
  event ("properties");
 end_trace("adapt", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 362); } return 0; } 
#line 5 "milk_crown.c"
#line 1 "two-phase.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/two-phase.h"
#line 13 "/home/miro3/Documents/Programming/basilisk/src/two-phase.h"
#line 1 "vof.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
#line 26 "/home/miro3/Documents/Programming/basilisk/src/vof.h"




#line 43 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
extern scalar * interfaces;
extern vector uf;
extern double dt;





static int defaults_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int defaults_0 (const int i, const double t, Event * _ev) { trace ("defaults_0", "/home/miro3/Documents/Programming/basilisk/src/vof.h", 51); 
{

  if (interfaces) for (scalar c = *interfaces, *_i90 = interfaces; ((scalar *)&c)->i >= 0; c = *++_i90)
    _attribute[c.i].refine = _attribute[c.i].prolongation = fraction_refine;

 end_trace("defaults_0", "/home/miro3/Documents/Programming/basilisk/src/vof.h", 57); } return 0; } 





static int stability_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int stability_0 (const int i, const double t, Event * _ev) { trace ("stability_0", "/home/miro3/Documents/Programming/basilisk/src/vof.h", 63);  {
  if (CFL > 0.5)
    CFL = 0.5;
 end_trace("stability_0", "/home/miro3/Documents/Programming/basilisk/src/vof.h", 66); } return 0; } 
#line 80 "/home/miro3/Documents/Programming/basilisk/src/vof.h"

#line 80

static void sweep_x (scalar c, scalar cc)
{
  vector n= new_vector("n");
  scalar alpha= new_scalar("alpha"), flux= new_scalar("flux");
  double cfl = 0.;
#line 94 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
  scalar * tracers = _attribute[c.i].tracers, * gfl = NULL, * tfluxl = NULL;
  if (tracers) {
    if (tracers) for (scalar t = *tracers, *_i91 = tracers; ((scalar *)&t)->i >= 0; t = *++_i91) {
      scalar gf = new_scalar("gf"), flux = new_scalar("flux");
      gfl = list_append (gfl, gf);
      tfluxl = list_append (tfluxl, flux);
    }







     { foreach(){

#line 108 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
 {
      scalar t, gf;
      scalar * _i6 = tracers; scalar * _i7 = gfl; if (tracers) for (t = *tracers, gf = *gfl; ((scalar *)&t)->i >= 0; t = *++_i6, gf = *++_i7) {
 double cl = val(c,-1,0,0), cc = val(c,0,0,0), cr = val(c,1,0,0);
 if (_attribute[t.i].inverse)
   cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;
 val(gf,0,0,0) = 0.;
 static const double cmin = 0.5;
 if (cc >= cmin) {
   if (cr >= cmin) {
     if (cl >= cmin) {
       if (_attribute[t.i].gradient)
  val(gf,0,0,0) = _attribute[t.i].gradient (val(t,-1,0,0)/cl, val(t,0,0,0)/cc, val(t,1,0,0)/cr)/Delta;
       else
  val(gf,0,0,0) = (val(t,1,0,0)/cr - val(t,-1,0,0)/cl)/(2.*Delta);
     }
     else
        val(gf,0,0,0) = (val(t,1,0,0)/cr - val(t,0,0,0)/cc)/Delta;
   }
   else if (cl >= cmin)
     val(gf,0,0,0) = (val(t,0,0,0)/cc - val(t,-1,0,0)/cl)/Delta;
 }
      }
    } } end_foreach(); }
    boundary (gfl);
  }






  reconstruction (c, n, alpha);

   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _cfl = cfl; 
#line 142

if (!is_constant(fm.x) && !is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 142
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 142
{

#line 142 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
 {






    double un = val(uf.x,0,0,0)*dt/(Delta*val_fm_x(fm.x,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_x(fm.x,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_x(fm.x,0,0,0)*s/val_cm(cm,0,0,0);
#line 169 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
    double cf = (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.) ? val(c,i,0,0) :
      rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), val(n.z,i,0,0)}, val(alpha,i,0,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.x,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; if (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = val(c,i,0,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,i,0,0)/ci + s*min(1., 1. - s*un)*val(gf,i,0,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.x,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_generic()
#line 197
 end_foreach_face(); }
if (is_constant(fm.x) && !is_constant(cm)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 142
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 142
{

#line 142 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
 {






    double un = val(uf.x,0,0,0)*dt/(Delta*val_fm_x(fm.x,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_x(fm.x,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_x(fm.x,0,0,0)*s/val_cm(cm,0,0,0);
#line 169 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
    double cf = (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.) ? val(c,i,0,0) :
      rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), val(n.z,i,0,0)}, val(alpha,i,0,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.x,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; if (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = val(c,i,0,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,i,0,0)/ci + s*min(1., 1. - s*un)*val(gf,i,0,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.x,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_generic()
#line 197
 end_foreach_face(); }
if (!is_constant(fm.x) && is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 142
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 142
{

#line 142 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
 {






    double un = val(uf.x,0,0,0)*dt/(Delta*val_fm_x(fm.x,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_x(fm.x,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_x(fm.x,0,0,0)*s/val_cm(cm,0,0,0);
#line 169 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
    double cf = (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.) ? val(c,i,0,0) :
      rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), val(n.z,i,0,0)}, val(alpha,i,0,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.x,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; if (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = val(c,i,0,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,i,0,0)/ci + s*min(1., 1. - s*un)*val(gf,i,0,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.x,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_generic()
#line 197
 end_foreach_face(); }
if (is_constant(fm.x) && is_constant(cm)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 142
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 142
{

#line 142 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
 {






    double un = val(uf.x,0,0,0)*dt/(Delta*val_fm_x(fm.x,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_x(fm.x,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_x(fm.x,0,0,0)*s/val_cm(cm,0,0,0);
#line 169 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
    double cf = (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.) ? val(c,i,0,0) :
      rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), val(n.z,i,0,0)}, val(alpha,i,0,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.x,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; if (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = val(c,i,0,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,i,0,0)/ci + s*min(1., 1. - s*un)*val(gf,i,0,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.x,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_generic()
#line 197
 end_foreach_face(); }OMP(omp critical) if (_cfl > cfl) cfl = _cfl;
mpi_all_reduce_double (cfl, MPI_MAX);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 197
 }
  delete (gfl); pfree (gfl,__func__,__FILE__,__LINE__);
#line 208 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
  scalar * fluxl = list_concat (NULL, tfluxl);
  fluxl = list_append (fluxl, flux);
  for (int l = depth() - 1; l >= 0; l--)
     { foreach_halo (prolongation, l){

#line 211 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
 {
#line 227 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
      if ((!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0))
 if (fluxl) for (scalar fl = *fluxl, *_i92 = fluxl; ((scalar *)&fl)->i >= 0; fl = *++_i92)
   val(fl,0,0,0) = (fine(fl,0,0,0) + fine(fl,0,1,0) +
    fine(fl,0,0,1) + fine(fl,0,1,1))/4.;
      if ((!is_leaf (neighbor(1,0,0)) && neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0))
 if (fluxl) for (scalar fl = *fluxl, *_i93 = fluxl; ((scalar *)&fl)->i >= 0; fl = *++_i93)
   val(fl,1,0,0) = (fine(fl,2,0,0) + fine(fl,2,1,0) +
     fine(fl,2,0,1) + fine(fl,2,1,1))/4.;

    } } end_foreach_halo(); }
  pfree (fluxl,__func__,__FILE__,__LINE__);





  if (cfl > 0.5 + 1e-6)
    fprintf (ferr,
      "WARNING: CFL must be <= 0.5 for VOF (cfl - 0.5 = %g)\n",
      cfl - 0.5), fflush (ferr);
#line 265 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
   { 
if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 265
foreach(){

#line 265 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
 {
    val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,1,0,0) + val(cc,0,0,0)*(val(uf.x,1,0,0) - val(uf.x,0,0,0)))/(val_cm(cm,0,0,0)*Delta);
    scalar t, tflux;
    scalar * _i11 = tracers; scalar * _i12 = tfluxl; if (tracers) for (t = *tracers, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i11, tflux = *++_i12)
      val(t,0,0,0) += dt*(val(tflux,0,0,0) - val(tflux,1,0,0))/(val_cm(cm,0,0,0)*Delta);
  } } end_foreach(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 265
foreach(){

#line 265 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
 {
    val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,1,0,0) + val(cc,0,0,0)*(val(uf.x,1,0,0) - val(uf.x,0,0,0)))/(val_cm(cm,0,0,0)*Delta);
    scalar t, tflux;
    scalar * _i11 = tracers; scalar * _i12 = tfluxl; if (tracers) for (t = *tracers, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i11, tflux = *++_i12)
      val(t,0,0,0) += dt*(val(tflux,0,0,0) - val(tflux,1,0,0))/(val_cm(cm,0,0,0)*Delta);
  } } end_foreach(); } }
  boundary (((scalar []){c,{-1}}));
  boundary (tracers);

  delete (tfluxl); pfree (tfluxl,__func__,__FILE__,__LINE__);
 delete (((scalar []){flux,alpha,n.x,n.y,n.z,{-1}})); }
#line 80

static void sweep_y (scalar c, scalar cc)
{
  vector n= new_vector("n");
  scalar alpha= new_scalar("alpha"), flux= new_scalar("flux");
  double cfl = 0.;
#line 94 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
  scalar * tracers = _attribute[c.i].tracers, * gfl = NULL, * tfluxl = NULL;
  if (tracers) {
    if (tracers) for (scalar t = *tracers, *_i91 = tracers; ((scalar *)&t)->i >= 0; t = *++_i91) {
      scalar gf = new_scalar("gf"), flux = new_scalar("flux");
      gfl = list_append (gfl, gf);
      tfluxl = list_append (tfluxl, flux);
    }







     { foreach(){

#line 108 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
 {
      scalar t, gf;
      scalar * _i6 = tracers; scalar * _i7 = gfl; if (tracers) for (t = *tracers, gf = *gfl; ((scalar *)&t)->i >= 0; t = *++_i6, gf = *++_i7) {
 double cl = val(c,0,-1,0), cc = val(c,0,0,0), cr = val(c,0,1,0);
 if (_attribute[t.i].inverse)
   cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;
 val(gf,0,0,0) = 0.;
 static const double cmin = 0.5;
 if (cc >= cmin) {
   if (cr >= cmin) {
     if (cl >= cmin) {
       if (_attribute[t.i].gradient)
  val(gf,0,0,0) = _attribute[t.i].gradient (val(t,0,-1,0)/cl, val(t,0,0,0)/cc, val(t,0,1,0)/cr)/Delta;
       else
  val(gf,0,0,0) = (val(t,0,1,0)/cr - val(t,0,-1,0)/cl)/(2.*Delta);
     }
     else
        val(gf,0,0,0) = (val(t,0,1,0)/cr - val(t,0,0,0)/cc)/Delta;
   }
   else if (cl >= cmin)
     val(gf,0,0,0) = (val(t,0,0,0)/cc - val(t,0,-1,0)/cl)/Delta;
 }
      }
    } } end_foreach(); }
    boundary (gfl);
  }






  reconstruction (c, n, alpha);

   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _cfl = cfl; 
#line 142

if (!is_constant(fm.y) && !is_constant(cm)) {
#undef val_fm_y
#define val_fm_y(a,k,i,j) val(a,k,i,j)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,k,i,j)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,k,i,j)
#undef val_fm_z
#define val_fm_z(a,k,i,j) val(a,k,i,j)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,k,i,j)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,k,i,j)
#undef val_fm_x
#define val_fm_x(a,k,i,j) val(a,k,i,j)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,k,i,j)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,k,i,j)
#undef val_cm
#define val_cm(a,k,i,j) val(a,k,i,j)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,k,i,j)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,k,i,j)
#line 142
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_y()) {
#line 142
{

#line 142 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
 {






    double un = val(uf.y,0,0,0)*dt/(Delta*val_fm_y(fm.y,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_y(fm.y,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_y(fm.y,0,0,0)*s/val_cm(cm,0,0,0);
#line 169 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
    double cf = (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.) ? val(c,0,i,0) :
      rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.z,0,i,0), val(n.x,0,i,0)}, val(alpha,0,i,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.y,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; if (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = val(c,0,i,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,i,0)/ci + s*min(1., 1. - s*un)*val(gf,0,i,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.y,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_generic()
#line 197
 end_foreach_face(); }
if (is_constant(fm.y) && !is_constant(cm)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.y.i -_NVARMAX], _constant[fm.z.i - _NVARMAX], _constant[fm.x.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_y
#define val_fm_y(a,k,i,j) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,k,i,j) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_fm_x
#define val_fm_x(a,k,i,j) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_cm
#define val_cm(a,k,i,j) val(a,k,i,j)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,k,i,j)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,k,i,j)
#line 142
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_y()) {
#line 142
{

#line 142 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
 {






    double un = val(uf.y,0,0,0)*dt/(Delta*val_fm_y(fm.y,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_y(fm.y,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_y(fm.y,0,0,0)*s/val_cm(cm,0,0,0);
#line 169 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
    double cf = (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.) ? val(c,0,i,0) :
      rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.z,0,i,0), val(n.x,0,i,0)}, val(alpha,0,i,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.y,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; if (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = val(c,0,i,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,i,0)/ci + s*min(1., 1. - s*un)*val(gf,0,i,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.y,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_generic()
#line 197
 end_foreach_face(); }
if (!is_constant(fm.y) && is_constant(cm)) {
#undef val_fm_y
#define val_fm_y(a,k,i,j) val(a,k,i,j)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,k,i,j)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,k,i,j)
#undef val_fm_z
#define val_fm_z(a,k,i,j) val(a,k,i,j)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,k,i,j)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,k,i,j)
#undef val_fm_x
#define val_fm_x(a,k,i,j) val(a,k,i,j)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,k,i,j)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,k,i,j)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,k,i,j) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 142
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_y()) {
#line 142
{

#line 142 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
 {






    double un = val(uf.y,0,0,0)*dt/(Delta*val_fm_y(fm.y,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_y(fm.y,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_y(fm.y,0,0,0)*s/val_cm(cm,0,0,0);
#line 169 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
    double cf = (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.) ? val(c,0,i,0) :
      rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.z,0,i,0), val(n.x,0,i,0)}, val(alpha,0,i,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.y,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; if (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = val(c,0,i,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,i,0)/ci + s*min(1., 1. - s*un)*val(gf,0,i,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.y,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_generic()
#line 197
 end_foreach_face(); }
if (is_constant(fm.y) && is_constant(cm)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.y.i -_NVARMAX], _constant[fm.z.i - _NVARMAX], _constant[fm.x.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_y
#define val_fm_y(a,k,i,j) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,k,i,j) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_fm_x
#define val_fm_x(a,k,i,j) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,k,i,j) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 142
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_y()) {
#line 142
{

#line 142 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
 {






    double un = val(uf.y,0,0,0)*dt/(Delta*val_fm_y(fm.y,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_y(fm.y,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_y(fm.y,0,0,0)*s/val_cm(cm,0,0,0);
#line 169 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
    double cf = (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.) ? val(c,0,i,0) :
      rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.z,0,i,0), val(n.x,0,i,0)}, val(alpha,0,i,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.y,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; if (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = val(c,0,i,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,i,0)/ci + s*min(1., 1. - s*un)*val(gf,0,i,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.y,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_generic()
#line 197
 end_foreach_face(); }OMP(omp critical) if (_cfl > cfl) cfl = _cfl;
mpi_all_reduce_double (cfl, MPI_MAX);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 197
 }
  delete (gfl); pfree (gfl,__func__,__FILE__,__LINE__);
#line 208 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
  scalar * fluxl = list_concat (NULL, tfluxl);
  fluxl = list_append (fluxl, flux);
  for (int l = depth() - 1; l >= 0; l--)
     { foreach_halo (prolongation, l){

#line 211 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
 {
#line 227 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
      if ((!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0))
 if (fluxl) for (scalar fl = *fluxl, *_i92 = fluxl; ((scalar *)&fl)->i >= 0; fl = *++_i92)
   val(fl,0,0,0) = (fine(fl,0,0,0) + fine(fl,0,0,1) +
    fine(fl,1,0,0) + fine(fl,1,0,1))/4.;
      if ((!is_leaf (neighbor(0,1,0)) && neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0))
 if (fluxl) for (scalar fl = *fluxl, *_i93 = fluxl; ((scalar *)&fl)->i >= 0; fl = *++_i93)
   val(fl,0,1,0) = (fine(fl,0,2,0) + fine(fl,0,2,1) +
     fine(fl,1,2,0) + fine(fl,1,2,1))/4.;

    } } end_foreach_halo(); }
  pfree (fluxl,__func__,__FILE__,__LINE__);





  if (cfl > 0.5 + 1e-6)
    fprintf (ferr,
      "WARNING: CFL must be <= 0.5 for VOF (cfl - 0.5 = %g)\n",
      cfl - 0.5), fflush (ferr);
#line 265 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
   { 
if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,k,i,j) val(a,k,i,j)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,k,i,j)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,k,i,j)
#line 265
foreach(){

#line 265 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
 {
    val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,0,1,0) + val(cc,0,0,0)*(val(uf.y,0,1,0) - val(uf.y,0,0,0)))/(val_cm(cm,0,0,0)*Delta);
    scalar t, tflux;
    scalar * _i11 = tracers; scalar * _i12 = tfluxl; if (tracers) for (t = *tracers, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i11, tflux = *++_i12)
      val(t,0,0,0) += dt*(val(tflux,0,0,0) - val(tflux,0,1,0))/(val_cm(cm,0,0,0)*Delta);
  } } end_foreach(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,k,i,j) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 265
foreach(){

#line 265 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
 {
    val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,0,1,0) + val(cc,0,0,0)*(val(uf.y,0,1,0) - val(uf.y,0,0,0)))/(val_cm(cm,0,0,0)*Delta);
    scalar t, tflux;
    scalar * _i11 = tracers; scalar * _i12 = tfluxl; if (tracers) for (t = *tracers, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i11, tflux = *++_i12)
      val(t,0,0,0) += dt*(val(tflux,0,0,0) - val(tflux,0,1,0))/(val_cm(cm,0,0,0)*Delta);
  } } end_foreach(); } }
  boundary (((scalar []){c,{-1}}));
  boundary (tracers);

  delete (tfluxl); pfree (tfluxl,__func__,__FILE__,__LINE__);
 delete (((scalar []){flux,alpha,n.x,n.y,n.z,{-1}})); }
#line 80

static void sweep_z (scalar c, scalar cc)
{
  vector n= new_vector("n");
  scalar alpha= new_scalar("alpha"), flux= new_scalar("flux");
  double cfl = 0.;
#line 94 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
  scalar * tracers = _attribute[c.i].tracers, * gfl = NULL, * tfluxl = NULL;
  if (tracers) {
    if (tracers) for (scalar t = *tracers, *_i91 = tracers; ((scalar *)&t)->i >= 0; t = *++_i91) {
      scalar gf = new_scalar("gf"), flux = new_scalar("flux");
      gfl = list_append (gfl, gf);
      tfluxl = list_append (tfluxl, flux);
    }







     { foreach(){

#line 108 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
 {
      scalar t, gf;
      scalar * _i6 = tracers; scalar * _i7 = gfl; if (tracers) for (t = *tracers, gf = *gfl; ((scalar *)&t)->i >= 0; t = *++_i6, gf = *++_i7) {
 double cl = val(c,0,0,-1), cc = val(c,0,0,0), cr = val(c,0,0,1);
 if (_attribute[t.i].inverse)
   cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;
 val(gf,0,0,0) = 0.;
 static const double cmin = 0.5;
 if (cc >= cmin) {
   if (cr >= cmin) {
     if (cl >= cmin) {
       if (_attribute[t.i].gradient)
  val(gf,0,0,0) = _attribute[t.i].gradient (val(t,0,0,-1)/cl, val(t,0,0,0)/cc, val(t,0,0,1)/cr)/Delta;
       else
  val(gf,0,0,0) = (val(t,0,0,1)/cr - val(t,0,0,-1)/cl)/(2.*Delta);
     }
     else
        val(gf,0,0,0) = (val(t,0,0,1)/cr - val(t,0,0,0)/cc)/Delta;
   }
   else if (cl >= cmin)
     val(gf,0,0,0) = (val(t,0,0,0)/cc - val(t,0,0,-1)/cl)/Delta;
 }
      }
    } } end_foreach(); }
    boundary (gfl);
  }






  reconstruction (c, n, alpha);

   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _cfl = cfl; 
#line 142

if (!is_constant(fm.z) && !is_constant(cm)) {
#undef val_fm_z
#define val_fm_z(a,j,k,i) val(a,j,k,i)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,j,k,i)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,j,k,i)
#undef val_fm_x
#define val_fm_x(a,j,k,i) val(a,j,k,i)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,j,k,i)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,j,k,i)
#undef val_fm_y
#define val_fm_y(a,j,k,i) val(a,j,k,i)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,j,k,i)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,j,k,i)
#undef val_cm
#define val_cm(a,j,k,i) val(a,j,k,i)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,j,k,i)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,j,k,i)
#line 142
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_z()) {
#line 142
{

#line 142 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
 {






    double un = val(uf.z,0,0,0)*dt/(Delta*val_fm_z(fm.z,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_z(fm.z,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_z(fm.z,0,0,0)*s/val_cm(cm,0,0,0);
#line 169 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
    double cf = (val(c,0,0,i) <= 0. || val(c,0,0,i) >= 1.) ? val(c,0,0,i) :
      rectangle_fraction ((coord){-s*val(n.z,0,0,i), val(n.x,0,0,i), val(n.y,0,0,i)}, val(alpha,0,0,i),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.z,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; if (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = val(c,0,0,i);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,0,i)/ci + s*min(1., 1. - s*un)*val(gf,0,0,i)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.z,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_generic()
#line 197
 end_foreach_face(); }
if (is_constant(fm.z) && !is_constant(cm)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.z.i -_NVARMAX], _constant[fm.x.i - _NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_z
#define val_fm_z(a,j,k,i) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_fm_x
#define val_fm_x(a,j,k,i) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,j,k,i) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_cm
#define val_cm(a,j,k,i) val(a,j,k,i)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,j,k,i)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,j,k,i)
#line 142
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_z()) {
#line 142
{

#line 142 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
 {






    double un = val(uf.z,0,0,0)*dt/(Delta*val_fm_z(fm.z,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_z(fm.z,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_z(fm.z,0,0,0)*s/val_cm(cm,0,0,0);
#line 169 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
    double cf = (val(c,0,0,i) <= 0. || val(c,0,0,i) >= 1.) ? val(c,0,0,i) :
      rectangle_fraction ((coord){-s*val(n.z,0,0,i), val(n.x,0,0,i), val(n.y,0,0,i)}, val(alpha,0,0,i),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.z,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; if (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = val(c,0,0,i);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,0,i)/ci + s*min(1., 1. - s*un)*val(gf,0,0,i)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.z,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_generic()
#line 197
 end_foreach_face(); }
if (!is_constant(fm.z) && is_constant(cm)) {
#undef val_fm_z
#define val_fm_z(a,j,k,i) val(a,j,k,i)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,j,k,i)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,j,k,i)
#undef val_fm_x
#define val_fm_x(a,j,k,i) val(a,j,k,i)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,j,k,i)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,j,k,i)
#undef val_fm_y
#define val_fm_y(a,j,k,i) val(a,j,k,i)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,j,k,i)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,j,k,i)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,j,k,i) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 142
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_z()) {
#line 142
{

#line 142 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
 {






    double un = val(uf.z,0,0,0)*dt/(Delta*val_fm_z(fm.z,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_z(fm.z,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_z(fm.z,0,0,0)*s/val_cm(cm,0,0,0);
#line 169 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
    double cf = (val(c,0,0,i) <= 0. || val(c,0,0,i) >= 1.) ? val(c,0,0,i) :
      rectangle_fraction ((coord){-s*val(n.z,0,0,i), val(n.x,0,0,i), val(n.y,0,0,i)}, val(alpha,0,0,i),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.z,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; if (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = val(c,0,0,i);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,0,i)/ci + s*min(1., 1. - s*un)*val(gf,0,0,i)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.z,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_generic()
#line 197
 end_foreach_face(); }
if (is_constant(fm.z) && is_constant(cm)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.z.i -_NVARMAX], _constant[fm.x.i - _NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_z
#define val_fm_z(a,j,k,i) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_fm_x
#define val_fm_x(a,j,k,i) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,j,k,i) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,j,k,i) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 142
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_z()) {
#line 142
{

#line 142 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
 {






    double un = val(uf.z,0,0,0)*dt/(Delta*val_fm_z(fm.z,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_z(fm.z,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_z(fm.z,0,0,0)*s/val_cm(cm,0,0,0);
#line 169 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
    double cf = (val(c,0,0,i) <= 0. || val(c,0,0,i) >= 1.) ? val(c,0,0,i) :
      rectangle_fraction ((coord){-s*val(n.z,0,0,i), val(n.x,0,0,i), val(n.y,0,0,i)}, val(alpha,0,0,i),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.z,0,0,0);






    scalar t, gf, tflux;
    scalar * _i8 = tracers; scalar * _i9 = gfl; scalar * _i10 = tfluxl; if (tracers) for (t = *tracers, gf = *gfl, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i8, gf = *++_i9, tflux = *++_i10) {
      double cf1 = cf, ci = val(c,0,0,i);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,0,i)/ci + s*min(1., 1. - s*un)*val(gf,0,0,i)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.z,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }
  } }  }}  end_foreach_face_generic()
#line 197
 end_foreach_face(); }OMP(omp critical) if (_cfl > cfl) cfl = _cfl;
mpi_all_reduce_double (cfl, MPI_MAX);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 197
 }
  delete (gfl); pfree (gfl,__func__,__FILE__,__LINE__);
#line 208 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
  scalar * fluxl = list_concat (NULL, tfluxl);
  fluxl = list_append (fluxl, flux);
  for (int l = depth() - 1; l >= 0; l--)
     { foreach_halo (prolongation, l){

#line 211 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
 {
#line 227 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
      if ((!is_leaf (neighbor(0,0,-1)) && neighbor(0,0,-1).neighbors && neighbor(0,0,-1).pid >= 0))
 if (fluxl) for (scalar fl = *fluxl, *_i92 = fluxl; ((scalar *)&fl)->i >= 0; fl = *++_i92)
   val(fl,0,0,0) = (fine(fl,0,0,0) + fine(fl,1,0,0) +
    fine(fl,0,1,0) + fine(fl,1,1,0))/4.;
      if ((!is_leaf (neighbor(0,0,1)) && neighbor(0,0,1).neighbors && neighbor(0,0,1).pid >= 0))
 if (fluxl) for (scalar fl = *fluxl, *_i93 = fluxl; ((scalar *)&fl)->i >= 0; fl = *++_i93)
   val(fl,0,0,1) = (fine(fl,0,0,2) + fine(fl,1,0,2) +
     fine(fl,0,1,2) + fine(fl,1,1,2))/4.;

    } } end_foreach_halo(); }
  pfree (fluxl,__func__,__FILE__,__LINE__);





  if (cfl > 0.5 + 1e-6)
    fprintf (ferr,
      "WARNING: CFL must be <= 0.5 for VOF (cfl - 0.5 = %g)\n",
      cfl - 0.5), fflush (ferr);
#line 265 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
   { 
if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,j,k,i) val(a,j,k,i)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,j,k,i)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,j,k,i)
#line 265
foreach(){

#line 265 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
 {
    val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,0,0,1) + val(cc,0,0,0)*(val(uf.z,0,0,1) - val(uf.z,0,0,0)))/(val_cm(cm,0,0,0)*Delta);
    scalar t, tflux;
    scalar * _i11 = tracers; scalar * _i12 = tfluxl; if (tracers) for (t = *tracers, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i11, tflux = *++_i12)
      val(t,0,0,0) += dt*(val(tflux,0,0,0) - val(tflux,0,0,1))/(val_cm(cm,0,0,0)*Delta);
  } } end_foreach(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,j,k,i) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 265
foreach(){

#line 265 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
 {
    val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,0,0,1) + val(cc,0,0,0)*(val(uf.z,0,0,1) - val(uf.z,0,0,0)))/(val_cm(cm,0,0,0)*Delta);
    scalar t, tflux;
    scalar * _i11 = tracers; scalar * _i12 = tfluxl; if (tracers) for (t = *tracers, tflux = *tfluxl; ((scalar *)&t)->i >= 0; t = *++_i11, tflux = *++_i12)
      val(t,0,0,0) += dt*(val(tflux,0,0,0) - val(tflux,0,0,1))/(val_cm(cm,0,0,0)*Delta);
  } } end_foreach(); } }
  boundary (((scalar []){c,{-1}}));
  boundary (tracers);

  delete (tfluxl); pfree (tfluxl,__func__,__FILE__,__LINE__);
 delete (((scalar []){flux,alpha,n.x,n.y,n.z,{-1}})); }






void vof_advection (scalar * interfaces, int i)
{
  if (interfaces) for (scalar c = *interfaces, *_i94 = interfaces; ((scalar *)&c)->i >= 0; c = *++_i94) {
#line 294 "/home/miro3/Documents/Programming/basilisk/src/vof.h"
    scalar cc= new_scalar("cc");
     { foreach(){

#line 295 "/home/miro3/Documents/Programming/basilisk/src/vof.h"

      val(cc,0,0,0) = (val(c,0,0,0) > 0.5); } end_foreach(); }






    void (* sweep[3]) (scalar, scalar);
    int d = 0;
    {
#line 305

      sweep[d++] = sweep_x;
#line 305

      sweep[d++] = sweep_y;
#line 305

      sweep[d++] = sweep_z;}
    boundary (((scalar []){c,{-1}}));
    for (d = 0; d < 3; d++)
      sweep[(i + d) % 3] (c, cc);
   delete (((scalar []){cc,{-1}})); }
}

static int vof_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int vof_0 (const int i, const double t, Event * _ev) { trace ("vof_0", "/home/miro3/Documents/Programming/basilisk/src/vof.h", 313); 
  vof_advection (interfaces, i); end_trace("vof_0", "/home/miro3/Documents/Programming/basilisk/src/vof.h", 314);  return 0; } 
#line 14 "/home/miro3/Documents/Programming/basilisk/src/two-phase.h"

scalar f= {11}, * interfaces = ((scalar []){{11},{-1}});
double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0.;





vector alphav= {{12},{13},{14}};
scalar rhov= {15};

static int defaults_1_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int defaults_1 (const int i, const double t, Event * _ev) { trace ("defaults_1", "/home/miro3/Documents/Programming/basilisk/src/two-phase.h", 25);  {
  alpha = alphav;
  rho = rhov;





  if (mu1 || mu2)
    mu = new_face_vector("mu");
 end_trace("defaults_1", "/home/miro3/Documents/Programming/basilisk/src/two-phase.h", 35); } return 0; } 
#line 59 "/home/miro3/Documents/Programming/basilisk/src/two-phase.h"
static int properties_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int properties_0 (const int i, const double t, Event * _ev) { trace ("properties_0", "/home/miro3/Documents/Programming/basilisk/src/two-phase.h", 59);  {
#line 84 "/home/miro3/Documents/Programming/basilisk/src/two-phase.h"
  _attribute[f.i].prolongation = refine_bilinear;
  boundary (((scalar []){f,{-1}}));


   { 
if (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 88
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 88
{

#line 88 "/home/miro3/Documents/Programming/basilisk/src/two-phase.h"
 {
    double ff = (val(f,0,0,0) + val(f,-1,0,0))/2.;
    val(alphav.x,0,0,0) = val_fm_x(fm.x,0,0,0)/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    if (mu1 || mu2) {
      vector muv = mu;
      val(muv.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 88
{

#line 88 "/home/miro3/Documents/Programming/basilisk/src/two-phase.h"
 {
    double ff = (val(f,0,0,0) + val(f,0,-1,0))/2.;
    val(alphav.y,0,0,0) = val_fm_y(fm.y,0,0,0)/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    if (mu1 || mu2) {
      vector muv = mu;
      val(muv.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 88
{

#line 88 "/home/miro3/Documents/Programming/basilisk/src/two-phase.h"
 {
    double ff = (val(f,0,0,0) + val(f,0,0,-1))/2.;
    val(alphav.z,0,0,0) = val_fm_z(fm.z,0,0,0)/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    if (mu1 || mu2) {
      vector muv = mu;
      val(muv.z,0,0,0) = val_fm_z(fm.z,0,0,0)*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  } }  }}  end_foreach_face_generic()
#line 95
 end_foreach_face(); }
if (is_constant(fm.x)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 88
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 88
{

#line 88 "/home/miro3/Documents/Programming/basilisk/src/two-phase.h"
 {
    double ff = (val(f,0,0,0) + val(f,-1,0,0))/2.;
    val(alphav.x,0,0,0) = val_fm_x(fm.x,0,0,0)/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    if (mu1 || mu2) {
      vector muv = mu;
      val(muv.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 88
{

#line 88 "/home/miro3/Documents/Programming/basilisk/src/two-phase.h"
 {
    double ff = (val(f,0,0,0) + val(f,0,-1,0))/2.;
    val(alphav.y,0,0,0) = val_fm_y(fm.y,0,0,0)/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    if (mu1 || mu2) {
      vector muv = mu;
      val(muv.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 88
{

#line 88 "/home/miro3/Documents/Programming/basilisk/src/two-phase.h"
 {
    double ff = (val(f,0,0,0) + val(f,0,0,-1))/2.;
    val(alphav.z,0,0,0) = val_fm_z(fm.z,0,0,0)/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    if (mu1 || mu2) {
      vector muv = mu;
      val(muv.z,0,0,0) = val_fm_z(fm.z,0,0,0)*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  } }  }}  end_foreach_face_generic()
#line 95
 end_foreach_face(); } }
   { 
if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 96
foreach(){

#line 96 "/home/miro3/Documents/Programming/basilisk/src/two-phase.h"

    val(rhov,0,0,0) = val_cm(cm,0,0,0)*(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2); } end_foreach(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 96
foreach(){

#line 96 "/home/miro3/Documents/Programming/basilisk/src/two-phase.h"

    val(rhov,0,0,0) = val_cm(cm,0,0,0)*(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2); } end_foreach(); } }


  _attribute[f.i].prolongation = fraction_refine;
  boundary (((scalar []){f,{-1}}));

 end_trace("properties_0", "/home/miro3/Documents/Programming/basilisk/src/two-phase.h", 103); } return 0; } 
#line 6 "milk_crown.c"
#line 1 "navier-stokes/conserving.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/conserving.h"
#line 14 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/conserving.h"
static void momentum_refine (Point point, scalar u) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 14 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/conserving.h"

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 14

  refine_bilinear (point, u);
  double rhou = 0.;
   { foreach_child()
    rhou += val_cm(cm,0,0,0)*(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2)*val(u,0,0,0); end_foreach_child(); }
  double du = val(u,0,0,0) - rhou/((1 << 3)*val_cm(cm,0,0,0)*(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2));
   { foreach_child()
    val(u,0,0,0) += du; end_foreach_child(); }
 }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 14

  refine_bilinear (point, u);
  double rhou = 0.;
   { foreach_child()
    rhou += val_cm(cm,0,0,0)*(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2)*val(u,0,0,0); end_foreach_child(); }
  double du = val(u,0,0,0) - rhou/((1 << 3)*val_cm(cm,0,0,0)*(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2));
   { foreach_child()
    val(u,0,0,0) += du; end_foreach_child(); }
 }}

static void momentum_restriction (Point point, scalar u)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 25 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/conserving.h"

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 25

  double rhou = 0.;
   { foreach_child()
    rhou += val_cm(cm,0,0,0)*(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2)*val(u,0,0,0); end_foreach_child(); }
  val(u,0,0,0) = rhou/((1 << 3)*val_cm(cm,0,0,0)*(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2));
 }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 25

  double rhou = 0.;
   { foreach_child()
    rhou += val_cm(cm,0,0,0)*(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2)*val(u,0,0,0); end_foreach_child(); }
  val(u,0,0,0) = rhou/((1 << 3)*val_cm(cm,0,0,0)*(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2));
 }}






static int defaults_2_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int defaults_2 (const int i, const double t, Event * _ev) { trace ("defaults_2", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/conserving.h", 37);  {
  stokes = true;
#line 48 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/conserving.h"
  assert (all[f.i].i == f.i);
  for (int i = f.i; i >= 2; i--)
    all[i] = all[i-1];
  all[1] = f;





  {
#line 57
 {
    _attribute[u.x.i].refine = _attribute[u.x.i].prolongation = momentum_refine;
    _attribute[u.x.i].restriction = momentum_restriction;
  }
#line 57
 {
    _attribute[u.y.i].refine = _attribute[u.y.i].prolongation = momentum_refine;
    _attribute[u.y.i].restriction = momentum_restriction;
  }
#line 57
 {
    _attribute[u.z.i].refine = _attribute[u.z.i].prolongation = momentum_refine;
    _attribute[u.z.i].restriction = momentum_restriction;
  }}

 end_trace("defaults_2", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/conserving.h", 62); } return 0; } 








#line 70

static double boundary_q1_x (Point neighbor, Point point, scalar q1)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 72 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/conserving.h"

  return clamp(val(f,0,0,0),0.,1.)*rho1*val(u.x,0,0,0);
}
#line 70

static double boundary_q1_y (Point neighbor, Point point, scalar q1)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 72 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/conserving.h"

  return clamp(val(f,0,0,0),0.,1.)*rho1*val(u.y,0,0,0);
}
#line 70

static double boundary_q1_z (Point neighbor, Point point, scalar q1)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 72 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/conserving.h"

  return clamp(val(f,0,0,0),0.,1.)*rho1*val(u.z,0,0,0);
}


#line 76

static double boundary_q2_x (Point neighbor, Point point, scalar q2)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 78 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/conserving.h"

  return (1. - clamp(val(f,0,0,0),0.,1.))*rho2*val(u.x,0,0,0);
}
#line 76

static double boundary_q2_y (Point neighbor, Point point, scalar q2)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 78 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/conserving.h"

  return (1. - clamp(val(f,0,0,0),0.,1.))*rho2*val(u.y,0,0,0);
}
#line 76

static double boundary_q2_z (Point neighbor, Point point, scalar q2)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 78 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/conserving.h"

  return (1. - clamp(val(f,0,0,0),0.,1.))*rho2*val(u.z,0,0,0);
}







#line 87

static void prolongation_q1_x (Point point, scalar q1) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 88 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/conserving.h"

   { foreach_child()
    val(q1,0,0,0) = clamp(val(f,0,0,0),0.,1.)*rho1*val(u.x,0,0,0); end_foreach_child(); }
}
#line 87

static void prolongation_q1_y (Point point, scalar q1) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 88 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/conserving.h"

   { foreach_child()
    val(q1,0,0,0) = clamp(val(f,0,0,0),0.,1.)*rho1*val(u.y,0,0,0); end_foreach_child(); }
}
#line 87

static void prolongation_q1_z (Point point, scalar q1) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 88 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/conserving.h"

   { foreach_child()
    val(q1,0,0,0) = clamp(val(f,0,0,0),0.,1.)*rho1*val(u.z,0,0,0); end_foreach_child(); }
}


#line 93

static void prolongation_q2_x (Point point, scalar q2) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 94 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/conserving.h"

   { foreach_child()
    val(q2,0,0,0) = (1. - clamp(val(f,0,0,0),0.,1.))*rho2*val(u.x,0,0,0); end_foreach_child(); }
}
#line 93

static void prolongation_q2_y (Point point, scalar q2) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 94 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/conserving.h"

   { foreach_child()
    val(q2,0,0,0) = (1. - clamp(val(f,0,0,0),0.,1.))*rho2*val(u.y,0,0,0); end_foreach_child(); }
}
#line 93

static void prolongation_q2_z (Point point, scalar q2) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 94 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/conserving.h"

   { foreach_child()
    val(q2,0,0,0) = (1. - clamp(val(f,0,0,0),0.,1.))*rho2*val(u.z,0,0,0); end_foreach_child(); }
}






static scalar * interfaces1 = NULL;

static int vof_1_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int vof_1 (const int i, const double t, Event * _ev) { trace ("vof_1", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/conserving.h", 106);  {






  vector q1= new_vector("q1"), q2= new_vector("q2");
  for (scalar s = *((scalar []){q1.x,q1.y,q1.z,q2.x,q2.y,q2.z,{-1}}), *_i95 = ((scalar []){q1.x,q1.y,q1.z,q2.x,q2.y,q2.z,{-1}}); ((scalar *)&s)->i >= 0; s = *++_i95)
    {
#line 115

      _attribute[s.i].v.x.i = -1;
#line 115

      _attribute[s.i].v.y.i = -1;
#line 115

      _attribute[s.i].v.z.i = -1;}
  for (int i = 0; i < nboundary; i++)
    {
#line 118
 {
      _attribute[q1.x.i].boundary[i] = boundary_q1_x;
      _attribute[q2.x.i].boundary[i] = boundary_q2_x;
    }
#line 118
 {
      _attribute[q1.y.i].boundary[i] = boundary_q1_y;
      _attribute[q2.y.i].boundary[i] = boundary_q2_y;
    }
#line 118
 {
      _attribute[q1.z.i].boundary[i] = boundary_q1_z;
      _attribute[q2.z.i].boundary[i] = boundary_q2_z;
    }}

  {
#line 123
 {
    _attribute[q1.x.i].prolongation = prolongation_q1_x;
    _attribute[q2.x.i].prolongation = prolongation_q2_x;
  }
#line 123
 {
    _attribute[q1.y.i].prolongation = prolongation_q1_y;
    _attribute[q2.y.i].prolongation = prolongation_q2_y;
  }
#line 123
 {
    _attribute[q1.z.i].prolongation = prolongation_q1_z;
    _attribute[q2.z.i].prolongation = prolongation_q2_z;
  }}






   { foreach(){

#line 133 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/conserving.h"

    {
#line 134
 {
      double fc = clamp(val(f,0,0,0),0,1);
      val(q1.x,0,0,0) = fc*rho1*val(u.x,0,0,0);
      val(q2.x,0,0,0) = (1. - fc)*rho2*val(u.x,0,0,0);
    }
#line 134
 {
      double fc = clamp(val(f,0,0,0),0,1);
      val(q1.y,0,0,0) = fc*rho1*val(u.y,0,0,0);
      val(q2.y,0,0,0) = (1. - fc)*rho2*val(u.y,0,0,0);
    }
#line 134
 {
      double fc = clamp(val(f,0,0,0),0,1);
      val(q1.z,0,0,0) = fc*rho1*val(u.z,0,0,0);
      val(q2.z,0,0,0) = (1. - fc)*rho2*val(u.z,0,0,0);
    }} } end_foreach(); }
  boundary ((scalar *)((vector []){{q1.x,q1.y,q1.z},{q2.x,q2.y,q2.z},{{-1},{-1},{-1}}}));






  {
#line 146
 {
    _attribute[q2.x.i].inverse = true;
    _attribute[q1.x.i].gradient = _attribute[q2.x.i].gradient = _attribute[u.x.i].gradient;
  }
#line 146
 {
    _attribute[q2.y.i].inverse = true;
    _attribute[q1.y.i].gradient = _attribute[q2.y.i].gradient = _attribute[u.y.i].gradient;
  }
#line 146
 {
    _attribute[q2.z.i].inverse = true;
    _attribute[q1.z.i].gradient = _attribute[q2.z.i].gradient = _attribute[u.z.i].gradient;
  }}





  _attribute[f.i].tracers = (scalar *)((vector []){{q1.x,q1.y,q1.z},{q2.x,q2.y,q2.z},{{-1},{-1},{-1}}});
  vof_advection (((scalar []){f,{-1}}), i);





   { foreach(){

#line 162 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/conserving.h"

    {
#line 163

      val(u.x,0,0,0) = (val(q1.x,0,0,0) + val(q2.x,0,0,0))/(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2);
#line 163

      val(u.y,0,0,0) = (val(q1.y,0,0,0) + val(q2.y,0,0,0))/(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2);
#line 163

      val(u.z,0,0,0) = (val(q1.z,0,0,0) + val(q2.z,0,0,0))/(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2);}; } end_foreach(); }
  boundary ((scalar *)((vector []){{u.x,u.y,u.z},{{-1},{-1},{-1}}}));





  interfaces1 = interfaces, interfaces = NULL;
 delete (((scalar []){q2.x,q2.y,q2.z,q1.x,q1.y,q1.z,{-1}}));  end_trace("vof_1", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/conserving.h", 172); } return 0; } 




static int tracer_advection_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int tracer_advection_0 (const int i, const double t, Event * _ev) { trace ("tracer_advection_0", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/conserving.h", 177);  {
  interfaces = interfaces1;
 end_trace("tracer_advection_0", "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/conserving.h", 179); } return 0; } 
#line 7 "milk_crown.c"
#line 1 "SGS.h"
#line 1 "./SGS.h"






#line 1 "tracer.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/tracer.h"
#line 15 "/home/miro3/Documents/Programming/basilisk/src/tracer.h"
extern scalar * tracers;
extern vector uf;
extern double dt;







static int defaults_3_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int defaults_3 (const int i, const double t, Event * _ev) { trace ("defaults_3", "/home/miro3/Documents/Programming/basilisk/src/tracer.h", 25);  {
  if (tracers) for (scalar s = *tracers, *_i96 = tracers; ((scalar *)&s)->i >= 0; s = *++_i96) {
    _attribute[s.i].refine = refine_linear;
    _attribute[s.i].restriction = restriction_volume_average;
  }
 end_trace("defaults_3", "/home/miro3/Documents/Programming/basilisk/src/tracer.h", 30); } return 0; } 





#line 1 "bcg.h"
#line 37 "/home/miro3/Documents/Programming/basilisk/src/tracer.h"

static int tracer_advection_1_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int tracer_advection_1 (const int i, const double t, Event * _ev) { trace ("tracer_advection_1", "/home/miro3/Documents/Programming/basilisk/src/tracer.h", 38);  {
  advection ((struct Advection){tracers, uf, dt});
 end_trace("tracer_advection_1", "/home/miro3/Documents/Programming/basilisk/src/tracer.h", 40); } return 0; } 




static int tracer_diffusion_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int tracer_diffusion_0 (const int i, const double t, Event * _ev) { trace ("tracer_diffusion_0", "/home/miro3/Documents/Programming/basilisk/src/tracer.h", 45); ; end_trace("tracer_diffusion_0", "/home/miro3/Documents/Programming/basilisk/src/tracer.h", 45);  return 0; } 
#line 8 "./SGS.h"

#line 1 "Vreman.h"
#line 1 "./Vreman.h"
#line 18 "./Vreman.h"
void eddyviscosity(double Cs, vector u, scalar Evis){
  double d1v1, d2v1, d3v1, d1v2, d2v2, d3v2, d1v3, d2v3, d3v3;
  double b11, b12, b13, b22, b23, b33;
  double abeta, bbeta;
   { foreach(){

#line 22 "./Vreman.h"
{
    d1v1=(val(u.x,1,0,0)-val(u.x,-1,0,0))/2/Delta;
    d2v1=(val(u.x,0,1,0)-val(u.x,0,-1,0))/2/Delta;
    d3v1=(val(u.x,0,0,1)-val(u.x,0,0,-1))/2/Delta;
    d1v2=(val(u.y,1,0,0)-val(u.y,-1,0,0))/2/Delta;
    d2v2=(val(u.y,0,1,0)-val(u.y,0,-1,0))/2/Delta;
    d3v2=(val(u.y,0,0,1)-val(u.y,0,0,-1))/2/Delta;
    d1v3=(val(u.z,1,0,0)-val(u.z,-1,0,0))/2/Delta;
    d2v3=(val(u.z,0,1,0)-val(u.z,0,-1,0))/2/Delta;
    d3v3=(val(u.z,0,0,1)-val(u.z,0,0,-1))/2/Delta;
    b11 = Delta*Delta*(d1v1*d1v1+d2v1*d2v1+d3v1*d3v1);
    b12 = Delta*Delta*(d1v1*d1v2+d2v1*d2v2+d3v1*d3v2);
    b13 = Delta*Delta*(d1v1*d1v3+d2v1*d2v3+d3v1*d3v3);
    b22 = Delta*Delta*(d1v2*d1v2+d2v2*d2v2+d3v2*d3v2);
    b23 = Delta*Delta*(d1v2*d1v3+d2v2*d2v3+d3v2*d3v3);
    b33 = Delta*Delta*(d1v3*d1v3+d2v3*d2v3+d3v3*d3v3);
    abeta = sq(d1v1)+sq(d2v1)+sq(d3v1)+
      sq(d1v2)+sq(d2v2)+sq(d3v2)+
      sq(d1v3)+sq(d2v3)+sq(d3v3);
    bbeta = b11*b22-sq(b12)+b11*b33-sq(b13)+b22*b33-sq(b23);

    double eddy = abeta>0 ? ( 2.5*Cs*Cs*sqrt(bbeta/(abeta))) : 0.;
    val(Evis,0,0,0) = max(eddy, 0);
  } } end_foreach(); }
}
#line 10 "./SGS.h"
#line 21 "./SGS.h"
vector mu_t= {{16},{17},{18}};
scalar Evis= {19};
double Csmag;
scalar * tracers;




static inline void Evisprol(Point point,scalar s){ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 29 "./SGS.h"

   { foreach_child()
    val(Evis,0,0,0)=bilinear(point,Evis)/4.; end_foreach_child(); }
}

static inline void Evisres(Point point,scalar s){ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 34 "./SGS.h"

  double sum = 0.;
   { foreach_child()
    sum += val(s,0,0,0); end_foreach_child(); }
  val(s,0,0,0) = sum/2.;
}




static int defaults_4_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i=0);   *ip = i; *tp = t;   return ret; } static int defaults_4 (const int i, const double t, Event * _ev) { trace ("defaults_4", "./SGS.h", 44); {
  if (3!=3)
    fprintf(qstdout(),"Warning %dD grid. The used formulations only make sense for 3D turbulence simulations\n",3);
  Csmag=0.12;
  _attribute[Evis.i].prolongation=Evisprol;
  _attribute[Evis.i].restriction=Evisres;






  _attribute[Evis.i].refine=no_restriction;
  _attribute[Evis.i].coarsen=no_restriction;
  {
#line 58
{
    _attribute[mu_t.x.i].coarsen=no_restriction;
  }
#line 58
{
    _attribute[mu_t.y.i].coarsen=no_restriction;
  }
#line 58
{
    _attribute[mu_t.z.i].coarsen=no_restriction;
  }}

 end_trace("defaults_4", "./SGS.h", 62); } return 0; } 
#line 86 "./SGS.h"
static int viscous_term_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int viscous_term_0 (const int i, const double t, Event * _ev) { trace ("viscous_term_0", "./SGS.h", 86); {

  correction (dt);

  eddyviscosity(Csmag,u,Evis);

  boundary(((scalar []){Evis,{-1}}));
   { 
if (!is_constant(rho)) {
#undef val_rho
#define val_rho(a,i,j,k) val(a,i,j,k)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,i,j,k)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,i,j,k)
#line 93
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 93
{

#line 93 "./SGS.h"
 {
    val(mu_t.x,0,0,0) = 0.5*(val_rho(rho,0,0,0)*val(Evis,0,0,0)+val_rho(rho,-1,0,0)*val(Evis,-1,0,0));
    } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 93
{

#line 93 "./SGS.h"
 {
    val(mu_t.y,0,0,0) = 0.5*(val_rho(rho,0,0,0)*val(Evis,0,0,0)+val_rho(rho,0,-1,0)*val(Evis,0,-1,0));
    } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 93
{

#line 93 "./SGS.h"
 {
    val(mu_t.z,0,0,0) = 0.5*(val_rho(rho,0,0,0)*val(Evis,0,0,0)+val_rho(rho,0,0,-1)*val(Evis,0,0,-1));
    } }  }}  end_foreach_face_generic()
#line 95
 end_foreach_face(); }
if (is_constant(rho)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,i,j,k) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
#line 93
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 93
{

#line 93 "./SGS.h"
 {
    val(mu_t.x,0,0,0) = 0.5*(val_rho(rho,0,0,0)*val(Evis,0,0,0)+val_rho(rho,-1,0,0)*val(Evis,-1,0,0));
    } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 93
{

#line 93 "./SGS.h"
 {
    val(mu_t.y,0,0,0) = 0.5*(val_rho(rho,0,0,0)*val(Evis,0,0,0)+val_rho(rho,0,-1,0)*val(Evis,0,-1,0));
    } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 93
{

#line 93 "./SGS.h"
 {
    val(mu_t.z,0,0,0) = 0.5*(val_rho(rho,0,0,0)*val(Evis,0,0,0)+val_rho(rho,0,0,-1)*val(Evis,0,0,-1));
    } }  }}  end_foreach_face_generic()
#line 95
 end_foreach_face(); } }
  correction (-dt);
  boundary_flux(((vector []){{mu_t.x,mu_t.y,mu_t.z},{{-1},{-1},{-1}}}));


  vector muv = mu;


   { 
if (!is_constant(mu.x)) {
#undef val_mu_x
#define val_mu_x(a,i,j,k) val(a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_y
#define val_mu_y(a,i,j,k) val(a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_z
#define val_mu_z(a,i,j,k) val(a,i,j,k)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,i,j,k)
#line 103
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 103
{

#line 103 "./SGS.h"

    val(muv.x,0,0,0) = val_mu_x(mu.x,0,0,0) + val(mu_t.x,0,0,0); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 103
{

#line 103 "./SGS.h"

    val(muv.y,0,0,0) = val_mu_y(mu.y,0,0,0) + val(mu_t.y,0,0,0); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 103
{

#line 103 "./SGS.h"

    val(muv.z,0,0,0) = val_mu_z(mu.z,0,0,0) + val(mu_t.z,0,0,0); }  }}  end_foreach_face_generic()
#line 104
 end_foreach_face(); }
if (is_constant(mu.x)) {
const struct { double x, y, z; } _const_mu = {_constant[mu.x.i -_NVARMAX], _constant[mu.y.i - _NVARMAX], _constant[mu.z.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_x
#define val_mu_x(a,i,j,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,i,j,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_z
#define val_mu_z(a,i,j,k) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#line 103
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 103
{

#line 103 "./SGS.h"

    val(muv.x,0,0,0) = val_mu_x(mu.x,0,0,0) + val(mu_t.x,0,0,0); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 103
{

#line 103 "./SGS.h"

    val(muv.y,0,0,0) = val_mu_y(mu.y,0,0,0) + val(mu_t.y,0,0,0); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 103
{

#line 103 "./SGS.h"

    val(muv.z,0,0,0) = val_mu_z(mu.z,0,0,0) + val(mu_t.z,0,0,0); }  }}  end_foreach_face_generic()
#line 104
 end_foreach_face(); } }
   end_trace("viscous_term_0", "./SGS.h", 105); } return 0; } 
#line 8 "milk_crown.c"



#line 1 "reduced.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/reduced.h"
#line 19 "/home/miro3/Documents/Programming/basilisk/src/reduced.h"
coord G = {0.,0.,0.}, Z = {0.,0.,0.};





#line 1 "iforce.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"
#line 20 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"










static int defaults_5_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int defaults_5 (const int i, const double t, Event * _ev) { trace ("defaults_5", "/home/miro3/Documents/Programming/basilisk/src/iforce.h", 30);  {
  if (is_constant(a.x)) {
    a = new_face_vector("a");
     { foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 33
{

#line 33 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"

      val(a.x,0,0,0) = 0.; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 33
{

#line 33 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"

      val(a.y,0,0,0) = 0.; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 33
{

#line 33 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"

      val(a.z,0,0,0) = 0.; }  }}  end_foreach_face_generic()
#line 34
 end_foreach_face(); }
    boundary ((scalar *)((vector []){{a.x,a.y,a.z},{{-1},{-1},{-1}}}));
  }
 end_trace("defaults_5", "/home/miro3/Documents/Programming/basilisk/src/iforce.h", 37); } return 0; } 






static int acceleration_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int acceleration_0 (const int i, const double t, Event * _ev) { trace ("acceleration_0", "/home/miro3/Documents/Programming/basilisk/src/iforce.h", 44); 
{





  scalar * list = NULL;
  if (interfaces) for (scalar f = *interfaces, *_i97 = interfaces; ((scalar *)&f)->i >= 0; f = *++_i97)
    if (_attribute[f.i].phi.i) {
      list = list_add (list, f);






       { foreach(){

#line 61 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"

 val(f,0,0,0) = clamp (val(f,0,0,0), 0., 1.); } end_foreach(); }
      boundary (((scalar []){f,{-1}}));
    }
#line 74 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"
  if (list) for (scalar f = *list, *_i98 = list; ((scalar *)&f)->i >= 0; f = *++_i98)
    _attribute[f.i].prolongation = _attribute[p.i].prolongation;
  boundary (list);
#line 87 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"
  vector ia = a;
   { 
if (!is_constant(alpha.x) && !is_constant(fm.x)) {
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 88
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 88
{

#line 88 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"

    if (list) for (scalar f = *list, *_i99 = list; ((scalar *)&f)->i >= 0; f = *++_i99)
      if (val(f,0,0,0) != val(f,-1,0,0)) {
#line 100 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < nodata && val(phi,-1,0,0) < nodata) ?
   (val(phi,0,0,0) + val(phi,-1,0,0))/2. :
   val(phi,0,0,0) < nodata ? val(phi,0,0,0) :
   val(phi,-1,0,0) < nodata ? val(phi,-1,0,0) :
   0.;

 val(ia.x,0,0,0) += val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*phif*(val(f,0,0,0) - val(f,-1,0,0))/Delta;
      } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 88
{

#line 88 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"

    if (list) for (scalar f = *list, *_i99 = list; ((scalar *)&f)->i >= 0; f = *++_i99)
      if (val(f,0,0,0) != val(f,0,-1,0)) {
#line 100 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < nodata && val(phi,0,-1,0) < nodata) ?
   (val(phi,0,0,0) + val(phi,0,-1,0))/2. :
   val(phi,0,0,0) < nodata ? val(phi,0,0,0) :
   val(phi,0,-1,0) < nodata ? val(phi,0,-1,0) :
   0.;

 val(ia.y,0,0,0) += val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*phif*(val(f,0,0,0) - val(f,0,-1,0))/Delta;
      } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 88
{

#line 88 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"

    if (list) for (scalar f = *list, *_i99 = list; ((scalar *)&f)->i >= 0; f = *++_i99)
      if (val(f,0,0,0) != val(f,0,0,-1)) {
#line 100 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < nodata && val(phi,0,0,-1) < nodata) ?
   (val(phi,0,0,0) + val(phi,0,0,-1))/2. :
   val(phi,0,0,0) < nodata ? val(phi,0,0,0) :
   val(phi,0,0,-1) < nodata ? val(phi,0,0,-1) :
   0.;

 val(ia.z,0,0,0) += val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*phif*(val(f,0,0,0) - val(f,0,0,-1))/Delta;
      } }  }}  end_foreach_face_generic()
#line 109
 end_foreach_face(); }
if (is_constant(alpha.x) && !is_constant(fm.x)) {
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 88
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 88
{

#line 88 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"

    if (list) for (scalar f = *list, *_i99 = list; ((scalar *)&f)->i >= 0; f = *++_i99)
      if (val(f,0,0,0) != val(f,-1,0,0)) {
#line 100 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < nodata && val(phi,-1,0,0) < nodata) ?
   (val(phi,0,0,0) + val(phi,-1,0,0))/2. :
   val(phi,0,0,0) < nodata ? val(phi,0,0,0) :
   val(phi,-1,0,0) < nodata ? val(phi,-1,0,0) :
   0.;

 val(ia.x,0,0,0) += val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*phif*(val(f,0,0,0) - val(f,-1,0,0))/Delta;
      } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 88
{

#line 88 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"

    if (list) for (scalar f = *list, *_i99 = list; ((scalar *)&f)->i >= 0; f = *++_i99)
      if (val(f,0,0,0) != val(f,0,-1,0)) {
#line 100 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < nodata && val(phi,0,-1,0) < nodata) ?
   (val(phi,0,0,0) + val(phi,0,-1,0))/2. :
   val(phi,0,0,0) < nodata ? val(phi,0,0,0) :
   val(phi,0,-1,0) < nodata ? val(phi,0,-1,0) :
   0.;

 val(ia.y,0,0,0) += val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*phif*(val(f,0,0,0) - val(f,0,-1,0))/Delta;
      } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 88
{

#line 88 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"

    if (list) for (scalar f = *list, *_i99 = list; ((scalar *)&f)->i >= 0; f = *++_i99)
      if (val(f,0,0,0) != val(f,0,0,-1)) {
#line 100 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < nodata && val(phi,0,0,-1) < nodata) ?
   (val(phi,0,0,0) + val(phi,0,0,-1))/2. :
   val(phi,0,0,0) < nodata ? val(phi,0,0,0) :
   val(phi,0,0,-1) < nodata ? val(phi,0,0,-1) :
   0.;

 val(ia.z,0,0,0) += val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*phif*(val(f,0,0,0) - val(f,0,0,-1))/Delta;
      } }  }}  end_foreach_face_generic()
#line 109
 end_foreach_face(); }
if (!is_constant(alpha.x) && is_constant(fm.x)) {
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 88
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 88
{

#line 88 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"

    if (list) for (scalar f = *list, *_i99 = list; ((scalar *)&f)->i >= 0; f = *++_i99)
      if (val(f,0,0,0) != val(f,-1,0,0)) {
#line 100 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < nodata && val(phi,-1,0,0) < nodata) ?
   (val(phi,0,0,0) + val(phi,-1,0,0))/2. :
   val(phi,0,0,0) < nodata ? val(phi,0,0,0) :
   val(phi,-1,0,0) < nodata ? val(phi,-1,0,0) :
   0.;

 val(ia.x,0,0,0) += val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*phif*(val(f,0,0,0) - val(f,-1,0,0))/Delta;
      } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 88
{

#line 88 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"

    if (list) for (scalar f = *list, *_i99 = list; ((scalar *)&f)->i >= 0; f = *++_i99)
      if (val(f,0,0,0) != val(f,0,-1,0)) {
#line 100 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < nodata && val(phi,0,-1,0) < nodata) ?
   (val(phi,0,0,0) + val(phi,0,-1,0))/2. :
   val(phi,0,0,0) < nodata ? val(phi,0,0,0) :
   val(phi,0,-1,0) < nodata ? val(phi,0,-1,0) :
   0.;

 val(ia.y,0,0,0) += val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*phif*(val(f,0,0,0) - val(f,0,-1,0))/Delta;
      } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 88
{

#line 88 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"

    if (list) for (scalar f = *list, *_i99 = list; ((scalar *)&f)->i >= 0; f = *++_i99)
      if (val(f,0,0,0) != val(f,0,0,-1)) {
#line 100 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < nodata && val(phi,0,0,-1) < nodata) ?
   (val(phi,0,0,0) + val(phi,0,0,-1))/2. :
   val(phi,0,0,0) < nodata ? val(phi,0,0,0) :
   val(phi,0,0,-1) < nodata ? val(phi,0,0,-1) :
   0.;

 val(ia.z,0,0,0) += val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*phif*(val(f,0,0,0) - val(f,0,0,-1))/Delta;
      } }  }}  end_foreach_face_generic()
#line 109
 end_foreach_face(); }
if (is_constant(alpha.x) && is_constant(fm.x)) {
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 88
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 88
{

#line 88 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"

    if (list) for (scalar f = *list, *_i99 = list; ((scalar *)&f)->i >= 0; f = *++_i99)
      if (val(f,0,0,0) != val(f,-1,0,0)) {
#line 100 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < nodata && val(phi,-1,0,0) < nodata) ?
   (val(phi,0,0,0) + val(phi,-1,0,0))/2. :
   val(phi,0,0,0) < nodata ? val(phi,0,0,0) :
   val(phi,-1,0,0) < nodata ? val(phi,-1,0,0) :
   0.;

 val(ia.x,0,0,0) += val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*phif*(val(f,0,0,0) - val(f,-1,0,0))/Delta;
      } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 88
{

#line 88 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"

    if (list) for (scalar f = *list, *_i99 = list; ((scalar *)&f)->i >= 0; f = *++_i99)
      if (val(f,0,0,0) != val(f,0,-1,0)) {
#line 100 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < nodata && val(phi,0,-1,0) < nodata) ?
   (val(phi,0,0,0) + val(phi,0,-1,0))/2. :
   val(phi,0,0,0) < nodata ? val(phi,0,0,0) :
   val(phi,0,-1,0) < nodata ? val(phi,0,-1,0) :
   0.;

 val(ia.y,0,0,0) += val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*phif*(val(f,0,0,0) - val(f,0,-1,0))/Delta;
      } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 88
{

#line 88 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"

    if (list) for (scalar f = *list, *_i99 = list; ((scalar *)&f)->i >= 0; f = *++_i99)
      if (val(f,0,0,0) != val(f,0,0,-1)) {
#line 100 "/home/miro3/Documents/Programming/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < nodata && val(phi,0,0,-1) < nodata) ?
   (val(phi,0,0,0) + val(phi,0,0,-1))/2. :
   val(phi,0,0,0) < nodata ? val(phi,0,0,0) :
   val(phi,0,0,-1) < nodata ? val(phi,0,0,-1) :
   0.;

 val(ia.z,0,0,0) += val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*phif*(val(f,0,0,0) - val(f,0,0,-1))/Delta;
      } }  }}  end_foreach_face_generic()
#line 109
 end_foreach_face(); } }






  if (list) for (scalar f = *list, *_i100 = list; ((scalar *)&f)->i >= 0; f = *++_i100)
    _attribute[f.i].prolongation = fraction_refine;
  boundary (list);






  if (list) for (scalar f = *list, *_i101 = list; ((scalar *)&f)->i >= 0; f = *++_i101) {
    scalar phi = _attribute[f.i].phi;
    delete (((scalar []){phi,{-1}}));
    _attribute[f.i].phi.i = 0;
  }
  pfree (list,__func__,__FILE__,__LINE__);
 end_trace("acceleration_0", "/home/miro3/Documents/Programming/basilisk/src/iforce.h", 131); } return 0; } 
#line 26 "/home/miro3/Documents/Programming/basilisk/src/reduced.h"
#line 1 "curvature.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"
#line 12 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"
static void curvature_restriction (Point point, scalar kappa)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 13 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"

  double k = 0., s = 0.;
   { foreach_child()
    if (val(kappa,0,0,0) != nodata)
      k += val(kappa,0,0,0), s++; end_foreach_child(); }
  val(kappa,0,0,0) = s ? k/s : nodata;
}







static void curvature_prolongation (Point point, scalar kappa)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 28 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"

   { foreach_child() {
    double sk = 0., s = 0.;
    for (int i = 0; i <= 1; i++)

      for (int j = 0; j <= 1; j++)


 for (int k = 0; k <= 1; k++)

   if (coarse(kappa,child.x*i,child.y*j,child.z*k) != nodata)
     sk += coarse(kappa,child.x*i,child.y*j,child.z*k), s++;
    val(kappa,0,0,0) = s ? sk/s : nodata;
  } end_foreach_child(); }
}
#line 62 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"
#line 1 "heights.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/heights.h"
#line 29 "/home/miro3/Documents/Programming/basilisk/src/heights.h"
static inline double height (double H) {
  return H > 20./2. ? H - 20. : H < -20./2. ? H + 20. : H;
}

static inline int orientation (double H) {
  return fabs(H) > 20./2.;
}
#line 49 "/home/miro3/Documents/Programming/basilisk/src/heights.h"
static void half_column (Point point, scalar c, vector h, vector cs, int j)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 50 "/home/miro3/Documents/Programming/basilisk/src/heights.h"







  const int complete = -1;

  {
#line 59
 {







    double S = val(c,0,0,0), H = S, ci, a;







    typedef struct { int s; double h; } HState;
    HState state = {0, 0};
    if (j == 1) {




      if (val(h.x,0,0,0) == 300.)
 state.s = complete, state.h = nodata;




      else {
 int s = (val(h.x,0,0,0) + 20./2.)/100.;
 state.h = val(h.x,0,0,0) - 100.*s;
 state.s = s - 1;
      }





      if (state.s != complete)
 S = state.s, H = state.h;
    }
#line 109 "/home/miro3/Documents/Programming/basilisk/src/heights.h"
    for (int i = 1; i <= 4; i++) {
      ci = i <= 2 ? val(c,i*j,0,0) : val(cs.x,(i - 2)*j,0,0);
      H += ci;




      if (S > 0. && S < 1.) {
 S = ci;
 if (ci <= 0. || ci >= 1.) {







   H -= i*ci;
   break;
 }
      }
#line 138 "/home/miro3/Documents/Programming/basilisk/src/heights.h"
      else if (S >= 1. && ci <= 0.) {
 H = (H - 0.5)*j + (j == -1)*20.;
 S = complete;
 break;
      }
      else if (S <= 0. && ci >= 1.) {
 H = (i + 0.5 - H)*j + (j == 1)*20.;
 S = complete;
 break;
      }
#line 156 "/home/miro3/Documents/Programming/basilisk/src/heights.h"
      else if (S == ci && modf(H, &a))
 break;
    }





    if (j == -1) {







      if (S != complete && ((val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) ||
       (S > 0. && S < 1.)))
 val(h.x,0,0,0) = 300.;
      else if (S == complete)
 val(h.x,0,0,0) = H;
      else





 val(h.x,0,0,0) = H + 100.*(1. + (S >= 1.));
    }
    else {
#line 195 "/home/miro3/Documents/Programming/basilisk/src/heights.h"
      if (state.s != complete ||
   (S == complete && fabs(height(H)) < fabs(height(state.h))))
 state.s = S, state.h = H;





      if (state.s != complete)
 val(h.x,0,0,0) = nodata;
      else
 val(h.x,0,0,0) = (state.h > 1e10 ? nodata : state.h);
    }
  }
#line 59
 {







    double S = val(c,0,0,0), H = S, ci, a;







    typedef struct { int s; double h; } HState;
    HState state = {0, 0};
    if (j == 1) {




      if (val(h.y,0,0,0) == 300.)
 state.s = complete, state.h = nodata;




      else {
 int s = (val(h.y,0,0,0) + 20./2.)/100.;
 state.h = val(h.y,0,0,0) - 100.*s;
 state.s = s - 1;
      }





      if (state.s != complete)
 S = state.s, H = state.h;
    }
#line 109 "/home/miro3/Documents/Programming/basilisk/src/heights.h"
    for (int i = 1; i <= 4; i++) {
      ci = i <= 2 ? val(c,0,i*j,0) : val(cs.y,0,(i - 2)*j,0);
      H += ci;




      if (S > 0. && S < 1.) {
 S = ci;
 if (ci <= 0. || ci >= 1.) {







   H -= i*ci;
   break;
 }
      }
#line 138 "/home/miro3/Documents/Programming/basilisk/src/heights.h"
      else if (S >= 1. && ci <= 0.) {
 H = (H - 0.5)*j + (j == -1)*20.;
 S = complete;
 break;
      }
      else if (S <= 0. && ci >= 1.) {
 H = (i + 0.5 - H)*j + (j == 1)*20.;
 S = complete;
 break;
      }
#line 156 "/home/miro3/Documents/Programming/basilisk/src/heights.h"
      else if (S == ci && modf(H, &a))
 break;
    }





    if (j == -1) {







      if (S != complete && ((val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) ||
       (S > 0. && S < 1.)))
 val(h.y,0,0,0) = 300.;
      else if (S == complete)
 val(h.y,0,0,0) = H;
      else





 val(h.y,0,0,0) = H + 100.*(1. + (S >= 1.));
    }
    else {
#line 195 "/home/miro3/Documents/Programming/basilisk/src/heights.h"
      if (state.s != complete ||
   (S == complete && fabs(height(H)) < fabs(height(state.h))))
 state.s = S, state.h = H;





      if (state.s != complete)
 val(h.y,0,0,0) = nodata;
      else
 val(h.y,0,0,0) = (state.h > 1e10 ? nodata : state.h);
    }
  }
#line 59
 {







    double S = val(c,0,0,0), H = S, ci, a;







    typedef struct { int s; double h; } HState;
    HState state = {0, 0};
    if (j == 1) {




      if (val(h.z,0,0,0) == 300.)
 state.s = complete, state.h = nodata;




      else {
 int s = (val(h.z,0,0,0) + 20./2.)/100.;
 state.h = val(h.z,0,0,0) - 100.*s;
 state.s = s - 1;
      }





      if (state.s != complete)
 S = state.s, H = state.h;
    }
#line 109 "/home/miro3/Documents/Programming/basilisk/src/heights.h"
    for (int i = 1; i <= 4; i++) {
      ci = i <= 2 ? val(c,0,0,i*j) : val(cs.z,0,0,(i - 2)*j);
      H += ci;




      if (S > 0. && S < 1.) {
 S = ci;
 if (ci <= 0. || ci >= 1.) {







   H -= i*ci;
   break;
 }
      }
#line 138 "/home/miro3/Documents/Programming/basilisk/src/heights.h"
      else if (S >= 1. && ci <= 0.) {
 H = (H - 0.5)*j + (j == -1)*20.;
 S = complete;
 break;
      }
      else if (S <= 0. && ci >= 1.) {
 H = (i + 0.5 - H)*j + (j == 1)*20.;
 S = complete;
 break;
      }
#line 156 "/home/miro3/Documents/Programming/basilisk/src/heights.h"
      else if (S == ci && modf(H, &a))
 break;
    }





    if (j == -1) {







      if (S != complete && ((val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) ||
       (S > 0. && S < 1.)))
 val(h.z,0,0,0) = 300.;
      else if (S == complete)
 val(h.z,0,0,0) = H;
      else





 val(h.z,0,0,0) = H + 100.*(1. + (S >= 1.));
    }
    else {
#line 195 "/home/miro3/Documents/Programming/basilisk/src/heights.h"
      if (state.s != complete ||
   (S == complete && fabs(height(H)) < fabs(height(state.h))))
 state.s = S, state.h = H;





      if (state.s != complete)
 val(h.z,0,0,0) = nodata;
      else
 val(h.z,0,0,0) = (state.h > 1e10 ? nodata : state.h);
    }
  }}
}
#line 222 "/home/miro3/Documents/Programming/basilisk/src/heights.h"
static void column_propagation (vector h)
{
   { foreach(){

#line 224 "/home/miro3/Documents/Programming/basilisk/src/heights.h"

    for (int i = -2; i <= 2; i++)
      {
#line 226

 if (fabs(height(val(h.x,i,0,0))) <= 3.5 &&
     fabs(height(val(h.x,i,0,0)) + i) < fabs(height(val(h.x,0,0,0))))
   val(h.x,0,0,0) = val(h.x,i,0,0) + i;
#line 226

 if (fabs(height(val(h.y,0,i,0))) <= 3.5 &&
     fabs(height(val(h.y,0,i,0)) + i) < fabs(height(val(h.y,0,0,0))))
   val(h.y,0,0,0) = val(h.y,0,i,0) + i;
#line 226

 if (fabs(height(val(h.z,0,0,i))) <= 3.5 &&
     fabs(height(val(h.z,0,0,i)) + i) < fabs(height(val(h.z,0,0,0))))
   val(h.z,0,0,0) = val(h.z,0,0,i) + i;}; } end_foreach(); }
  boundary ((scalar *)((vector []){{h.x,h.y,h.z},{{-1},{-1},{-1}}}));
}
#line 291 "/home/miro3/Documents/Programming/basilisk/src/heights.h"

#line 291

static void refine_h_x (Point point, scalar h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 293 "/home/miro3/Documents/Programming/basilisk/src/heights.h"





  bool complete = true;
   { foreach_child() {
    for (int i = -2; i <= 2; i++)
      if (allocated(i,0,0) &&
   !(!is_leaf(neighbor(i,0,0)) && !neighbor(i,0,0).neighbors && neighbor(i,0,0).pid >= 0) && !(neighbor(i,0,0).pid < 0) &&
   fabs(height(val(h,i,0,0))) <= 3.5 &&
   fabs(height(val(h,i,0,0)) + i) < fabs(height(val(h,0,0,0))))
 val(h,0,0,0) = val(h,i,0,0) + i;
    if (val(h,0,0,0) == nodata)
      complete = false;
  } end_foreach_child(); }
  if (complete)
    return;
#line 319 "/home/miro3/Documents/Programming/basilisk/src/heights.h"
  int ori = orientation(val(h,0,0,0));
#line 332 "/home/miro3/Documents/Programming/basilisk/src/heights.h"
  double H[3][3], H0 = height(val(h,0,0,0));
  for (int i = -1; i <= 1; i++)
    for (int j = -1; j <= 1; j++)
      if (val(h,0,i,j) == nodata || orientation(val(h,0,i,j)) != ori)
 return;
      else
 H[i+1][j+1] = height(val(h,0,i,j)) - H0;

  double h0 =
    2.*H0 + (H[2][2] + H[2][0] + H[0][0] + H[0][2] +
      30.*(H[2][1] + H[0][1] + H[1][0] + H[1][2]))/512.
    + 20.*ori;
  double h1 = (H[2][2] + H[2][0] - H[0][0] - H[0][2] +
        30.*(H[2][1] - H[0][1]))/128.;
  double h2 = (H[2][2] - H[2][0] - H[0][0] + H[0][2] +
        30.*(H[1][2] - H[1][0]))/128.;
  double h3 = (H[0][0] + H[2][2] - H[0][2] - H[2][0])/32.;
   { foreach_child()
    if (val(h,0,0,0) == nodata)
      val(h,0,0,0) = h0 + h1*child.y + h2*child.z + h3*child.y*child.z - child.x/2.; end_foreach_child(); }

}
#line 291

static void refine_h_y (Point point, scalar h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 293 "/home/miro3/Documents/Programming/basilisk/src/heights.h"





  bool complete = true;
   { foreach_child() {
    for (int i = -2; i <= 2; i++)
      if (allocated(0,i,0) &&
   !(!is_leaf(neighbor(0,i,0)) && !neighbor(0,i,0).neighbors && neighbor(0,i,0).pid >= 0) && !(neighbor(0,i,0).pid < 0) &&
   fabs(height(val(h,0,i,0))) <= 3.5 &&
   fabs(height(val(h,0,i,0)) + i) < fabs(height(val(h,0,0,0))))
 val(h,0,0,0) = val(h,0,i,0) + i;
    if (val(h,0,0,0) == nodata)
      complete = false;
  } end_foreach_child(); }
  if (complete)
    return;
#line 319 "/home/miro3/Documents/Programming/basilisk/src/heights.h"
  int ori = orientation(val(h,0,0,0));
#line 332 "/home/miro3/Documents/Programming/basilisk/src/heights.h"
  double H[3][3], H0 = height(val(h,0,0,0));
  for (int i = -1; i <= 1; i++)
    for (int j = -1; j <= 1; j++)
      if (val(h,j,0,i) == nodata || orientation(val(h,j,0,i)) != ori)
 return;
      else
 H[i+1][j+1] = height(val(h,j,0,i)) - H0;

  double h0 =
    2.*H0 + (H[2][2] + H[2][0] + H[0][0] + H[0][2] +
      30.*(H[2][1] + H[0][1] + H[1][0] + H[1][2]))/512.
    + 20.*ori;
  double h1 = (H[2][2] + H[2][0] - H[0][0] - H[0][2] +
        30.*(H[2][1] - H[0][1]))/128.;
  double h2 = (H[2][2] - H[2][0] - H[0][0] + H[0][2] +
        30.*(H[1][2] - H[1][0]))/128.;
  double h3 = (H[0][0] + H[2][2] - H[0][2] - H[2][0])/32.;
   { foreach_child()
    if (val(h,0,0,0) == nodata)
      val(h,0,0,0) = h0 + h1*child.z + h2*child.x + h3*child.z*child.x - child.y/2.; end_foreach_child(); }

}
#line 291

static void refine_h_z (Point point, scalar h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 293 "/home/miro3/Documents/Programming/basilisk/src/heights.h"





  bool complete = true;
   { foreach_child() {
    for (int i = -2; i <= 2; i++)
      if (allocated(0,0,i) &&
   !(!is_leaf(neighbor(0,0,i)) && !neighbor(0,0,i).neighbors && neighbor(0,0,i).pid >= 0) && !(neighbor(0,0,i).pid < 0) &&
   fabs(height(val(h,0,0,i))) <= 3.5 &&
   fabs(height(val(h,0,0,i)) + i) < fabs(height(val(h,0,0,0))))
 val(h,0,0,0) = val(h,0,0,i) + i;
    if (val(h,0,0,0) == nodata)
      complete = false;
  } end_foreach_child(); }
  if (complete)
    return;
#line 319 "/home/miro3/Documents/Programming/basilisk/src/heights.h"
  int ori = orientation(val(h,0,0,0));
#line 332 "/home/miro3/Documents/Programming/basilisk/src/heights.h"
  double H[3][3], H0 = height(val(h,0,0,0));
  for (int i = -1; i <= 1; i++)
    for (int j = -1; j <= 1; j++)
      if (val(h,i,j,0) == nodata || orientation(val(h,i,j,0)) != ori)
 return;
      else
 H[i+1][j+1] = height(val(h,i,j,0)) - H0;

  double h0 =
    2.*H0 + (H[2][2] + H[2][0] + H[0][0] + H[0][2] +
      30.*(H[2][1] + H[0][1] + H[1][0] + H[1][2]))/512.
    + 20.*ori;
  double h1 = (H[2][2] + H[2][0] - H[0][0] - H[0][2] +
        30.*(H[2][1] - H[0][1]))/128.;
  double h2 = (H[2][2] - H[2][0] - H[0][0] + H[0][2] +
        30.*(H[1][2] - H[1][0]))/128.;
  double h3 = (H[0][0] + H[2][2] - H[0][2] - H[2][0])/32.;
   { foreach_child()
    if (val(h,0,0,0) == nodata)
      val(h,0,0,0) = h0 + h1*child.x + h2*child.y + h3*child.x*child.y - child.z/2.; end_foreach_child(); }

}







void heights (scalar c, vector h)
{ trace ("heights", "/home/miro3/Documents/Programming/basilisk/src/heights.h", 362);
  vector cs= new_vector("cs");
  {
#line 364

    for (int i = 0; i < nboundary; i++)
      _attribute[cs.x.i].boundary[i] = _attribute[c.i].boundary[i];
#line 364

    for (int i = 0; i < nboundary; i++)
      _attribute[cs.y.i].boundary[i] = _attribute[c.i].boundary[i];
#line 364

    for (int i = 0; i < nboundary; i++)
      _attribute[cs.z.i].boundary[i] = _attribute[c.i].boundary[i];}





  restriction (((scalar []){c,{-1}}));
  for (int j = -1; j <= 1; j += 2) {





     { foreach_level(0){

#line 379 "/home/miro3/Documents/Programming/basilisk/src/heights.h"

      {
#line 380

        val(h.x,0,0,0) = nodata;
#line 380

        val(h.y,0,0,0) = nodata;
#line 380

        val(h.z,0,0,0) = nodata;}; } end_foreach_level(); }

    for (int l = 1; l <= depth(); l++) {




       { foreach_level (l){

#line 388 "/home/miro3/Documents/Programming/basilisk/src/heights.h"

 {
#line 389

   val(cs.x,0,0,0) = val(c,2*j,0,0);
#line 389

   val(cs.y,0,0,0) = val(c,0,2*j,0);
#line 389

   val(cs.z,0,0,0) = val(c,0,0,2*j);}; } end_foreach_level(); }
#line 400 "/home/miro3/Documents/Programming/basilisk/src/heights.h"
       { foreach_level (l - 1){

#line 400 "/home/miro3/Documents/Programming/basilisk/src/heights.h"

 {
#line 401
 {
   val(cs.x,0,0,0) = val(c,j,0,0);
   val(cs.x,j,0,0) = val(c,2*j,0,0);
        }
#line 401
 {
   val(cs.y,0,0,0) = val(c,0,j,0);
   val(cs.y,0,j,0) = val(c,0,2*j,0);
        }
#line 401
 {
   val(cs.z,0,0,0) = val(c,0,0,j);
   val(cs.z,0,0,j) = val(c,0,0,2*j);
        }} } end_foreach_level(); }






       { foreach_halo (prolongation, l - 1){

#line 411 "/home/miro3/Documents/Programming/basilisk/src/heights.h"

 {
#line 412

   _attribute[c.i].prolongation (point, cs.x);
#line 412

   _attribute[c.i].prolongation (point, cs.y);
#line 412

   _attribute[c.i].prolongation (point, cs.z);}; } end_foreach_halo(); }
      { Boundary ** _i = boundaries, * _b; while ((_b = *_i++)) if (_b->level) _b->level (_b, (scalar *)((vector []){{cs.x,cs.y,cs.z},{{-1},{-1},{-1}}}), l); };





       { foreach_level (l){

#line 420 "/home/miro3/Documents/Programming/basilisk/src/heights.h"

        half_column (point, c, h, cs, j); } end_foreach_level(); }
    }
  }






  {
#line 430
 {
    _attribute[h.x.i].prolongation = no_data;
    _attribute[h.x.i].restriction = no_restriction;
  }
#line 430
 {
    _attribute[h.y.i].prolongation = no_data;
    _attribute[h.y.i].restriction = no_restriction;
  }
#line 430
 {
    _attribute[h.z.i].prolongation = no_data;
    _attribute[h.z.i].restriction = no_restriction;
  }}
  boundary ((scalar *)((vector []){{h.x,h.y,h.z},{{-1},{-1},{-1}}}));






  {
#line 441

    _attribute[h.x.i].prolongation = refine_h_x;
#line 441

    _attribute[h.y.i].prolongation = refine_h_y;
#line 441

    _attribute[h.z.i].prolongation = refine_h_z;}




  column_propagation (h);
 delete (((scalar []){cs.x,cs.y,cs.z,{-1}}));  end_trace("heights", "/home/miro3/Documents/Programming/basilisk/src/heights.h", 448); }
#line 63 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"
#line 76 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"

#line 76

static double kappa_z (Point point, vector h) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 77 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"

  int ori = orientation(val(h.z,0,0,0));
  for (int i = -1; i <= 1; i++)
    for (int j = -1; j <= 1; j++)
      if (val(h.z,i,j,0) == nodata || orientation(val(h.z,i,j,0)) != ori)
 return nodata;
  double hx = (val(h.z,1,0,0) - val(h.z,-1,0,0))/2.;
  double hy = (val(h.z,0,1,0) - val(h.z,0,-1,0))/2.;
  double hxx = (val(h.z,1,0,0) + val(h.z,-1,0,0) - 2.*val(h.z,0,0,0))/Delta;
  double hyy = (val(h.z,0,1,0) + val(h.z,0,-1,0) - 2.*val(h.z,0,0,0))/Delta;
  double hxy = (val(h.z,1,1,0) + val(h.z,-1,-1,0) - val(h.z,1,-1,0) - val(h.z,-1,1,0))/(4.*Delta);
  return (hxx*(1. + sq(hy)) + hyy*(1. + sq(hx)) - 2.*hxy*hx*hy)/
    pow(1. + sq(hx) + sq(hy), 3/2.);
}
#line 76

static double kappa_x (Point point, vector h) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 77 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"

  int ori = orientation(val(h.x,0,0,0));
  for (int i = -1; i <= 1; i++)
    for (int j = -1; j <= 1; j++)
      if (val(h.x,0,i,j) == nodata || orientation(val(h.x,0,i,j)) != ori)
 return nodata;
  double hx = (val(h.x,0,1,0) - val(h.x,0,-1,0))/2.;
  double hy = (val(h.x,0,0,1) - val(h.x,0,0,-1))/2.;
  double hxx = (val(h.x,0,1,0) + val(h.x,0,-1,0) - 2.*val(h.x,0,0,0))/Delta;
  double hyy = (val(h.x,0,0,1) + val(h.x,0,0,-1) - 2.*val(h.x,0,0,0))/Delta;
  double hxy = (val(h.x,0,1,1) + val(h.x,0,-1,-1) - val(h.x,0,1,-1) - val(h.x,0,-1,1))/(4.*Delta);
  return (hxx*(1. + sq(hy)) + hyy*(1. + sq(hx)) - 2.*hxy*hx*hy)/
    pow(1. + sq(hx) + sq(hy), 3/2.);
}
#line 76

static double kappa_y (Point point, vector h) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 77 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"

  int ori = orientation(val(h.y,0,0,0));
  for (int i = -1; i <= 1; i++)
    for (int j = -1; j <= 1; j++)
      if (val(h.y,j,0,i) == nodata || orientation(val(h.y,j,0,i)) != ori)
 return nodata;
  double hx = (val(h.y,0,0,1) - val(h.y,0,0,-1))/2.;
  double hy = (val(h.y,1,0,0) - val(h.y,-1,0,0))/2.;
  double hxx = (val(h.y,0,0,1) + val(h.y,0,0,-1) - 2.*val(h.y,0,0,0))/Delta;
  double hyy = (val(h.y,1,0,0) + val(h.y,-1,0,0) - 2.*val(h.y,0,0,0))/Delta;
  double hxy = (val(h.y,1,0,1) + val(h.y,-1,0,-1) - val(h.y,-1,0,1) - val(h.y,1,0,-1))/(4.*Delta);
  return (hxx*(1. + sq(hy)) + hyy*(1. + sq(hx)) - 2.*hxy*hx*hy)/
    pow(1. + sq(hx) + sq(hy), 3/2.);
}
#line 99 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"
static double height_curvature (Point point, scalar c, vector h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 100 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"







  typedef struct {
    double n;
    double (* kappa) (Point, vector);
  } NormKappa;
  struct { NormKappa x, y, z; } n;
  {
#line 112

    n.x.n = val(c,1,0,0) - val(c,-1,0,0), n.x.kappa = kappa_x;
#line 112

    n.y.n = val(c,0,1,0) - val(c,0,-1,0), n.y.kappa = kappa_y;
#line 112

    n.z.n = val(c,0,0,1) - val(c,0,0,-1), n.z.kappa = kappa_z;}
  double (* kappaf) (Point, vector) = NULL; NOT_UNUSED (kappaf);




  if (fabs(n.x.n) < fabs(n.y.n))
    swap (NormKappa, n.x, n.y);

  if (fabs(n.x.n) < fabs(n.z.n))
    swap (NormKappa, n.x, n.z);
  if (fabs(n.y.n) < fabs(n.z.n))
    swap (NormKappa, n.y, n.z);





  double kappa = nodata;
  {
#line 132

    if (kappa == nodata) {
      kappa = n.x.kappa (point, h);
      if (kappa != nodata) {
 kappaf = n.x.kappa;
 if (n.x.n < 0.)
   kappa = - kappa;
      }
    }
#line 132

    if (kappa == nodata) {
      kappa = n.y.kappa (point, h);
      if (kappa != nodata) {
 kappaf = n.y.kappa;
 if (n.y.n < 0.)
   kappa = - kappa;
      }
    }
#line 132

    if (kappa == nodata) {
      kappa = n.z.kappa (point, h);
      if (kappa != nodata) {
 kappaf = n.z.kappa;
 if (n.z.n < 0.)
   kappa = - kappa;
      }
    }}

  if (kappa != nodata) {




    if (fabs(kappa) > 1./Delta)
      kappa = sign(kappa)/Delta;
#line 167 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"
  }

  return kappa;
}
#line 184 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"
#line 1 "parabola.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/parabola.h"
#line 1 "utils.h"
#line 2 "/home/miro3/Documents/Programming/basilisk/src/parabola.h"






typedef struct {
  coord o;




  double t[3][3];



  double ** M, rhs[6], a[6];


} ParabolaFit;

static void parabola_fit_init (ParabolaFit * p, coord o, coord m)
{
  {
#line 25

    p->o.x = o.x;
#line 25

    p->o.y = o.y;
#line 25

    p->o.z = o.z;}






  double max;
  coord nx = {0., 0., 0.}, ny, nz;
  int d = 0;

  {
#line 37

    nz.x = m.x;
#line 37

    nz.y = m.y;
#line 37

    nz.z = m.z;}
  normalize (&nz);
  max = sq(nz.x);

  if (sq(nz.y) > max) { max = sq(nz.y); d = 1; }
  if (sq(nz.z) > max) d = 2;
  switch (d) {
  case 0: nx.x = - nz.z/nz.x; nx.z = 1.0; break;
  case 1: nx.y = - nz.z/nz.y; nx.z = 1.0; break;
  case 2: nx.z = - nz.x/nz.z; nx.x = 1.0; break;
  }
  normalize (&nx);


  {
#line 52

    ny.x = nz.y*nx.z - nz.z*nx.y;
#line 52

    ny.y = nz.z*nx.x - nz.x*nx.z;
#line 52

    ny.z = nz.x*nx.y - nz.y*nx.x;}


  p->t[0][0] = nx.x; p->t[0][1] = nx.y; p->t[0][2] = nx.z;
  p->t[1][0] = ny.x; p->t[1][1] = ny.y; p->t[1][2] = ny.z;
  p->t[2][0] = nz.x; p->t[2][1] = nz.y; p->t[2][2] = nz.z;



  int n = 6;


  p->M = (double **) matrix_new (n, n, sizeof(double));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      p->M[i][j] = 0.;
    p->rhs[i] = 0.;
  }
}

static void parabola_fit_add (ParabolaFit * p, coord m, double w)
{
#line 85 "/home/miro3/Documents/Programming/basilisk/src/parabola.h"
  double x1 = m.x - p->o.x, y1 = m.y - p->o.y, z1 = m.z - p->o.z;
  double x = p->t[0][0]*x1 + p->t[0][1]*y1 + p->t[0][2]*z1;
  double y = p->t[1][0]*x1 + p->t[1][1]*y1 + p->t[1][2]*z1;
  double z = p->t[2][0]*x1 + p->t[2][1]*y1 + p->t[2][2]*z1;
#line 98 "/home/miro3/Documents/Programming/basilisk/src/parabola.h"
  double x2 = x*x, x3 = x2*x, x4 = x3*x;
  double y2 = y*y, y3 = y2*y, y4 = y3*y;
  p->M[0][0] += w*x4; p->M[1][1] += w*y4; p->M[2][2] += w*x2*y2;
  p->M[3][3] += w*x2; p->M[4][4] += w*y2; p->M[5][5] += w;
  p->M[0][2] += w*x3*y; p->M[0][3] += w*x3; p->M[0][4] += w*x2*y;
  p->M[1][2] += w*x*y3; p->M[1][3] += w*x*y2; p->M[1][4] += w*y3;
  p->M[2][5] += w*x*y;
  p->M[3][5] += w*x;
  p->M[4][5] += w*y;
  p->rhs[0] += w*x2*z; p->rhs[1] += w*y2*z; p->rhs[2] += w*x*y*z;
  p->rhs[3] += w*x*z; p->rhs[4] += w*y*z; p->rhs[5] += w*z;


}

static double parabola_fit_solve (ParabolaFit * p)
{
#line 139 "/home/miro3/Documents/Programming/basilisk/src/parabola.h"
  p->M[0][1] = p->M[2][2]; p->M[0][5] = p->M[3][3];
  p->M[1][5] = p->M[4][4];
  p->M[2][3] = p->M[0][4]; p->M[2][4] = p->M[1][3];
  p->M[3][4] = p->M[2][5];
  for (int i = 1; i < 6; i++)
    for (int j = 0; j < i; j++)
      p->M[i][j] = p->M[j][i];
  double pivmin = matrix_inverse (p->M, 6, 1e-10);
  if (pivmin)
    for (int i = 0; i < 6; i++) {
      p->a[i] = 0.;
      for (int j = 0; j < 6; j++)
 p->a[i] += p->M[i][j]*p->rhs[j];
    }
  else
    for (int i = 0; i < 6; i++)
      p->a[i] = 0.;


  matrix_free (p->M);
  return pivmin;
}

static double parabola_fit_curvature (ParabolaFit * p,
          double kappamax, double * kmax)
{
  double kappa;
#line 176 "/home/miro3/Documents/Programming/basilisk/src/parabola.h"
  double hxx = 2.*p->a[0], hyy = 2.*p->a[1], hxy = p->a[2];
  double hx = p->a[3], hy = p->a[4];

  double dnm = 1. + sq(hx) + sq(hy);
  kappa = - (hxx*(1. + sq(hy)) + hyy*(1. + sq(hx)) - 2.*hxy*hx*hy)
    /sqrt (dnm*dnm*dnm);
  if (kmax) {
    double kg = (hxx*hyy - hxy*hxy)/(dnm*dnm);
    double a = kappa*kappa/4. - kg;
    *kmax = fabs (kappa/2.);
    if (a >= 0.)
      *kmax += sqrt (a);
  }

  if (fabs (kappa) > kappamax) {
    if (kmax)
      *kmax = kappamax;
    return kappa > 0. ? kappamax : - kappamax;
  }
  return kappa;
}
#line 185 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"






static int independents (coord * p, int n)
{
  if (n < 2)
    return n;
  int ni = 1;
  for (int j = 1; j < n; j++) {
    bool depends = false;
    for (int i = 0; i < j && !depends; i++) {
      double d2 = 0.;
      {
#line 200

 d2 += sq(p[i].x - p[j].x);
#line 200

 d2 += sq(p[i].y - p[j].y);
#line 200

 d2 += sq(p[i].z - p[j].z);}
      depends = (d2 < sq(0.5));
    }
    ni += !depends;
  }
  return ni;
}






static double height_curvature_fit (Point point, scalar c, vector h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 215 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"






  coord ip[3 == 2 ? 6 : 27];
  int n = 0;




  {
#line 227
 {





    int n1 = 0, n2 = 0;






    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
 if (val(h.z,i,j,0) != nodata) {
   if (orientation(val(h.z,i,j,0))) n1++; else n2++;
 }

    int ori = (n1 > n2);
#line 258 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"
    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
 if (val(h.z,i,j,0) != nodata && orientation(val(h.z,i,j,0)) == ori)
   ip[n].x = i, ip[n].y = j, ip[n++].z = height(val(h.z,i,j,0));

  }
#line 227
 {





    int n1 = 0, n2 = 0;






    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
 if (val(h.x,0,i,j) != nodata) {
   if (orientation(val(h.x,0,i,j))) n1++; else n2++;
 }

    int ori = (n1 > n2);
#line 258 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"
    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
 if (val(h.x,0,i,j) != nodata && orientation(val(h.x,0,i,j)) == ori)
   ip[n].y = i, ip[n].z = j, ip[n++].x = height(val(h.x,0,i,j));

  }
#line 227
 {





    int n1 = 0, n2 = 0;






    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
 if (val(h.y,j,0,i) != nodata) {
   if (orientation(val(h.y,j,0,i))) n1++; else n2++;
 }

    int ori = (n1 > n2);
#line 258 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"
    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
 if (val(h.y,j,0,i) != nodata && orientation(val(h.y,j,0,i)) == ori)
   ip[n].z = i, ip[n].x = j, ip[n++].y = height(val(h.y,j,0,i));

  }}





  if (independents (ip, n) < (3 == 2 ? 3 : 9))
    return nodata;





  coord m = mycs (point, c), fc;
  double alpha = plane_alpha (val(c,0,0,0), m);
  double area = plane_area_center (m, alpha, &fc);
  ParabolaFit fit;
  parabola_fit_init (&fit, fc, m);




  parabola_fit_add (&fit, fc, area*100.);






  for (int i = 0; i < n; i++)
    parabola_fit_add (&fit, ip[i], 1.);
  parabola_fit_solve (&fit);
  double kappa = parabola_fit_curvature (&fit, 2., NULL)/Delta;



  return kappa;
}






static double centroids_curvature_fit (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 308 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"






  coord m = mycs (point, c), fc;
  double alpha = plane_alpha (val(c,0,0,0), m);
  plane_area_center (m, alpha, &fc);
  ParabolaFit fit;
  parabola_fit_init (&fit, fc, m);





  coord r = {x,y,z};
   { foreach_neighbor(1)
    if (val(c,0,0,0) > 0. && val(c,0,0,0) < 1.) {
      coord m = mycs (point, c), fc;
      double alpha = plane_alpha (val(c,0,0,0), m);
      double area = plane_area_center (m, alpha, &fc);
      coord rn = {x,y,z};
      {
#line 331

 fc.x += (rn.x - r.x)/Delta;
#line 331

 fc.y += (rn.y - r.y)/Delta;
#line 331

 fc.z += (rn.z - r.z)/Delta;}
      parabola_fit_add (&fit, fc, area);
    } end_foreach_neighbor(); }
  parabola_fit_solve (&fit);
  double kappa = parabola_fit_curvature (&fit, 2., NULL)/Delta;



  return kappa;
}
#line 354 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"
static inline bool interfacial (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 355 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"

  if (val(c,0,0,0) >= 1.) {
    for (int i = -1; i <= 1; i += 2)
      {
#line 358

 if (val(c,i,0,0) <= 0.)
   return true;
#line 358

 if (val(c,0,i,0) <= 0.)
   return true;
#line 358

 if (val(c,0,0,i) <= 0.)
   return true;}
  }
  else if (val(c,0,0,0) <= 0.) {
    for (int i = -1; i <= 1; i += 2)
      {
#line 364

 if (val(c,i,0,0) >= 1.)
   return true;
#line 364

 if (val(c,0,i,0) >= 1.)
   return true;
#line 364

 if (val(c,0,0,i) >= 1.)
   return true;}
  }
  else
    return true;
  return false;
}
#line 384 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"
typedef struct {
  int h;
  int f;
  int a;
  int c;
} cstats;

struct Curvature {
  scalar c, kappa;
  double sigma;
  bool add;
};


cstats curvature (struct Curvature p)
{ trace ("curvature", "/home/miro3/Documents/Programming/basilisk/src/curvature.h", 399);
  scalar c = p.c, kappa = p.kappa;
  double sigma = p.sigma ? p.sigma : 1.;
  cstats s = {0,0,0,0};
  vector h= new_vector("h");
  heights (c, h);






  _attribute[kappa.i].refine = _attribute[kappa.i].prolongation = curvature_prolongation;
  _attribute[kappa.i].restriction = curvature_restriction;






  scalar k= new_scalar("k");
  scalar_clone (k, kappa);

   { foreach(){

#line 422 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"
 {




    if (!interfacial (point, c))
      val(k,0,0,0) = nodata;





    else if ((val(k,0,0,0) = height_curvature (point, c, h)) != nodata)
      s.h++;
    else if ((val(k,0,0,0) = height_curvature_fit (point, c, h)) != nodata)
      s.f++;
  } } end_foreach(); }
  boundary (((scalar []){k,{-1}}));

   { foreach(){

#line 441 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"
 {





    double kf;
    if (val(k,0,0,0) < nodata)
      kf = val(k,0,0,0);
    else if (interfacial (point, c)) {





      double sk = 0., a = 0.;
       { foreach_neighbor(1)
 if (val(k,0,0,0) < nodata)
   sk += val(k,0,0,0), a++; end_foreach_neighbor(); }
      if (a > 0.)
 kf = sigma*sk/a, s.a++;
      else




 kf = sigma*centroids_curvature_fit (point, c), s.c++;
    }
    else
      kf = nodata;




    if (kf == nodata)
      val(kappa,0,0,0) = nodata;
    else if (p.add)
      val(kappa,0,0,0) += sigma*kf;
    else
      val(kappa,0,0,0) = sigma*kf;
  } } end_foreach(); }
  boundary (((scalar []){kappa,{-1}}));

  { cstats _ret =  s; delete (((scalar []){k,h.x,h.y,h.z,{-1}}));  end_trace("curvature", "/home/miro3/Documents/Programming/basilisk/src/curvature.h", 484);  return _ret; }
 delete (((scalar []){k,h.x,h.y,h.z,{-1}}));  end_trace("curvature", "/home/miro3/Documents/Programming/basilisk/src/curvature.h", 485); }
#line 504 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"

#line 504

static double pos_x (Point point, vector h, coord * G, coord * Z) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 505 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"

  if (val(h.x,0,0,0) == nodata)
    return nodata;
  coord o = {x, y, z};
  o.x += height(val(h.x,0,0,0))*Delta;
  double pos = 0.;
  {
#line 511

    pos += (o.x - Z->x)*G->x;
#line 511

    pos += (o.y - Z->y)*G->y;
#line 511

    pos += (o.z - Z->z)*G->z;}
  return pos;
}
#line 504

static double pos_y (Point point, vector h, coord * G, coord * Z) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 505 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"

  if (val(h.y,0,0,0) == nodata)
    return nodata;
  coord o = {x, y, z};
  o.y += height(val(h.y,0,0,0))*Delta;
  double pos = 0.;
  {
#line 511

    pos += (o.y - Z->y)*G->y;
#line 511

    pos += (o.z - Z->z)*G->z;
#line 511

    pos += (o.x - Z->x)*G->x;}
  return pos;
}
#line 504

static double pos_z (Point point, vector h, coord * G, coord * Z) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 505 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"

  if (val(h.z,0,0,0) == nodata)
    return nodata;
  coord o = {x, y, z};
  o.z += height(val(h.z,0,0,0))*Delta;
  double pos = 0.;
  {
#line 511

    pos += (o.z - Z->z)*G->z;
#line 511

    pos += (o.x - Z->x)*G->x;
#line 511

    pos += (o.y - Z->y)*G->y;}
  return pos;
}







static double height_position (Point point, scalar f, vector h,
          coord * G, coord * Z)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 524 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"







  typedef struct {
    double n;
    double (* pos) (Point, vector, coord *, coord *);
  } NormPos;
  struct { NormPos x, y, z; } n;
  {
#line 536

    n.x.n = val(f,1,0,0) - val(f,-1,0,0), n.x.pos = pos_x;
#line 536

    n.y.n = val(f,0,1,0) - val(f,0,-1,0), n.y.pos = pos_y;
#line 536

    n.z.n = val(f,0,0,1) - val(f,0,0,-1), n.z.pos = pos_z;}




  if (fabs(n.x.n) < fabs(n.y.n))
    swap (NormPos, n.x, n.y);

  if (fabs(n.x.n) < fabs(n.z.n))
    swap (NormPos, n.x, n.z);
  if (fabs(n.y.n) < fabs(n.z.n))
    swap (NormPos, n.y, n.z);





  double pos = nodata;
  {
#line 555

    if (pos == nodata)
      pos = n.x.pos (point, h, G, Z);
#line 555

    if (pos == nodata)
      pos = n.y.pos (point, h, G, Z);
#line 555

    if (pos == nodata)
      pos = n.z.pos (point, h, G, Z);}

  return pos;
}
#line 571 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"
struct Position {
  scalar f, pos;
  coord G, Z;
  bool add;
};

void position (struct Position p)
{
  scalar f = p.f, pos = p.pos;
  coord * G = &p.G, * Z = &p.Z;






  _attribute[pos.i].refine = _attribute[pos.i].prolongation = curvature_prolongation;
  _attribute[pos.i].restriction = curvature_restriction;


  vector h= new_vector("h");
  heights (f, h);
   { foreach(){

#line 593 "/home/miro3/Documents/Programming/basilisk/src/curvature.h"
 {
    if (interfacial (point, f)) {
      double hp = height_position (point, f, h, G, Z);
      if (hp == nodata) {





 coord n = mycs (point, f), o = {x,y,z}, c;
 double alpha = plane_alpha (val(f,0,0,0), n);
 plane_area_center (n, alpha, &c);
 hp = 0.;
 {
#line 606

   hp += (o.x + Delta*c.x - Z->x)*G->x;
#line 606

   hp += (o.y + Delta*c.y - Z->y)*G->y;
#line 606

   hp += (o.z + Delta*c.z - Z->z)*G->z;}
      }
      if (p.add)
 val(pos,0,0,0) += hp;
      else
 val(pos,0,0,0) = hp;
    }
    else
      val(pos,0,0,0) = nodata;
  } } end_foreach(); }
  boundary (((scalar []){pos,{-1}}));
 delete (((scalar []){h.x,h.y,h.z,{-1}})); }
#line 27 "/home/miro3/Documents/Programming/basilisk/src/reduced.h"
#line 36 "/home/miro3/Documents/Programming/basilisk/src/reduced.h"
static int acceleration_1_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int acceleration_1 (const int i, const double t, Event * _ev) { trace ("acceleration_1", "/home/miro3/Documents/Programming/basilisk/src/reduced.h", 36); 
{
  scalar phi = _attribute[f.i].phi;
  coord G1;
  {
#line 40

    G1.x = (rho2 - rho1)*G.x;
#line 40

    G1.y = (rho2 - rho1)*G.y;
#line 40

    G1.z = (rho2 - rho1)*G.z;}

  if (phi.i)
    position ((struct Position){f, phi, G1, Z, .add = true});
  else {
    phi = new_scalar("phi");
    position ((struct Position){f, phi, G1, Z, .add = false});
    _attribute[f.i].phi = phi;
  }
 end_trace("acceleration_1", "/home/miro3/Documents/Programming/basilisk/src/reduced.h", 50); } return 0; } 
#line 12 "milk_crown.c"

#line 1 "tension.h"
#line 1 "/home/miro3/Documents/Programming/basilisk/src/tension.h"
#line 21 "/home/miro3/Documents/Programming/basilisk/src/tension.h"



#line 36 "/home/miro3/Documents/Programming/basilisk/src/tension.h"
static int stability_1_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int stability_1 (const int i, const double t, Event * _ev) { trace ("stability_1", "/home/miro3/Documents/Programming/basilisk/src/tension.h", 36);  {





  double amin = HUGE, amax = -HUGE, dmin = HUGE;
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _amin = amin; double _amax = amax; double _dmin = dmin; 
#line 43

if (!is_constant(alpha.x) && !is_constant(fm.x)) {
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 43
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 43
{

#line 43 "/home/miro3/Documents/Programming/basilisk/src/tension.h"
 {
    if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) > _amax) _amax = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
    if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) < _amin) _amin = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 43
{

#line 43 "/home/miro3/Documents/Programming/basilisk/src/tension.h"
 {
    if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) > _amax) _amax = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
    if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) < _amin) _amin = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 43
{

#line 43 "/home/miro3/Documents/Programming/basilisk/src/tension.h"
 {
    if (val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0) > _amax) _amax = val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0);
    if (val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0) < _amin) _amin = val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  end_foreach_face_generic()
#line 47
 end_foreach_face(); }
if (is_constant(alpha.x) && !is_constant(fm.x)) {
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 43
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 43
{

#line 43 "/home/miro3/Documents/Programming/basilisk/src/tension.h"
 {
    if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) > _amax) _amax = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
    if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) < _amin) _amin = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 43
{

#line 43 "/home/miro3/Documents/Programming/basilisk/src/tension.h"
 {
    if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) > _amax) _amax = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
    if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) < _amin) _amin = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 43
{

#line 43 "/home/miro3/Documents/Programming/basilisk/src/tension.h"
 {
    if (val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0) > _amax) _amax = val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0);
    if (val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0) < _amin) _amin = val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  end_foreach_face_generic()
#line 47
 end_foreach_face(); }
if (!is_constant(alpha.x) && is_constant(fm.x)) {
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 43
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 43
{

#line 43 "/home/miro3/Documents/Programming/basilisk/src/tension.h"
 {
    if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) > _amax) _amax = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
    if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) < _amin) _amin = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 43
{

#line 43 "/home/miro3/Documents/Programming/basilisk/src/tension.h"
 {
    if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) > _amax) _amax = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
    if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) < _amin) _amin = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 43
{

#line 43 "/home/miro3/Documents/Programming/basilisk/src/tension.h"
 {
    if (val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0) > _amax) _amax = val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0);
    if (val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0) < _amin) _amin = val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  end_foreach_face_generic()
#line 47
 end_foreach_face(); }
if (is_constant(alpha.x) && is_constant(fm.x)) {
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 43
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 43
{

#line 43 "/home/miro3/Documents/Programming/basilisk/src/tension.h"
 {
    if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) > _amax) _amax = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
    if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) < _amin) _amin = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 43
{

#line 43 "/home/miro3/Documents/Programming/basilisk/src/tension.h"
 {
    if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) > _amax) _amax = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
    if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) < _amin) _amin = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 43
{

#line 43 "/home/miro3/Documents/Programming/basilisk/src/tension.h"
 {
    if (val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0) > _amax) _amax = val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0);
    if (val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0) < _amin) _amin = val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  end_foreach_face_generic()
#line 47
 end_foreach_face(); }OMP(omp critical) if (_amin < amin) amin = _amin;
mpi_all_reduce_double (amin, MPI_MIN);
OMP(omp critical) if (_amax > amax) amax = _amax;
mpi_all_reduce_double (amax, MPI_MAX);
OMP(omp critical) if (_dmin < dmin) dmin = _dmin;
mpi_all_reduce_double (dmin, MPI_MIN);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 47
 }
  double rhom = (1./amin + 1./amax)/2.;





  if (interfaces) for (scalar c = *interfaces, *_i102 = interfaces; ((scalar *)&c)->i >= 0; c = *++_i102)
    if (_attribute[c.i].sigma) {
      double dt = sqrt (rhom*cube(dmin)/(pi*_attribute[c.i].sigma));
      if (dt < dtmax)
 dtmax = dt;
    }
 end_trace("stability_1", "/home/miro3/Documents/Programming/basilisk/src/tension.h", 60); } return 0; } 







static int acceleration_2_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int acceleration_2 (const int i, const double t, Event * _ev) { trace ("acceleration_2", "/home/miro3/Documents/Programming/basilisk/src/tension.h", 68); 
{




  if (interfaces) for (scalar f = *interfaces, *_i103 = interfaces; ((scalar *)&f)->i >= 0; f = *++_i103)
    if (_attribute[f.i].sigma) {





      scalar phi = _attribute[f.i].phi;
      if (phi.i)
 curvature ((struct Curvature){f, phi, _attribute[f.i].sigma, .add = true});
      else {
 phi = new_scalar("phi");
 curvature ((struct Curvature){f, phi, _attribute[f.i].sigma, .add = false});
 _attribute[f.i].phi = phi;
      }
    }
 end_trace("acceleration_2", "/home/miro3/Documents/Programming/basilisk/src/tension.h", 90); } return 0; } 
#line 14 "milk_crown.c"
#line 1 "output_vtu_foreach.h"
#line 1 "./output_vtu_foreach.h"





void output_pvtu_ascii (scalar * list, vector * vlist, int n, FILE * fp, char * subname)
{
    fputs ("<?xml version=\"1.0\"?>\n"
    "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
    fputs ("\t <PUnstructuredGrid GhostLevel=\"0\">\n", fp);
    fputs ("\t\t\t <PCellData Scalars=\"scalars\">\n", fp);
    if (list) for (scalar s = *list, *_i104 = list; ((scalar *)&s)->i >= 0; s = *++_i104) {
      fprintf (fp,"\t\t\t\t <PDataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", _attribute[s.i].name);
      fputs ("\t\t\t\t </PDataArray>\n", fp);
    }
    if (vlist) for (vector v = *vlist, *_i105 = vlist; ((scalar *)&v)->i >= 0; v = *++_i105) {
      fprintf (fp,"\t\t\t\t <PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"%s\" format=\"ascii\">\n", _attribute[v.x.i].name);
      fputs ("\t\t\t\t </PDataArray>\n", fp);
    }
    fputs ("\t\t\t </PCellData>\n", fp);
    fputs ("\t\t\t <PPoints>\n", fp);
    fputs ("\t\t\t\t <PDataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);
    fputs ("\t\t\t\t </PDataArray>\n", fp);
    fputs ("\t\t\t </PPoints>\n", fp);

    for (int i = 0; i < npe(); i++)
      fprintf (fp, "<Piece Source=\"%s_n%3.3d.vtu\"/> \n", subname, i);

    fputs ("\t </PUnstructuredGrid>\n", fp);
    fputs ("</VTKFile>\n", fp);
}
#line 41 "./output_vtu_foreach.h"
void output_vtu_ascii_foreach (scalar * list, vector * vlist, int n, FILE * fp, bool linear)
{


#ifdef _OPENMP
  int num_omp = omp_get_max_threads();
  omp_set_num_threads(1);
#endif

  scalar marker= new_vertex_scalar("marker");
  int no_points = 0, no_cells=0 ;
   { foreach_vertex(){

#line 52 "./output_vtu_foreach.h"
{
    val(marker,0,0,0) = _k;
    no_points += 1;
  } } end_foreach_vertex(); }
   { foreach(){

#line 56 "./output_vtu_foreach.h"
{
    no_cells += 1;
  } } end_foreach(); }

  fputs ("<?xml version=\"1.0\"?>\n"
  "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
  fputs ("\t <UnstructuredGrid>\n", fp);
  fprintf (fp,"\t\t <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", no_points, no_cells);
  fputs ("\t\t\t <CellData Scalars=\"scalars\">\n", fp);
  if (list) for (scalar s = *list, *_i106 = list; ((scalar *)&s)->i >= 0; s = *++_i106) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", _attribute[s.i].name);
     { foreach(){

#line 67 "./output_vtu_foreach.h"
{
      fprintf (fp, "%g\n", val(s,0,0,0));
    } } end_foreach(); }
    fputs ("\t\t\t\t </DataArray>\n", fp);
  }
  if (vlist) for (vector v = *vlist, *_i107 = vlist; ((scalar *)&v)->i >= 0; v = *++_i107) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"%s\" format=\"ascii\">\n", _attribute[v.x.i].name);
     { foreach(){

#line 74 "./output_vtu_foreach.h"
{




      fprintf (fp, "%g %g %g\n", val(v.x,0,0,0), val(v.y,0,0,0), val(v.z,0,0,0));

    } } end_foreach(); }
    fputs ("\t\t\t\t </DataArray>\n", fp);
  }
  fputs ("\t\t\t </CellData>\n", fp);
  fputs ("\t\t\t <Points>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);
   { foreach_vertex(){

#line 87 "./output_vtu_foreach.h"
{




    fprintf (fp, "%g %g %g\n", x, y, z);

  } } end_foreach_vertex(); }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t </Points>\n", fp);
  fputs ("\t\t\t <Cells>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n", fp);
   { foreach(){

#line 99 "./output_vtu_foreach.h"
{
#line 114 "./output_vtu_foreach.h"
    int ape1 = val(marker,0,0,0);
    int ape2 = val(marker,1,0,0);
    int ape3 = val(marker,1,1,0);
    int ape4 = val(marker,0,1,0);
    int ape5 = val(marker,0,0,1);
    int ape6 = val(marker,1,0,1);
    int ape7 = val(marker,1,1,1);
    int ape8 = val(marker,0,1,1);
    fprintf (fp, "%u %u %u %u %u %u %u %u\n", ape1, ape2, ape3, ape4, ape5, ape6, ape7, ape8);


  } } end_foreach(); }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n", fp);

  for (int i = 1; i < no_cells+1; i++){




    fprintf (fp, "%d \n", i*8);

  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n", fp);
   { foreach(){

#line 139 "./output_vtu_foreach.h"
{




    fputs ("12 \n", fp);

  } } end_foreach(); }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t </Cells>\n", fp);
  fputs ("\t\t </Piece>\n", fp);
  fputs ("\t </UnstructuredGrid>\n", fp);
  fputs ("</VTKFile>\n", fp);
  fflush (fp);



#ifdef _OPENMP
  omp_set_num_threads(num_omp);
#endif
 delete (((scalar []){marker,{-1}})); }






void output_pvtu_bin (scalar * list, vector * vlist, int n, FILE * fp, char * subname)
{
    fputs ("<?xml version=\"1.0\"?>\n"
    "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
    fputs ("\t <PUnstructuredGrid GhostLevel=\"0\">\n", fp);
    fputs ("\t\t\t <PCellData Scalars=\"scalars\">\n", fp);
    if (list) for (scalar s = *list, *_i108 = list; ((scalar *)&s)->i >= 0; s = *++_i108) {
      fprintf (fp,"\t\t\t\t <PDataArray type=\"Float64\" Name=\"%s\" format=\"appended\">\n", _attribute[s.i].name);
      fputs ("\t\t\t\t </PDataArray>\n", fp);
    }
    if (vlist) for (vector v = *vlist, *_i109 = vlist; ((scalar *)&v)->i >= 0; v = *++_i109) {
      fprintf (fp,"\t\t\t\t <PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"%s\" format=\"appended\">\n", _attribute[v.x.i].name);
      fputs ("\t\t\t\t </PDataArray>\n", fp);
    }
    fputs ("\t\t\t </PCellData>\n", fp);
    fputs ("\t\t\t <PPoints>\n", fp);
    fputs ("\t\t\t\t <PDataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);
    fputs ("\t\t\t\t </PDataArray>\n", fp);
    fputs ("\t\t\t </PPoints>\n", fp);

    for (int i = 0; i < npe(); i++)
      fprintf (fp, "<Piece Source=\"%s_n%3.3d.vtu\"/> \n", subname, i);

    fputs ("\t </PUnstructuredGrid>\n", fp);
    fputs ("</VTKFile>\n", fp);
}
#line 201 "./output_vtu_foreach.h"
void output_vtu_bin_foreach (scalar * list, vector * vlist, int n, FILE * fp, bool linear)
{
#ifdef _OPENMP
  int num_omp = omp_get_max_threads();
  omp_set_num_threads(1);
#endif

  scalar marker= new_vertex_scalar("marker");
  int no_points = 0, no_cells=0 ;
   { foreach_vertex(){

#line 210 "./output_vtu_foreach.h"
{
    val(marker,0,0,0) = _k;
    no_points += 1;
  } } end_foreach_vertex(); }
   { foreach(){

#line 214 "./output_vtu_foreach.h"
{
    no_cells += 1;
  } } end_foreach(); }

  fputs ("<?xml version=\"1.0\"?>\n"
  "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
  fputs ("\t <UnstructuredGrid>\n", fp);
  fprintf (fp,"\t\t <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", no_points, no_cells);
  fputs ("\t\t\t <CellData Scalars=\"scalars\">\n", fp);

  int count = 0;
  if (list) for (scalar s = *list, *_i110 = list; ((scalar *)&s)->i >= 0; s = *++_i110) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"appended\" offset=\"%d\">\n", _attribute[s.i].name,count);
    count += ((no_cells)+1)*8;
    fputs ("\t\t\t\t </DataArray>\n", fp);
  }
  if (vlist) for (vector v = *vlist, *_i111 = vlist; ((scalar *)&v)->i >= 0; v = *++_i111) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"3\"  format=\"appended\" offset=\"%d\">\n", _attribute[v.x.i].name,count);
    count += ((no_cells*3)+1)*8;
    fputs ("\t\t\t\t </DataArray>\n", fp);
  }
  fputs ("\t\t\t </CellData>\n", fp);
  fputs ("\t\t\t <Points>\n", fp);
  fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\"  format=\"appended\" offset=\"%d\">\n",count);
  count += ((no_points*3)+1)*8;
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t </Points>\n", fp);
  fputs ("\t\t\t <Cells>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n", fp);
   { foreach(){

#line 243 "./output_vtu_foreach.h"
{
#line 258 "./output_vtu_foreach.h"
    int ape1 = val(marker,0,0,0);
    int ape2 = val(marker,1,0,0);
    int ape3 = val(marker,1,1,0);
    int ape4 = val(marker,0,1,0);
    int ape5 = val(marker,0,0,1);
    int ape6 = val(marker,1,0,1);
    int ape7 = val(marker,1,1,1);
    int ape8 = val(marker,0,1,1);
    fprintf (fp, "%u %u %u %u %u %u %u %u\n", ape1, ape2, ape3, ape4, ape5, ape6, ape7, ape8);


  } } end_foreach(); }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n", fp);
  for (int i = 1; i < no_cells+1; i++){




    fprintf (fp, "%d \n", i*8);

  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n", fp);
   { foreach(){

#line 282 "./output_vtu_foreach.h"
{




    fputs ("12 \n", fp);

  } } end_foreach(); }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t </Cells>\n", fp);
  fputs ("\t\t </Piece>\n", fp);
  fputs ("\t </UnstructuredGrid>\n", fp);
  fputs ("\t <AppendedData encoding=\"raw\">\n", fp);
  fputs ("_", fp);
  unsigned long long block_len=no_cells*8;



  if (list) for (scalar s = *list, *_i112 = list; ((scalar *)&s)->i >= 0; s = *++_i112) {
    fwrite (&block_len, sizeof (unsigned long long), 1, fp);
     { foreach(){

#line 302 "./output_vtu_foreach.h"

      fwrite (&val(s,0,0,0), sizeof (double), 1, fp); } end_foreach(); }
  }
  block_len=no_cells*8*3;
  if (vlist) for (vector v = *vlist, *_i113 = vlist; ((scalar *)&v)->i >= 0; v = *++_i113) {
    fwrite (&block_len, sizeof (unsigned long long), 1, fp);
     { foreach(){

#line 308 "./output_vtu_foreach.h"
{
      fwrite (&val(v.x,0,0,0), sizeof (double), 1, fp);
      fwrite (&val(v.y,0,0,0), sizeof (double), 1, fp);




      fwrite (&val(v.z,0,0,0), sizeof (double), 1, fp);

    } } end_foreach(); }
  }
  block_len=no_points*8*3;
  fwrite (&block_len, sizeof (unsigned long long), 1, fp);
   { foreach_vertex(){

#line 321 "./output_vtu_foreach.h"
{
    fwrite (&x, sizeof (double), 1, fp);
    fwrite (&y, sizeof (double), 1, fp);
    fwrite (&z, sizeof (double), 1, fp);
  } } end_foreach_vertex(); }
  fputs ("\t\n", fp);
  fputs ("\t </AppendedData>\n", fp);
  fputs ("</VTKFile>\n", fp);
  fflush (fp);
#ifdef _OPENMP
  omp_set_num_threads(num_omp);
#endif
 delete (((scalar []){marker,{-1}})); }
#line 15 "milk_crown.c"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#line 40 "milk_crown.c"
scalar f0= {20};



double z_shift, x_shift, x_left, y_left, z_left, x_right, y_right, z_right;

double drop_dx, drop_dy, release_height, u0_left, u0_right;

char parent_dir[256];

double RE_left, RE_right, WE_left, WE_right;


int main(int argc, char* argv[]) { _init_solver();
#line 68 "milk_crown.c"
  if(argc == 6){
    drop_dx = atof(argv[1]);
    drop_dy = atof(argv[2]);
    release_height = atof(argv[3]);
    u0_left = atof(argv[4]);
    u0_right = atof(argv[5]);
  }
  else{
    assert(argc == 4);

    drop_dx = atof(argv[1]);
    drop_dy = atof(argv[2]);
    release_height = atof(argv[3]);

    u0_left = 0.0;
    u0_right = 0.0;
  }

  assert(drop_dx > 0);
  assert(release_height > 0);
  assert(u0_left > -1E-15);
  assert(u0_right > -1E-15);


  assert(0.0035 < 0.053);
  assert(2*0.0035/2 + drop_dx < 0.053);

  assert(drop_dx > 0.0035);

  assert(release_height+0.0014001 +0.0035/2 < 0.053);
  assert(release_height-0.0035/2 > 0);

  if(drop_dy > 0){
    assert(release_height+0.0014001 +0.0035/2 + drop_dy < 0.053);
  }
  else{
    assert(release_height+drop_dy-0.0035/2 > 0.0014001);
  }


  z_shift = 0;
  x_shift = 0.5*(drop_dx - 0.0035);

  x_left = 0.053/2. - x_shift;
  y_left = release_height+0.0014001;
  z_left = 0.053/2.;

  x_right = 0.053/2. + x_shift;
  y_right = y_left + drop_dy;
  z_right = 0.053/2.;


  rho1 = 1025.8;
  rho2 = 1.225;
  mu1 = 0.0021;
  mu2 = 1.983E-5;
  const double sigma = 0.0467;
  _attribute[f.i].sigma = sigma;

    G.y = -9.81;






  DT = 0.1;
  TOLERANCE = 1e-2;
  CFL = 0.4;




  size (0.053);
  origin ((struct _origin){0., 0., 0.});
  init_grid (1 << 5);


  sprintf(parent_dir,
          "results_dx%g_dy%g_H%g_U0l%g_U0r%g",
          drop_dx, drop_dy, release_height, u0_left, u0_right);

  struct stat st = {0};
  if(stat(parent_dir, &st) == -1 && pid() == 0){
      mkdir(parent_dir, 0700);
  }


  const double t_imp_left = (-u0_left + sqrt(pow(u0_left, 2) + 2*9.81*(y_left-0.0014001)))/9.81;
  const double U_left = u0_left + t_imp_left*9.81;

  const double t_imp_right = (-u0_right + sqrt(pow(u0_right, 2) + 2*9.81*(y_right-0.0014001)))/9.81;
  const double U_right = u0_right + t_imp_right*9.81;

  RE_left = U_left*rho1*0.0035/mu1;
  RE_right = U_right*rho1*0.0035/mu1;

  WE_left = pow(U_left, 2)*rho1*0.0035/sigma;
  WE_right = pow(U_right, 2)*rho1*0.0035/sigma;


  run();
 free_solver(); }


static int init_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int init_0 (const int i, const double t, Event * _ev) { trace ("init_0", "milk_crown.c", 173);  {


  Csmag=0.12;

  do { int refined; do { refined = 0; ((Tree *)grid)->refined.n = 0;  { foreach_leaf(){

#line 178 "milk_crown.c"
 if (((fabs(1.-sqrt(pow(x-x_left, 2)+pow(y-y_left, 2)+pow(z-z_left, 2))/(0.0035/2.)) < 0.2) && (level < 10)) || ((fabs(1.-sqrt(pow(x-x_right, 2)+pow(y-y_right, 2)+pow(z-z_right, 2))/(0.0035/2.)) < 0.2) && (level < 10))) { refine_cell (point, all, 0, &((Tree *)grid)->refined); refined++; continue; } } end_foreach_leaf(); } mpi_all_reduce (refined, MPI_INT, MPI_SUM); if (refined) { mpi_boundary_refine (all); mpi_boundary_update (all); } } while (refined); } while(0)


   ;
  do { scalar phi= new_vertex_scalar("phi");  { foreach_vertex(){

#line 182 "milk_crown.c"
 val(phi,0,0,0) = max( max(0.0035/2.-sqrt(pow(x-x_left, 2)+pow(y-y_left, 2)+pow(z-z_left, 2)), 0.0035/2.-sqrt(pow(x-x_right, 2)+pow(y-y_right, 2)+pow(z-z_right, 2))), 0.0014001 -y); } end_foreach_vertex(); } fractions ((struct Fractions){phi, f0});  delete (((scalar []){phi,{-1}})); } while(0)


                ;

  _attribute[f0.i].refine = _attribute[f0.i].prolongation = fraction_refine;
  restriction (((scalar []){f0,{-1}}));

   { foreach(){

#line 190 "milk_crown.c"

    val(f,0,0,0) = val(f0,0,0,0); } end_foreach(); }
  boundary(((scalar []){f,f0,{-1}}));

   { foreach(){

#line 194 "milk_crown.c"
 {
      val(u.x,0,0,0) = 0.;
      val(u.y,0,0,0) = 0.;
      val(u.z,0,0,0) = 0.;

      if (sqrt(pow(x-x_left, 2) + pow(y-y_left, 2) + pow(z-z_left, 2)) < 1.1*0.0035/2.) {
        val(u.y,0,0,0) = -u0_left;
      }

      if (sqrt(pow(x-x_right, 2) + pow(y-y_right, 2) + pow(z-z_right, 2)) < 1.1*0.0035/2.) {
        val(u.y,0,0,0) = -u0_right;
      }
  } } end_foreach(); }
 end_trace("init_0", "milk_crown.c", 207); } return 0; } 


static int print_progress_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int print_progress (const int i, const double t, Event * _ev) { trace ("print_progress", "milk_crown.c", 210);  {
  if (pid()==0) {
    printf("Iteration %d, t = %g\n", i, t);
    }
 end_trace("print_progress", "milk_crown.c", 214); } return 0; } 


static int adapt_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int adapt_0 (const int i, const double t, Event * _ev) { trace ("adapt_0", "milk_crown.c", 217);  {

  adapt_wavelet ((struct Adapt){((scalar []){f,u.x,u.y,u.z,{-1}}), (double[]){0.025, 0.1, 0.1, 0.1}, 13});
 end_trace("adapt_0", "milk_crown.c", 220); } return 0; } 
#line 232 "milk_crown.c"
static int export_vtk_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (t=0.0001);   *ip = i; *tp = t;   return ret; } static int export_vtk_expr1 (int * ip, double * tp, Event * _ev) {   int i = *ip; double t = *tp;   int ret = ( t += 0.001);   *ip = i; *tp = t;   return ret; } static int export_vtk (const int i, const double t, Event * _ev) { trace ("export_vtk", "milk_crown.c", 232); 
{
  char vtk_dir[512];
  sprintf(vtk_dir, "%s/%s", parent_dir, "vtk");

  struct stat st = {0};
  if (stat(vtk_dir, &st) == -1 && pid() == 0){
    mkdir(vtk_dir, 0700);
  }


  scalar eddy_visc= new_scalar("eddy_visc");
   { 
if (!is_constant(mu.x)) {
#undef val_mu_x
#define val_mu_x(a,i,j,k) val(a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_y
#define val_mu_y(a,i,j,k) val(a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_z
#define val_mu_z(a,i,j,k) val(a,i,j,k)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,i,j,k)
#line 244
foreach(){

#line 244 "milk_crown.c"
 {
    val(eddy_visc,0,0,0) = 0. ;
    {
#line 246
 {
      val(eddy_visc,0,0,0) += (0.5/3)*(val_mu_x(mu.x,0,0,0) + val_mu_x(mu.x,1,0,0));
    }
#line 246
 {
      val(eddy_visc,0,0,0) += (0.5/3)*(val_mu_y(mu.y,0,0,0) + val_mu_y(mu.y,0,1,0));
    }
#line 246
 {
      val(eddy_visc,0,0,0) += (0.5/3)*(val_mu_z(mu.z,0,0,0) + val_mu_z(mu.z,0,0,1));
    }}
  } } end_foreach(); }
if (is_constant(mu.x)) {
const struct { double x, y, z; } _const_mu = {_constant[mu.x.i -_NVARMAX], _constant[mu.y.i - _NVARMAX], _constant[mu.z.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_x
#define val_mu_x(a,i,j,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,i,j,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_z
#define val_mu_z(a,i,j,k) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#line 244
foreach(){

#line 244 "milk_crown.c"
 {
    val(eddy_visc,0,0,0) = 0. ;
    {
#line 246
 {
      val(eddy_visc,0,0,0) += (0.5/3)*(val_mu_x(mu.x,0,0,0) + val_mu_x(mu.x,1,0,0));
    }
#line 246
 {
      val(eddy_visc,0,0,0) += (0.5/3)*(val_mu_y(mu.y,0,0,0) + val_mu_y(mu.y,0,1,0));
    }
#line 246
 {
      val(eddy_visc,0,0,0) += (0.5/3)*(val_mu_z(mu.z,0,0,0) + val_mu_z(mu.z,0,0,1));
    }}
  } } end_foreach(); } }

  static int j = 0;
  char cmd[1024];
  sprintf (cmd, "%s/snapshot_p%d_%.6i.vtu", vtk_dir, pid(), j++);
  FILE * fp = fopen(cmd, "w");


  output_vtu_bin_foreach ((scalar *) ((scalar []){f,p,{-1}}), (vector *) ((vector []){{u.x,u.y,u.z},{{-1},{-1},{-1}}}), N, fp, false);
  fclose (fp);
 delete (((scalar []){eddy_visc,{-1}}));  end_trace("export_vtk", "milk_crown.c", 259); } return 0; } 




static int output_interface_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (t=0.0001);   *ip = i; *tp = t;   return ret; } static int output_interface_expr1 (int * ip, double * tp, Event * _ev) {   int i = *ip; double t = *tp;   int ret = ( t += 0.0005);   *ip = i; *tp = t;   return ret; } static int output_interface (const int i, const double t, Event * _ev) { trace ("output_interface", "milk_crown.c", 264);  {
  printf("(left) RE = %g, WE= %g | (right) RE=%g, WE = %g", RE_left, WE_left, RE_right, WE_right);

  char iface_dir[512];
  sprintf(iface_dir, "%s/%s", parent_dir, "interface");

  struct stat st = {0};
  if(stat(iface_dir, &st) == -1){
    mkdir(iface_dir, 0700);
  }

  static int j = 0;
  char name[1024];
  sprintf (name, "%s/interface_p%d_%.6i.vtu", iface_dir, pid(), j++);
  FILE * fp = fopen (name, "w");
  output_vtu ((struct OutputFacets){f, fp});

  fclose (fp);
 end_trace("output_interface", "milk_crown.c", 282); } return 0; } 


static int end_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (t = 0.002);   *ip = i; *tp = t;   return ret; } static int end (const int i, const double t, Event * _ev) { trace ("end", "milk_crown.c", 285);  {
  printf ("Simulation ended");
  printf("(left) RE = %g, WE= %g | (right) RE=%g, WE = %g", RE_left, WE_left, RE_right, WE_right);
 end_trace("end", "milk_crown.c", 288); } return 0; } 
#line 77 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
static double _boundary0 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 76 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 77
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 77
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 77
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 77
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 77
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 77
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 77
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 77
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); } return 0.; } static double _boundary0_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 76 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 77
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 77
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 77
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 77
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 77
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 77
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 77
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 77
return  neumann_homogeneous(); } return 0.; }
#line 78 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
static double _boundary1 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 77 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 78
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 78
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 78
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 78
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 78
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 78
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 78
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 78
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); } return 0.; } static double _boundary1_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 77 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 78
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 78
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 78
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 78
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 78
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 78
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 78
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 78
return  neumann_homogeneous(); } return 0.; }
#line 84 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
static double _boundary2 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 83 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 84
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 84
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 84
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 84
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 84
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 84
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 84
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 84
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); } return 0.; } static double _boundary2_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 83 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 84
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 84
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 84
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 84
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 84
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 84
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 84
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 84
return  neumann_homogeneous(); } return 0.; }
#line 85 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
static double _boundary3 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 84 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 85
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 85
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 85
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 85
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 85
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 85
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 85
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 85
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); } return 0.; } static double _boundary3_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 84 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 85
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 85
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 85
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 85
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 85
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 85
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 85
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 85
return  neumann_homogeneous(); } return 0.; }
#line 88 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
static double _boundary4 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 87 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 88
return  neumann(val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 88
return  neumann(val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 88
return  neumann(val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 88
return  neumann(val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 88
return  neumann(val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 88
return  neumann(val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 88
return  neumann(val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 88
return  neumann(val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); } return 0.; } static double _boundary4_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 87 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 88
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 88
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 88
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 88
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 88
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 88
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 88
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 88
return  neumann_homogeneous(); } return 0.; }
#line 89 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"
static double _boundary5 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 88 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann(-val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0)); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann(-val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0)); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann(-val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0)); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann(-val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0)); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 89
return  neumann(-val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0)); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 89
return  neumann(-val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0)); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 89
return  neumann(-val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0)); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 89
return  neumann(-val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0)); } return 0.; } static double _boundary5_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 88 "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 89
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 89
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 89
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 89
return  neumann_homogeneous(); } return 0.; }
size_t datasize = 21*sizeof (double);
static int defaults (const int i, const double t, Event * _ev);
static int defaults_expr0 (int * ip, double * tp, Event * _ev);
static int init (const int i, const double t, Event * _ev);
static int init_expr0 (int * ip, double * tp, Event * _ev);
static int set_dtmax (const int i, const double t, Event * _ev);
static int set_dtmax_expr0 (int * ip, double * tp, Event * _ev);
static int stability (const int i, const double t, Event * _ev);
static int stability_expr0 (int * ip, double * tp, Event * _ev);
static int vof (const int i, const double t, Event * _ev);
static int vof_expr0 (int * ip, double * tp, Event * _ev);
static int tracer_advection (const int i, const double t, Event * _ev);
static int tracer_advection_expr0 (int * ip, double * tp, Event * _ev);
static int tracer_diffusion (const int i, const double t, Event * _ev);
static int tracer_diffusion_expr0 (int * ip, double * tp, Event * _ev);
static int properties (const int i, const double t, Event * _ev);
static int properties_expr0 (int * ip, double * tp, Event * _ev);
static int advection_term (const int i, const double t, Event * _ev);
static int advection_term_expr0 (int * ip, double * tp, Event * _ev);
static int viscous_term (const int i, const double t, Event * _ev);
static int viscous_term_expr0 (int * ip, double * tp, Event * _ev);
static int acceleration (const int i, const double t, Event * _ev);
static int acceleration_expr0 (int * ip, double * tp, Event * _ev);
static int projection (const int i, const double t, Event * _ev);
static int projection_expr0 (int * ip, double * tp, Event * _ev);
static int adapt (const int i, const double t, Event * _ev);
static int adapt_expr0 (int * ip, double * tp, Event * _ev);
static int defaults_0 (const int i, const double t, Event * _ev);
static int defaults_0_expr0 (int * ip, double * tp, Event * _ev);
static int stability_0 (const int i, const double t, Event * _ev);
static int stability_0_expr0 (int * ip, double * tp, Event * _ev);
static int vof_0 (const int i, const double t, Event * _ev);
static int vof_0_expr0 (int * ip, double * tp, Event * _ev);
static int defaults_1 (const int i, const double t, Event * _ev);
static int defaults_1_expr0 (int * ip, double * tp, Event * _ev);
static int properties_0 (const int i, const double t, Event * _ev);
static int properties_0_expr0 (int * ip, double * tp, Event * _ev);
static int defaults_2 (const int i, const double t, Event * _ev);
static int defaults_2_expr0 (int * ip, double * tp, Event * _ev);
static int vof_1 (const int i, const double t, Event * _ev);
static int vof_1_expr0 (int * ip, double * tp, Event * _ev);
static int tracer_advection_0 (const int i, const double t, Event * _ev);
static int tracer_advection_0_expr0 (int * ip, double * tp, Event * _ev);
static int defaults_3 (const int i, const double t, Event * _ev);
static int defaults_3_expr0 (int * ip, double * tp, Event * _ev);
static int tracer_advection_1 (const int i, const double t, Event * _ev);
static int tracer_advection_1_expr0 (int * ip, double * tp, Event * _ev);
static int tracer_diffusion_0 (const int i, const double t, Event * _ev);
static int tracer_diffusion_0_expr0 (int * ip, double * tp, Event * _ev);
static int defaults_4 (const int i, const double t, Event * _ev);
static int defaults_4_expr0 (int * ip, double * tp, Event * _ev);
static int viscous_term_0 (const int i, const double t, Event * _ev);
static int viscous_term_0_expr0 (int * ip, double * tp, Event * _ev);
static int defaults_5 (const int i, const double t, Event * _ev);
static int defaults_5_expr0 (int * ip, double * tp, Event * _ev);
static int acceleration_0 (const int i, const double t, Event * _ev);
static int acceleration_0_expr0 (int * ip, double * tp, Event * _ev);
static int acceleration_1 (const int i, const double t, Event * _ev);
static int acceleration_1_expr0 (int * ip, double * tp, Event * _ev);
static int stability_1 (const int i, const double t, Event * _ev);
static int stability_1_expr0 (int * ip, double * tp, Event * _ev);
static int acceleration_2 (const int i, const double t, Event * _ev);
static int acceleration_2_expr0 (int * ip, double * tp, Event * _ev);
static int init_0 (const int i, const double t, Event * _ev);
static int init_0_expr0 (int * ip, double * tp, Event * _ev);
static int print_progress (const int i, const double t, Event * _ev);
static int print_progress_expr0 (int * ip, double * tp, Event * _ev);
static int adapt_0 (const int i, const double t, Event * _ev);
static int adapt_0_expr0 (int * ip, double * tp, Event * _ev);
static int export_vtk (const int i, const double t, Event * _ev);
static int export_vtk_expr0 (int * ip, double * tp, Event * _ev);
static int export_vtk_expr1 (int * ip, double * tp, Event * _ev);
static int output_interface (const int i, const double t, Event * _ev);
static int output_interface_expr0 (int * ip, double * tp, Event * _ev);
static int output_interface_expr1 (int * ip, double * tp, Event * _ev);
static int end (const int i, const double t, Event * _ev);
static int end_expr0 (int * ip, double * tp, Event * _ev);
static void _set_boundary0 (void);
static void _set_boundary1 (void);
static void _set_boundary2 (void);
static void _set_boundary3 (void);
static void _set_boundary4 (void);
static void _set_boundary5 (void);
void _init_solver (void) {
  void init_solver();
  init_solver();
  Events = (Event *) pmalloc (sizeof (Event), __func__, __FILE__, __LINE__);
  Events[0].last = 1;
  event_register ((Event){ 0, 1, defaults, {defaults_expr0}, ((int *)0), ((double *)0),
    "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 96, "defaults"});
  event_register ((Event){ 0, 1, defaults_0, {defaults_0_expr0}, ((int *)0), ((double *)0),
    "/home/miro3/Documents/Programming/basilisk/src/vof.h", 51, "defaults"});
  event_register ((Event){ 0, 1, defaults_1, {defaults_1_expr0}, ((int *)0), ((double *)0),
    "/home/miro3/Documents/Programming/basilisk/src/two-phase.h", 25, "defaults"});
  event_register ((Event){ 0, 1, defaults_2, {defaults_2_expr0}, ((int *)0), ((double *)0),
    "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/conserving.h", 37, "defaults"});
  event_register ((Event){ 0, 1, defaults_3, {defaults_3_expr0}, ((int *)0), ((double *)0),
    "/home/miro3/Documents/Programming/basilisk/src/tracer.h", 25, "defaults"});
  event_register ((Event){ 0, 1, defaults_4, {defaults_4_expr0}, ((int *)0), ((double *)0),
    "./SGS.h", 44, "defaults"});
  event_register ((Event){ 0, 1, defaults_5, {defaults_5_expr0}, ((int *)0), ((double *)0),
    "/home/miro3/Documents/Programming/basilisk/src/iforce.h", 30, "defaults"});
  event_register ((Event){ 0, 1, init, {init_expr0}, ((int *)0), ((double *)0),
    "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 135, "init"});
  event_register ((Event){ 0, 1, init_0, {init_0_expr0}, ((int *)0), ((double *)0),
    "milk_crown.c", 173, "init"});
  event_register ((Event){ 0, 1, print_progress, {print_progress_expr0}, ((int *)0), ((double *)0),
    "milk_crown.c", 210, "print_progress"});
  event_register ((Event){ 0, 2, export_vtk, {export_vtk_expr0, export_vtk_expr1}, ((int *)0), ((double *)0),
    "milk_crown.c", 232, "export_vtk"});
  event_register ((Event){ 0, 2, output_interface, {output_interface_expr0, output_interface_expr1}, ((int *)0), ((double *)0),
    "milk_crown.c", 264, "output_interface"});
  event_register ((Event){ 0, 1, end, {end_expr0}, ((int *)0), ((double *)0),
    "milk_crown.c", 285, "end"});
  event_register ((Event){ 0, 1, set_dtmax, {set_dtmax_expr0}, ((int *)0), ((double *)0),
    "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 163, "set_dtmax"});
  event_register ((Event){ 0, 1, stability, {stability_expr0}, ((int *)0), ((double *)0),
    "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 165, "stability"});
  event_register ((Event){ 0, 1, stability_0, {stability_0_expr0}, ((int *)0), ((double *)0),
    "/home/miro3/Documents/Programming/basilisk/src/vof.h", 63, "stability"});
  event_register ((Event){ 0, 1, stability_1, {stability_1_expr0}, ((int *)0), ((double *)0),
    "/home/miro3/Documents/Programming/basilisk/src/tension.h", 36, "stability"});
  event_register ((Event){ 0, 1, vof, {vof_expr0}, ((int *)0), ((double *)0),
    "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 175, "vof"});
  event_register ((Event){ 0, 1, vof_0, {vof_0_expr0}, ((int *)0), ((double *)0),
    "/home/miro3/Documents/Programming/basilisk/src/vof.h", 313, "vof"});
  event_register ((Event){ 0, 1, vof_1, {vof_1_expr0}, ((int *)0), ((double *)0),
    "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/conserving.h", 106, "vof"});
  event_register ((Event){ 0, 1, tracer_advection, {tracer_advection_expr0}, ((int *)0), ((double *)0),
    "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 176, "tracer_advection"});
  event_register ((Event){ 0, 1, tracer_advection_0, {tracer_advection_0_expr0}, ((int *)0), ((double *)0),
    "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/conserving.h", 177, "tracer_advection"});
  event_register ((Event){ 0, 1, tracer_advection_1, {tracer_advection_1_expr0}, ((int *)0), ((double *)0),
    "/home/miro3/Documents/Programming/basilisk/src/tracer.h", 38, "tracer_advection"});
  event_register ((Event){ 0, 1, tracer_diffusion, {tracer_diffusion_expr0}, ((int *)0), ((double *)0),
    "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 177, "tracer_diffusion"});
  event_register ((Event){ 0, 1, tracer_diffusion_0, {tracer_diffusion_0_expr0}, ((int *)0), ((double *)0),
    "/home/miro3/Documents/Programming/basilisk/src/tracer.h", 45, "tracer_diffusion"});
  event_register ((Event){ 0, 1, properties, {properties_expr0}, ((int *)0), ((double *)0),
    "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 184, "properties"});
  event_register ((Event){ 0, 1, properties_0, {properties_0_expr0}, ((int *)0), ((double *)0),
    "/home/miro3/Documents/Programming/basilisk/src/two-phase.h", 59, "properties"});
  event_register ((Event){ 0, 1, advection_term, {advection_term_expr0}, ((int *)0), ((double *)0),
    "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 245, "advection_term"});
  event_register ((Event){ 0, 1, viscous_term, {viscous_term_expr0}, ((int *)0), ((double *)0),
    "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 275, "viscous_term"});
  event_register ((Event){ 0, 1, viscous_term_0, {viscous_term_0_expr0}, ((int *)0), ((double *)0),
    "./SGS.h", 86, "viscous_term"});
  event_register ((Event){ 0, 1, acceleration, {acceleration_expr0}, ((int *)0), ((double *)0),
    "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 310, "acceleration"});
  event_register ((Event){ 0, 1, acceleration_0, {acceleration_0_expr0}, ((int *)0), ((double *)0),
    "/home/miro3/Documents/Programming/basilisk/src/iforce.h", 44, "acceleration"});
  event_register ((Event){ 0, 1, acceleration_1, {acceleration_1_expr0}, ((int *)0), ((double *)0),
    "/home/miro3/Documents/Programming/basilisk/src/reduced.h", 36, "acceleration"});
  event_register ((Event){ 0, 1, acceleration_2, {acceleration_2_expr0}, ((int *)0), ((double *)0),
    "/home/miro3/Documents/Programming/basilisk/src/tension.h", 68, "acceleration"});
  event_register ((Event){ 0, 1, projection, {projection_expr0}, ((int *)0), ((double *)0),
    "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 325, "projection"});
  event_register ((Event){ 0, 1, adapt, {adapt_expr0}, ((int *)0), ((double *)0),
    "/home/miro3/Documents/Programming/basilisk/src/navier-stokes/centered.h", 360, "adapt"});
  event_register ((Event){ 0, 1, adapt_0, {adapt_0_expr0}, ((int *)0), ((double *)0),
    "milk_crown.c", 217, "adapt"});
  _attribute = (_Attributes *) pcalloc (datasize/sizeof(double), sizeof (_Attributes), __func__, __FILE__, __LINE__);
  all = (scalar *) pmalloc (sizeof (scalar)*22,__func__, __FILE__, __LINE__);
  for (int i = 0; i < 21; i++)
    all[i].i = i;
  all[21].i = -1;
  set_fpe();
  octree_methods();
  init_scalar ((scalar){20}, "f0");
  init_scalar ((scalar){19}, "Evis");
  init_face_vector ((vector){{16},{17},{18}}, "mu_t");
  init_scalar ((scalar){15}, "rhov");
  init_face_vector ((vector){{12},{13},{14}}, "alphav");
  init_scalar ((scalar){11}, "f");
  init_face_vector ((vector){{8},{9},{10}}, "uf");
  init_scalar ((scalar){7}, "pf");
  init_vector ((vector){{4},{5},{6}}, "g");
  init_vector ((vector){{1},{2},{3}}, "u");
  init_scalar ((scalar){0}, "p");
  init_const_scalar ((scalar){_NVARMAX+7}, "zeroc",  0.);
  init_const_scalar ((scalar){_NVARMAX+6}, "unity",  1.);
  init_const_vector ((vector){{_NVARMAX+3},{_NVARMAX+4},{_NVARMAX+5}}, "unityf", (double []) {1.,1.,1.});
  init_const_vector ((vector){{_NVARMAX+0},{_NVARMAX+1},{_NVARMAX+2}}, "zerof", (double []) {0.,0.,0.});
  _set_boundary0();
  _set_boundary1();
  _set_boundary2();
  _set_boundary3();
  _set_boundary4();
  _set_boundary5();
}
