# R interface to some OpenMP runtime controls

# number of processors (hardware hyperthreads):
Romp_get_num_procs <- function() .C("Romp_get_num_procs",integer(1))[[1]]

# number of threads available, as specified by
# OMP_NUM_THREADS if set, otherwise omp_get_num_procs()
#
Romp_get_max_threads <- function() .C("Romp_get_max_threads",integer(1))[[1]]

# set max_threads to be used by subsequent parallel sections
#
Romp_set_num_threads <- function(nthreads)
  .C("Romp_set_num_threads",as.integer(nthreads))[[1]]

Romp_schedules <- c("static","dynamic","guided","auto")

# Return the kind of schedule and chunksize
#
Romp_get_schedule <- function()
{
  ret <- .C("Romp_get_schedule",
    kind=integer(1),chunk_size=integer(1))
  ret$kind <- Romp_schedules[ret$kind]
  ret
}

Romp_set_schedule <- function(kind,chunk_size=0)
{
  k <- pmatch(kind,Romp_schedules)
  if(is.na(k)) stop("unknown kind of schedule")
  ret <- .C("Romp_set_schedule",
    kind=as.integer(k),chunk_size=as.integer(chunk_size))
  ret$kind <- Romp_schedules[ret$kind]
  ret
}
