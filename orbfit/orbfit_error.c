/* Error reporting for the orbfit library.  Errors are recorded so the
 * Python caller can read them with orbfit_last_error(); orbfit_fail()
 * also unwinds to the guard set by the current entry point instead of
 * terminating the host process. */
#include <setjmp.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "orbfit.h"
#include "orbfit_api.h"

static jmp_buf *active_guard = NULL;
static char last_error[1024] = "";

static void
record(const char *fmt, va_list ap)
{
  vsnprintf(last_error, sizeof last_error, fmt, ap);
  fprintf(stderr, "%s\n", last_error);
}

void
orbfit_error(const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  record(fmt, ap);
  va_end(ap);
}

void
orbfit_rethrow(void)
{
  if (last_error[0] == 0) orbfit_error("orbfit: unspecified error");
  if (active_guard != NULL) longjmp(*active_guard, 1);
  exit(1);
}

void
orbfit_fail(const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  record(fmt, ap);
  va_end(ap);
  orbfit_rethrow();
}

void
orbfit_begin(jmp_buf *guard)
{
  last_error[0] = 0;
  active_guard = guard;
}

void
orbfit_end(void)
{
  active_guard = NULL;
}

const char *
orbfit_last_error(void)
{
  return last_error;
}
