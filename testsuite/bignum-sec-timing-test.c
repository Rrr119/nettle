#include "testutils.h"

#if WITH_HOGWEED

#include "bignum.h"

#include <time.h>

#define DIST_COUNT 3
#define ROUNDS 41
#define TRIM 6
#define TRIMMED_ROUNDS (ROUNDS - 2 * TRIM)
#define POWM_REPS 2
#define CMP_REPS 50000

/* A two-sided trimmed paired t-test with |t| <= 6, plus a 5%
   relative-difference guard to reject small but reproducible timing
   separations. Pairing and trimming keep frequency drift and scheduler
   outliers from making the acceptance test itself flaky. */
#define T_SQUARED_MAX 36.0
#define RELATIVE_DIFF_MAX 0.05

/* The compatibility threshold is measured against GMP's existing
   mpz_powm_sec entry point. The wrapper only adds a function call. */
#define POWM_PERF_MAX_RATIO 1.5

/* mini-gmp has no prior constant-time entry point. Keep the fixed-work
   implementation within a small factor of its variable-time mpz_powm. */
#define MINI_POWM_PERF_MAX_RATIO 4.0

struct timing_samples
{
  double op[ROUNDS];
};

static volatile unsigned long timing_sink;

static double
timestamp(void)
{
#if HAVE_CLOCK_GETTIME && defined CLOCK_PROCESS_CPUTIME_ID
  struct timespec ts;

  /* CPU time excludes scheduler preemption while retaining the
     instruction-count differences that the test is looking for. */
  if (clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts) != 0)
    die("clock_gettime failed\n");

  return ts.tv_sec + 1e-9 * ts.tv_nsec;
#elif HAVE_CLOCK_GETTIME && defined CLOCK_MONOTONIC
  struct timespec ts;

  if (clock_gettime(CLOCK_MONOTONIC, &ts) != 0)
    die("clock_gettime failed\n");

  return ts.tv_sec + 1e-9 * ts.tv_nsec;
#else
  return (double) clock() / CLOCKS_PER_SEC;
#endif
}

static int
compare_double(const void *ap, const void *bp)
{
  double a = *(const double *) ap;
  double b = *(const double *) bp;

  return (a > b) - (a < b);
}

static void
sample_stats(const struct timing_samples *samples,
	     double *mean, double *variance)
{
  double values[ROUNDS];
  double sum = 0;
  double sum_sq = 0;
  unsigned i;

  for (i = 0; i < ROUNDS; i++)
    values[i] = samples->op[i];
  qsort(values, ROUNDS, sizeof(*values), compare_double);

  for (i = TRIM; i < ROUNDS - TRIM; i++)
    {
      sum += values[i];
      sum_sq += values[i] * values[i];
    }

  *mean = sum / TRIMMED_ROUNDS;
  *variance = (sum_sq - sum * sum / TRIMMED_ROUNDS)
    / (TRIMMED_ROUNDS - 1);
  if (*variance < 0)
    *variance = 0;
}

static void
paired_stats(const struct timing_samples *left,
	     const struct timing_samples *right,
	     double *mean, double *variance)
{
  double differences[ROUNDS];
  double sum = 0;
  double sum_sq = 0;
  unsigned i;

  for (i = 0; i < ROUNDS; i++)
    differences[i] = left->op[i] - right->op[i];
  qsort(differences, ROUNDS, sizeof(*differences), compare_double);

  for (i = TRIM; i < ROUNDS - TRIM; i++)
    {
      sum += differences[i];
      sum_sq += differences[i] * differences[i];
    }

  *mean = sum / TRIMMED_ROUNDS;
  *variance = (sum_sq - sum * sum / TRIMMED_ROUNDS)
    / (TRIMMED_ROUNDS - 1);
  if (*variance < 0)
    *variance = 0;
}

static void
check_timing_pair(const char *operation,
		  const char *left_name, const char *right_name,
		  const struct timing_samples *left,
		  const struct timing_samples *right)
{
  double left_mean, left_var, right_mean, right_var;
  double diff, diff_var, abs_diff, average, t_squared;

  sample_stats(left, &left_mean, &left_var);
  sample_stats(right, &right_mean, &right_var);
  paired_stats(left, right, &diff, &diff_var);

  abs_diff = diff < 0 ? -diff : diff;
  average = (left_mean + right_mean) / 2;
  t_squared = diff * diff
    / (diff_var / TRIMMED_ROUNDS);

  printf("%s: %s vs %s: %.3f vs %.3f us, t^2=%.2f, "
	 "relative difference %.3f\n",
	 operation, left_name, right_name,
	 1e6 * left_mean, 1e6 * right_mean, t_squared,
	 average > 0 ? abs_diff / average : 0);

  if (t_squared > T_SQUARED_MAX
      || (average > 0 && abs_diff > RELATIVE_DIFF_MAX * average))
    {
      fprintf(stderr, "Timing separation detected for %s (%s vs %s).\n",
	      operation, left_name, right_name);
      FAIL();
    }
}

static void
init_modulus(mpz_t m)
{
  mpz_set_ui(m, 1);
  mpz_mul_2exp(m, m, 1023);
  mpz_add_ui(m, m, 0x1b3);
}

static void
init_exponents(gmp_randstate_t rands, mpz_t e[DIST_COUNT])
{
  unsigned i;

  for (i = 0; i < DIST_COUNT; i++)
    mpz_init(e[i]);

  mpz_urandomb(e[0], rands, 384);
  mpz_setbit(e[0], 383);

  /* Low Hamming weight, with the same bit length. */
  mpz_set_ui(e[1], 1);
  mpz_mul_2exp(e[1], e[1], 383);
  mpz_setbit(e[1], 1);

  /* High Hamming weight, with the same bit length. */
  mpz_set_ui(e[2], 1);
  mpz_mul_2exp(e[2], e[2], 384);
  mpz_sub_ui(e[2], e[2], 1);
}

static void
init_comparisons(gmp_randstate_t rands, mpz_t a, mpz_t b[DIST_COUNT])
{
  unsigned i;

  mpz_init(a);
  for (i = 0; i < DIST_COUNT; i++)
    mpz_init(b[i]);

  mpz_urandomb(a, rands, 1024);
  mpz_setbit(a, 1023);

  mpz_set(b[0], a);

  mpz_set(b[1], a);
  mpz_combit(b[1], 0);

  mpz_set(b[2], a);
  mpz_combit(b[2], 768);
}

static double
measure_powm(const mpz_t b, const mpz_t e, const mpz_t m,
	     mpz_t result, unsigned reps)
{
  double start = timestamp();
  unsigned i;

  for (i = 0; i < reps; i++)
    {
      nettle_mpz_powm_sec(result, b, e, m);
      timing_sink ^= mpz_getlimbn(result, 0);
    }

  return (timestamp() - start) / reps;
}

static double
measure_baseline_powm(const mpz_t b, const mpz_t e, const mpz_t m,
		      mpz_t result, unsigned reps)
{
  double start = timestamp();
  unsigned i;

  for (i = 0; i < reps; i++)
    {
#if NETTLE_USE_MINI_GMP
      mpz_powm(result, b, e, m);
#else
      mpz_powm_sec(result, b, e, m);
#endif
      timing_sink ^= mpz_getlimbn(result, 0);
    }

  return (timestamp() - start) / reps;
}

static double
measure_cmp(const mpz_t a, const mpz_t b, unsigned reps)
{
  double start = timestamp();
  unsigned i;

  for (i = 0; i < reps; i++)
    timing_sink += (unsigned long) nettle_mpz_cmp_sec(a, b);

  return (timestamp() - start) / reps;
}

static void
run_powm_timing(const mpz_t b, const mpz_t m, mpz_t e[DIST_COUNT])
{
  static const char *const names[DIST_COUNT] =
    { "uniform", "low-hamming", "high-hamming" };
  struct timing_samples samples[DIST_COUNT];
  struct timing_samples reference;
  mpz_t result;
  unsigned round;

  mpz_init(result);

  for (round = 0; round < ROUNDS; round++)
    {
      unsigned slot;

      for (slot = 0; slot < DIST_COUNT; slot++)
	{
	  unsigned dist = (round + slot) % DIST_COUNT;
	  samples[dist].op[round] = measure_powm(b, e[dist], m, result,
						 POWM_REPS);

	  if (dist == 0)
	    reference.op[round] = measure_baseline_powm(b, e[dist], m, result,
							POWM_REPS);
	}
    }

  check_timing_pair("powm", names[0], names[1],
		    &samples[0], &samples[1]);
  check_timing_pair("powm", names[0], names[2],
		    &samples[0], &samples[2]);
  check_timing_pair("powm", names[1], names[2],
		    &samples[1], &samples[2]);

  {
    double ct_mean, ct_var, ref_mean, ref_var, ratio;
#if NETTLE_USE_MINI_GMP
    const char *baseline = "mpz_powm";
    const double max_ratio = MINI_POWM_PERF_MAX_RATIO;
#else
    const char *baseline = "gmp mpz_powm_sec";
    const double max_ratio = POWM_PERF_MAX_RATIO;
#endif

    sample_stats(&samples[0], &ct_mean, &ct_var);
    sample_stats(&reference, &ref_mean, &ref_var);
    ratio = ref_mean > 0 ? ct_mean / ref_mean : 1;

    printf("powm performance: nettle %.3f us, %s %.3f us, "
	   "ratio %.3f (limit %.1f)\n",
	   1e6 * ct_mean, baseline, 1e6 * ref_mean, ratio, max_ratio);

    if (ratio > max_ratio)
      {
	fprintf(stderr, "Constant-time powm performance regression: %.3f\n",
		ratio);
	FAIL();
      }
  }

  mpz_clear(result);
}

static void
run_cmp_timing(const mpz_t a, mpz_t b[DIST_COUNT])
{
  static const char *const names[DIST_COUNT] =
    { "equal", "low-bit-difference", "high-bit-difference" };
  struct timing_samples samples[DIST_COUNT];
  unsigned round;

  for (round = 0; round < ROUNDS; round++)
    {
      unsigned slot;

      for (slot = 0; slot < DIST_COUNT; slot++)
	{
	  unsigned dist = (round + slot) % DIST_COUNT;
	  samples[dist].op[round] = measure_cmp(a, b[dist], CMP_REPS);
	}
    }

  check_timing_pair("cmp", names[0], names[1],
		    &samples[0], &samples[1]);
  check_timing_pair("cmp", names[0], names[2],
		    &samples[0], &samples[2]);
  check_timing_pair("cmp", names[1], names[2],
		    &samples[1], &samples[2]);
}

void
test_main(void)
{
  gmp_randstate_t rands;
  mpz_t base;
  mpz_t modulus;
  mpz_t exponents[DIST_COUNT];
  mpz_t a;
  mpz_t b[DIST_COUNT];
  unsigned i;

  test_randinit(rands);
  mpz_init(base);
  mpz_init(modulus);

  init_modulus(modulus);
  mpz_urandomb(base, rands, 1024);
  mpz_fdiv_r(base, base, modulus);

  init_exponents(rands, exponents);
  init_comparisons(rands, a, b);

  run_powm_timing(base, modulus, exponents);
  run_cmp_timing(a, b);

  for (i = 0; i < DIST_COUNT; i++)
    {
      mpz_clear(b[i]);
      mpz_clear(exponents[i]);
    }
  mpz_clear(a);
  mpz_clear(modulus);
  mpz_clear(base);
}

#else /* !WITH_HOGWEED */

void
test_main(void)
{
  SKIP();
}

#endif /* !WITH_HOGWEED */
