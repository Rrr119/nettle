#include "testutils.h"

#if WITH_HOGWEED

#include "bignum.h"

static void
check_powm(const mpz_t b, const mpz_t e, const mpz_t m)
{
  mpz_t expected;
  mpz_t result;

  mpz_init(expected);
  mpz_init(result);

  mpz_powm(expected, b, e, m);
  nettle_mpz_powm_sec(result, b, e, m);
  ASSERT(mpz_cmp(result, expected) == 0);

  /* Check in-place operation, which is supported by the mpz_* API. */
  mpz_set(result, b);
  nettle_mpz_powm_sec(result, result, e, m);
  ASSERT(mpz_cmp(result, expected) == 0);

  mpz_clear(result);
  mpz_clear(expected);
}

static void
check_cmp(const mpz_t a, const mpz_t b)
{
  ASSERT(nettle_mpz_cmp_sec(a, b) == mpz_cmp(a, b));
  ASSERT(mpz_cmp_sec(a, b) == mpz_cmp(a, b));
  ASSERT(nettle_mpz_cmp_sec(b, a) == mpz_cmp(b, a));
  ASSERT(mpz_cmp_sec(b, a) == mpz_cmp(b, a));
}

void
test_main(void)
{
  gmp_randstate_t rands;
  mpz_t a;
  mpz_t b;
  mpz_t e;
  mpz_t m;
  mpz_t t;
  unsigned bits;

  test_randinit(rands);
  mpz_init(a);
  mpz_init(b);
  mpz_init(e);
  mpz_init(m);
  mpz_init(t);

  mpz_set_ui(b, 2);
  mpz_set_ui(e, 0);
  mpz_set_ui(m, 5);
  check_powm(b, e, m);

  mpz_set_ui(m, 1);
  check_powm(b, e, m);
  mpz_set_ui(e, 1);
  check_powm(b, e, m);

  mpz_set_si(b, -2);
  mpz_set_ui(e, 3);
  mpz_set_ui(m, 5);
  check_powm(b, e, m);

  mpz_set_ui(b, 11);
  mpz_set_ui(e, 3);
  mpz_set_ui(m, 5);
  check_powm(b, e, m);

  for (bits = 17; bits <= 521; bits += 24)
    {
      mpz_urandomb(b, rands, bits + 32);
      mpz_urandomb(e, rands, bits);
      mpz_urandomb(m, rands, bits);

      /* Make m positive and odd, with the requested bit size. */
      mpz_setbit(m, bits - 1);
      mpz_setbit(m, 0);

      check_powm(b, e, m);

      /* Exercise negative bases and bases larger than the modulus. */
      mpz_neg(t, b);
      check_powm(t, e, m);
      mpz_add(t, b, m);
      check_powm(t, e, m);
    }

  mpz_set_ui(a, 0);
  mpz_set_ui(b, 0);
  check_cmp(a, b);

  mpz_set_si(a, -1);
  check_cmp(a, b);
  mpz_set_ui(b, 1);
  check_cmp(a, b);
  check_cmp(b, a);

  mpz_set_si(a, -1234567);
  mpz_set_si(b, -1234568);
  check_cmp(a, b);
  check_cmp(b, a);

  mpz_set_str(a,
	      "ffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
	      "ffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
	      16);
  mpz_add_ui(b, a, 1);
  check_cmp(a, b);
  check_cmp(b, a);

  for (bits = 9; bits <= 521; bits += 23)
    {
      mpz_urandomb(a, rands, bits);
      mpz_setbit(a, bits - 1);
      mpz_set(b, a);
      check_cmp(a, b);

      mpz_combit(b, 0);
      check_cmp(a, b);
      check_cmp(b, a);

      mpz_set(b, a);
      mpz_combit(b, bits / 2);
      check_cmp(a, b);
      check_cmp(b, a);

      mpz_set(b, a);
      mpz_combit(b, bits - 2);
      check_cmp(a, b);
      check_cmp(b, a);
    }

  mpz_clear(t);
  mpz_clear(m);
  mpz_clear(e);
  mpz_clear(b);
  mpz_clear(a);
}

#else /* !WITH_HOGWEED */

void
test_main(void)
{
  SKIP();
}

#endif /* !WITH_HOGWEED */
