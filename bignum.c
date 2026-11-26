/* bignum.c

   Bignum operations that are missing from gmp.

   Copyright (C) 2001 Niels Möller

   This file is part of GNU Nettle.

   GNU Nettle is free software: you can redistribute it and/or
   modify it under the terms of either:

     * the GNU Lesser General Public License as published by the Free
       Software Foundation; either version 3 of the License, or (at your
       option) any later version.

   or

     * the GNU General Public License as published by the Free
       Software Foundation; either version 2 of the License, or (at your
       option) any later version.

   or both in parallel, as here.

   GNU Nettle is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received copies of the GNU General Public License and
   the GNU Lesser General Public License along with this program.  If
   not, see http://www.gnu.org/licenses/.
*/

#if HAVE_CONFIG_H
# include "config.h"
#endif

#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

#include "bignum.h"

/* Two's complement negation means that -x = ~x + 1, ~x = -(x+1),
 * and we use that x = ~~x = ~(-x-1).
 *
 * Examples:
 *
 *   x  ~x = -x+1     ~~x = x
 *  -1          0          ff
 *  -2          1          fe
 * -7f         7e          81
 * -80         7f          80
 * -81         80        ff7f
 */

/* Including extra sign bit, if needed. Also one byte for zero. */
size_t
nettle_mpz_sizeinbase_256_s(const mpz_t x)
{
  if (mpz_sgn(x) >= 0)
    return 1 + mpz_sizeinbase(x, 2) / 8;
  else
    {
      /* We'll output ~~x, so we need as many bits as for ~x */
      size_t size;
      mpz_t c;

      mpz_init(c);
      mpz_com(c, x); /* Same as c = - x - 1 = |x| + 1 */
      size = 1 + mpz_sizeinbase(c,2) / 8;
      mpz_clear(c);

      return size;
    }
}

size_t
nettle_mpz_sizeinbase_256_u(const mpz_t x)
{
  return (mpz_sizeinbase(x,2) + 7) / 8;
}

static void
nettle_mpz_to_octets(size_t length, uint8_t *s,
		     const mpz_t x, uint8_t sign)
{
  uint8_t *dst = s + length - 1;
  size_t size = mpz_size(x);
  size_t i;
  
  for (i = 0; i<size; i++)
    {
      mp_limb_t limb = mpz_getlimbn(x, i);
      size_t j;

      for (j = 0; length && j < sizeof(mp_limb_t); j++)
        {
          *dst-- = sign ^ (limb & 0xff);
          limb >>= 8;
          length--;
	}
    }
  
  if (length)
    memset(s, sign, length);
}

void
nettle_mpz_get_str_256(size_t length, uint8_t *s, const mpz_t x)
{
  if (!length)
    {
      /* x must be zero */
      assert(!mpz_sgn(x));
      return;
    }

  if (mpz_sgn(x) >= 0)
    {
      assert(nettle_mpz_sizeinbase_256_u(x) <= length);
      nettle_mpz_to_octets(length, s, x, 0);
    }
  else
    {
      mpz_t c;
      mpz_init(c);
      mpz_com(c, x);

      assert(nettle_mpz_sizeinbase_256_u(c) <= length);
      nettle_mpz_to_octets(length, s, c, 0xff);

      mpz_clear(c);
    }
}

/* Converting from strings */

/* mpz_import was introduced in GMP-4.1 */
#define nettle_mpz_from_octets(x, length, s) \
   mpz_import((x), (length), 1, 1, 0, 0, (s))

void
nettle_mpz_set_str_256_u(mpz_t x,
			 size_t length, const uint8_t *s)
{
  nettle_mpz_from_octets(x, length, s);
}

void
nettle_mpz_init_set_str_256_u(mpz_t x,
			      size_t length, const uint8_t *s)
{
  mpz_init(x);
  nettle_mpz_from_octets(x, length, s);
}

void
nettle_mpz_set_str_256_s(mpz_t x,
			 size_t length, const uint8_t *s)
{
  if (!length)
    {
      mpz_set_ui(x, 0);
      return;
    }
  
  nettle_mpz_from_octets(x, length, s);

  if (s[0] & 0x80)
    {
      mpz_t t;

      mpz_init_set_ui(t, 1);
      mpz_mul_2exp(t, t, length*8);
      mpz_sub(x, x, t);
      mpz_clear(t);
    }
}

void
nettle_mpz_init_set_str_256_s(mpz_t x,
			      size_t length, const uint8_t *s)
{
  mpz_init(x);
  nettle_mpz_set_str_256_s(x, length, s);
}

int
nettle_mpz_cmp_sec(const mpz_t a, const mpz_t b)
{
  mp_size_t an = mpz_size(a);
  mp_size_t bn = mpz_size(b);
  mp_size_t n = an > bn ? an : bn;
  mp_limb_t greater = 0;
  mp_limb_t less = 0;
  mp_size_t i;

  for (i = n; i-- > 0;)
    {
      mp_limb_t av = mpz_getlimbn(a, i);
      mp_limb_t bv = mpz_getlimbn(b, i);
      mp_limb_t undecided = (greater | less) ^ 1;

      greater |= undecided & (mp_limb_t) (av > bv);
      less |= undecided & (mp_limb_t) (av < bv);
    }

  {
    int sa = (a[0]._mp_size > 0) - (a[0]._mp_size < 0);
    int sb = (b[0]._mp_size > 0) - (b[0]._mp_size < 0);
    int sign_diff = (sa > sb) - (sa < sb);
    int same_sign_mask = - (int) (sa == sb);

    return (same_sign_mask & (sa * ((int) greater - (int) less)))
      | (~same_sign_mask & sign_diff);
  }
}

#if NETTLE_USE_MINI_GMP

#define CT_LIMB_BITS ((unsigned) (sizeof(mp_limb_t) * CHAR_BIT))

static void *
ct_alloc(size_t n)
{
  void *p = malloc(n);
  if (!p)
    abort();

  return p;
}

static void
ct_limbs_zero(mp_limb_t *rp, mp_size_t n)
{
  mp_size_t i;

  for (i = 0; i < n; i++)
    rp[i] = 0;
}

/* r = (2*x + bit) mod m. r and x may be the same array. */
static void
ct_mod_double_add(mp_limb_t *rp, const mp_limb_t *xp, mp_limb_t bit,
		  const mp_limb_t *mp, mp_size_t n, mp_limb_t *tp)
{
  mp_limb_t carry = bit;
  mp_limb_t borrow;
  mp_limb_t use_sub;
  mp_limb_t mask;
  mp_size_t i;

  for (i = 0; i < n; i++)
    {
      mp_limb_t x = xp[i];
      rp[i] = (x << 1) | carry;
      carry = x >> (CT_LIMB_BITS - 1);
    }

  borrow = mpn_sub_n(tp, rp, mp, n);
  use_sub = carry | (borrow ^ 1);
  mask = -use_sub;

  for (i = 0; i < n; i++)
    rp[i] = (rp[i] & ~mask) | (tp[i] & mask);
}

/* r = x mod m, processing a fixed number of input limbs. */
static void
ct_mod_reduce(mp_limb_t *rp, const mp_limb_t *xp, mp_size_t xn,
	      const mp_limb_t *mp, mp_size_t n, mp_limb_t *tp)
{
  mp_size_t i;
  unsigned j;

  ct_limbs_zero(rp, n);

  for (i = xn; i-- > 0;)
    for (j = CT_LIMB_BITS; j-- > 0;)
      ct_mod_double_add(rp, rp, (xp[i] >> j) & 1, mp, n, tp);
}

/* r = x * 2^(n * GMP limb bits) mod m. */
static void
ct_to_mont(mp_limb_t *rp, const mp_limb_t *mp, mp_size_t n, mp_limb_t *tp)
{
  mp_size_t i;
  unsigned j;

  for (i = 0; i < n; i++)
    for (j = 0; j < CT_LIMB_BITS; j++)
      ct_mod_double_add(rp, rp, 0, mp, n, tp);
}

static mp_limb_t
ct_limb_inverse(mp_limb_t x)
{
  mp_limb_t inv = x;
  unsigned i;

  /* Newton iteration doubles the number of correct low bits. */
  for (i = 0; i < 6; i++)
    inv *= 2 - x * inv;

  return inv;
}

/* Montgomery multiplication. All loops have lengths determined by n. */
static void
ct_mont_mul(mp_limb_t *rp, const mp_limb_t *ap, const mp_limb_t *bp,
	    const mp_limb_t *mp, mp_limb_t m0inv, mp_size_t n,
	    mp_limb_t *tp, mp_limb_t *sub)
{
  mp_limb_t borrow;
  mp_limb_t use_sub;
  mp_limb_t mask;
  mp_size_t i;
  mp_size_t j;

  mpn_mul_n(tp, ap, bp, n);
  tp[2 * n] = 0;

  for (i = 0; i < n; i++)
    {
      mp_limb_t u = tp[i] * m0inv;
      mp_limb_t carry = mpn_addmul_1(tp + i, mp, n, u);

      for (j = i + n; j <= 2 * n; j++)
	{
	  mp_limb_t v = tp[j] + carry;
	  carry = v < carry;
	  tp[j] = v;
	}
    }

  borrow = mpn_sub_n(sub, tp + n, mp, n);
  use_sub = (tp[2 * n] != 0) | (borrow ^ 1);
  mask = -use_sub;

  for (i = 0; i < n; i++)
    rp[i] = (tp[n + i] & ~mask) | (sub[i] & mask);
}

static void
ct_powm_mini(mpz_t r, const mpz_t b, const mpz_t e, const mpz_t m)
{
  mp_size_t n = mpz_size(m);
  mp_size_t bn = mpz_size(b);
  mp_size_t en = mpz_size(e);
  mp_size_t i;
  mp_size_t k;
  unsigned j;
  mp_limb_t *mp;
  mp_limb_t *bp;
  mp_limb_t *base;
  mp_limb_t *one;
  mp_limb_t *acc;
  mp_limb_t *prod;
  mp_limb_t *ep;
  mp_limb_t *tp;
  mp_limb_t *sub;
  mp_limb_t m0inv;
  mp_limb_t negative;
  mp_limb_t nonzero;
  mp_limb_t mask;

  if (n == 0)
    abort();

  mp = ct_alloc(n * sizeof(*mp));
  bp = ct_alloc((bn ? bn : 1) * sizeof(*bp));
  base = ct_alloc(n * sizeof(*base));
  one = ct_alloc(n * sizeof(*one));
  acc = ct_alloc(n * sizeof(*acc));
  prod = ct_alloc(n * sizeof(*prod));
  ep = ct_alloc((en ? en : 1) * sizeof(*ep));
  tp = ct_alloc((2 * n + 1) * sizeof(*tp));
  sub = ct_alloc(n * sizeof(*sub));

  for (i = 0; i < n; i++)
    mp[i] = mpz_getlimbn(m, i);
  for (i = 0; i < bn; i++)
    bp[i] = mpz_getlimbn(b, i);
  for (i = 0; i < en; i++)
    ep[i] = mpz_getlimbn(e, i);

  ct_mod_reduce(base, bp, bn, mp, n, sub);

  /* Map a negative input to the corresponding non-negative residue. */
  negative = (mp_limb_t) (b[0]._mp_size < 0);
  nonzero = 0;
  for (i = 0; i < n; i++)
    nonzero |= base[i];
  mask = - (negative & (mp_limb_t) (nonzero != 0));
  mpn_sub_n(sub, mp, base, n);
  for (i = 0; i < n; i++)
    base[i] = (base[i] & ~mask) | (sub[i] & mask);

  ct_limbs_zero(one, n);
  one[0] = 1;
  ct_to_mont(one, mp, n, sub);

  ct_to_mont(base, mp, n, sub);
  for (i = 0; i < n; i++)
    acc[i] = one[i];

  m0inv = -ct_limb_inverse(mp[0]);

  for (i = en; i-- > 0;)
    for (j = CT_LIMB_BITS; j-- > 0;)
      {
	mp_limb_t bit = (ep[i] >> j) & 1;
	mp_limb_t bit_mask;

	ct_mont_mul(acc, acc, acc, mp, m0inv, n, tp, sub);
	ct_mont_mul(prod, acc, base, mp, m0inv, n, tp, sub);

	bit_mask = -bit;
	for (k = 0; k < n; k++)
	  acc[k] = (acc[k] & ~bit_mask) | (prod[k] & bit_mask);
      }

  ct_limbs_zero(one, n);
  one[0] = 1;
  ct_mont_mul(prod, acc, one, mp, m0inv, n, tp, sub);

  mpz_import(r, n, -1, sizeof(mp_limb_t), 0, 0, prod);

  free(sub);
  free(tp);
  free(ep);
  free(prod);
  free(acc);
  free(one);
  free(base);
  free(bp);
  free(mp);
}

void
nettle_mpz_powm_sec(mpz_t r, const mpz_t b, const mpz_t e, const mpz_t m)
{
  if (mpz_sgn(m) <= 0 || mpz_even_p(m) || mpz_sgn(e) < 0)
    {
      /* Match mpz_powm's behaviour for inputs outside the
	 constant-time domain. */
      mpz_powm(r, b, e, m);
      return;
    }

  ct_powm_mini(r, b, e, m);
}

#else /* !NETTLE_USE_MINI_GMP */

void
nettle_mpz_powm_sec(mpz_t r, const mpz_t b, const mpz_t e, const mpz_t m)
{
  mpz_powm_sec(r, b, e, m);
}

#endif /* !NETTLE_USE_MINI_GMP */
