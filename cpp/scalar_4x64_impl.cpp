#include <cstdint>
#include <boost/multiprecision/cpp_int.hpp>
#include "scalar_4x64.h"

using uint8_t = std::uint8_t;
using uint64_t = std::uint64_t;
using uint128_t = boost::multiprecision::uint128_t;

/* Limbs of the secp256k1 order. */
static constexpr uint64_t SECP256K1_N_0 = 0xBFD25E8CD0364141ULL;
static constexpr uint64_t SECP256K1_N_1 = 0xBAAEDCE6AF48A03BULL;
static constexpr uint64_t SECP256K1_N_2 = 0xFFFFFFFFFFFFFFFEULL;
static constexpr uint64_t SECP256K1_N_3 = 0xFFFFFFFFFFFFFFFFULL;

/* Limbs of 2^256 minus the secp256k1 order. */
static constexpr uint64_t SECP256K1_N_C_0 = ~SECP256K1_N_0 + 1;
static constexpr uint64_t SECP256K1_N_C_1 = ~SECP256K1_N_1;
static constexpr uint64_t SECP256K1_N_C_2 = 1;

static int secp256k1_scalar_check_overflow(const secp256k1_scalar_t *a) {
	int yes = 0;
	int no = 0;
	no |= (a->d[3] < SECP256K1_N_3); /* No need for a > check. */
	no |= (a->d[2] < SECP256K1_N_2);
	yes |= (a->d[2] > SECP256K1_N_2) & ~no;
	no |= (a->d[1] < SECP256K1_N_1);
	yes |= (a->d[1] > SECP256K1_N_1) & ~no;
	yes |= (a->d[0] >= SECP256K1_N_0) & ~no;
	return yes;
}

extern "C"
void secp256k1_scalar_add_bit_cpp(secp256k1_scalar_t *r, uint8_t bit)
{
	uint128_t t;
	t = (uint128_t)r->d[0] + (((uint64_t)((bit >> 6) == 0)) << (bit & 0x3F));
	r->d[0] = static_cast<uint64_t>(t & 0xFFFFFFFFFFFFFFFFULL); t >>= 64;
	t += (uint128_t)r->d[1] + (((uint64_t)((bit >> 6) == 1)) << (bit & 0x3F));
	r->d[1] = static_cast<uint64_t>(t & 0xFFFFFFFFFFFFFFFFULL); t >>= 64;
	t += (uint128_t)r->d[2] + (((uint64_t)((bit >> 6) == 2)) << (bit & 0x3F));
	r->d[2] = static_cast<uint64_t>(t & 0xFFFFFFFFFFFFFFFFULL); t >>= 64;
	t += (uint128_t)r->d[3] + (((uint64_t)((bit >> 6) == 3)) << (bit & 0x3F));
	r->d[3] = static_cast<uint64_t>(t & 0xFFFFFFFFFFFFFFFFULL);
}

extern "C"
int secp256k1_scalar_reduce_cpp(secp256k1_scalar_t *r, unsigned int overflow)
{
	uint128_t t;
	t = (uint128_t)r->d[0] + overflow * SECP256K1_N_C_0;
	r->d[0] = static_cast<uint64_t>(t & 0xFFFFFFFFFFFFFFFFULL); t >>= 64;
	t += (uint128_t)r->d[1] + overflow * SECP256K1_N_C_1;
	r->d[1] = static_cast<uint64_t>(t & 0xFFFFFFFFFFFFFFFFULL); t >>= 64;
	t += (uint128_t)r->d[2] + overflow * SECP256K1_N_C_2;
	r->d[2] = static_cast<uint64_t>(t & 0xFFFFFFFFFFFFFFFFULL); t >>= 64;
	t += (uint64_t)r->d[3];
	r->d[3] = static_cast<uint64_t>(t & 0xFFFFFFFFFFFFFFFFULL);
	return overflow;
}

extern "C"
int secp256k1_scalar_add_cpp(secp256k1_scalar_t *r, const secp256k1_scalar_t *a, const secp256k1_scalar_t *b)
{
	int overflow;
	uint128_t t = (uint128_t)a->d[0] + b->d[0];
	r->d[0] = static_cast<uint64_t>(t & 0xFFFFFFFFFFFFFFFFULL); t >>= 64;
	t += (uint128_t)a->d[1] + b->d[1];
	r->d[1] = static_cast<uint64_t>(t & 0xFFFFFFFFFFFFFFFFULL); t >>= 64;
	t += (uint128_t)a->d[2] + b->d[2];
	r->d[2] = static_cast<uint64_t>(t & 0xFFFFFFFFFFFFFFFFULL); t >>= 64;
	t += (uint128_t)a->d[3] + b->d[3];
	r->d[3] = static_cast<uint64_t>(t & 0xFFFFFFFFFFFFFFFFULL); t >>= 64;
	overflow = static_cast<int>(t) + secp256k1_scalar_check_overflow(r);
	secp256k1_scalar_reduce_cpp(r, overflow);
	return overflow;
}

static int secp256k1_scalar_is_zero(const secp256k1_scalar_t *a) {
	return (a->d[0] | a->d[1] | a->d[2] | a->d[3]) == 0;
}

extern "C"
int secp256k1_scalar_wnaf_force_odd_cpp(secp256k1_scalar_t *r)
{
	/* If we are odd, mask = 0 and this is a no-op;
	* if we are even, mask = 11...11 and this is identical to secp256k1_scalar_negate */
	uint64_t mask = (r->d[0] & 1) - 1;
	uint64_t nonzero = (secp256k1_scalar_is_zero(r) != 0) - 1;
	uint128_t t = (uint128_t)(r->d[0] ^ mask) + ((SECP256K1_N_0 + 1) & mask);
	r->d[0] = static_cast<uint64_t>(t & nonzero); t >>= 64;
	t += (uint128_t)(r->d[1] ^ mask) + (SECP256K1_N_1 & mask);
	r->d[1] = static_cast<uint64_t>(t & nonzero); t >>= 64;
	t += (uint128_t)(r->d[2] ^ mask) + (SECP256K1_N_2 & mask);
	r->d[2] = static_cast<uint64_t>(t & nonzero); t >>= 64;
	t += (uint128_t)(r->d[3] ^ mask) + (SECP256K1_N_3 & mask);
	r->d[3] = static_cast<uint64_t>(t & nonzero);
	return 2 * (mask == 0) - 1;
}

extern "C"
void secp256k1_scalar_negate_cpp(secp256k1_scalar_t *r, const secp256k1_scalar_t *a)
{
	uint64_t nonzero = 0xFFFFFFFFFFFFFFFFULL * (secp256k1_scalar_is_zero(a) == 0);
	uint128_t t = (uint128_t)(~a->d[0]) + SECP256K1_N_0 + 1;
	r->d[0] = static_cast<uint64_t>(t & nonzero); t >>= 64;
	t += (uint128_t)(~a->d[1]) + SECP256K1_N_1;
	r->d[1] = static_cast<uint64_t>(t & nonzero); t >>= 64;
	t += (uint128_t)(~a->d[2]) + SECP256K1_N_2;
	r->d[2] = static_cast<uint64_t>(t & nonzero); t >>= 64;
	t += (uint128_t)(~a->d[3]) + SECP256K1_N_3;
	r->d[3] = static_cast<uint64_t>(t & nonzero);
}

/** Add a*b to the number defined by (c0,c1). c1 must never overflow. */
#define muladd_fast(a,b) { \
    uint32_t tl, th; \
    { \
        uint64_t t = (uint64_t)a * b; \
        th = t >> 32;         /* at most 0xFFFFFFFE */ \
        tl = t; \
    } \
    c0 += tl;                 /* overflow is handled on the next line */ \
    th += (c0 < tl) ? 1 : 0;  /* at most 0xFFFFFFFF */ \
    c1 += th;                 /* never overflows by contract (verified in the next line) */ \
}

/** Extract the lowest 64 bits of (c0,c1,c2) into n, and left shift the number 64 bits. c2 is required to be zero. */
#define extract_fast(n) { \
    (n) = c0; \
    c0 = c1; \
    c1 = 0; \
}

/** Add a to the number defined by (c0,c1). c1 must never overflow, c2 must be zero. */
#define sumadd_fast(a) { \
    c0 += (a);                 /* overflow is handled on the next line */ \
    c1 += (c0 < (a)) ? 1 : 0;  /* never overflows by contract (verified the next line) */ \
}

/** Add a*b to the number defined by (c0,c1,c2). c2 must never overflow. */
#define muladd(a,b) { \
    uint64_t tl, th; \
    { \
        uint128_t t = (uint128_t)a * b; \
        th = static_cast<uint64_t>(t >> 64);         /* at most 0xFFFFFFFFFFFFFFFE */ \
        tl = static_cast<uint64_t>(t); \
    } \
    c0 += tl;                 /* overflow is handled on the next line */ \
    th += (c0 < tl) ? 1 : 0;  /* at most 0xFFFFFFFFFFFFFFFF */ \
    c1 += th;                 /* overflow is handled on the next line */ \
    c2 += (c1 < th) ? 1 : 0;  /* never overflows by contract (verified in the next line) */ \
}

/** Extract the lowest 64 bits of (c0,c1,c2) into n, and left shift the number 64 bits. */
#define extract(n) { \
    (n) = c0; \
    c0 = c1; \
    c1 = c2; \
    c2 = 0; \
}

/** Add a to the number defined by (c0,c1,c2). c2 must never overflow. */
#define sumadd(a) { \
    unsigned int over; \
    c0 += (a);                  /* overflow is handled on the next line */ \
    over = (c0 < (a)) ? 1 : 0; \
    c1 += over;                 /* overflow is handled on the next line */ \
    c2 += (c1 < over) ? 1 : 0;  /* never overflows by contract */ \
}

#if 0
extern "C"
void secp256k1_scalar_reduce_512_cpp(secp256k1_scalar_t *r, const uint64_t *l)
{
	std::cout << "Reduce " << r->d[0] << " " << r->d[1] << " " << r->d[2] << " " << r->d[3] << "\n";
	std::cout << "l " << l[0] << " " << l[1] << " " << l[2] << " " << l[3] << "\n" <<
		             l[4] << " " << l[5] << " " << l[6] << " " << l[7] << "\n";
	uint128_t c;
	uint64_t c0, c1, c2;
	uint64_t n0 = l[4], n1 = l[5], n2 = l[6], n3 = l[7];
	uint64_t m0, m1, m2, m3, m4, m5;
	uint32_t m6;
	uint64_t p0, p1, p2, p3;
	uint32_t p4;

	/* Reduce 512 bits into 385. */
	/* m[0..6] = l[0..3] + n[0..3] * SECP256K1_N_C. */
	c0 = l[0]; c1 = 0; c2 = 0;
	muladd_fast(n0, SECP256K1_N_C_0);
	extract_fast(m0);
	sumadd_fast(l[1]);
	muladd(n1, SECP256K1_N_C_0);
	muladd(n0, SECP256K1_N_C_1);
	extract(m1);
	sumadd(l[2]);
	muladd(n2, SECP256K1_N_C_0);
	muladd(n1, SECP256K1_N_C_1);
	sumadd(n0);
	extract(m2);
	sumadd(l[3]);
	muladd(n3, SECP256K1_N_C_0);
	muladd(n2, SECP256K1_N_C_1);
	sumadd(n1);
	extract(m3);
	muladd(n3, SECP256K1_N_C_1);
	sumadd(n2);
	extract(m4);
	sumadd_fast(n3);
	extract_fast(m5);
	m6 = c0;

	/* Reduce 385 bits into 258. */
	/* p[0..4] = m[0..3] + m[4..6] * SECP256K1_N_C. */
	c0 = m0; c1 = 0; c2 = 0;
	muladd_fast(m4, SECP256K1_N_C_0);
	extract_fast(p0);
	sumadd_fast(m1);
	muladd(m5, SECP256K1_N_C_0);
	muladd(m4, SECP256K1_N_C_1);
	extract(p1);
	sumadd(m2);
	muladd(m6, SECP256K1_N_C_0);
	muladd(m5, SECP256K1_N_C_1);
	sumadd(m4);
	extract(p2);
	sumadd_fast(m3);
	muladd_fast(m6, SECP256K1_N_C_1);
	sumadd_fast(m5);
	extract_fast(p3);
	p4 = c0 + m6;

	/* Reduce 258 bits into 256. */
	/* r[0..3] = p[0..3] + p[4] * SECP256K1_N_C. */
	c = p0 + (uint128_t)SECP256K1_N_C_0 * p4;
	r->d[0] = static_cast<uint64_t>(c & 0xFFFFFFFFFFFFFFFFULL); c >>= 64;
	c += p1 + (uint128_t)SECP256K1_N_C_1 * p4;
	r->d[1] = static_cast<uint64_t>(c & 0xFFFFFFFFFFFFFFFFULL); c >>= 64;
	c += p2 + (uint128_t)p4;
	r->d[2] = static_cast<uint64_t>(c & 0xFFFFFFFFFFFFFFFFULL); c >>= 64;
	c += p3;
	r->d[3] = static_cast<uint64_t>(c & 0xFFFFFFFFFFFFFFFFULL); c >>= 64;

	/* Final reduction of r. */
	secp256k1_scalar_reduce_cpp(r, static_cast<unsigned int>(c) + secp256k1_scalar_check_overflow(r));
	std::cout << "Final Reduce " << r->d[0] << " " << r->d[1] << " " << r->d[2] << " " << r->d[3] << "\n";
}
#endif
