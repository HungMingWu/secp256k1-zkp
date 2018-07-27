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
