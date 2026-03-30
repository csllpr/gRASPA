#include "data_struct.h"

#include <cstdlib>
#include <cstdint>
#include <random>

namespace {
bool g_use_fast_host_rng = false;
uint64_t g_fast_host_rng_state[4] = {
  0x243f6a8885a308d3ULL,
  0x13198a2e03707344ULL,
  0xa4093822299f31d0ULL,
  0x082efa98ec4e6c89ULL,
};

inline uint64_t SplitMix64Step(uint64_t& state)
{
  uint64_t z = (state += 0x9e3779b97f4a7c15ULL);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

inline uint64_t RotL64(uint64_t x, int k)
{
  return (x << k) | (x >> (64 - k));
}

inline uint64_t NextFastHostRNG()
{
  const uint64_t result = RotL64(g_fast_host_rng_state[1] * 5ULL, 7) * 9ULL;
  const uint64_t t = g_fast_host_rng_state[1] << 17;

  g_fast_host_rng_state[2] ^= g_fast_host_rng_state[0];
  g_fast_host_rng_state[3] ^= g_fast_host_rng_state[1];
  g_fast_host_rng_state[1] ^= g_fast_host_rng_state[2];
  g_fast_host_rng_state[0] ^= g_fast_host_rng_state[3];
  g_fast_host_rng_state[2] ^= t;
  g_fast_host_rng_state[3] = RotL64(g_fast_host_rng_state[3], 45);

  return result;
}
}

void ConfigureHostUniformRNG(int seed, bool use_fast)
{
  g_use_fast_host_rng = use_fast;
  if(!g_use_fast_host_rng)
  {
    std::srand(seed);
    return;
  }

  uint64_t splitmix_state = static_cast<uint64_t>(static_cast<uint32_t>(seed)) + 0x9e3779b97f4a7c15ULL;
  for(size_t i = 0; i < 4; i++)
    g_fast_host_rng_state[i] = SplitMix64Step(splitmix_state);
}

double Get_Uniform_Random()
{
  if(g_use_fast_host_rng)
    return static_cast<double>(NextFastHostRNG() >> 11) * (1.0 / 9007199254740992.0);
  return static_cast<double>(std::rand()) / RAND_MAX;
}

// The polar form of the Box-Muller transformation, See Knuth v2, 3rd ed, p122, adapted from RASPA2
double Get_Gaussian_Random(void)
{
  double ran1,ran2,r2 = 0.0;

  while((r2>1.0)||(r2==0.0))
  {
    ran1=2.0*Get_Uniform_Random()-1.0;
    ran2=2.0*Get_Uniform_Random()-1.0;
    r2=pow(ran1,2)+pow(ran2,2);
  }
  return ran2*sqrt(-2.0*std::log(r2)/r2);
}
