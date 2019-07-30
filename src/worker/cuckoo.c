#include <stdint.h>    // for types uint32_t,uint64_t
#include "portable_endian.h"    // for htole32/64

#define rotl(x, b) ((x) << (b)) | ((x) >> (64 - (b)))
#define EBIT 15
#define CLEN 12

#define EN 1 << EBIT
#define M EN << 1
#define MASK (1 << EBIT) - 1

// set siphash keys from 32 byte char array
#define setkeys() \
  k0 = htole64(((uint64_t *)mesg)[0]); \
  k1 = htole64(((uint64_t *)mesg)[1]); \
  k2 = htole64(((uint64_t *)mesg)[2]); \
  k3 = htole64(((uint64_t *)mesg)[3]);


#define sip_round() \
  v0 += v1; v2 += v3; v1 = rotl(v1,13); \
  v3 = rotl(v3,16); v1 ^= v0; v3 ^= v2; \
  v0 = rotl(v0,32); v2 += v1; v0 += v3; \
  v1 = rotl(v1,17); v3 = rotl(v3,21); \
  v1 ^= v2; v3 ^= v0; v2 = rotl(v2,32); 

#define siphash24( nonce ) ({\
  v0 = k0; v1 = k1; v2 = k2; v3 = k3; \
  v3 ^= (nonce); \
  sip_round(); sip_round(); \
  v0 ^= (nonce); \
  v2 ^= 0xff; \
  sip_round(); sip_round(); sip_round(); sip_round(); \
  (v0 ^ v1 ^ v2  ^ v3); \
})

uint32_t c_solve(uint32_t *prof, const uint8_t *mesg) {
  int graph[M];
  int V[EN], U[EN];
  int path[CLEN];

  uint64_t k0, k1, k2, k3;
  uint64_t v0, v1, v2, v3;
  
  setkeys();
 
  for(int i=0; i<M; ++i) {
      graph[i] = -1;
  }
  
  #pragma ivdep
  for(uint64_t i=0; i<EN; ++i) {
      U[i] = ( siphash24((i << 1)) & MASK) << 1;
      V[i] = (((siphash24(((i<<1)+1))) & MASK) << 1) + 1;
  }
  
  for(uint64_t i=0; i<EN; ++i) {
      int u = U[i];
      int v = V[i];

      int pre = -1;
      int cur = u;
      int next;
      while(cur != -1) {
          next = graph[cur];
          graph[cur] = pre;
          pre = cur;
          cur = next;
      }

      int m = 0;
      cur = v;
      while(graph[cur] != -1 && m < CLEN) {
          cur = graph[cur];
          ++m;
      }

      if(cur != u) {
        graph[u] = v;
      } else if(m == CLEN-1) {
            int j;
            
            cur = v;
            for(j=0; j<=m; ++j) {
                path[j] = cur;
                cur = graph[cur];
            }

            for(j=0; j<M; ++j) {
                graph[j] = -1;
            }
            
            for(j=1; j<=m; ++j) {
                graph[path[j]] = path[j-1];
            }

            int k = 0;
            int b = CLEN -1;
            for(j=0; k < b; ++j) {
                int u = U[j];
                int v = V[j];

                if(graph[u] == v || graph[v] == u) {
                    prof[k] = j;
                    ++k;
                }
            }
            prof[k] = i;
            return 1;
      }

  }

  return 0;
}