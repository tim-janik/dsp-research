#include <cstring>

typedef unsigned int uint;

////////////// start: code based on log2 code from Anklang/ASE by Tim Janik

/** Union to compartmentalize an IEEE-754 float.
 * IEEE 754 single precision floating point layout:
 * ```
 *        31 30           23 22            0
 * +--------+---------------+---------------+
 * | s 1bit | e[30:23] 8bit | f[22:0] 23bit |
 * +--------+---------------+---------------+
 * B0------------------->B1------->B2-->B3-->
 * ```
 */
union FloatIEEE754 {
  float         v_float;
  struct {
#if   __BYTE_ORDER == __LITTLE_ENDIAN
    uint mantissa : 23, biased_exponent : 8, sign : 1;
#elif __BYTE_ORDER == __BIG_ENDIAN
    uint sign : 1, biased_exponent : 8, mantissa : 23;
#endif
  } mpn;
  static constexpr const int   BIAS = 127;                       ///< Exponent bias.
};

/** Fast approximation of logarithm to base 2.
 * The parameter `x` is the exponent within `[1.1e-38…2^127]`.
 * Within `1e-7…+1`, the error stays below 3.8e-6 which corresponds to a sample
 * precision of 18 bit. When `x` is an exact power of 2, the error approaches
 * zero. With FMA instructions and `-ffast-math enabled`, execution times should
 * be below 10ns on 3GHz machines.
 */

static inline float
fast_log2 (float value)
{
  const int EXPONENT_MASK = 0x7F800000;
  int iv;
  memcpy (&iv, &value, sizeof (float));                 // iv = *(int *) &values[k]
  int fexp = (iv >> 23) - FloatIEEE754::BIAS;            // extract exponent without bias (rely on sign bit == 0)
  iv = (iv & ~EXPONENT_MASK) | FloatIEEE754::BIAS << 23; // reset exponent to 2^0 so v_float is mantissa in [1..2]
  float r, x;
  memcpy (&x, &iv, sizeof (float));                      // x = *(float *) &iv
  x -= 1;
  // x=[0..1]; r = log2 (x + 1);
  // h=0.0113916; // offset to reduce error at origin
  // f=(1/log(2)) * log(x+1); dom=[0-h;1+h]; p=remez(f, 6, dom, 1);
  // p = p - p(0); // discard non-0 offset
  // err=p-f; plot(err,[0;1]); plot(f,p,dom); // result in sollya
  r = x *  -0.0259366993544709205147977455165000143561553284592936f;
  r = x * (+0.122047857676447181074792747820717519424533931189428f + r);
  r = x * (-0.27814297685064327713977752916286528359628147166014f + r);
  r = x * (+0.45764712300320092992105460899527194244236573556309f + r);
  r = x * (-0.71816105664624015087225994551041120290062342459945f + r);
  r = x * (+1.44254540258782520489769598315182363877204824648687f + r);
  return fexp + r; // log2 (i) + log2 (x)
}


