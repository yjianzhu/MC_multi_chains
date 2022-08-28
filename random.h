// random.h
//============================================================
// The first set of subroutines are used to generate random
// numbers with a uniform distribution in the range [0,1),
// using the MERSENNE TWISTER RANDOM NUMBER GENERATOR method.
// The final subroutine converts this to a normal distribution.
//------------------------------------------------------------
#define BB 624
#define CC 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UMASK 0x80000000UL /* most significant w-r bits */
#define LMASK 0x7fffffffUL /* least significant r bits */
#define MIXBITS(u,v) ( ((u) & UMASK) | ((v) & LMASK) )
#define TWIST(u,v) ((MIXBITS(u,v) >> 1) ^ ((v)&1UL ? MATRIX_A : 0UL))
static unsigned long state[BB]; /* the array for the state vector  */
static int left = 1;
static int initf = 0;
static unsigned long *next;
void init_genrand(unsigned long s);
static void next_state(void);
unsigned long genrand_int32(void);
double genrand_real2(void);

//============================================================
// initializes state[BB] with a seed
//------------------------------------------------------------
void init_genrand(unsigned long s)
{
	int j;
	state[0] = s & 0xffffffffUL;
	for (j = 1; j<BB; j++) {
		state[j] = (1812433253UL * (state[j - 1] ^ (state[j - 1] >> 30)) + j);
		/* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
		/* In the previous versions, MSBs of the seed affect   */
		/* only MSBs of the array state[].                     */
		/* 2002/01/09 modified by Makoto Matsumoto             */
		state[j] &= 0xffffffffUL;  /* for >32 bit machines */
	}
	left = 1; initf = 1;
}
//============================================================
static void next_state(void)
{
	unsigned long *p = state;
	int j;

	/* if init_genrand() has not been called, */
	/* a default initial seed is used         */
	if (initf == 0) init_genrand(5489UL);

	left = BB;
	next = state;

	for (j = BB - CC + 1; --j; p++)
		*p = p[CC] ^ TWIST(p[0], p[1]);

	for (j = CC; --j; p++)
		*p = p[CC - BB] ^ TWIST(p[0], p[1]);

	*p = p[CC - BB] ^ TWIST(p[0], state[0]);
}
//============================================================
// generates a random number on [0,0xffffffff]-interval
//------------------------------------------------------------
unsigned long genrand_int32(void)
{
	unsigned long y;

	if (--left == 0) next_state();
	y = *next++;

	/* Tempering */
	y ^= (y >> 11);
	y ^= (y << 7) & 0x9d2c5680UL;
	y ^= (y << 15) & 0xefc60000UL;
	y ^= (y >> 18);

	return y;
}
//============================================================
// generates a random number on [0,1)-real-interval
//------------------------------------------------------------
double genrand_real2(void)
{
	unsigned long y;

	if (--left == 0) next_state();
	y = *next++;

	/* Tempering */
	y ^= (y >> 11);
	y ^= (y << 7) & 0x9d2c5680UL;
	y ^= (y << 15) & 0xefc60000UL;
	y ^= (y >> 18);

	return (double)y * (1.0 / 4294967296.0);
	/* divided by 2^32 */
}
//============================================================
// generates a random number on [0,1) with 53-bit resolution
//------------------------------------------------------------
double genrand_res53(void)
{
        unsigned long a = genrand_int32() >> 5, b = genrand_int32() >> 6;
        return(a*67108864.0 + b)*(1.0 / 9007199254740992.0);
}
