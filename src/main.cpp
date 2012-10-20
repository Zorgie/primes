#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

using namespace std;

mpz_t* fermatFactors(mpz_t N);

int main(int argc, char* argv[])
{
	if(argc < 2)
	{
		printf("Please supply a number to factor.\n");
		return 1;
	}

	mpz_t N;
	mpz_init(N);
	mpz_init_set_str(N, argv[1], 10);
	
	mpz_t* factors = fermatFactors(N);
	int numFactors = mpz_get_ui(factors[0]);

	for(int i=0; i<numFactors; i++)
	{
		printf("Factor %d: %s\n", i+1, mpz_get_str(NULL, 10, factors[i+1]));
	}
	delete[] factors;
	return 0;
}


mpz_t* fermatFactors(mpz_t N)
{
//	printf("Factoring %s\n", mpz_get_str(NULL, 10, N));
	// N is already a prime number.
	if(mpz_probab_prime_p(N, 10))
	{
		mpz_t* factors = (mpz_t*)malloc(2*sizeof(mpz_t));
		mpz_init(factors[0]);
		mpz_init(factors[1]);
		mpz_set_ui(factors[0], 1);
		mpz_set(factors[1],N);
		return factors;
	}

//	printf("%s is not a prime number.\n", mpz_get_str(NULL, 10, N));
	mpz_t a, b, b2, one;

	mpz_init(a);
	mpz_init(b);
	mpz_init(b2);
	mpz_init_set_ui(one, 1);
	
	
	mpz_sqrt(a, N);
	mpz_add(a, a, one);
	mpz_mul(b2, a, a);
	mpz_sub(b2, b2, N);

	while(!mpz_perfect_square_p(b2) && mpz_cmp(a, N) < 0)
	{
	//	printf("%s is not a square, trying the next number.\n", mpz_get_str(NULL, 10, b2));
		mpz_add(b2, b2, a);
		mpz_add(b2, b2, a);
		mpz_add(b2, b2, one);
		mpz_add(a, a, one);
	}	
	mpz_sqrt(b, b2);
	mpz_t fac1, fac2;
	mpz_init(fac1);
	mpz_init(fac2);

	mpz_add(fac1, a, b);
	mpz_sub(fac2, a, b);

	mpz_t* first;
	mpz_t* second;
	first = fermatFactors(fac1);
	second = fermatFactors(fac2);
	int quant1, quant2;
	
	quant1 = mpz_get_ui(first[0]);
	quant2 = mpz_get_ui(second[0]);

	mpz_t* factors = new mpz_t[quant1 + quant2 + 1];
	for(int i=0; i<quant1 + quant2 + 1; i++)
	{	
		mpz_init(factors[i]);
	}
	mpz_set_ui(factors[0], quant1+quant2);
	for(int i=0; i<quant1; i++)
	{
		mpz_set(factors[i+1], first[i+1]);
	}
	for(int i=0; i<quant2; i++)
	{
		mpz_set(factors[i+1+quant1], second[i+1]);
	}
	delete[] first;
	delete[] second;
	return factors;
}

