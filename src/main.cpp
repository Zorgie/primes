#include <stdio.h>
#include <stdlib.h>
#include <list>
#include <gmp.h>
#include <time.h>
#include <string>

using namespace std;

static const int TIMEOUT = 1000000000;
static const int INPUT = 10;

int solve(char*);
int main2();
int fermatFactors(const mpz_t&, mpz_t*);
mpz_t* getFactors(mpz_t &N, mpz_t*);
mpz_t* fermatFactors(mpz_t N);

int main()
{
	bool first = true;
	for(int i=0; i<INPUT; i++)
	{
		char* s = new char[100];
		scanf("%s", s);
		if(first)
		{
			first = false;
		}
		else
			printf("\n");
		int count = 0;
		for(int j=0; j<100; j++)
		{
			count++;
			if(s[j] == 0)
				j = 100;
		}
		if(count <= 12)
		{
			solve(s);
		}
		else
		{
			printf("fail\n");
		}
	}
}

int main2()
{
	char** input = new char*[INPUT];
	char* s;
	for(int i=0; i<INPUT; i++)
	{
		input[i] = new char[INPUT];
		scanf("%s", input[i]);
		if(input[i][0] == ' ')
		{
			break;
		}
	}
	
	for(int i=0; i<INPUT; i++)
	{
		printf("\n");
		int count = 0;
		for(int j=0; j<100; j++)
		{
			count++;
			if(input[i][j] == 0)
			{
				j = 100;
			}
		}
		if(count <= 12)
		{
			solve(input[i]);
		}
		else
		{
			printf("fail\n");
		}
	}
	return 0;
}

int solve(char* argv)
{
	mpz_t N;
	mpz_init(N);
	mpz_init_set_str(N, argv, 10);
	
	mpz_t* factors = new mpz_t[200];
	mpz_t* start = factors;
	for(int i=0; i<200; i++)
	{
		mpz_init(factors[i]);
	}
	if(fermatFactors(N, factors) == -1)
	{
		printf("fail");
		return 1;
	}
	
	while(mpz_cmp_ui(*start, 0) != 0)
	{
		printf("%s\n", mpz_get_str(NULL, 10, *start));
		start++;
	}

	
	//mpz_t* factors = fermatFactors(N);
//	int numFactors = mpz_get_ui(factors[0]);

//	for(int i=0; i<numFactors; i++)
//	{
//		printf("Factor %d: %s\n", i+1, mpz_get_str(NULL, 10, factors[i+1]));
//	}
//	delete[] factors;
	return 0;
}

int fermatFactors(const mpz_t& _N, mpz_t* factors)
{
	mpz_t N;
	mpz_init(N);
	mpz_set(N, _N);
	while(mpz_cmp_ui(N, 1) != 0 && !mpz_probab_prime_p(N, 10))
	{
		//printf("%s is not a prime, factorizing it.\n", mpz_get_str(NULL, 10, N));
		mpz_t oldN;
		mpz_init(oldN);
		mpz_set(oldN, N);
		factors = getFactors(N, factors);
		if(mpz_cmp(oldN, N) == 0)
			return -1;
	}
	return 0;
}

mpz_t* getFactors(mpz_t &N, mpz_t* factors)
{
//	printf("Getting factors for %s\n", mpz_get_str(NULL, 10, N));
	if(mpz_cmp_ui(N, 1) == 0)
		return factors;
	mpz_t one, cmp;
	mpz_init(one);
	mpz_init(cmp);
	mpz_set_ui(one, 1);
	mpz_and(cmp, N, one);
	while(mpz_cmp(cmp, one) != 0)
	{
		// Special case, N is even.
		mpz_t two;
		mpz_init(two);
		mpz_set_ui(two, 2);
		mpz_set(*factors, two);
		factors++;
		mpz_div(N, N, two);
		mpz_and(cmp, N, one);
	}
	mpz_t a, b, b2;
	mpz_init(a);
	mpz_init(b);
	mpz_init(b2);
	mpz_sqrt(a, N);
	while(mpz_cmp(a, N) < 0)
	{
		mpz_add_ui(a, a, 1);
		mpz_mul(b2, a, a);
		mpz_sub(b2, b2, N);
		if(mpz_perfect_square_p(b2))
		{
			mpz_sqrt(b, b2);
			mpz_t fact1, fact2;
			mpz_init(fact1);
			mpz_init(fact2);
			mpz_add(fact1, a, b);
			mpz_sub(fact2, a, b);
			if(mpz_cmp_ui(fact1, 1) != 0)
			{
				mpz_set(*factors, fact1);
				factors++;
			}
			if(mpz_cmp_ui(fact2, 1) != 0)
			{
				mpz_set(*factors, fact2);
				factors++;
			}
			mpz_div(N, N, fact1);
			mpz_div(N, N, fact2);
			//printf("Factors %s and %s found, N is now %s.\n", mpz_get_str(NULL, 10, fact1), mpz_get_str(NULL, 10, fact2), mpz_get_str(NULL, 10, N));
			return factors;
		}
	}
	return factors;	
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

