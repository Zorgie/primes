#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <vector>
#include <iostream>

using namespace std;

void trialDivision(mpz_t number,vector<mpz_t*>& v);

int main (int argc, char* argv[])
{

	

	if(argc < 2)
	{
		printf("Please supply a number to factor.\n");
		return 0;
	}
	
	mpz_t N;
	mpz_init(N);
	mpz_init_set_str(N, argv[1], 10);

	
	vector<mpz_t*> vec;
	trialDivision(N,vec);

	int end = vec.size();

	for(int i=0; i < end; i++)
	{
		cout << "Fac: "<< mpz_get_str(NULL, 10, *(vec.at(i))) << endl;
	}
		
	cout << "WTF!" << endl;
	
	return 0;


}
// seems relevant to the report to have a naive implementation just to see how it would peform compared to the other algorithms
void trialDivision(mpz_t number,vector<mpz_t*>& v)
{
	mpz_t tmp;
	mpz_init(tmp);
	mpz_set_ui(tmp,(unsigned long int)2);

	if(mpz_divisible_p(number,tmp))
	{
		v.push_back();
		mpz_t isLeft;
		mpz_init(isLeft);
		mpz_divexact(isLeft, number, ny);
		cout << mpz_get_str(NULL,10,isLeft) << endl;
		trialDivision(isLeft, v);
	}

	return;

	mpz_t rot;
	mpz_init(rot);
	mpz_sqrt (rot, number);
	
	mpz_set_ui(tmp,(unsigned long int)3);


	/*while(mpz_cmp(rot, tmp)>0)
	{
			cout << "hej" << endl;	
		if(mpz_divisible_p(number,tmp))
		{
			mpz_t ny;
			mpz_init(ny);
			mpz_set(ny,tmp);
			v.push_back(&ny);
			mpz_t isLeft;
			mpz_init(isLeft);
			mpz_divexact(isLeft, number, ny);
			trialDivision(isLeft, v);
		}
		mpz_add_ui(tmp,tmp,(unsigned long int)2);
		
	}
	return;
*/
	
}



