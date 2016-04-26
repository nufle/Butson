/* main.cpp

Main program used to search butson matrices

Arnaud Mallen mallen.arnaud@gmail.com
Vincent Popie vincent.popie@gmail.com

 */
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <ctime>

#include "ObjectPool.h"
#include "matrix.h"
#include "completer.h"

using std::cout;
using std::endl;
using std::cin;
using std::ostringstream;
using std::clock_t;


int main (int argc, char ** argv)
{
	// TO BE REPLACED by bash script
	constexpr int lconst = 2;
	constexpr int Nconst = NVALUE;

	const int sizeT = lconst * Nconst;
	Matrix<int,sizeT,sizeT> MN;
	MN.addZeroValuedRow();
	MN.addZeroValuedRow();
	for(int valeurs=0;valeurs<(lconst);valeurs++){
		for( int position=0;position<(Nconst);position++)
			MN(1,valeurs*Nconst+position)=static_cast<int>(valeurs);
	}	
	Completer<int,Nconst,lconst> cp;

	vector<int> L(lconst*lconst,0);
	
	clock_t begin = clock();
	std::cout<<"Start completing..."<<std::endl;

	cp.complete(&MN,2,0,&L);

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

	std::cout<<"Elapsed time : "<<elapsed_secs<<" s"<<std::endl;

	cp.checkSolutions();

	return 0;
}







