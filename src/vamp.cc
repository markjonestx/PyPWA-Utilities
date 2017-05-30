
#include <complex>
#include <iostream>
#include <cstdlib>
#include <unistd.h>

using namespace std;


void printUsage(char* prog) {
	cerr << "usage:" << endl;
	cerr << "  " << prog << " [-h] [-m maxEvents] [-B] < inputfile" << endl;
	cerr << "-B\tconvert ascii to binary mode" << endl;

}


int main(int argc, char** argv) {

	int maxEvents = 0;
	int nRead = 0;
	int binaryMode = 0;

	

    extern char* optarg;
    int c;
    while ( (c = getopt(argc,argv, "m:h:B")) != -1 )
        switch(c) {
            case 'm':
                maxEvents = atoi(optarg);
                break;
		case 'B':
			binaryMode = 1;
			break;
            case 'h':
            case '?':
                printUsage(argv[0]);
				exit(0);
	       
        }


	complex<double> a;

	if (binaryMode) {
	while ((cin >> a) && (maxEvents?nRead<maxEvents:1) ) {
		nRead++;
		cout.write((char*) &a,sizeof(a));
	}
	}
	else {

	while( (cin.read((char*) &a,sizeof(a))) && (maxEvents?nRead<maxEvents:1) ) {
		nRead++;
		cout << a << endl;
	}
}

	return 0;

}

