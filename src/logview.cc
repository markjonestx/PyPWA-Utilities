#include <fitlog.h>

char *progname;
int lineno = 1;

int main(int argc, char** argv) {

	int debug = 0;
	int printWaves = 0;
	int phase_mode = 0;
	int split = 0;
	int nevents_mode = 1;
	int first_wave_set = 1;
	int second_wave_set = 0;
	ifstream logFile;
	char* progfile = NULL;
	char* logfile = NULL;
	fitlog log,log2;
	hash_map<const char*, parameter> params;
	hash_map<const char*, integral> ints;
	list<fitlog> bins;

	progname = argv[0];
	extern char* optarg;
	extern int optind;
	int c;
	while ( (c = getopt(argc,argv, "sn:l:wpdh")) != -1 )
		switch(c) {
			case 'd':
				debug = atoi(optarg);
				break;
			case 'n':
				progfile = optarg;
				break;
			case 'l':
				logfile = optarg;
				break;
			case 's':
				split = 1;
				break;
			case 'p':
				phase_mode = 1;
				nevents_mode = 0;
				break;
			case 'w':
				printWaves = 1;
				break;
			case 'h':
			case '?':
				cerr << "useage: " << argv[0] << " [-d] -n nevents.camp -l fit.pars [-w] [-p] [-s]" << endl;
				exit(0);
		}

	logFile.open(logfile);
	while( !(logFile >> log).eof() ) {
	    bins.push_back(log);
	}
	// cout << bins.size() << " bins read" << endl;
	// list<fitlog>::iterator bin_i = bins.begin();
	// while(bin_i != bins.end()) {
	//     bin_i->print();
	//     bin_i++;
	// }

	if(split) {
		char bin_name[1024];
		list<fitlog>::iterator bin_i = bins.begin();
		while(bin_i != bins.end()) {
			// make filename
			sprintf(bin_name,"%f.log",bin_i->binCenter());
			ofstream ofs(bin_name);
			ofs << *bin_i;
			bin_i++;
		}
		exit (0);
	}
	if(printWaves) {
		list<string> waves;
		waves = log.waveNames();
		list<string>::iterator wv = waves.begin();
		while(wv != waves.end()) {
			cout << *wv << " ";
			wv++;
		}
		cout << endl;
		exit(0);
	}

	// list<fitlog>::iterator bin_i = bins.begin();
	// while(bin_i != bins.end()) {
	//     bin_i->print();
	//     bin_i++;
	// }
	list<string> wvs;
	list<string> wvs1;
	list<string> wvs2;
	matrix<double> deriv;

	first_wave_set = 1;
	second_wave_set = 0;
	for(int iarg = optind; iarg < argc; iarg++) {
	    if(!strcmp(argv[iarg],"minus")) {
		first_wave_set = 0;
	    }
	    else if(first_wave_set) {
		wvs.push_back(argv[iarg]);
		wvs1.push_back(argv[iarg]);
	    }
	    else {
		second_wave_set = 1;
		wvs.push_back(argv[iarg]);
		wvs2.push_back(argv[iarg]);
	    }
	}

	if(nevents_mode) {
		// all fitlog's share the parser: it is like a static member
		log.progFile(progfile);
		log.progParse();

		list<fitlog>::iterator bin_i = bins.begin();
		while(bin_i != bins.end()) {

		    bin_i->setIntegrals();
		    bin_i->setParams(wvs);
		    bin_i->setNevents();
		    //		    deriv = bin_i->ddp( (double (fitlog::*)())&fitlog::nevents );
		    deriv = bin_i->ddp(&call_fitlog_nevents);
		    cout << bin_i->binCenter() << " " << bin_i->binWidth()/2 << " ";
		    cout << bin_i->nevents() << " " <<
		      bin_i->delta(&call_fitlog_nevents) << endl;
		    // bin_i->print();
		    bin_i++;
		}
	} else if(phase_mode) {
		log.progFile(progfile);
		log.progParse();

		list<fitlog>::iterator bin_i = bins.begin();
		while(bin_i != bins.end()) {

		    bin_i->setParams(wvs);
		    bin_i->setParams(wvs1,1);
		    bin_i->setParams(wvs2,2);
		    // bin_i->setNevents();
		    // deriv = bin_i->ddp( (double (fitlog::*)())&(fitlog::phase) );
		    cout << bin_i->binCenter() << " " << bin_i->binWidth()/2 << " ";
		    cout << bin_i->phase() << " " <<
		      bin_i->delta(&call_fitlog_phase) << endl;
		    // bin_i->print();
		    bin_i++;
		}
	}

	return 0;
}

// dummy function needed since we don't link to but not minuit.o
int mnparm(int, string, double, double, double, double) {
    cerr << "this is impossible" << endl;
    throw "aFit";
    return 0;
}
