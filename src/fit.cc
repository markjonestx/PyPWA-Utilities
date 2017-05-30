#include <iostream>
#include <fstream>
#include <string>
#include <complex>
#include <stdio.h>
#include <signal.h>
#include <setjmp.h>
#include <unistd.h>
#include <fit.h>
#include <parseutils.h>
#include <minuit.h>
#include <fitlog.h>

extern FILE* yyin;

void PrintUsage(char *processName);

char *progname;
char *ofile_name = "fit.log";
int lineno = 1;

float low_mass = 0.0;
float high_mass = 0.0;
int quiet = 0;


// Added by P Eugenio
int open_mnresult(string fname){
  //int ierr,unit=1;
  const char* name = fname.c_str();
  fortopen_(name,strlen(name));
}
int open_mninput(string fname){
  
  //int ierr,unit=1;
  const char* name = fname.c_str();
  fortopeninput_(name,strlen(name));

}
//---------------



void fcn(int *npar, double *grad, double *fval, double *xval, int *iflag)
{

    static int fcnCallno = 0;
    static double minus_log_like = 0;

    switch (*iflag) {
    case 1:
	// parse file
	if(!quiet) cout << "iflag = " << *iflag << endl;
	if(!quiet) cout << "parsing file." << endl;
	init_event_code();
	yyparse();
	read_amps();
	read_rmat();
	set_params();

	cout << "calling  mnintr(fcn)" << endl;
	mnintr(fcn);


	break;
    case 2:
	if(!quiet) cout << "iflag = " << *iflag << endl;
	if(!quiet) cout << "no gradient calculation yet." << endl;
	break;
    case 3:
	{

	    if(!quiet) cout << "iflag = " << *iflag << endl;
	    if(!quiet) cout << "outputing values" << endl;
	    ofstream outfile(ofile_name);
	    
	    // Rounding off # of events
	    int nevents = static_cast<int>(floor(real(lookup("nevents")->val) + 0.5));

	    hash_map <const char *, parameter> p = get_params();
	    hash_map <const char *, integral, hash<const char*>, eqstr> ints = get_integrals();
	    matrix<double> emat = minEmat(nMinpars());
	    fitlog fl(nevents, low_mass, high_mass, minus_log_like, p, ints, emat);
	    outfile << fl;
	}
	break;

    default:
	break;
    }

    lookup("fcn")->val = 0.0;

    // cerr << "fcn(initial): " << lookup("fcn")->val << endl;

    init_event_code();
    rewind_amps();
    rewind_rmat();
    update_params(xval);
    int event;
    int iset = 1;
    event = 1;
    while (update_symbols()) {

	execute(event_prog);
	double d = lookup("fcn")->val.real();
	if (!finite(d) && iset) {
	    cerr << "fcn not finite at event " << event << endl;
	    cerr << "symbol table on event " << event << ":" << endl;
	    print_stable();
	    iset = 0;
	}
	event++;
    }
    complex < double >fcn_unnormalized;
    complex < double >fcn_normalized;
    complex < double >nevents;

    if(!quiet) cout << "fcn call: " << fcnCallno++ << endl;
    fcn_unnormalized = lookup("fcn")->val;
    if(!quiet) cout << "fcn(unnormalized): " << fcn_unnormalized << endl;

    init_norm_code();
    execute(norm_prog);

    fcn_normalized = real(lookup("fcn")->val);
    nevents = real(lookup("nevents")->val);
    if(!quiet) cout << "fcn(normalized): " << fcn_normalized << endl;
    if(!quiet) cout << "nevents: " << nevents << endl;
    if(!quiet) cout << "normalization: " << (fcn_normalized -
				  fcn_unnormalized) / real(nevents) << endl;

    minus_log_like = real(fcn_normalized);
    *fval = minus_log_like;
    if(!quiet) cout << endl;

  // Added by P Eugenio  -----------------------------------
    if(*iflag ==3){
      // **
      // ** Now write the loglikehood
      // **    fcn_normalized = real(lookup("fcn")->val);
      if(!quiet) cerr << "my_loglikelihood:  " << real(fcn_unnormalized) <<endl;
      
      double args[10];
      int narg = 0;
      string command = "SAVE";  // save to "fit.result"
      args[0] = 0;
      mnexcm(fcn,command,args,narg);
      
    }
    //---------------------------------------------
}



int main(int argc, char **argv)
{
    int status, interactive=0;
    progname = argv[0];
    srand48((long)time(NULL));          //////PME
    extern char *optarg;
    // extern int optind;
    int c;
    //--- Added by Paul
    int fix_parm[8]={0,0,0,0,0,0,0,0};
    int nfixed=0;
    int UseRandomStarts=0;
    string mninput_fname("minuit.input");
    string mnresult_fname("fit.result");

    while ((c = getopt(argc, argv, "p:r:i:l:u:Io:qh")) != -1)
	switch (c) {
	  //----------- Added by P Eugenio 
	case 'p': 
	  if(nfixed>7){
	    cerr <<"You can only fix 8 parms\n";exit(0);
	  }
	  fix_parm[nfixed] = atoi(optarg);
	  cerr<<"Fixing minuet parm: "<<fix_parm[nfixed++];
	  break;
	case 'R':
	  UseRandomStarts = 1;
	  break;
	case 'r':
	  mnresult_fname = optarg;
	case 'i':
	  mninput_fname = optarg;
	  cerr<<"Using Minuit input file: "<<mninput_fname<<endl;
	  break;
	  //--------------------------
	case 'I':
	    interactive = 1;
	    break;
	case 'l':
	    low_mass = atof(optarg);
	    break;
	case 'u':
	    high_mass = atof(optarg);
	    break;
	case 'o':
	    ofile_name = optarg;
	    break;
	case 'q':
	    quiet = 1;
	    break;
	case 'h':
	    cerr << "usage:" << endl;
	    cerr << progname << " [-i -l # -u # -h -o file] < inputfile" << endl;
	    cerr << "    where:" << endl;
	    //-------------Added by P Eugenio
	    cerr << "        -p #  fix a minuet parameter" << endl;
	    cerr << "        -i <minuit_input_file> \n";
	    cerr << "        -r <fit_result_file> \n";
	    cerr << "        -R        use random starting parameter values \n";
	    ///-----------------------
	    cerr << "        -I:       interactive minuit session" << endl;
	    cerr << "        -l #:     lower edge of bin" << endl;
	    cerr << "        -u #:     upper edge of bin" << endl;
	    cerr << "        -o file:  set output file" << endl;
	    cerr << "        -q:       run quietly" << endl;
	    cerr << "        -h:       print help" << endl;
	    exit(0);
	    break;
	}

	yyin  = fopen(argv[optind],"r");

    init();

    string command, title;
    double args[10];
    int narg = 0;

    //  was 
    //    mninit(5, 6, 7);
    //    fortopen_();
    // 
    //------Added by Paul
    mninit(1,6,7);
    open_mnresult(mnresult_fname);
    open_mninput(mninput_fname);
    //-------------------


    title = "pwa fit";
    mnseti(title);

    command = "CALL FCN";
    args[0] = 1;
    narg = 1;
    mnexcm(fcn, command, args, narg);

    command = "SET ERR";
    args[0] = 0.5;
    narg = 1;
    mnexcm(fcn, command, args, narg);

    if(interactive) {
        mnintr(fcn);
	exit(0);
    }

  // Added by Paul--------------------------


    while(nfixed--){
      command = "SET PAR";
      args[0] = fix_parm[nfixed];
      args[1] = 0.0;
      narg = 2;
      mnexcm(fcn,command,args,narg);
      
      command = "FIX";
      args[0] = fix_parm[nfixed];
      narg = 1;
      mnexcm(fcn,command,args,narg);
    }
    //---------------------------
    if(UseRandomStarts)
      for(int r=0;r<50;r++){ //randon starts
	command = "SET PAR";
	args[0] = r+1;
	args[1] = drand48();
	narg = 2;
	mnexcm(fcn,command,args,narg);
      }
    
    
    
    command = "MINI";
    args[0] = 50000; //----- modified by paul
    narg = 0;
    status = mnexcm(fcn, command, args, narg);

    command = "CALL FCN";
    args[0] = 3;
    narg = 1;
    status = mnexcm(fcn, command, args, narg);

//      if(status == 0) {
//              command = "SAVE";
//              args[0] = 0;
//              narg = 0;
//              mnexcm(fcn,command,args,narg);
//      }

    return 0;
}

