#include <integral.h>

    

integral::integral() {
    this->_nwaves = 0;
    this->_nevents = 0;
    this->_maxEvents = 0;
}


integral::integral(char** files) {
    this->_nwaves = 0;
    this->_nevents = 0;
    this->_maxEvents = 0;
    this->files(files);
    this->_sum = matrix<complex<double> >(this->_nwaves, this->_nwaves);    
}


integral::integral(const integral& ni) {
    this->_nwaves = ni._nwaves;
    this->_nevents = ni._nevents;
    this->_maxEvents = ni._maxEvents;
    this->_index = ni._index;
    this->_sum = ni._sum;
}


integral::~integral() {
    ;
}


integral& integral::operator=(const integral& ni) {
    this->_nwaves = ni._nwaves;
    this->_nevents = ni._nevents;
    this->_maxEvents = ni._maxEvents;
    this->_index = ni._index;
    this->_sum = ni._sum;
    return *this;
}


    

integral& integral::integrate() {
    if (this->_nwaves == 0) throw "no waves";

    int nRead = 0;
    ifstream *ampfile = new ifstream[this->_nwaves];
    complex<double> *amps = new complex<double>[this->_nwaves];

    

    string filename;
    int fileindex;

    map<string, int>::iterator ampName = this->_index.begin();
    while( ampName != this->_index.end() ) {
        filename = (*ampName).first;
        fileindex = (*ampName).second;
        ampfile[fileindex].open((filename).c_str());
        if(!ampfile[fileindex]) {
            cerr << "error: cannot open " << filename << endl;
            throw "file not found";
        }
        ampName++;
    }



    int eof = 0;
    while (!eof && _maxEvents?nRead<_maxEvents:1) {
        

    int index;
    map<string, int>::iterator ampName = this->_index.begin();
    while( ampName != this->_index.end() ) {
        index = (*ampName).second;
        ampfile[index].read((char*) &amps[index], sizeof(complex<double>));
        if ( (eof = ampfile[index].eof()) ) break;
        ampName++;
    }
    if (eof) break;
    this->_nevents++;
	nRead++;

    if ( !(nRead % 100) ) cerr << nRead << "\r" << flush;


        

    int i, j;
    map<string, int>::iterator ampName_i = this->_index.begin();
    while( ampName_i != this->_index.end() ) {
        i = (*ampName_i).second;
        map<string, int>::iterator ampName_j = this->_index.begin();
        while( ampName_j != this->_index.end() ) {
            j = (*ampName_j).second;
            this->_sum.el(i,j) += amps[i]*conj(amps[j]);
            ampName_j++;
        }
        ampName_i++;
    }


    }

    delete [] ampfile;
    delete [] amps;
    return *this;
}


integral& integral::files(list<string> files) {
    list<string>::iterator file = files.begin();
    while(file != files.end()) {
        this->_index[*file] = this->_nwaves++;
        file++;
    }
    this->_sum = matrix<complex<double> >(this->_nwaves, this->_nwaves);    
    return *this;
}

integral& integral::files(char** files) {
    list<string> slist;
    while(*files) {
        slist.push_back(*files);
        files++;
    }
    this->files(slist);
    return *this;
}


    integral& integral::renormalize(int n) {
        _sum = ((complex<double>) ((double) n/ (double) _nevents))*_sum;
        _nevents = n;
        return *this;
    }


    integral& integral::max(int m) { 
        this->_maxEvents = m; 
        return *this; 
    }

    integral& integral::events(int n) { 
        this->_nevents = n; 
        return *this; 
    }


complex<double>& integral::el(string iname, string jname) {
    return (_sum.el(_index[iname],_index[jname]));
}


    

int integral::nevents() const {
	return _nevents;
}


list<string> integral::files() const {
    list<string> fileList;
    map<string, int>::const_iterator ampName = this->_index.begin();
    while( ampName != this->_index.end() ) {
        fileList.push_back((*ampName).first);
        ampName++;
    }
    return fileList;
}

char** integral::files_c_str() const {
    string fname;
    int i = 0;
    char** fileList = (char**) malloc((_nwaves+1)*sizeof(char*));

    map<string, int>::const_iterator ampName = this->_index.begin();
    while( ampName != this->_index.end() ) {
        fname = (*ampName).first;
        fileList[i] = (char*) malloc ((fname.size()+1)*sizeof(char));
        strcpy(fileList[i],fname.c_str());
        ampName++; i++;
    }
    fileList[i] = NULL;
    return fileList;
}


complex<double> integral::val(string iname, string jname) {

    if( _index.find(iname) == _index.end() ) {
        cerr << "error: " << iname << " not in integral" << endl;
        throw "bad wave access";
    }
    if( _index.find(jname) == _index.end() ) {
        cerr << "error: " << jname << " not in integral" << endl;
        throw "bad wave access";
    }

    return (this->el(iname,jname)/((double) _nevents));
}


integral integral::get(list<string> flist) {
    // need to check that all requested files are in list
    list<string>::iterator iname;
    list<string>::iterator jname;
    integral ret;

    iname = flist.begin();
    while (iname != flist.end()) {
        if( _index.find(*iname) == _index.end() ) {
            cerr << "error: " << *iname << " not in integral" << endl;
            throw "bad wave access";
        }
        iname++;
    }

    ret.files(flist);
    ret.events(_nevents);
    iname = flist.begin();
    while (iname != flist.end()) {
        jname = flist.begin();
        while (jname != flist.end()) {
            ret.el(*iname,*jname) = this->el(*iname,*jname);
            jname++;
        }
        iname++;
    }
    return ret;
}

integral integral::get(char** flist) {
    list<string> slist;
    while(*flist) {
        slist.push_back(*flist);
        flist++;
    }
    return this->get(slist);
}


    int integral::index(string s) {
        return _index[s];
    }

    int integral::index(char* s) {
        return _index[s];
    }


    matrix<complex<double> > integral::mat() {
        return ((complex<double>) (1.0/((double) _nevents)))*_sum;
    }


    

const integral& integral::print(ostream& os) const {
    
    os << _nwaves << endl;
    os << _nevents << endl;
    os << _sum;
    os << _index.size() << endl;
    map<string, int>::const_iterator ampName = this->_index.begin();
    while( ampName != _index.end() ) {
        os << (*ampName).first << " " << (*ampName).second << endl;
        ampName++;
    }

    return *this;
}


const integral& integral::print_events(ostream& os) const {

    os << _nevents << endl;
    return *this;
}


integral& integral::scan(istream& is) {
    int indexSize = 0, index = 0;
    string name;

    is >> _nwaves;
    is >> _nevents;
    is >> _sum;
    is >> indexSize;
    while(indexSize--) {
        is >> name >> index;
        _index[name] = index;
    }
    return *this;
}


