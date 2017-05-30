#include <keyfile.h>

extern particleDataTable table;
extern FILE* keyin;
extern char* fname;
extern "C++" int keylex();
extern "C++" int keyparse();


keyfile::keyfile() {
	_filename = "";
}

keyfile::keyfile(string filename) {
	_filename = filename;
	fname = (char*) malloc( ((filename.length())+1)*sizeof(char));
	strcpy(fname,filename.c_str());
}

keyfile::~keyfile() {
	;
}

keyfile::keyfile(const keyfile& kf) {
	this->_atomicWaves = kf._atomicWaves;
	this->_filename = kf._filename;
	this->_file = kf._file;
}

keyfile keyfile::operator=(const keyfile& kf) {
	this->_atomicWaves = kf._atomicWaves;
	this->_filename = kf._filename;
	this->_file = kf._file;
	return *this;
}

keyfile keyfile::open(const string filename) {
	fname = (char*) malloc( ((filename.length())+1)*sizeof(char));
	strcpy(fname,filename.c_str());
	this->_file = fopen(fname,"r");
	return *this;
}


keyfile keyfile::run() {
	keyin = this->_file;
	keyparse();
	return *this;
}

event e;
keyfile keyfile::run(event& ev) {
	e = ev; 
	keyin = this->_file;
	keyparse();
	return *this;
}

void keyfile::rewind() {
	fseek(this->_file,0L,SEEK_SET);
}

void keyfile::close() {
	;
}

FILE* keyfile::file() {
	return this->_file;
}
		
wave keyfile::addWave(const wave w) {
	this->_atomicWaves.push_back(w);
	return w;
}

wave keyfile::operator[](int index) const {
	int i = 0;
	list<wave>::const_iterator w;
	w = this->_atomicWaves.begin();
	while(w != this->_atomicWaves.end() ) {
		if ( i == index ) {
			return *w;
		}
		w++; i++;
	}
	return *w;
}

int keyfile::nWaves() const {
	return this->_atomicWaves.size();
}

