
#line 8 "keyfile.nw"


#include <cstdio>
#include <list>
#include <string>
#include <gamp/wave.h>


class keyfile {
private:
    list<wave> _atomicWaves;
    string _filename;
    FILE *_file;
public:
    keyfile();

    keyfile(string filename);

    ~keyfile();

    keyfile(const keyfile &);

    keyfile operator=(const keyfile &);

    keyfile open(const string filename);

    keyfile run(event &e);

    keyfile run();

    void rewind();

    void close();

    FILE *file();

    wave addWave(const wave w);

    wave operator[](int index) const;

    int nWaves() const;
};

keyfile readKey(keyfile);

