#include <vamp/Arguments.hpp>


using namespace std;


void outputBinary(int maxEvents) {
    int nRead = 0;
    complex<double> amplitude;

    while ((cin >> amplitude) && maxEvents != nRead) {
        cout.write((char *) &amplitude, sizeof(complex<double>));
        nRead++;
    }
}

void outputAscii(int maxEvents) {
    int nRead = 0;
    complex<double> amplitude;

    while ((cin.read((char *) &amplitude, sizeof(complex<double>))) && maxEvents != nRead) {
        cout << amplitude << endl;
        nRead++;
    }
}


int main(int argc, char **argv) {
    vamp::Arguments args(argc, argv);

    if (args.askedForHelp()) {
        args.printHelp();
        return 0;
    }

    if (args.getBinary()) {
        outputBinary(args.getMaxEvents());
    } else {
        outputAscii(args.getMaxEvents());
    }

    return 0;

}

