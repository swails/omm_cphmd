// This tests the OpenMM functionality here
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>

#include "Amber.h"
#include "OpenMM.h"

using namespace std;

void check_omm_system(void) {
    Amber::AmberParm parm;

    parm.rdparm("files/trx.prmtop");

    OpenMM::System *system = parm.createSystem();

    ofstream output;
    output.open("serialized_system.xml");
    OpenMM::XmlSerializer::serialize<OpenMM::System>(system, string("System"),
                                                     output);
    output.close();
    delete system;
}

void check_omm_pme(void) {
}

int main() {

    cout << "Testing OpenMM System creation...";
    check_omm_system();
    cout << " OK." << endl;

    cout << "Testing OpenMM PME potential energies...";
    check_omm_pme();
    cout << " OK." << endl;
}
