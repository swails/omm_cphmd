// This tests the OpenMM functionality here
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>

#include "amberparm.h"
#include "OpenMM.h"

using namespace std;

void check_omm_system(void) {
    Amber::AmberParm parm;

    parm.rdparm("trx.prmtop");

    OpenMM::System *system = parm.createSystem();

    ofstream output;
    output.open("serialized_system.xml");
    OpenMM::XmlSerializer::serialize<OpenMM::System>(system, string("System"),
                                                     output);
    output.close();
    delete system;
}

int main() {

    cout << "Testing OpenMM System creation...";
    check_omm_system();
    cout << " OK." << endl;
}
