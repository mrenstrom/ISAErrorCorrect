//
//  readFasta.cpp
//  mapris
//
//  Created by mark enstrom on 3/26/18.
//  Copyright © 2018 mark enstrom. All rights reserved.
//

#include "readFasta.hpp"
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>
#include <iomanip>
#include <numeric>
#include <ctime>
#include "NWAlign.hpp"
using namespace std;
//
// vector sort classes
//
struct myclassS {
    bool operator() (Capture *pa1,Capture *pa2) { return (pa1->_count > pa2->_count);}
} smallCompObj;


struct myclassL {
    bool operator() (Capture *pa1,Capture *pa2) { return (pa1->_count < pa2->_count);}
} largeCompObj;


struct myclassLen {
    bool operator() (Capture *pa1,Capture *pa2) { return (pa1->_seq.length() > pa2->_seq.length());}
} smallLenCompObj;


struct myclassString {
    bool operator() (Capture *pa1,Capture *pa2) { return (pa1->_seq > pa2->_seq);}
} smallStrCompObj;

//  readFile.cpp
//  Created by mark enstrom on 12/3/17.
//  Copyright © 2017 mark enstrom. All rights reserved.
//
/*--------------------------------------------------------------------------------------------
 *
 * like strncmp - exact string comp
 * compare a dna(char) sequence and return # differences
 *
 *
 *--------------------------------------------------------------------------------------------*/
int dnancmp(char *p1, char *p2,size_t n) {
    int error = 0;
    for (int i = 0; i < n; i++) {
        if (p1[i] != p2[i]) {
            error += 1;
            if (error > 4) return error;
        }
    }
    return error;
}

/*--------------------------------------------------------------------------------------------
 *
 * Read for fasta
 *
 *
 *>BMTW6:02212:01662
 *AGTAGTATGTGGCCCGCTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTT
 *TAGTCAGTGCAAGGCAGAGATATACTGGCCTTCCTGCTAGGGGAGAGTGTGAGCCACACT
 *CACGCTCTCTCGTTCTCCACTGGAGGGGACACCCGCACACAAGGAGCTTGCAGCCCAGGA
 *GGACCAGGAAACCGG
 *>BMTW6:03334:00712
 *CATAGTGGTGA
 *
 *
 *--------------------------------------------------------------------------------------------*/
bool
ReadFasta::readFastA(string filename)
{
    std::ifstream inFile;
    inFile.open(filename);
    if (!inFile) return false;
    char hdr[1024];
    char seq[1024];
    //
    // read first ID line
    //
    int state = 0;   // 0 = ID, 1 = SEQ
    string id = "";
    do {
        if (state == 0) {
            inFile.getline(hdr, 1024);
            if (!inFile) {
                break;
            }
            id = hdr;
            id = id.substr(1);
            state = 1;
        } else if (state == 1) {
            inFile.getline(seq, 1024);
            if (!inFile) {
                break;
            }
            string s = seq;
            //
            // id should be unique
            //
            _isaMap[id] = s;
            state = 0;
        }
    } while (true);
    cout << "read " <<  _isaMap.size() << " Fasta sequences\n";
    if (_isaMap.size() > 0)
        return true;
    return false;
}


