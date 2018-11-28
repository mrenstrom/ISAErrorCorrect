//
//  main.cpp
//  mapris
//
//  Created by mark enstrom on 3/26/18.
//  Copyright © 2018 mark enstrom. All rights reserved.
//

#include "readFasta.hpp"
#include "NWAlign.hpp"
#include "ReadAlign.hpp"
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>
#include <iomanip>
#include <unistd.h>
using namespace std;



/*--------------------------------------------------------------------------------------------
 * main routine for filtering fasta for ISA processing
 *    fasta file
 *    align file
 *    output
 *--------------------------------------------------------------------------------------------*/
int main(int argc, const char * argv[]) {
    //
    // initial read of source
    //
    std::cout << "Read FASTA\n";
    if (argc < 4) {
        cout << "Usage: ISAErrorCorrent fasta align output log\n";
        exit(0);
    }
    string fasta = argv[1];
    string align = argv[2];
    string output = argv[3];
    string logf   = argv[4];
    
    cout << "Read fasta file " << fasta << endl;
    cout << "Read   aln file " << align << endl;
    
    ReadFasta rf = ReadFasta();
    bool bRet = rf.readFastA(fasta);
    if (!bRet) {
        cout << "error reading file " << fasta << endl;
        exit(0);
    }
    ReadAlign rd = ReadAlign(rf._isaMap);
    bRet = rd.read(align);
    if (!bRet) {
        cout << "error reading file " << align << endl;
        exit(0);
    }
    //
    // search inidiviual alignments for exact matches
    //
    rd.buildAlignmentGroups();
    //
    // assign sequence to groups
    //
    rd.assignSequence();
    //
    // build consensus sequence for each group
    //
    rd.buildConsensusSequence(output);
    //
    // attempt to errror-correct unassigned sequeneces
    //     -- those that failed alignment by blat
    //
    rd.errorCorrectSequence();
    //
    // do any seq on different align match
    //
    rd.errorCorrectAlignments(logf);
    //
    // output results
    //
    rd.writeOutput(output);
    return 0;
}
