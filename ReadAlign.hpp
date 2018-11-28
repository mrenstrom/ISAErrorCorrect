//
//  ReadAlign.hpp
//  ISAErrorCorrect
//
//  Created by mark enstrom on 10/10/18.
//  Copyright © 2018 mark enstrom. All rights reserved.
//

#ifndef ReadAlign_hpp
#define ReadAlign_hpp

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <vector>
#include <map>
#include <unordered_map>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <list>
#include <stdlib.h>
#include <assert.h>
#include <ctime>
#include <unistd.h>
#include <algorithm>

#include <memory>
#include <stdexcept>
#include <array>
#include <iomanip>
#include <unistd.h>



#include "readFasta.hpp"
#include "NWAlign.hpp"



using namespace std;


class Align {
public:
    string _id;
    int    _count;
    string _chr;
    int    _startPos;
    string _dir;
    double _r;
    double _q;
};

typedef Align* PAlign;
typedef vector<PAlign> AlignVector;
typedef AlignVector* PAlignVector;

typedef std::unordered_map<std::string,PAlignVector> AlignMap;
typedef std::unordered_map<std::string,int> SeqMap;


typedef vector<pair<string,PAlignVector>> PAM;



class AlignGroup {
public:
    string _chr;
    int _startPos;
    string _dir;
    int _count;
    bool _valid;
    string _parentID;
    // alignment ID,<PAlign>
    AlignMap _alignMap;
    // exact alignmment match ID/sequences
    SeqMap _seqMap;
    // unaligned error-corrected ID/sequences
    SeqMap _unalignedSequences;

    string _conSeq = "";
    AlignGroup() {
        _alignMap = AlignMap();
        _seqMap = SeqMap();
        _unalignedSequences = SeqMap();
        _valid = true;
    }
};
typedef AlignGroup *PAlignGroup;
typedef vector<PAlignGroup> AGVec;
typedef AGVec *PAGVec;



class ReadAlign {
public:
    AlignMap _alignMap;
    PAGVec _pagVec;
    ISAMap _availableSeq;
    
    
    ReadAlign(ISAMap &iMap) {
        _alignMap = AlignMap();
        _pagVec = new AGVec();
        for (auto pr:iMap) {
            _availableSeq[pr.first] = pr.second;
        }
    }
    
    
    bool read(string filename);
    void buildAlignmentGroups();
    void assignSequence();
    void buildConsensusSequence(string fname);
    void errorCorrectSequence();
    void errorCorrectAlignments(string logf);
    void writeOutput(string filename);
};









#endif /* ReadAlign_hpp */
