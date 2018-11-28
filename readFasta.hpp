//
//  readFasta.hpp
//  mapris
//
//  Created by mark enstrom on 3/26/18.
//  Copyright © 2018 mark enstrom. All rights reserved.
//

#ifndef readFasta_hpp
#define readFasta_hpp

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
using namespace std;





#define KMER 8

class Capture {
public:
    string  _seq;
    string  _id;
    size_t  _count;
    unordered_map<std::string,int> _famSeq;
    unordered_map<std::string,size_t> _seqMap;
    bool    _valid;
    int     _span;
    Capture(string seq,string id,int i) {
        _seq = seq;
        _id = id;
        _count = i;
        _valid=true;
        _famSeq = unordered_map<std::string,int>();
        _seqMap = unordered_map<std::string,size_t>();
        //
        // !!! duplicate just to keep counts accurate
        // !!! can delete later
        //
        _famSeq[seq] = i;
        //
        // build seq map (8) out to max of 40
        //
        size_t l = seq.length() - KMER;
        if (l > (40-KMER)) l = (40-KMER);
        for (int i = 0; i < (l-KMER);i++) {
            string s = seq.substr(i,KMER);
            ++_seqMap[s];
        }
    };
    void addSeq(string seq) {
        ++_famSeq[seq];
        _count++;
    }
    void merge(Capture *p) {
        for (auto pr:p->_famSeq) {
            string s = pr.first;
            int count = pr.second;
            _famSeq[s] += count;
        }
    }
    
    bool similar(string s2) {
        for (int i = 0; i < (s2.length()-KMER); i += KMER) {
            if  (_seqMap.find(s2.substr(i,KMER)) != _seqMap.end()) {
                return true;
            }
        }
        return false;
    }
    
    void dump(string s) {
        for (int i = 0; i < (s.length()-8); i += 8) {
            cout << "comp " << s.substr(i,8) << endl;
        }
        
        for (auto pr:_seqMap) {
            cout << pr.first << "\t" << pr.second << endl;
        }
    }
};


typedef Capture *PCapture;
typedef vector<PCapture> CapVec;
typedef CapVec *PCapVec;


typedef std::unordered_map<std::string,std::string> ISAMap;


class ReadFasta {
public:
    vector<string> _sourceFiles;
    string _dstDir;
    string _test;
    string _fileID;
    //
    // read stats
    //
    CapVec _capture;

    ISAMap _isaMap;
    
    ReadFasta() {
        _isaMap = ISAMap();
        _capture = CapVec();
    }
    bool readFastA(string filename);
    bool findSequence(string seq, string tag, string target);

    bool readQ(string filename,string outName);
    bool addSequence(string s);

};

#endif /* readFasta_hpp */
