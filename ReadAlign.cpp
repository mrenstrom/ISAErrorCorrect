//
//  ReadAlign.cpp
//  ISAErrorCorrect
//
//  Created by mark enstrom on 10/10/18.
//  Copyright © 2018 mark enstrom. All rights reserved.
//

#include "ReadAlign.hpp"


struct StringSortL {
    bool operator() (string s1,string s2) { return (s1.length() > s2.length());}
} seqStrComp;


struct PagSort {
    bool operator() (PAlignGroup p1,PAlignGroup p2) { return (p1->_count > p2->_count);}
} pagSortSmall;


struct PAMSort {
    bool operator() (pair<string,PAlignVector> p1,pair<string,PAlignVector> p2) { return (p1.second->at(0)->_count > p2.second->at(0)->_count);}
} pamSortSmall;

std::string exec(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");
    while (!feof(pipe.get())) {
        if (fgets(buffer.data(), 128, pipe.get()) != nullptr)
            result += buffer.data();
    }
    return result;
}

/*--------------------------------------------------------------------------------------------
 *
 * find good fasta sequences:
 *      good primer and ltr
 *      no vector
 *
 *
 *
 *--------------------------------------------------------------------------------------------*/

int countBases(char *p) {
    int count = 0;
    while (*p++ != '\0') {
        if (
            (*p == 'A') ||
            (*p == 'C') ||
            (*p == 'G') ||
            (*p == 'T')
            )
            count++;
    }
    return count;
}
/*--------------------------------------------------------------------------------------------
 *
 *  read output FASTA from clustal
 *
 *
 *--------------------------------------------------------------------------------------------*/
bool readTmpFasta(string fname,vector<string> &alStr) {
    std::ifstream inFile;
    inFile.open(fname);
    char seq[512];
    //
    // read first header line
    //
    inFile.getline(seq, 512);
    if (!inFile) {
        return(false);
    }
    bool bSequenceComplete = false;
    string s = "";
    do {
        inFile.getline(seq, 512);
        if (!inFile) {
            break;
        }
        
        if (seq[0] == '>') {
            bSequenceComplete = true;
        } else {
            s += string(seq);
        }
        
        if (bSequenceComplete) {
            alStr.push_back(s);
            s = "";
            bSequenceComplete = false;
        }
    } while (true);
    //
    // save final
    //
    if (s != "") {
        alStr.push_back(s);
    }
    //cout << "# seq = " << alStr.size() << endl;
    
    return (true);
}
/*--------------------------------------------------------------------------------------------
 * FindMajorSequence
 * given vector of aligned sequences, find max at each position. If max is "-" delete.
 *
 *
 *--------------------------------------------------------------------------------------------*/
string findMajorSequence(vector<string> &alStr) {
    if (alStr.size() < 2) return "";
    //    for (auto s:alStr) {
    //        cout << "\t\t" << s <<endl;
    //    }
    //
    // get seq length
    //
    string s1 = alStr[0];
    size_t sl = s1.length();
    char cMajor[1024] = {0};
    int targetIndex = 0;
    for (int i = 0; i < sl; i++) {
        int iA = 0;
        int iC = 0;
        int iG = 0;
        int iT = 0;
        int iD = 0;
        for (auto s:alStr) {
            if (s[i] == 'A') {
                iA++;
            } else if (s[i] == 'C') {
                iC++;
            } else if (s[i] == 'G') {
                iG++;
            } else if (s[i] == 'T') {
                iT++;
            } else if (s[i] == '-') {
                iD++;
            } else {
                cout << "Error ins eq " << s << endl;
            }
        }
        if ((iA >= iC) && (iA >= iG) && (iA >= iT) && (iA >= iD)) {
            cMajor[targetIndex++] = 'A';
        } else if ((iC >= iA) && (iC >= iG) && (iC >= iT) && (iC >= iD)) {
            cMajor[targetIndex++] = 'C';
        } else if ((iG >= iA) && (iG >= iC) && (iG >= iT) && (iG >= iD)) {
            cMajor[targetIndex++] = 'G';
        } else if ((iT >= iA) && (iT >= iC) && (iT >= iG) && (iT >= iD)) {
            cMajor[targetIndex++] = 'T';
        } else if ((iD > iA) && (iD > iC) && (iD > iG) && (iD > iT)) {
            // nothing added
        } else {
            cout << "tie!!!\n";
        }
        
    }
    
    string s = string(cMajor);
    //cout << "..........................\n";
    //cout << "\t\t" << s << endl;
    return s;
}


std::string trim(const std::string& str,
                const std::string& whitespace = " \t")
{
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content
    
    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;
    
    return str.substr(strBegin, strRange);
}


/*--------------------------------------------------------------------------------------------
 *
 * build groups from exact match alignments
 * id    count    alignStr    class    Align    chromo    start    dir    span    Alignments    r    refseq
 * MAPRIS:2:17326:1       1    chr8                114403515    -    1.00    169.00
 * MAPRIS:2:17326:1       1    chr6                 72417397    +    1.00    169.00
 * MAPRIS:2:17326:1       1    chr11                  523799    -    1.00    169.00*
 * MAPRIS:2:17327:1       .......................
 *
 *
 * one sequence may show on several lines...these are multi-align. Only the first count is included
 *
 *
 *--------------------------------------------------------------------------------------------*/

bool
ReadAlign::read(string filename)
{
    std::ifstream inFile;
    inFile.open(filename);
    if (!inFile) {
        return false;
    } 
    //
    // read first header line
    //
    string s = "";
    string s_id;
    string s_count;
    string s_chr;
    string s_startPos;
    string s_dir;
    string s_r;
    string s_q;
    //
    // skip header
    //
    getline(inFile, s_q);
    do {
        getline(inFile, s_id, '\t');
        if (!inFile) {
            break;
        }
        if (s_id.substr(0,1) == "#") {
            // comment, read rest of line
            getline(inFile, s_id);
            continue;
        }
        
        getline(inFile, s_count, '\t');
        getline(inFile, s_chr,'\t');
        getline(inFile, s_startPos,'\t');
        getline(inFile, s_dir,'\t');
        getline(inFile, s_r, '\t');
        getline(inFile, s_q);
        
        s_count = trim(s_count);
        s_chr = trim(s_chr);
        s_startPos = trim(s_startPos);
        s_dir = trim(s_dir);
        
        Align *pa = new Align();
        pa->_id = s_id;
        pa->_count = stoi(s_count);
        pa->_chr = s_chr;
        pa->_startPos = stoi(s_startPos);
        pa->_dir = s_dir;
        pa->_r = stod(s_r);
        pa->_q = stod(s_q);
        //
        // new or add to existing...all pa have some [count] but different alignment
        //
        try {
            PAlignVector pav = _alignMap.at(s_id);
            //
            // adding multi-alignemnts to one sequenceID
            //
            pav->push_back(pa);
            
        } catch(std::out_of_range) {
            //
            // add new sequenceID with one alignment
            //
            PAlignVector pval = new AlignVector();
            pval->push_back(pa);
            _alignMap[s_id] = pval;
        }
        
    } while (true);
    
    cout << "read " << _alignMap.size() << " Alignments " << endl;
    
    return true;
}

/*--------------------------------------------------------------------------------------------
 *
 * build groups from exact match alignments
 *
 *
 *
 *
 *--------------------------------------------------------------------------------------------*/
void ReadAlign::buildAlignmentGroups() {
    //
    // _alignMap is a map of id:alignmentVector read from input file.
    // _alignMap contains all alignments
    //
    // convert this map to a vector and sort vector on count
    //
    PAM pam = PAM();
    for (auto pr:_alignMap) {
        pam.push_back(pr);
    }
    std::sort(pam.begin(),pam.end(),pamSortSmall);
    //
    // start with empty group vector
    //
    _pagVec->clear();
    //
    // group together sequences with matching alignments
    //
    // pam is a vector of (ID,AlignmentVector)
    //
    for (auto pr:pam) {
        string id = pr.first;
        //cout << id << " " << _pagVec->size() << endl;
        PAlignVector pav = pr.second;
        // get top alignment
//        PAlign pa = pav->at(0);
//        bool found = false;
//        for (auto pg: *_pagVec) {
//            if ((pg->_chr == pa->_chr) &&
//                (pg->_dir == pa->_dir) &&
//                (abs(pg->_startPos - pa->_startPos) < 100)
//                ) {
//                pg->_alignMap[id] = pav;
//                pg->_count += pa->_count;
//                found = true;
//                break;
//            }
//        }
        //
        // multi
        //
        bool found = false;
        //
        // for each alignment in vector of alignments for this ID
        //
        int pa1Count = 0;
        for (PAlign pa1 : *pav) {
            //
            // for all alignments in each ID of the group
            //
            for (auto pg: *_pagVec) {
                for (auto pr : pg->_alignMap) {
                    PAlignVector pav2 = pr.second;
                    int pa2Count = 0;
                    for (auto pa2: *pav2) {
                        if ((pa1->_chr == pa2->_chr) &&
                            (pa1->_dir == pa1->_dir) &&
                            (abs(pa1->_startPos - pa2->_startPos) < 10)
                         ) {
                            pg->_alignMap[id] = pav;
                            pg->_count += pa1->_count;
                            found = true;
                            //cout << "Found " << endl;
                            //cout << "    " << pa1->_chr << " " << pa1->_startPos << " " << pa1Count << endl;
                            //cout << "    " << pa2->_chr << " " << pa2->_startPos << " " << pa2Count << endl;
                            //if (abs(pa1->_startPos - pa2->_startPos) > 10) {
                            //    cout << abs(pa1->_startPos - pa2->_startPos) << endl;
                            //}
                            break;
                        }
                        pa2Count++;
                    }
                    if (found) break;
                }
                if (found) break;
            }
            if (found) break;
            pa1Count++;
        }
        //
        // if not found start new group
        //
        PAlign pa = pav->at(0);
    
        if (!found) {
            //cout << "New " << pa->_chr << " " << pa->_startPos << endl;
            AlignGroup *pg = new AlignGroup();
            pg->_parentID = pa->_id;
            pg->_chr = pa->_chr;
            pg->_dir = pa->_dir;
            pg->_startPos = pa->_startPos;
            pg->_count = pa->_count;
            pg->_alignMap[id] = pav;
            _pagVec->push_back(pg);
        }
    }
    cout << "Found " << _pagVec->size() << " alignment groups " << endl;
}

/*--------------------------------------------------------------------------------------------
 *
 * assign sequences to alignment groups
 *
 *
 *
 *
 *--------------------------------------------------------------------------------------------*/
void ReadAlign::assignSequence()
{
    for (auto pg:*_pagVec) {
        //
        // table lookup sequence for each id in the group
        //
        for (auto pr:pg->_alignMap) {
            string s_id = pr.first;
            //
            // look for seq s
            //
            string seq = "";
            try {
                seq = _availableSeq.at(s_id);
                _availableSeq.erase(s_id);
                //cout << s << " " << seq << endl;
            } catch(std::out_of_range) {
                cout << "ERROR: can't find seq for ID:" << s_id << endl;
            }
            pg->_seqMap[seq] += 1;
        }
        //cout << "SeqMap size = " << pg->_seqMap.size() << endl;
    }
    cout << "\n" << "Unaligned sequences after group assigment = " << _availableSeq.size() << endl;
}

/*--------------------------------------------------------------------------------------------
 *
 * build consensus seq for each alignment group
 *
 *
 *
 *
 *--------------------------------------------------------------------------------------------*/
void ReadAlign::buildConsensusSequence(string fname)
{
    cout << "Build consensus sequences\n";
    std::string sTmpFasta = fname + "_tmp.fasta";
    std::string sTmpOut  =  fname + "_tmp.aln";
    for (auto pg:*_pagVec) {
        if (pg->_seqMap.size() > 1) {
            //
            // convert to vector
            //
            vector<string> sVec = vector<string>();
            for (auto pr:pg->_seqMap) {
                string s = pr.first;
                sVec.push_back(s);
            }
            //
            // sort for length
            //
            std::sort(sVec.begin(),sVec.end(),seqStrComp);
            //
            // find min for first 100
            //
            int minS = 1000;
            int si = 1;
            for (auto s: sVec) {
                if (s.length() < minS) minS = (int)s.length();
                if (++si > 100) break;
            }
            //cout << "Find consensus for " << pg->_seqMap.size() << endl;
            std::ofstream tmpFasta;
            tmpFasta.open(sTmpFasta);
            si = 1;
            for (auto s:sVec) {
                tmpFasta << ">MAPRIS" << to_string(si++) << "\n";
                string s2 = s.substr(0,minS);
                tmpFasta << s2 << endl;
                //cout << s << endl;
                if (si > 100) break;
            }
            tmpFasta.close();
            string cmdLine = "clustalo -i " + sTmpFasta + " --full --outfmt=a2m --outfile=" + sTmpOut + " --force";
            //string cmdLine = "/Users/mark_enstrom/Barcode/RIS/clustal/clustalo -i " + sTmpFasta + " --full --outfmt=a2m --outfile=" + sTmpOut + " --force";
            exec(cmdLine.c_str());
            vector<string> alStr = vector<string>();
            readTmpFasta(sTmpOut,alStr);
            string s = findMajorSequence(alStr);
            //cout << s << endl;
            pg->_conSeq = s;
            if (pg->_conSeq.length() < 30) {
               cout << "Error in consenses : length < 30\n";
               cout << pg->_conSeq << endl;
               for (auto as: alStr) {
                 cout << "    " << as << endl;
               }
               si = 0;
               for (auto s:sVec) {
                 tmpFasta << ">MAPRIS" << to_string(si++) << "\n";
                 string s2 = s.substr(0,minS);
                 tmpFasta << s2 << endl;
                 cout << s << endl;
                 if (si > 100) break;
               }
            }
        } else {
            //
            // assign consenus to only seq
            //
            pg->_conSeq = pg->_seqMap.begin()->first;
            if (pg->_conSeq.length() < 30) {
               cout << "consensus seq = 0 from default\n";
            }
        }
    }
    //
    // cleanup
    //
    std::remove(sTmpFasta.c_str());
    std::remove(sTmpOut.c_str());
}

/*--------------------------------------------------------------------------------------------
 *
 * attempt to error-correct unaligned sequences
 *
 *
 *
 *
 *--------------------------------------------------------------------------------------------*/
void ReadAlign::errorCorrectSequence() {
    cout << "\n" << "Attempt to match unaligned to aligned sequences= " << _availableSeq.size() << endl;
    int numMatch = 0;
    for (auto pr:_availableSeq) {
        string UID = pr.first;
        string s1 = pr.second;
        for (auto pg:*_pagVec) {
            NWAlign nwAlign = NWAlign();
            string s2 = pg->_conSeq;
            if (s1.length() > s2.length()) {
                s1 = s1.substr(0,s2.length());
            } else {
                s2 = s2.substr(0,s1.length());
            }
            int maxError = (int)s1.length() / 8;
            int editDist = nwAlign.alignWithLeadingGap(s1,s2);
            if (editDist < maxError) {
                //cout << "match " << s1.length() << " " << editDist << endl;
                //cout << "    " << s1 << endl;
                //cout << "    " << s2 << endl;
                //
                // add sequence to alignment
                //
                pg->_unalignedSequences[pr.first] += 1;
                //
                // get count from ID    MAPRIS:file:seq:count
                //
                size_t lc = UID.find_last_of(":");
                string sCount = UID.substr(lc+1);
                int c = stoi(sCount);
                pg->_count += c;
                numMatch++;
                break;
            }
        }
    }
    cout << "Error corrected unaligned sequences =  " << numMatch << endl;
}

/*--------------------------------------------------------------------------------------------
 *
 * attempt to error-correct DIFFERENT alignments with similar consensus sequence
 *
 *   only error-correct Multi alignments, not unique ???
 *
 *
 *--------------------------------------------------------------------------------------------*/

void ReadAlign::errorCorrectAlignments(string logf) {
    ofstream myLog;
    myLog.open(logf);
    int numMatch = 0;
    cout << "Search align groups" << endl;
    for (int i = 0; i < _pagVec->size(); i++) {
        auto pi = _pagVec->at(i);
        if (!pi->_valid) continue;
        for (int j = (int)_pagVec->size()-1; j > i; j--) {
            auto pj = _pagVec->at(j);
            if (!pj->_valid) continue;
            string s1 = pi->_conSeq;
            string s2 = pj->_conSeq;
            if (s1.length() > s2.length()) {
                s1 = s1.substr(0,s2.length());
            } else {
                s2 = s2.substr(0,s1.length());
            }
            if ((s1.length() == 0) || (s2.length() == 0)) continue;
            NWAlign nwAlign = NWAlign();
            int maxError = (int)s1.length() / 8;
            int editDist = nwAlign.alignWithLeadingGap(s1,s2);
            if (editDist <= maxError) {
                
                
                
                myLog << "Match: distance = " << editDist << endl;
                myLog << pi->_parentID << endl;
                myLog << s1 << endl;
                myLog << pj->_parentID << endl;
                myLog << s2 << endl;
                
                myLog << setw(22) << pi->_chr << " " << setw(8) << pi->_startPos << " " << pi->_count << endl;
                myLog << setw(22) << pj->_chr << " " << setw(8) << pj->_startPos << " " << pj->_count << endl;
                
                pi->_count+= pj->_count;
                pj->_valid = false;
                numMatch++;
            }
        }
    }
    cout << "Num Align Match = " << numMatch << endl;
    myLog.close();
}
/*--------------------------------------------------------------------------------------------
 *
 * write output
 *
 *
 *
 *--------------------------------------------------------------------------------------------*/
void ReadAlign::writeOutput(string filename)
{
    //
    //
    //
    ofstream myStats;
    myStats.open(filename);
    //
    // header
    //
    myStats << "id\t";
    myStats << "count\t";
    myStats << "alignStr\t";
    myStats << "chromo\t";
    myStats << "start\t";
    myStats << "dir\t";
    myStats << "span\t";
    myStats << "rUnique\t";
    myStats << "cUniqye\t";
    myStats << "refseq\n";
    //
    //
    //
    std::sort(_pagVec->begin(),_pagVec->end(),pagSortSmall);
    
    for (auto pg:*_pagVec) {
        // id of first alignment
        PAlign paMax = nullptr;
        for (auto pr:pg->_alignMap) {
            PAlign pNew = pr.second->at(0);
            if ((paMax == nullptr) || (pNew->_count > paMax->_count)) {
                paMax = pNew;
            }
        }
        myStats << paMax->_id << "\t";
        myStats << pg->_count << "\t";
        myStats << pg->_chr << "_" << pg->_startPos << "_" << pg->_dir << "\t";
        myStats << pg->_chr << "\t";
        myStats << pg->_startPos << "\t";
        myStats << pg->_dir << "\t";
        myStats << 0 << "\t";
        myStats << paMax->_r << "\t";
        myStats << paMax->_q << "\t";
        myStats << pg->_conSeq << "\n";
    }
    myStats.close();
}

