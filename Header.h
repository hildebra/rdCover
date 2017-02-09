#pragma once

#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>

#include <map>
#include <unordered_map>
#include <algorithm>
#include <time.h>

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
#define _gziprea//d
#else
#define _gzipread
#endif

#define __De //bug

using namespace std;

#ifdef _gzipread
#include "gzstream.h"
#endif


typedef unsigned int uint;
typedef unsigned long ulong;

class Chromo{
public:
	Chromo(string na);
	void addGene(int s, int e, bool c) { start.push_back(s); end.push_back(e); complete.push_back(c); }
	int getSt() { return start[idx]; }
	int getEn() { return end[idx]; }
	void idxIncr() { idx++; }
	int size() { return start.size(); }
	string getID() { return ID; }
	float geneSize() { return float(end[idx] - start[idx]); }
	vector<int> BpGeneCovered();
private:
	string ID;
	int idx; int ChrSize;
	vector<int> start, end;
	vector<bool> complete;
};

typedef std::unordered_map<string, Chromo*>::iterator CHRidmapsIT;
typedef std::unordered_map<string, Chromo*> CHRidmaps;

class GFFcont{
public:
	GFFcont(const string &, float );
	~GFFcont() {for (auto itr = track.begin(); itr != track.end(); ++itr) {		delete itr->second;	}}
	void read_coverage(string, GFFcont*);
private:
	void writeGeneCnts(Chromo*, const vector<int>&, ofstream&, ofstream&, float = 0.01);
	void makeWindowCnts(const vector<int>&, ofstream&, const string);
	CHRidmaps track;
	float readL;
	int totGenes, geneBP, uncovBP, geneComplBP;
	int genesCompl, gene5pCompl, gene3pCompl, geneIncompl;
};


inline bool fileExists(const std::string& name) {
	struct stat buffer;
	return (stat(name.c_str(), &buffer) == 0);
}