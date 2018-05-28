
#include "Header.h"

const char* rar_ver = "0.55 alpha";

int main(int argc, char* argv[]) {
	if (argc <= 3) { cerr << "Coverage per gene estimator.\nNot enough args (3). Aborting.\n(c) Falk Hildebrand"; exit(88); }
	string covF = argv[1];
	string gffF = argv[2];
	string readLs = argv[3];
	float readL = stof(readLs);
	GFFcont * gff = new GFFcont(gffF, readL);
	clock_t tStart = clock();
	cout << "Analysing " << covF << " and " << gffF << endl;
#ifdef _Debug
	cerr << "DEBUG mode\n";
#endif


	gff->read_coverage(covF, gff);

	printf("Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

	cout << "Finished\n";
}

void GFFcont::read_coverage(string covF, GFFcont* gff) {
	istream *opt = NULL;
	string line;
	string goF = covF + ".pergene", coF = covF + ".percontig", coMF = covF + ".median.percontig" , goMF = covF + ".median.pergene";
	string geneStatsF = covF + ".geneStats", wiF = covF + ".window";
	string countF = covF + ".count_pergene";
	float geneCovCutoff = 0.1f;
	float ctgCovCutoff = 0.1f;
	if (!fileExists(covF)) {
		cerr << "Can't find infile " << covF << endl;
		exit(63);
	}

	ofstream genesOut, ctgOut, winOut, CountsOut, genesMedOut, ctgMedOut;
	//first write stats on genes
	genesOut.open(geneStatsF.c_str());
	genesOut << "GeneNumber\tAvgGeneLength\tAvgComplGeneLength\tBpGenes\tBpNotGenes\t";
	genesOut << "Gcomplete\tG5pComplete\tG3pComplete\tGincomplete\n";
	if (totGenes == 0 || genesCompl == 0) {
		cerr << "TotGenes 0\n";
	}
	else {
		genesOut << totGenes << "\t" << geneBP / totGenes << "\t" << geneComplBP / genesCompl << "\t" << geneBP << "\t" << uncovBP << "\t";
		genesOut << genesCompl << "\t" << gene5pCompl << "\t" << gene3pCompl << "\t" << geneIncompl << endl;
		genesOut.close();
	}
	//now real readCoverage
	string substGZ = covF.substr(covF.length() - 3);
	if (substGZ == ".gz") {
#ifdef _Debug
		cerr << "Read g";
#endif

#ifdef _gzipread	
		opt = new igzstream(covF.c_str(), ios::in);
#else
		cerr << "gzip not supported in your build\n"; exit(33);
#endif
#ifdef _Debug
		cerr << "zip\n";
#endif
	}
	else {
		opt = new ifstream(covF.c_str(), ios::in);
	}
	if (!*opt) { cerr << "Can not open coverage file " << covF << endl; exit(33); }
	genesOut.open(goF.c_str());
	//set formatting
	genesOut.setf(std::ios_base::fixed, std::ios_base::floatfield);
	genesOut.precision(2);
	if (!genesOut) { cerr << "Cant open out1 file " << goF << endl; exit(34); }

	
	genesMedOut.open(goMF.c_str());
	//set formatting
	genesMedOut.setf(std::ios_base::fixed, std::ios_base::floatfield);
	genesMedOut.precision(2);
	if (!genesMedOut) { cerr << "Cant open out1b file " << goMF << endl; exit(34); }

	
	CountsOut.open(countF.c_str());
	//set formatting
	CountsOut.setf(std::ios_base::fixed, std::ios_base::floatfield);
	CountsOut.precision(2);
	if (!CountsOut) { cerr << "Cant open out3 file " << countF << endl; exit(34); }

	ctgOut.open(coF.c_str());
	ctgOut.setf(std::ios_base::fixed, std::ios_base::floatfield); ctgOut.precision(2);
	//median contig coverage.. might be better to get rid of random deviations
	ctgMedOut.open(coMF.c_str());
	ctgMedOut.setf(std::ios_base::fixed, std::ios_base::floatfield); ctgOut.precision(2);

	

	winOut.open(wiF.c_str());
	winOut.setf(std::ios_base::fixed, std::ios_base::floatfield);
	winOut.precision(2);
	string currentChrStr("");
	bool dynami = false;
	CHRidmapsIT curChr;// int cst(0), cen(0); 
	vector<int> cov(0); bool notCount(true);
	//chromosome specific counts
	int len(0); int chrCnt(0);
	while (std::getline(*opt, line, '\n')) {
		if (line.substr(0, 1) == "#") { continue; }
#ifdef _Debug
		cerr << line;
#endif

		string segs;
		stringstream ss;
		ss << line;
		std::getline(ss, segs, '\t');//0
		if (segs != currentChrStr) {//new chromosome container
			if (!notCount) {//counted something, make connection
				writeGeneCnts(curChr->second, cov, genesOut, genesMedOut,CountsOut, geneCovCutoff);
				makeWindowCnts(cov, winOut, curChr->second->getID());
			}
			//first write old counts out
			if (currentChrStr != "") {
				float ctgCov((float)chrCnt / (float)len);
				writeCtgCov(ctgCov, ctgCovCutoff, currentChrStr, ctgOut);
				//and write median
				float medCtCov = (float)CalcMHWScore(cov);
				writeCtgCov(medCtCov, ctgCovCutoff, currentChrStr, ctgMedOut);

			}
			curChr = track.find(segs);
			if (curChr == track.end()) {
				cout << "Couldn't find chromosome " << segs << endl;
				notCount = true;
			}
			else {
				notCount = false;
			}


			//reset to new chromosome
			currentChrStr = segs;
			//size_t pos1 = currentChrStr.find("_length_");
			size_t pos1 = currentChrStr.find("_L=");
			size_t pos2 = currentChrStr.find_first_of(";=", pos1 + 3);
			if (pos1 == string::npos || pos2 == string::npos) {
				dynami = true;
				cov.clear();
			}
			else {
				string tmp = currentChrStr.substr(pos1 + 3, pos2 - pos1 - 3);
				len = atoi(tmp.c_str());
				len++;
				cov.clear();
				cov.resize(len, 0);
			}

			chrCnt = 0;
		}

		std::getline(ss, segs, '\t');//1
		int st = atoi(segs.c_str());
		std::getline(ss, segs, '\t');//2
		int en = atoi(segs.c_str());
		std::getline(ss, segs, '\t');//3
		int ccov = atoi(segs.c_str());
		chrCnt += ccov* (en - st);
		if (dynami) {
			cov.resize(en, 0);
			for (; st < en; st++) {
				cov[st] = ccov;
			}
		}
		else {
			if (en >= len) { en = len; }
			if (notCount) { continue; }//no use in position specific counting
			for (; st < en; st++) {
				cov[st] = ccov;
			}
		}
	}
#ifdef _Debug
	cerr << "File read\n";
#endif
	if (!notCount) {//counted something, make connection
		writeGeneCnts(curChr->second, cov, genesOut, genesMedOut, CountsOut, geneCovCutoff);
		makeWindowCnts(cov, winOut, curChr->second->getID());
		
		//write coverage
		float ctgCov((float)chrCnt / (float)len);
		if (dynami) {
			ctgCov = (float)chrCnt / (float)cov.size();
		}
		writeCtgCov(ctgCov, ctgCovCutoff, currentChrStr, ctgOut);
		//and write median
		float medCtCov = (float)CalcMHWScore(cov);
		writeCtgCov(medCtCov, ctgCovCutoff, currentChrStr, ctgMedOut);

	}

	delete opt;
	genesOut.close(); ctgOut.close(); winOut.close();
	ctgMedOut.close(); genesMedOut.close();
}

void GFFcont::writeCtgCov(float ctgCov, float ctgCovCutoff, string currentChrStr, ofstream& ctgOut){
	if (ctgCov >= ctgCovCutoff) {
		ctgOut << currentChrStr << "\t" << ctgCov << endl;
	}
	else {
		ctgOut << currentChrStr << "\t0\n";
	}
}

void GFFcont::writeGeneCnts(Chromo* chr, const vector<int>& cov, ofstream& of, ofstream& ofMed, ofstream& cntsOF, float geneCovCutoff) {
	string chid = chr->getID();
	int chsize = chr->size();
	int cvsiz = (int)cov.size();
	for (int i = 0; i < chsize; i++) {
		float curcov(0); int x(chr->getSt());
		uint sta (x);
		for ( ; x < chr->getEn(); x++) {
			if (x >= cvsiz) { break; }
			curcov += (float)cov[x];
			
		}
		//create subvector, to calc median of
		vector<int> newVec(1, 0); // empty ini
		if (sta < cov.size()) {
			vector<int>::const_iterator first = cov.begin() + sta;
			vector<int>::const_iterator last = cov.begin() + x;
			newVec = vector<int>(first, last);
		}
		double medCov = CalcMHWScore(newVec);

		float curCnt = curcov;
		curcov /= chr->geneSize(); chr->idxIncr();
		if (curcov >= geneCovCutoff) {
			of << chid + "_" << i + 1 << "\t" << curcov << endl;
		}
		if (medCov >= geneCovCutoff) {
			ofMed << chid + "_" << i + 1 << "\t" << medCov << endl;
		}
		curCnt /= readL;
		if (curCnt >= geneCovCutoff) {
			cntsOF << chid + "_" << i + 1 << "\t" << curCnt << endl;
		}
	}
}
//write coverage in floating window
void GFFcont::makeWindowCnts(const vector<int>& cov, ofstream& of, const string chrN) {
	int AQS = 0;
	//int TotC = 0;
	uint W = 100; uint X = 50;
	unsigned int CS = cov.size();//static_cast<unsigned int> (Qual.size());
	if (W >= CS) { W = CS; } // too short

							 //1st loop to ini window
	for (unsigned int i = 0; i<(unsigned int)W; i++) {
		AQS += cov[i];
	}
	//TotC = AQS;
	//1st window
	of << chrN << "\t" << CS << endl;
	int W2 = W / 2;
	of << W2 << "\t" << (AQS / W) << endl;
	uint subCnt(0);
	for (unsigned int i = W; i<(unsigned int)CS; i++) {
		AQS += cov[i] - cov[i - W]; subCnt++;
		//TotC += cov[i];
		if (subCnt >= X) {
			subCnt = 0;
			of << (i - W2) + 1 << "\t" << (AQS / W) << endl;
		}
	}

}

GFFcont::GFFcont(const string & inf, float RL) :
	readL(RL), totGenes(0), geneBP(0), uncovBP(0), geneComplBP(0), genesCompl(0), gene5pCompl(0), gene3pCompl(0)
{
	ifstream opt;
	string line;
	opt.open(inf.c_str(), ios::in);
	if (!opt) { cerr << "Can't open " << inf << endl; exit(33); }
	string currentChrStr(""); //int curLength = 0;
	CHRidmapsIT curChr; int cst(-1), cen(-1);
	while (std::getline(opt, line, '\n')) {
		if (line.substr(0, 1) == "#") { continue; }
		string segs;
		stringstream ss;
		ss << line;
		std::getline(ss, segs, '\t');//0
		if (segs != currentChrStr) {//new chromosome container
			if (cst >= 0) {//already read a chromosome in
				vector<int> tmp = curChr->second->BpGeneCovered();
				geneBP += tmp[0];
				uncovBP += tmp[1];
				geneComplBP += tmp[2];
				totGenes += curChr->second->size();
			}
			track[segs] = new Chromo(segs);
			curChr = track.find(segs);
			currentChrStr = segs;
		}
		std::getline(ss, segs, '\t');//1
		std::getline(ss, segs, '\t');//2
		std::getline(ss, segs, '\t');//3
		cst = atoi(segs.c_str());
		std::getline(ss, segs, '\t');//4
		cen = atoi(segs.c_str());
		std::getline(ss, segs, '\t');//5
		std::getline(ss, segs, '\t');//6
		std::getline(ss, segs, '\t');//7
		std::getline(ss, segs, '\t');//8
		size_t pos(segs.find("partial="));
		bool cpl(false);
		string geneEnds = segs.substr(pos + 8, segs.find(";", pos + 8) - pos - 8);
		if (geneEnds == "00") {
			genesCompl++;
			cpl = true;
		}
		else if (geneEnds == "10") { //5 prime correct
			gene5pCompl++;
		}
		else if (geneEnds == "01") {
			gene3pCompl++;
		}
		else {
			geneIncompl++;
		}

		curChr->second->addGene(cst, cen, cpl);
	}
	opt.close();
	if (cst >= 0) {//already read a chromosome in
		vector<int> tmp = curChr->second->BpGeneCovered();
		geneBP += tmp[0];
		uncovBP += tmp[1];
		geneComplBP += tmp[2];
		totGenes += curChr->second->size();
	}

	cout << "Read GFF\n";

}

////////////////////////////////////////////////////////////
//                    Chromo

Chromo::Chromo(string na) :ID(na), idx(0), start(0), end(0) {
	size_t pos = ID.find("_L=");
	size_t pos2 = ID.find_first_of(";=", pos + 1);
	if (pos == string::npos || pos2 == string::npos) {
		ChrSize = -1;
	}
	else {
		string x = ID.substr(pos + 3, pos2 - pos - 3);
		ChrSize = atoi(ID.substr(pos + 3, pos2 - pos - 3).c_str());
	}
}
vector<int> Chromo::BpGeneCovered() {
	vector<int> BPcov(3, 0);//0=covered by gene; 1= not covered by gene; 2=covered by complete gene
	for (size_t i = 0; i < start.size(); i++) {
		//not covered by genes
		if (i == 0) {
			BPcov[1] += start[i];
		}
		if (i == (start.size() - 1)) {
			if (ChrSize > -1) {
				BPcov[1] += ChrSize - end[i];
			}
		}
		else {
			BPcov[1] += start[i + 1] - end[i];
		}
		//covered by genes
		BPcov[0] += end[i] - start[i];
		if (complete[i]) {
			BPcov[2] += end[i] - start[i];
		}
	}
	return BPcov;
}
