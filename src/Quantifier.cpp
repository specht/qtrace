/*
Copyright (c) 2007-2008 Michael Specht

This file is part of SimQuant.

SimQuant is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SimQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with SimQuant.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "Quantifier.h"
#include <QtCore>
#include <QtSvg>
#include <math.h> 


k_Quantifier::k_Quantifier(r_ScanType::Enumeration ae_ScanType,
						   QList<tk_IntPair> ak_MsLevels,
						   int ai_IsotopeCount, int ai_MinCharge, int ai_MaxCharge, 
						   double ad_MinSnr, double ad_MassAccuracy, QString as_SvgOutPath, 
						   QIODevice* ak_TextOutDevice_, QIODevice* ak_YamlOutDevice_,
						   bool ab_PrintStatistics)
	: k_ScanIterator(ae_ScanType, ak_MsLevels)
	, mi_WatchIsotopesCount(ai_IsotopeCount)
	, mk_TextOutStream(ak_TextOutDevice_)
	, mk_YamlOutStream(ak_YamlOutDevice_)
	, mi_MinCharge(ai_MinCharge)
	, mi_MaxCharge(ai_MaxCharge)
	, md_MinSnr(ad_MinSnr)
	, md_MassAccuracy(ad_MassAccuracy)
	, md_ElutionProfilePeakWidth(0.1)
	, ms_SvgOutPath(as_SvgOutPath)
	, mb_PrintStatistics(ab_PrintStatistics)
{
	// test parameter sanity
	if (mi_MinCharge < 2 || mi_MinCharge > 10)
	{
		printf("Error: minimum charge (%d) out of range.\n", mi_MinCharge);
		exit(1);
	}
	if (mi_MaxCharge < 2 || mi_MaxCharge > 10)
	{
		printf("Error: maximum charge (%d) out of range.\n", mi_MaxCharge);
		exit(1);
	}
	if (mi_MaxCharge < mi_MinCharge)
	{
		printf("Error: minimum charge (%d) bigger than maximum charge(%d).\n", mi_MinCharge, mi_MaxCharge);
		exit(1);
	}
	if (md_MassAccuracy < 0.0)
	{
		printf("Error: mass accuracy must not be less than zero.\n");
		exit(1);
	}
	
	if (!ms_SvgOutPath.isEmpty())
	{
		QStringList lk_Files = QDir(ms_SvgOutPath).entryList(QDir::NoDotAndDotDot);
		if (!lk_Files.empty())
		{
			printf("Error: There are files at %s. Please remove them or specify an SVG output path which contains no files.\n", ms_SvgOutPath.toStdString().c_str());
			exit(1);
		}
	}
	
	Q_INIT_RESOURCE(qtrace);
	
	QFile lk_File(":res/AminoAcids.csv");
	lk_File.open(QFile::ReadOnly);
	QTextStream lk_TextStream(&lk_File);

	while (!lk_TextStream.atEnd())
	{
		QString ls_Line = lk_TextStream.readLine().trimmed();
		if (ls_Line.isEmpty() || ls_Line.startsWith("#"))
			continue;

		QStringList lk_List = ls_Line.split(QChar(';'));
		mk_AminoAcidWeight[lk_List[3][0].toAscii()] = lk_List[4].toDouble();
	}
	lk_File.close();
}


k_Quantifier::~k_Quantifier()
{
}


bool compareScanQuantitationBySnr(const r_ScanQuantitationResult& a, const r_ScanQuantitationResult& b)
{
	return a.md_Snr > b.md_Snr;
}


bool compareScanQuantitationByRetentionTime(const r_ScanQuantitationResult& a, const r_ScanQuantitationResult& b)
{
	return a.md_RetentionTime < b.md_RetentionTime;
}


bool sortByMz(const QPair<double, QString>& a, const QPair<double, QString>& b)
{
	return a.first < b.first;
}


QString dtos(double ad_Value)
{
	static char lc_String_[2048];
	sprintf(lc_String_, "%1.2f", ad_Value);
	return QString(lc_String_);
}


void k_Quantifier::quantify(QStringList ak_SpectraFiles, QStringList ak_Peptides)
{
	mk_Peptides = ak_Peptides;
	
	printf("Searching for %d peptide%s in %d file%s, trying charge states %d to %d.\n",
		ak_Peptides.size(), ak_Peptides.size() != 1 ? "s" : "",
		ak_SpectraFiles.size(), ak_SpectraFiles.size() != 1 ? "s" : "",
		mi_MinCharge, mi_MaxCharge);
	QStringList lk_ScanTypes;
	if (me_ScanType & r_ScanType::MS1)
		lk_ScanTypes << "Full";
	if (me_ScanType & r_ScanType::SIM)
		lk_ScanTypes << "SIM";
	
	printf("Looking in %s scans, expecting %d isotope peaks at a mass accuracy of %1.2f ppm.\n",
		lk_ScanTypes.join("/").toStdString().c_str(), mi_WatchIsotopesCount, md_MassAccuracy);
		
	// determine all target m/z values
	typedef QPair<double, QString> tk_DoubleStringPair;
	QList<tk_DoubleStringPair> lk_TempList;
	for (int li_PeptideIndex = 0; li_PeptideIndex < mk_Peptides.size(); ++li_PeptideIndex)
	{
		QString ls_Peptide = mk_Peptides[li_PeptideIndex];
		for (int li_Charge = mi_MinCharge; li_Charge <= mi_MaxCharge; ++li_Charge)
		{
			double ld_PeptideMz = this->calculatePeptideMass(ls_Peptide, li_Charge);
			double ld_ModMz = HEAVY_ARGININE * ls_Peptide.count("R");

			double ld_Mz;
			QString ls_Key;
			
			for (int i = 0; i < mi_WatchIsotopesCount; ++i)
			{
				
				// save unlabeled mass
				ld_Mz = ld_PeptideMz + i * NEUTRON / li_Charge;
				ls_Key = QString("%1-%2-unlabeled-%3").arg(ls_Peptide).arg(li_Charge).arg(i);
				lk_TempList.push_back(tk_DoubleStringPair(ld_Mz, ls_Key));
				
				// save labeled mass
				ld_Mz = ld_PeptideMz + i * NEUTRON / li_Charge + ld_ModMz / li_Charge;
				ls_Key = QString("%1-%2-labeled-%3").arg(ls_Peptide).arg(li_Charge).arg(i);
				lk_TempList.push_back(tk_DoubleStringPair(ld_Mz, ls_Key));
			}
			// save forbidden peak (one to the left from the unlabeled A+0 peak)
			ld_Mz = ld_PeptideMz - NEUTRON / li_Charge;
			ls_Key = QString("%1-%2-forbidden").arg(ls_Peptide).arg(li_Charge);
			lk_TempList.push_back(tk_DoubleStringPair(ld_Mz, ls_Key));
		}
	}

	qSort(lk_TempList.begin(), lk_TempList.end(), sortByMz);

	mk_AllTargetMasses = QList<double>();
	mk_TargetMzIndex = QHash<QString, int>();

	foreach (tk_DoubleStringPair lk_Pair, lk_TempList)
	{
		mk_TargetMzIndex[lk_Pair.second] = mk_AllTargetMasses.size();
		mk_AllTargetMasses.push_back(lk_Pair.first);
	}
	
	if (mk_YamlOutStream.device())
		mk_YamlOutStream << "results:" << endl;
		
	// parse all bands
	foreach (QString ls_Path, ak_SpectraFiles)
	{
		mk_ScanFailureReason.clear();
		mk_ScanHash = QHash<QString, r_Scan>();
		ms_CurrentSpot = QFileInfo(ls_Path).baseName();
		mk_SpotResults = tk_SpotResults();
		// parse spot
		this->parseFile(ls_Path);
		
/*		
		foreach(QString ls_Peptide, mk_Peptides)
		{
			QFile lk_File("elution profile/" + ls_Peptide + ".csv");
			lk_File.open(QIODevice::WriteOnly);
			QTextStream lk_Stream(&lk_File);
			lk_Stream << "retention time;unlabeled;labeled;ratio" << endl;
			for (int i = 0; i < mk_ElutionProfile[ls_Peptide][0].size(); ++i)
			{
				lk_Stream << mk_ElutionProfile[ls_Peptide][0][i] << ";"
					<< mk_ElutionProfile[ls_Peptide][1][i] << ";"
					<< mk_ElutionProfile[ls_Peptide][2][i] << ";";
				if (mk_ElutionProfile[ls_Peptide][1][i] < 0.000001 || mk_ElutionProfile[ls_Peptide][2][i] < 0.000001)
					lk_Stream << "0.0;";
				else
					lk_Stream << mk_ElutionProfile[ls_Peptide][1][i] / mk_ElutionProfile[ls_Peptide][2][i] << ";";
				lk_Stream << endl;
			}
			lk_Stream.flush();
			lk_File.close();
		}
*/
		
		printf(" done.\n");
		
		if (mk_YamlOutStream.device())
			mk_YamlOutStream << "  \"" << ms_CurrentSpot << "\":" << endl;
			
		foreach (QString ls_Peptide, mk_SpotResults.keys())
		{
			// put all scan quantitation results in a list
			QList<r_ScanQuantitationResult> lk_ScanResults;
			foreach (r_ScanQuantitationResult lr_ScanQuantiationResult, mk_SpotResults[ls_Peptide])
				lk_ScanResults.push_back(lr_ScanQuantiationResult);
				
			// ...and sort them by retention time
			qSort(lk_ScanResults.begin(), lk_ScanResults.end(), compareScanQuantitationByRetentionTime);
			
/*			
			// find the scan quantitation result with the best SNR
			int li_BestIndex = 0;
			for (int i = 1; i < lk_ScanResults.size(); ++i)
				if (lk_ScanResults[i].md_Snr > lk_ScanResults[li_BestIndex].md_Snr)
					li_BestIndex = i;
					
			// now starting from the best scan quantitation result, and expand in both
			// directions while the time difference is less or equal to the maximum
			// allowed time difference
			int li_LeftIndex = li_BestIndex;
			int li_RightIndex = li_BestIndex;
			
			while (li_LeftIndex > 0 && fabs(lk_ScanResults[li_LeftIndex].md_RetentionTime - lk_ScanResults[li_LeftIndex - 1].md_RetentionTime) < md_MaxTimeDifference / 60.0)
				--li_LeftIndex;
			while (li_RightIndex < lk_ScanResults.size() - 1 && fabs(lk_ScanResults[li_RightIndex].md_RetentionTime - lk_ScanResults[li_RightIndex + 1].md_RetentionTime) < md_MaxTimeDifference / 60.0)
				++li_RightIndex;
				
			// keep all results that are good regarding the time difference
			lk_ScanResults = lk_ScanResults.mid(li_LeftIndex, li_RightIndex - li_LeftIndex + 1);
*/			

/*
			// ATTENTION ATTENTION
			// This crop upper by snr block has been removed because I think it is 
			// crap and that it's better to use as many scans as you can get.
			// This saves us one voodoo parameter and if you don't like some SNR
			// change the SNR threshold.
			
			// ...sort remaining results by snr to crop the upper 25% (user adjustable value)
			qSort(lk_ScanResults.begin(), lk_ScanResults.end(), compareScanQuantitationBySnr);
			double ld_SnrCutOff = lk_ScanResults.first().md_Snr * (1.0 - md_CropUpper);
			QList<r_ScanQuantitationResult>::iterator lk_Iter = lk_ScanResults.begin();
			
			// always include the best scan into the merged result
			++lk_Iter;
			
			// include remaining scans
			while (lk_Iter != lk_ScanResults.end() && lk_Iter->md_Snr >= ld_SnrCutOff)
				++lk_Iter;
				
			// erase remaining scans
			lk_ScanResults.erase(lk_Iter, lk_ScanResults.end());
*/			
			
			if (mk_YamlOutStream.device())
			{
 				mk_YamlOutStream << "      " << ls_Peptide << ":" << endl;
				foreach (r_ScanQuantitationResult lr_QuantitationResult, lk_ScanResults)
				{
					r_Scan& lr_Scan = mk_ScanHash[lr_QuantitationResult.ms_ScanHashKey];
					mk_YamlOutStream << "        - { id: \"" << lr_Scan.ms_Id << "\""
						<< ", retentionTime: " << lr_Scan.md_RetentionTime 
						<< ", filterLine: \"" << lr_Scan.ms_FilterLine << "\""
						<< ", svg: \"" << lr_QuantitationResult.ms_ScanHashKey << "\""
						<< ", snr: " << lr_QuantitationResult.md_Snr
						<< ", amountUnlabeled: " << lr_QuantitationResult.md_AmountUnlabeled
						<< ", amountLabeled: " << lr_QuantitationResult.md_AmountLabeled
						<< ", ratio: " << lr_QuantitationResult.md_Ratio 
						<< ", charge: " << lr_QuantitationResult.mi_Charge
						<< " }" << endl;
				}
			}
		}
		if (mb_PrintStatistics)
		{
			// print failure evaluation
			QHash<r_QuantitationFailureReason::Enumeration, int> lk_Reasons;
			foreach (QString ls_Key, mk_ScanFailureReason.keys())
			{
				r_QuantitationFailureReason::Enumeration le_Reason = mk_ScanFailureReason[ls_Key];
				if (!lk_Reasons.contains(le_Reason))
					lk_Reasons[le_Reason] = 0;
				++lk_Reasons[le_Reason];
			}
			for (int i = 0; i < r_QuantitationFailureReason::Size; ++i)
			{
				QString ls_Description;
				r_QuantitationFailureReason::Enumeration le_Reason = (r_QuantitationFailureReason::Enumeration)i;
				switch (le_Reason)
				{
					case r_QuantitationFailureReason::NoMatchedTargetMass:
						ls_Description = "No target mass could be matched within the specified mass accuracy.";
						break;
					case r_QuantitationFailureReason::IsotopePeaksMissing:
						ls_Description = "One or more isotope peaks were missing.";
						break;
					case r_QuantitationFailureReason::UnlabeledIsotopePeaksNotDescending:
						ls_Description = "The unlabeled isotope peaks were not descending.";
						break;
					case r_QuantitationFailureReason::ForbiddenPeakPresent:
						ls_Description = "A forbidden peak (left from the unlabeled ion series) was present.";
						break;
					case r_QuantitationFailureReason::LowSnr:
						ls_Description = "The signal to noise ratio was too low.";
						break;
					case r_QuantitationFailureReason::Success:
						ls_Description = "A peptide was successfully quantified.";
						break;
				}
				int li_Count = 0;
				if (lk_Reasons.contains(le_Reason))
					li_Count = lk_Reasons[le_Reason];
				printf("%6d %s\n", li_Count, ls_Description.toStdString().c_str());
			}
		}
	}
}


void k_Quantifier::handleScan(r_Scan& ar_Scan)
{
	if (ar_Scan.mr_Spectrum.mi_PeaksCount == 0)
	{
		printf("Warning: Empty spectrum (scan #%s @ %1.2f minutes)!\n", ar_Scan.ms_Id.toStdString().c_str(), ar_Scan.md_RetentionTime);
		return;
	}
	
	updateFailureReason(ar_Scan.ms_Id, r_QuantitationFailureReason::NoMatchedTargetMass);
	
	QHash<QString, int> lk_TargetPeakStartIndex;
	QHash<QString, int> lk_TargetPeakEndIndex;
	
	// determine each peptide's elution profile contribution in this scan
	
	/*
	// start with finding start/end pairs of peaks that correspond to each target m/z value
	for (int i = 0; i < ar_Scan.mr_Spectrum.mi_PeaksCount; ++i)
	{
		double ld_Mz = ar_Scan.mr_Spectrum.md_MzValues_[i];
		foreach (QString ls_Peptide, mk_Peptides)
		{
			for (int li_Labeled = 0; li_Labeled < 2; ++li_Labeled)
			{
				for (int li_Charge = mi_MinCharge; li_Charge <= mi_MaxCharge; ++li_Charge)
				{
					QString ls_Key = QString("%1-%2-%3-%4").arg(ls_Peptide).arg(li_Charge).arg(li_Labeled == 0 ? "unlabeled": "labeled").arg(0);
					if (ld_Mz < mk_AllTargetMasses[mk_TargetMzIndex[ls_Key]] - md_ElutionProfilePeakWidth)
						lk_TargetPeakStartIndex[ls_Key] = i;
						
					if (ld_Mz > mk_AllTargetMasses[mk_TargetMzIndex[ls_Key]] + md_ElutionProfilePeakWidth)
					{
						if (!lk_TargetPeakEndIndex.contains(ls_Key))
							lk_TargetPeakEndIndex[ls_Key] = i;
					}
				}
			}
		}
	}
	
	foreach (QString ls_Peptide, mk_Peptides)
	{
		for (int li_Labeled = 0; li_Labeled < 2; ++li_Labeled)
		{
			for (int li_Charge = mi_MinCharge; li_Charge <= mi_MaxCharge; ++li_Charge)
			{
				QString ls_Key = QString("%1-%2-%3-%4").arg(ls_Peptide).arg(li_Charge).arg(li_Labeled == 0 ? "unlabeled": "labeled").arg(0);
				// we must have a start and and end point,
				// and the start and end indicies must be different
				if (lk_TargetPeakStartIndex.contains(ls_Key) && lk_TargetPeakEndIndex.contains(ls_Key) && lk_TargetPeakStartIndex[ls_Key] < lk_TargetPeakEndIndex[ls_Key])
				{
					// determine product of signal and target m/z filter triangle
					// TODO: continue here!!
				}
			}
		}
	}
	*/
	
	// find all peaks in this spectrum
	QList<r_Peak> lk_AllPeaks = this->findAllPeaks(ar_Scan.mr_Spectrum);
	
	// we need at least a few peaks
	if (lk_AllPeaks.size() < mi_WatchIsotopesCount * 2)
		return;
	
	// match all target m/z values simultaneously to this spectrum's peaks
	// create root bucket
	r_Bucket lr_RootBucket(0, lk_AllPeaks.size());
	for (int i = 0; i < mk_AllTargetMasses.size(); ++i)
		lr_RootBucket.mk_Entries.push_back(i);
	
	QList<r_Bucket> lk_Buckets;
	lk_Buckets.push_back(lr_RootBucket);
	
	// after parallel searching, this list will contain for 
	// each target m/z value the peak which is closest to it
 	QHash<int, int> lk_PeakForTargetMz;
	
	// repeat until no more buckets are left to be processed
	while (!lk_Buckets.empty())
	{
		// put all target m/z values into the appropriate bucket
		QList<r_Bucket> lk_NewBuckets;
		foreach (r_Bucket lr_Bucket, lk_Buckets)
		{
			if (lr_Bucket.mi_Length <= 1)
			{
				for (int i = 0; i < lr_Bucket.mk_Entries.size(); ++i)
				{
					/*
					printf("%1.6f is close to %1.6f\n", 
						mk_AllTargetMasses[lr_Bucket.mk_Entries[i]],
						ar_Scan.mr_Spectrum.md_MzValues_[lr_Bucket.mi_Start]);
					*/
					lk_PeakForTargetMz[lr_Bucket.mk_Entries[i]] = lr_Bucket.mi_Start;
				}
			}
			else
			{
				// split this bucket and create left and right children
				int li_HalfSize = lr_Bucket.mi_Length / 2;
				r_Bucket lr_LeftChild(lr_Bucket.mi_Start, li_HalfSize);
				r_Bucket lr_RightChild(lr_Bucket.mi_Start + li_HalfSize, lr_Bucket.mi_Length - li_HalfSize);
				double ld_LeftBorder = lk_AllPeaks[lr_RightChild.mi_Start - 1].md_PeakMz;
				double ld_RightBorder = lk_AllPeaks[lr_RightChild.mi_Start].md_PeakMz;
				double ld_Razor = (ld_LeftBorder + ld_RightBorder) * 0.5;
				
				// now determine which target m/z entries go into the 
				// left child and which go into the right child
				
				// TODO: this could be sped up because everything is sorted
				for (int i = 0; i < lr_Bucket.mk_Entries.size(); ++i)
				{
					if (mk_AllTargetMasses[lr_Bucket.mk_Entries[i]] > ld_Razor)
						lr_RightChild.mk_Entries.push_back(lr_Bucket.mk_Entries[i]);
					else
						lr_LeftChild.mk_Entries.push_back(lr_Bucket.mk_Entries[i]);
				}
				
				if (!lr_LeftChild.mk_Entries.empty())
					lk_NewBuckets.push_back(lr_LeftChild);
				if (!lr_RightChild.mk_Entries.empty())
					lk_NewBuckets.push_back(lr_RightChild);
			}
		}
		lk_Buckets = lk_NewBuckets;
	}
	
	// discard all target m/z matches that are bogus because
	// they are not within the specified mass accuracy
	QHash<int, int> lk_GoodPeakForTargetMz;
	foreach (int li_TargetMzIndex, lk_PeakForTargetMz.keys())
	{
		double ld_TargetMz = mk_AllTargetMasses[li_TargetMzIndex];
		r_Peak lr_Peak = lk_AllPeaks[lk_PeakForTargetMz[li_TargetMzIndex]];
		double ld_PeakMz = lr_Peak.md_PeakMz;
		double ld_MzError = fabs(ld_TargetMz - ld_PeakMz);
		double ld_MaxMzError = ld_TargetMz * md_MassAccuracy / 1000000.0;
		if (ld_MzError <= ld_MaxMzError)
		{
			lk_GoodPeakForTargetMz[li_TargetMzIndex] = lk_PeakForTargetMz[li_TargetMzIndex];
			updateFailureReason(ar_Scan.ms_Id, r_QuantitationFailureReason::IsotopePeaksMissing);
		}
	}
	
	lk_PeakForTargetMz = lk_GoodPeakForTargetMz;
	
	// now check for each peptide/charge whether all peaks are there
	foreach (QString ls_Peptide, mk_Peptides)
	{
		for (int li_Charge = mi_MinCharge; li_Charge <= mi_MaxCharge; ++li_Charge)
		{
			QList<r_Peak> lk_Peaks;
			QList<double> lk_TargetMz;
			r_Peak* lr_ForbiddenPeak_ = NULL;
			for (int li_Labeled = 0; li_Labeled < 2; ++li_Labeled)
			{
				for (int li_Isotope = 0; li_Isotope < mi_WatchIsotopesCount; ++li_Isotope)
				{
					QString ls_Key = QString("%1-%2-%3-%4").arg(ls_Peptide).
						arg(li_Charge).arg(li_Labeled == 0 ? "unlabeled" : "labeled").arg(li_Isotope);
					int li_TargetMzIndex = mk_TargetMzIndex[ls_Key];
					if (lk_PeakForTargetMz.contains(li_TargetMzIndex))
					{
						lk_Peaks.push_back(lk_AllPeaks[lk_PeakForTargetMz[li_TargetMzIndex]]);
						lk_TargetMz.push_back(mk_AllTargetMasses[li_TargetMzIndex]);
					}
				}
			}
			QString ls_Key = QString("%1-%2-forbidden").arg(ls_Peptide).arg(li_Charge);
			int li_TargetMzIndex = mk_TargetMzIndex[ls_Key];
			if (lk_PeakForTargetMz.contains(li_TargetMzIndex))
				lr_ForbiddenPeak_ = &lk_AllPeaks[lk_PeakForTargetMz[li_TargetMzIndex]];

			if (lk_Peaks.size() == mi_WatchIsotopesCount * 2)
			{
				// yay! we found a peptide with a certain charge state!
				updateFailureReason(ar_Scan.ms_Id, r_QuantitationFailureReason::UnlabeledIsotopePeaksNotDescending);
				r_ScanQuantitationResult lr_ScanResult = this->checkResult(lk_Peaks, lr_ForbiddenPeak_, ar_Scan, ls_Peptide, li_Charge, lk_TargetMz);
				if (lr_ScanResult.mb_IsGood)
				{
					lr_ScanResult.ms_ScanHashKey = QString("%1.%2.%3.%4").arg(ms_CurrentSpot).arg(ar_Scan.ms_Id).arg(ls_Peptide).arg(li_Charge);
					lr_ScanResult.md_RetentionTime = ar_Scan.md_RetentionTime;
					lr_ScanResult.mi_Charge = li_Charge;
					lr_ScanResult.ms_ScanId = ar_Scan.ms_Id;
					mk_ScanHash.insert(lr_ScanResult.ms_ScanHashKey, r_Scan(ar_Scan));
					if (!mk_SpotResults.contains(ls_Peptide))
						mk_SpotResults[ls_Peptide] = QList<r_ScanQuantitationResult>();
					mk_SpotResults[ls_Peptide].push_back(lr_ScanResult);
					
					if (!ms_SvgOutPath.isEmpty())
					{
						QFile lk_File(QString("%1/%2.svg").arg(ms_SvgOutPath).arg(lr_ScanResult.ms_ScanHashKey));
						lk_File.open(QIODevice::WriteOnly);
						QTextStream lk_Stream(&lk_File);
						lk_Stream << this->renderScanAsSvg(ar_Scan, lr_ScanResult);
						lk_Stream.flush();
						lk_File.close();
					}
				}
			}
		}
	}
	
/*	
	foreach (QString ls_Peptide, mk_Peptides)
	{
		for (int li_Charge = mi_MinCharge; li_Charge <= mi_MaxCharge; ++li_Charge)
		{
			// check this scan for peptide with charge
			QPair<QString, int> lk_PeptideAndCharge(ls_Peptide, li_Charge);
			r_ScanQuantitationResult lr_Result = this->searchPeptide(ar_Scan.mr_Spectrum, mk_ChargedPeptideMasses[lk_PeptideAndCharge], li_Charge, mk_ChargedPeptideModificationMasses[lk_PeptideAndCharge]);
			if (lr_Result.mb_IsGood)
			{
				lr_Result.ms_ScanHashKey = QString("%1.%2.%3.%4").arg(ms_CurrentSpot).arg(ar_Scan.ms_Id).arg(ls_Peptide).arg(li_Charge);
				lr_Result.md_RetentionTime = ar_Scan.md_RetentionTime;
				lr_Result.mi_Charge = li_Charge;
				lr_Result.ms_ScanId = ar_Scan.ms_Id;
				mk_ScanHash.insert(lr_Result.ms_ScanHashKey, r_Scan(ar_Scan));
				
				if (!mk_SpotResults.contains(ls_Peptide))
					mk_SpotResults[ls_Peptide] = QList<r_ScanQuantitationResult>();
				if (!ms_SvgOutPath.isEmpty())
				{
					QFile lk_File(QString("%1/%2.svg").arg(ms_SvgOutPath).arg(lr_Result.ms_ScanHashKey));
					lk_File.open(QIODevice::WriteOnly);
					QTextStream lk_Stream(&lk_File);
					lk_Stream << this->renderScanAsSvg(ar_Scan, lr_Result);
					lk_Stream.flush();
					lk_File.close();
				}
				mk_SpotResults[ls_Peptide].push_back(lr_Result);
			}
		}
	}
	*/
}


void k_Quantifier::progressFunction(QString as_ScanId, bool)
{
	printf("\r%s: scan #%s...", ms_CurrentSpot.toStdString().c_str(), as_ScanId.toStdString().c_str());
}


double k_Quantifier::calculatePeptideMass(QString as_Peptide, int ai_Charge)
{
	double ld_Mass = WATER;
	for (int i = 0; i < as_Peptide.length(); ++i)
		ld_Mass += mk_AminoAcidWeight[as_Peptide.at(i).toAscii()];
	// UNSURE ABOUT THIS VALUE BUT OK WITH BIANCAS EXAMPLES, should be 1.0078250
	return (ld_Mass + HYDROGEN * ai_Charge) / ai_Charge;
}


r_ScanQuantitationResult k_Quantifier::searchPeptide(r_Spectrum& ar_Spectrum, double ad_PeptideMz, int ai_Charge, double ad_ModificationMass)
{
	r_ScanQuantitationResult lr_Result;
	lr_Result.mb_IsGood = false;
	/*
	
	// find the m/z value that is closest to the target peptide m/z value
	// TODO: maybe this search could be sped up with a binary search variation
	
	// determine unlabeled and labeled target masses
	QList<double> lk_TargetMz;
	for (int i = 0; i < mi_WatchIsotopesCount; ++i)
		lk_TargetMz.push_back(ad_PeptideMz + i * NEUTRON / ai_Charge);
	for (int i = 0; i < mi_WatchIsotopesCount; ++i)
		lk_TargetMz.push_back(ad_PeptideMz + i * NEUTRON / ai_Charge + ad_ModificationMass / ai_Charge);
		
	// find peaks
	QList<r_Peak> lk_Peaks = this->findPeaks(ar_Spectrum, lk_TargetMz);
	
	// check whether peaks are too much off-center
	for (int i = 0; i < mi_WatchIsotopesCount * 2; ++i)
	{
		double ld_LeftMz = ar_Spectrum.md_MzValues_[lk_Peaks[i].mi_LeftBorderIndex];
		double ld_RightMz = ar_Spectrum.md_MzValues_[lk_Peaks[i].mi_RightBorderIndex];
		double ld_PeakCenter = (ld_LeftMz + ld_RightMz) * 0.5;
		double ld_PeakWidth2 = fabs(ld_RightMz - ld_LeftMz) * 0.5;
		double ld_OffCenter = fabs(lk_TargetMz[i] - ld_PeakCenter) / ld_PeakWidth2;
		if (ld_OffCenter > md_MaxOffCenter)
			return lr_Result;
	}
	
	// check whether peaks are non-overlapping
	for (int i = 0; i < mi_WatchIsotopesCount * 2 - 1; ++i)
	{
		if (lk_Peaks[i].mi_RightBorderIndex >= lk_Peaks[i + 1].mi_LeftBorderIndex)
			return lr_Result;
	}
	
	// check whether unlabeled peaks are descending
	for (int i = 1; i < mi_WatchIsotopesCount; ++i)
	{
		if (lk_Peaks[i - 1].md_PeakIntensity <= lk_Peaks[i].md_PeakIntensity)
			return lr_Result;
	}
	
	QList<double> lk_PeakSnr;
	// determine SNR for each peak
	for (int i = 0; i < mi_WatchIsotopesCount * 2; ++i)
		lk_PeakSnr << lk_Peaks[i].md_PeakIntensity / lk_Peaks[i].md_OutsideBorderMaxIntensity;
		
	// determine gap maximum
	double ld_GapMaximum = 0.0;
	for (int i = 0; i < mi_WatchIsotopesCount * 2; ++i)
		ld_GapMaximum = std::max<double>(ld_GapMaximum, lk_Peaks[i].md_OutsideBorderMaxIntensity);
		
	// check one isotope peak to the left
	QList<double> lk_LeftTargetMz;
	// unlabeled ion
	lk_LeftTargetMz.push_back(ad_PeptideMz - NEUTRON / ai_Charge);
	// labeled ion
	//lk_LeftTargetMz.push_back(ad_PeptideMz - NEUTRON / ai_Charge + ad_ModificationMass / ai_Charge);
	
	QList<r_Peak> lk_LeftPeaks = this->findPeaks(ar_Spectrum, lk_LeftTargetMz);
	
	// make sure that there is no left isotope peak 
	foreach (r_Peak lr_Peak, lk_LeftPeaks)
	{
		if (lr_Peak.md_PeakIntensity > ld_GapMaximum)
			return lr_Result;
	}
		
	double ld_MinSnr = lk_PeakSnr[0];
	for (int i = 1; i < mi_WatchIsotopesCount * 2; ++i)
		ld_MinSnr = std::min<double>(ld_MinSnr, lk_PeakSnr[i]);
		
	if (ld_MinSnr < md_MinSnr)
		return lr_Result;
		
	lr_Result.md_Snr = ld_MinSnr;
	
	double ld_UnlabeledSum = 0.0;
	double ld_LabeledSum = 0.0;
	for (int i = 0; i < mi_WatchIsotopesCount; ++i)
	{
		ld_UnlabeledSum += lk_Peaks[i].md_PeakIntensity;
		ld_LabeledSum += lk_Peaks[i + mi_WatchIsotopesCount].md_PeakIntensity;
	}
	lr_Result.md_Ratio = ld_UnlabeledSum / ld_LabeledSum;
	
	for (int i = 0; i < mi_WatchIsotopesCount * 2; ++i)
		lr_Result.mk_Highlight.push_back(tk_IntPair(lk_Peaks[i].mi_LeftBorderIndex, lk_Peaks[i].mi_RightBorderIndex));
		
	for (int i = 0; i < mi_WatchIsotopesCount * 2; ++i)
		lr_Result.mk_TargetMz.push_back(lk_TargetMz[i]);
		
	lr_Result.md_MinMz = lk_Peaks.first().md_PeakMz - 2.0 / ai_Charge - 0.5;
	lr_Result.md_MaxMz = lk_Peaks.last().md_PeakMz + 0.5;
	
	lr_Result.mb_IsGood = true;
	lr_Result.mk_Peaks = lk_Peaks;
	*/
	
	return lr_Result;
}


QList<r_Peak> k_Quantifier::findPeaks(r_Spectrum& ar_Spectrum, QList<double> ak_TargetMasses)
{
	QList<int> lk_BestMatch;
	QList<double> lk_Error;
	for (int i = 0; i < ak_TargetMasses.size(); ++i)
	{
		lk_BestMatch.push_back(0);
		lk_Error.push_back(fabs(ar_Spectrum.md_MzValues_[0] - ak_TargetMasses[i]));
	}
	
	for (int i = 1; i < ar_Spectrum.mi_PeaksCount; ++i)
	{
		for (int k = 0; k < ak_TargetMasses.size(); ++k)
		{
			double ld_Error = fabs(ar_Spectrum.md_MzValues_[i] - ak_TargetMasses[k]);
			if (ld_Error < lk_Error[k])
			{
				lk_BestMatch[k] = i;
				lk_Error[k] = ld_Error;
			}
		}
	}
	
	QList<r_Peak> lk_Result;
	for (int i = 0; i < ak_TargetMasses.size(); ++i)
		lk_Result.push_back(this->findPeak(ar_Spectrum, lk_BestMatch[i]));
	
	return lk_Result;
}


QList<r_Peak> k_Quantifier::findAllPeaks(r_Spectrum& ar_Spectrum)
{
	QList<r_Peak> lk_Results;
	double ld_LastIntensity = ar_Spectrum.md_IntensityValues_[0];
	int li_LastDirection = -100;
	int li_PeakIndex = -1;
	int li_ValleyIndex = -1;
	for (int i = 1; i < ar_Spectrum.mi_PeaksCount; ++i)
	{
		double ld_ThisIntensity = ar_Spectrum.md_IntensityValues_[i];
		double ld_Slope = ld_ThisIntensity - ld_LastIntensity;
		int li_Direction = 0;
		if (ld_Slope > 0)
			li_Direction = 1;
		else if (ld_Slope < 0)
			li_Direction = -1;
			
		if (li_LastDirection >= -1)
		{
			if (li_Direction == -1 && li_LastDirection != -1)
			{
				// from ascend or equal to descend (peak!)
				li_PeakIndex = i - 1;
			}
			if (li_Direction != -1 && li_LastDirection == -1)
			{
				// end of peak
				if (li_PeakIndex >= 0 && li_ValleyIndex >= 0)
				{
					// we now have a good peak
					r_Peak lr_Peak;
					lr_Peak.mi_PeakIndex = li_PeakIndex;
					lr_Peak.mi_LeftBorderIndex = li_ValleyIndex;
					lr_Peak.mi_RightBorderIndex = i - 1;
					lr_Peak.md_OutsideBorderMaxIntensity = 1e-15;
					if (lr_Peak.mi_LeftBorderIndex > 0)
						lr_Peak.md_OutsideBorderMaxIntensity = std::max<double>(lr_Peak.md_OutsideBorderMaxIntensity, ar_Spectrum.md_IntensityValues_[lr_Peak.mi_LeftBorderIndex - 1]);
					if (lr_Peak.mi_RightBorderIndex < ar_Spectrum.mi_PeaksCount - 1)
						lr_Peak.md_OutsideBorderMaxIntensity = std::max<double>(lr_Peak.md_OutsideBorderMaxIntensity, ar_Spectrum.md_IntensityValues_[lr_Peak.mi_RightBorderIndex + 1]);
					lr_Peak.md_Snr = lr_Peak.md_PeakIntensity / lr_Peak.md_OutsideBorderMaxIntensity;
					fitGaussian(&lr_Peak.md_GaussA, &lr_Peak.md_GaussB, &lr_Peak.md_GaussC,
								ar_Spectrum.md_MzValues_[lr_Peak.mi_PeakIndex - 1],
								ar_Spectrum.md_IntensityValues_[lr_Peak.mi_PeakIndex - 1],
								ar_Spectrum.md_MzValues_[lr_Peak.mi_PeakIndex],
								ar_Spectrum.md_IntensityValues_[lr_Peak.mi_PeakIndex],
								ar_Spectrum.md_MzValues_[lr_Peak.mi_PeakIndex + 1],
								ar_Spectrum.md_IntensityValues_[lr_Peak.mi_PeakIndex + 1]);
					lr_Peak.md_PeakMz = lr_Peak.md_GaussB;
					lr_Peak.md_PeakIntensity = lr_Peak.md_GaussA;
					lr_Peak.md_PeakArea = lr_Peak.md_GaussA * lr_Peak.md_GaussC;
					lk_Results.push_back(lr_Peak);
				}
				li_ValleyIndex = i - 1;
			}
			if (li_Direction == 1 && li_LastDirection != 1)
			{
				// start of peak
				li_ValleyIndex = i - 1;
			}
		}
		li_LastDirection = li_Direction;
		ld_LastIntensity = ld_ThisIntensity;
	}
	
	return lk_Results;
}


r_Peak k_Quantifier::findPeak(r_Spectrum& ar_Spectrum, int ai_Index)
{
	r_Peak lr_Peak;
	
	lr_Peak.mi_PeakIndex = ai_Index;
	
	// ascend to the top...
	while (true)
	{
		double ld_Intensity = ar_Spectrum.md_IntensityValues_[lr_Peak.mi_PeakIndex];
		if ((lr_Peak.mi_PeakIndex > 0) && (ar_Spectrum.md_IntensityValues_[lr_Peak.mi_PeakIndex - 1] > ld_Intensity))
		{
			--lr_Peak.mi_PeakIndex;
		}
		else if ((lr_Peak.mi_PeakIndex < ar_Spectrum.mi_PeaksCount - 1) && (ar_Spectrum.md_IntensityValues_[lr_Peak.mi_PeakIndex + 1] > ld_Intensity))
		{
			++lr_Peak.mi_PeakIndex;
		}
		else
			break;
	}
	// mi_PeakIndex is at a local maximum now
	
	// descend to the left...
	lr_Peak.mi_LeftBorderIndex = lr_Peak.mi_PeakIndex;
	while ((lr_Peak.mi_LeftBorderIndex > 0) && (ar_Spectrum.md_IntensityValues_[lr_Peak.mi_LeftBorderIndex - 1] < ar_Spectrum.md_IntensityValues_[lr_Peak.mi_LeftBorderIndex]))
		--lr_Peak.mi_LeftBorderIndex;
		
	// descent to the right...
	lr_Peak.mi_RightBorderIndex = lr_Peak.mi_PeakIndex;
	while ((lr_Peak.mi_RightBorderIndex < ar_Spectrum.mi_PeaksCount - 1) && (ar_Spectrum.md_IntensityValues_[lr_Peak.mi_RightBorderIndex + 1] < ar_Spectrum.md_IntensityValues_[lr_Peak.mi_RightBorderIndex]))
		++lr_Peak.mi_RightBorderIndex;
		
	lr_Peak.md_PeakIntensity = ar_Spectrum.md_IntensityValues_[lr_Peak.mi_PeakIndex];
	lr_Peak.md_PeakMz = ar_Spectrum.md_MzValues_[lr_Peak.mi_PeakIndex];
	
	lr_Peak.md_OutsideBorderMaxIntensity = 1e-15;
	if (lr_Peak.mi_LeftBorderIndex > 0)
		lr_Peak.md_OutsideBorderMaxIntensity = 
			std::max<double>(lr_Peak.md_OutsideBorderMaxIntensity,
							 ar_Spectrum.md_IntensityValues_[lr_Peak.mi_LeftBorderIndex - 1]);
	if (lr_Peak.mi_RightBorderIndex < ar_Spectrum.mi_PeaksCount - 1)
		lr_Peak.md_OutsideBorderMaxIntensity = 
			std::max<double>(lr_Peak.md_OutsideBorderMaxIntensity,
							 ar_Spectrum.md_IntensityValues_[lr_Peak.mi_RightBorderIndex + 1]);
							 
	lr_Peak.md_Snr = lr_Peak.md_PeakIntensity / lr_Peak.md_OutsideBorderMaxIntensity;
	
	return lr_Peak;
}


QString k_Quantifier::renderScanAsSvg(r_Scan& ar_Scan, r_ScanQuantitationResult ar_QuantitationResult)
{
	double ld_Ratio = 4.0;
	double ld_Width = 1024.0;
	double ld_Height = ld_Width / ld_Ratio;
	double ld_Border = 4.0;
	
	if (ar_QuantitationResult.md_MinMz == 0.0)
		ar_QuantitationResult.md_MinMz = ar_Scan.mr_Spectrum.md_MzValues_[0];
	if (ar_QuantitationResult.md_MaxMz == 0.0)
		ar_QuantitationResult.md_MaxMz = ar_Scan.mr_Spectrum.md_MzValues_[ar_Scan.mr_Spectrum.mi_PeaksCount - 1];
		
	int li_Start = 0;
	int li_End = ar_Scan.mr_Spectrum.mi_PeaksCount - 1;
	
	// adjust start and end values
	while (li_Start < ar_Scan.mr_Spectrum.mi_PeaksCount - 1 && ar_Scan.mr_Spectrum.md_MzValues_[li_Start] < ar_QuantitationResult.md_MinMz)
		++li_Start;
	while (li_End > 0 && ar_Scan.mr_Spectrum.md_MzValues_[li_End] > ar_QuantitationResult.md_MaxMz)
		--li_End;
		
	if (li_Start > li_End)
	{
		// oops! start and end got crossed.
		li_Start = 0;
		li_End = 0;
	}
	
	// determine x and y extensions
	double xmin = ar_Scan.mr_Spectrum.md_MzValues_[li_Start];
	double xmax = ar_Scan.mr_Spectrum.md_MzValues_[li_End];
	xmin = ar_QuantitationResult.md_MinMz;
	xmax = ar_QuantitationResult.md_MaxMz;
	double ymax = ar_Scan.mr_Spectrum.md_IntensityValues_[li_Start];
	for (int i = li_Start + 1; i <= li_End; ++i)
		ymax = std::max<double>(ymax, ar_Scan.mr_Spectrum.md_IntensityValues_[i]);
	// ymax remains unscaled!
	double ymaxScaled = this->scale(1.0);
	
	QBuffer lk_Buffer;
	lk_Buffer.open(QBuffer::WriteOnly);
	
	double ld_OnePercent = 0.0;
	double ld_TenPercent = 0.0;
	double ld_OneHundredPercent = 0.0;
	
	double x0 = ld_Border;
	double y0 = ld_Height - ld_Border;
	double dx = (ld_Width - 2.0 * ld_Border) / (xmax - xmin);
	double dy = -(ld_Height - 2.0 * ld_Border);
	
	{
		QSvgGenerator lk_Generator;
		lk_Generator.setSize(QSize((int)ld_Width, (int)ld_Height));
		lk_Generator.setOutputDevice(&lk_Buffer);
		
		QPainter lk_Painter(&lk_Generator);
		
		// fill background
		lk_Painter.fillRect(QRectF(0.0, 0.0, ld_Width, ld_Height), QBrush(Qt::white));
		
		// draw frame
		QPen lk_Pen;
		lk_Pen.setWidthF(1.0);
		lk_Pen.setJoinStyle(Qt::BevelJoin);
		lk_Pen.setColor(QColor(192, 192, 192));
		lk_Painter.setPen(lk_Pen);
		
		lk_Painter.drawLine(QPointF(0.0, y0), QPointF(ld_Width, y0));
		//lk_Painter.drawLine(QPointF(0.0, y0), QPointF(ld_Width, y0));
		//lk_Painter.drawLine(QPointF(x0, y0), QPointF(x0, y0 + (ymax - ymin) * dy));
		
		// draw spectrum
		lk_Pen.setWidthF(1.0);
		lk_Pen.setJoinStyle(Qt::RoundJoin);
		lk_Pen.setColor(QColor(192, 192, 192));
		lk_Pen.setStyle(Qt::DashLine);
		lk_Painter.setPen(lk_Pen);
		
		// draw target m/z lines
		foreach (double ld_Mz, ar_QuantitationResult.mk_TargetMz)
		{
			double x = (ld_Mz - xmin) * dx + x0;
			lk_Painter.drawLine(QPointF(x, y0), QPointF(x, y0 + dy));
		}
		
		double y;
		y = this->scale(1.0 / 100.0) / ymaxScaled;
		ld_OnePercent = y * dy + y0;
		lk_Painter.drawLine(QPointF(x0, y * dy + y0), QPointF(x0 + dx * (xmax - xmin), y * dy + y0));
		y = this->scale(10.0 / 100.0) / ymaxScaled;
		ld_TenPercent = y * dy + y0;
		lk_Painter.drawLine(QPointF(x0, y * dy + y0), QPointF(x0 + dx * (xmax - xmin), y * dy + y0));
		y = this->scale(100.0 / 100.0) / ymaxScaled;
		ld_OneHundredPercent = y * dy + y0;
		lk_Painter.drawLine(QPointF(x0, y * dy + y0), QPointF(x0 + dx * (xmax - xmin), y * dy + y0));
		
		lk_Pen.setWidthF(1.0);
		lk_Pen.setJoinStyle(Qt::RoundJoin);
		lk_Pen.setColor(QColor(128, 128, 128));
		lk_Pen.setStyle(Qt::SolidLine);
		lk_Painter.setPen(lk_Pen);
		
		QVector<QPointF> lk_Points;
		for (int i = li_Start; i <= li_End; ++i)
			lk_Points.append(QPointF((ar_Scan.mr_Spectrum.md_MzValues_[i] - xmin) * dx + x0, 
									this->scale(ar_Scan.mr_Spectrum.md_IntensityValues_[i] / ymax) / ymaxScaled * dy + y0));
		lk_Painter.drawPolyline(QPolygonF(lk_Points));
		
		/*
		typedef tk_IntPair tk_Pair;
		foreach (tk_Pair lk_Pair, ar_QuantitationResult.mk_Highlight)
		{
			// draw highlight
			lk_Pen.setWidthF(1.0);
			lk_Pen.setJoinStyle(Qt::RoundJoin);
			lk_Pen.setColor(QColor(192, 0, 0));
			lk_Painter.setPen(lk_Pen);
			
			QVector<QPointF> lk_Points;
			for (int i = lk_Pair.first; i <= lk_Pair.second; ++i)
				lk_Points.append(QPointF((ar_Scan.mr_Spectrum.md_MzValues_[i] - xmin) * dx + x0, 
										this->scale(ar_Scan.mr_Spectrum.md_IntensityValues_[i] / ymax) / ymaxScaled * dy + y0));
			lk_Painter.drawPolyline(QPolygonF(lk_Points));
		}
		*/
		
		foreach (r_Peak lr_Peak, ar_QuantitationResult.mk_Peaks)
		{
			// draw Gaussian
			lk_Pen.setWidthF(1.0);
			lk_Pen.setJoinStyle(Qt::RoundJoin);
			lk_Pen.setColor(QColor(0, 128, 255, 255));
			lk_Pen.setStyle(Qt::SolidLine);
			lk_Painter.setPen(lk_Pen);
			
			QVector<QPointF> lk_Points;
			double w = sqrt(lr_Peak.md_GaussC) / 2.0;
			for (double x = lr_Peak.md_PeakMz - w; x <= lr_Peak.md_PeakMz + w; x += w * 0.05)
				lk_Points.append(QPointF((x - xmin) * dx + x0, 
										this->scale(gaussian(x, lr_Peak.md_GaussA, lr_Peak.md_GaussB, lr_Peak.md_GaussC) / ymax) / ymaxScaled * dy + y0));
			lk_Painter.drawPolyline(QPolygonF(lk_Points));
		}
	}
	lk_Buffer.close();
	lk_Buffer.open(QBuffer::ReadOnly);
	lk_Buffer.seek(0);
	QString ls_Result = QString(lk_Buffer.readAll());
	lk_Buffer.close();
	QString ls_Labels;
	ls_Labels += QString("<text x='4.0' y='%1' style='font-size:10px; font-family:Verdana;'>1%</text>\n").arg(ld_OnePercent + 13.0);
	ls_Labels += QString("<text x='4.0' y='%1' style='font-size:10px; font-family:Verdana;'>10%</text>\n").arg(ld_TenPercent + 13.0);
	ls_Labels += QString("<text x='4.0' y='%1' style='font-size:10px; font-family:Verdana;'>100%</text>\n").arg(ld_OneHundredPercent + 13.0);
	foreach (r_Peak lr_Peak, ar_QuantitationResult.mk_Peaks)
	{
		double ld_X, ld_Y;
		ld_X = (lr_Peak.md_PeakMz - xmin) * dx + x0;
		ld_Y = this->scale(lr_Peak.md_PeakIntensity / ymax) / ymaxScaled * dy + y0;
		char lc_Number_[1024];
		sprintf(lc_Number_, "%.2g", lr_Peak.md_PeakIntensity);
		ls_Labels += QString("<text x='%1' y='%2' style='font-size:10px; font-family:Verdana;'>%3</text>\n").arg(ld_X + 3.0).arg(ld_Y + 10.0).arg(lc_Number_);
	}
	ls_Result.replace("</svg>", ls_Labels + "</svg>");
	return ls_Result;
}


double k_Quantifier::scale(const double ad_Value) const
{
	return log(ad_Value * 100.0 + 1.0);
}


void k_Quantifier::calculateMeanAndStandardDeviation(QList<double> ak_Values, double* ad_Mean_, double* ad_StandardDeviation_)
{
	double ld_Mean = 0.0;
	foreach (double ld_Value, ak_Values)
		ld_Mean += ld_Value;
		
	ld_Mean /= ak_Values.size();
	
	double ld_StandardDeviation = 0.0;
	foreach (double ld_Value, ak_Values)
		ld_StandardDeviation += pow(ld_Value - ld_Mean, 2.0);
		
	ld_StandardDeviation /= ak_Values.size();
	ld_StandardDeviation = sqrt(ld_StandardDeviation);
	
	*ad_Mean_ = ld_Mean;
	*ad_StandardDeviation_ = ld_StandardDeviation;
}


r_ScanQuantitationResult 
k_Quantifier::checkResult(QList<r_Peak> ak_Peaks, r_Peak* ar_ForbiddenPeak_,
						  r_Scan& ar_Scan, QString as_Peptide, int ai_Charge,
						  QList<double> ak_TargetMz)
{
	r_ScanQuantitationResult lr_Result;
	lr_Result.mb_IsGood = false;
	
	// check whether unlabeled peaks are descending
	for (int i = 1; i < mi_WatchIsotopesCount; ++i)
	{
		if (ak_Peaks[i - 1].md_PeakIntensity <= ak_Peaks[i].md_PeakIntensity)
			return lr_Result;
	}
	
	updateFailureReason(ar_Scan.ms_Id, r_QuantitationFailureReason::ForbiddenPeakPresent);
	
	// determine gap maximum
	double ld_GapMaximum = 0.0;
	foreach (r_Peak lr_Peak, ak_Peaks)
		ld_GapMaximum = std::max<double>(ld_GapMaximum, lr_Peak.md_OutsideBorderMaxIntensity);

	// make sure that there is no forbidden left isotope peak 
	//if (ar_ForbiddenPeak_ && ar_ForbiddenPeak_->md_PeakIntensity > ld_GapMaximum)
	if (ar_ForbiddenPeak_)
		return lr_Result;
	
	updateFailureReason(ar_Scan.ms_Id, r_QuantitationFailureReason::LowSnr);
		
	double ld_MinSnr = ak_Peaks.first().md_Snr;
	for (int i = 1; i < ak_Peaks.size(); ++i)
		ld_MinSnr = std::min<double>(ld_MinSnr, ak_Peaks[i].md_Snr);
		
	if (ld_MinSnr < md_MinSnr)
		return lr_Result;
	
	updateFailureReason(ar_Scan.ms_Id, r_QuantitationFailureReason::Success);
		
	lr_Result.md_Snr = ld_MinSnr;
	
	double ld_UnlabeledSum = 0.0;
	double ld_LabeledSum = 0.0;
	for (int i = 0; i < mi_WatchIsotopesCount; ++i)
	{
		ld_UnlabeledSum += ak_Peaks[i].md_PeakArea;
		ld_LabeledSum += ak_Peaks[i + mi_WatchIsotopesCount].md_PeakArea;
	}
	lr_Result.md_AmountUnlabeled = ld_UnlabeledSum;
	lr_Result.md_AmountLabeled = ld_LabeledSum;
	lr_Result.md_Ratio = ld_UnlabeledSum / ld_LabeledSum;
	
	for (int i = 0; i < mi_WatchIsotopesCount * 2; ++i)
		lr_Result.mk_Highlight.push_back(tk_IntPair(ak_Peaks[i].mi_LeftBorderIndex, ak_Peaks[i].mi_RightBorderIndex));
		
	lr_Result.mk_TargetMz = ak_TargetMz;
		
	lr_Result.md_MinMz = ak_Peaks.first().md_PeakMz - 2.0 / ai_Charge - 0.5;
	lr_Result.md_MaxMz = ak_Peaks.last().md_PeakMz + 0.5;
	
	lr_Result.mb_IsGood = true;
	lr_Result.mk_Peaks = ak_Peaks;
	
	return lr_Result;
}


void k_Quantifier::fitGaussian(double* a_, double* b_, double* c_, double x0, double y0, 
							   double x1, double y1, double x2, double y2)
{
	double A = 2.0 * (x1 - x0);
	double B = x0 * x0 - x1 * x1;
	double C = 2.0 * (x2 - x0);
	double D = x0 * x0 - x2 * x2;
	double lny0 = log(y0);
	double lny1 = log(y1);
	double lny2 = log(y2);
	double F = lny1 - lny0;
	double G = lny2 - lny0;
	double b = (F * D - B * G) / (A * G - F * C);
	double x0b = x0 - b;
	double x1b = x1 - b;
	double d = ((x0b * x0b) - (x1b * x1b)) / (lny1 - lny0);
	double a = y0 / (exp(-((x0b * x0b) / d)));
	double c = sqrt(d * 0.5);
	*a_ = a;
	*b_ = b;
	*c_ = c;
}
							   

double k_Quantifier::gaussian(double x, double a, double b, double c)
{
	return a * exp(-(pow(x - b, 2.0) / (2 * c * c)));
}


void k_Quantifier::updateFailureReason(const QString& as_Key, r_QuantitationFailureReason::Enumeration ae_Value)
{
	if (!mk_ScanFailureReason.contains(as_Key))
		mk_ScanFailureReason[as_Key] = ae_Value;
	mk_ScanFailureReason[as_Key] = (r_QuantitationFailureReason::Enumeration)(std::max<int>((int)mk_ScanFailureReason[as_Key], (int)ae_Value));
}
