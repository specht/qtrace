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
#include <math.h> 
#include <limits>


k_Quantifier::k_Quantifier(r_LabelType::Enumeration ae_LabelType,
						   r_ScanType::Enumeration ae_ScanType,
						   QList<tk_IntPair> ak_MsLevels,
						   int ai_IsotopeCount, int ai_MinCharge, int ai_MaxCharge, 
						   double ad_MinSnr, double ad_MassAccuracy, 
						   double ad_ExcludeMassAccuracy,
						   QIODevice* ak_CsvOutDevice_, QIODevice* ak_XhtmlOutDevice_,
						   bool ab_PrintStatistics)
	: k_ScanIterator(ae_ScanType, ak_MsLevels)
	, me_LabelType(ae_LabelType)
	, mi_WatchIsotopesCount(ai_IsotopeCount)
	, mk_CsvOutStream(ak_CsvOutDevice_)
	, mk_XhtmlOutStream(ak_XhtmlOutDevice_)
	, mi_MinCharge(ai_MinCharge)
	, mi_MaxCharge(ai_MaxCharge)
	, md_MinSnr(ad_MinSnr)
	, md_MassAccuracy(ad_MassAccuracy)
	, md_ExcludeMassAccuracy(ad_ExcludeMassAccuracy)
	, md_ElutionProfilePeakWidth(0.1)
	, mb_PrintStatistics(ab_PrintStatistics)
	, mui_QuantitationResultCount(0)
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
	if (md_ExcludeMassAccuracy < 0.0)
	{
		printf("Error: exclude mass accuracy must not be less than zero.\n");
		exit(1);
	}
	
	Q_INIT_RESOURCE(qtrace);
	
	QFile lk_File(":res/AminoAcids.csv");
	lk_File.open(QFile::ReadOnly);
	QTextStream lk_TextStream(&lk_File);

	lk_TextStream.readLine();
	while (!lk_TextStream.atEnd())
	{
		QString ls_Line = lk_TextStream.readLine().trimmed();
		if (ls_Line.isEmpty() || ls_Line.startsWith("#"))
			continue;

		QStringList lk_List;
		foreach (QString ls_Entry, ls_Line.split(QChar(';')))
		{
			if (ls_Entry.at(0) == QChar('"') && ls_Entry.at(ls_Entry.length() - 1) == QChar('"'))
				ls_Entry = ls_Entry.mid(1, ls_Entry.length() - 2);
			lk_List << ls_Entry;
		}
		char lc_AminoAcid = lk_List[3][0].toAscii();
		mk_AminoAcidWeight[lc_AminoAcid] = lk_List[4].toDouble();
		QString ls_Composition = lk_List[6];
		ls_Composition.remove("\"");
		QHash<char, int> lk_ElementCount;
		while (!ls_Composition.isEmpty())
		{
			char lc_Element = ls_Composition.at(0).toAscii();
			lk_ElementCount[lc_Element] = 1;
			ls_Composition = ls_Composition.right(ls_Composition.length() - 1);
			if (ls_Composition.isEmpty())
				break;
			QString ls_Number;
			char lc_NextChar = ls_Composition.at(0).toAscii();
			while (lc_NextChar >= '0' && lc_NextChar <= '9')
			{
				ls_Number += QChar(lc_NextChar);
				ls_Composition = ls_Composition.right(ls_Composition.length() - 1);
				if (ls_Composition.isEmpty())
					break;
				lc_NextChar = ls_Composition.at(0).toAscii();
			}
			if (!ls_Number.isEmpty())
				lk_ElementCount[lc_Element] = ls_Number.toInt();
		}
		int li_NitrogenCount = 0;
		if (lk_ElementCount.contains('N'))
			li_NitrogenCount = lk_ElementCount['N'];
		mk_AminoAcidNitrogenCount[lc_AminoAcid] = li_NitrogenCount;
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
	
	printf("Looking in %s scans, expecting %d isotope peaks at a presence mass accuracy of %1.2f ppm and an absence mass accuracy of %1.2f ppm.\n",
		lk_ScanTypes.join("/").toStdString().c_str(), mi_WatchIsotopesCount, md_MassAccuracy, md_ExcludeMassAccuracy);
		
	// determine all target m/z values
	typedef QPair<double, QString> tk_DoubleStringPair;
	QList<tk_DoubleStringPair> lk_TempList;
	for (int li_PeptideIndex = 0; li_PeptideIndex < mk_Peptides.size(); ++li_PeptideIndex)
	{
		QString ls_Peptide = mk_Peptides[li_PeptideIndex];
		switch (me_LabelType)
		{
			case r_LabelType::HeavyArginine:
				mk_LabeledEnvelopeCountForPeptide[ls_Peptide] = 1;
				break;
			case r_LabelType::HeavyArginineAndProline:
				mk_LabeledEnvelopeCountForPeptide[ls_Peptide] = ls_Peptide.count("P") + 1;
				break;
			case r_LabelType::N15Labeling:
				mk_LabeledEnvelopeCountForPeptide[ls_Peptide] = 1;
				break;
			default:
				printf("Invalid label type specified!\n");
				break;
		}
		for (int li_Charge = mi_MinCharge; li_Charge <= mi_MaxCharge; ++li_Charge)
		{
			double ld_PeptideMz = this->calculatePeptideMass(ls_Peptide, li_Charge);
			double ld_ModMz = 0.0;
			if (me_LabelType == r_LabelType::N15Labeling)
			{
				int li_NitrogenCount = 0;
				for (int i = 0; i < ls_Peptide.length(); ++i)
					li_NitrogenCount += mk_AminoAcidNitrogenCount[ls_Peptide.at(i).toAscii()];
// 				printf("%d N in %s\n", li_NitrogenCount, ls_Peptide.toStdString().c_str());
				ld_ModMz = HEAVY_NITROGEN * li_NitrogenCount;
			}
			else
				ld_ModMz = HEAVY_ARGININE * ls_Peptide.count("R");

			double ld_Mz;
			QString ls_Key;
			
			for (int i = 0; i < mi_WatchIsotopesCount; ++i)
			{
				// save unlabeled mass
				ld_Mz = ld_PeptideMz + i * NEUTRON / li_Charge;
// 				printf("%s (%d+): %1.6f\n", ls_Peptide.toStdString().c_str(), li_Charge, ld_Mz);
				ls_Key = QString("%1-%2-unlabeled-%3").arg(ls_Peptide).arg(li_Charge).arg(i);
				lk_TempList.push_back(tk_DoubleStringPair(ld_Mz, ls_Key));
				
				// save labeled mass
				for (int k = 0; k < mk_LabeledEnvelopeCountForPeptide[ls_Peptide]; ++k)
				{
					ld_Mz = ld_PeptideMz + i * NEUTRON / li_Charge + (ld_ModMz + k * HEAVY_PROLINE) / li_Charge;
// 					printf("%s* (%d+): %1.6f\n", ls_Peptide.toStdString().c_str(), li_Charge, ld_Mz);
					ls_Key = QString("%1-%2-labeled-%3-%4").arg(ls_Peptide).arg(li_Charge).arg(k).arg(i);
					lk_TempList.push_back(tk_DoubleStringPair(ld_Mz, ls_Key));
				}
			}
			// save forbidden peak (one to the left from the unlabeled A+0 peak)
			ld_Mz = ld_PeptideMz - NEUTRON / li_Charge;
			ls_Key = QString("%1-%2-forbidden").arg(ls_Peptide).arg(li_Charge);
			lk_TempList.push_back(tk_DoubleStringPair(ld_Mz, ls_Key));
			// save labeled forbidden peak (one to the left from the labeled A+0 peak)
			for (int k = 0; k < mk_LabeledEnvelopeCountForPeptide[ls_Peptide]; ++k)
			{
				ld_Mz = ld_PeptideMz - NEUTRON / li_Charge + (ld_ModMz + k * HEAVY_PROLINE) / li_Charge;
				ls_Key = QString("%1-%2-forbidden-labeled-%3").arg(ls_Peptide).arg(li_Charge).arg(k);
				lk_TempList.push_back(tk_DoubleStringPair(ld_Mz, ls_Key));
			}
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
	
	if (mk_CsvOutStream.device())
		mk_CsvOutStream << "id,filename,scan id,peptide,amount light,amount heavy,retention time,charge,filter line,snr" << endl;
	
	if (mk_XhtmlOutStream.device())
	{
		QFile lk_File(":res/qtrace-xhtml-header.xhtml.part");
		lk_File.open(QIODevice::ReadOnly);
		QByteArray lk_Content = lk_File.readAll();
		mk_XhtmlOutStream << QString(lk_Content);
		lk_File.close();
	}
		
	// parse all bands
	foreach (QString ls_Path, ak_SpectraFiles)
	{
		mk_ScanFailureReason.clear();
		ms_CurrentSpot = QFileInfo(ls_Path).baseName();
		
		// parse spot
		this->parseFile(ls_Path);
		
		printf(" done.\n");
		
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

	if (mk_XhtmlOutStream.device())
	{
		QFile lk_File(":res/qtrace-xhtml-footer.xhtml.part");
		lk_File.open(QIODevice::ReadOnly);
		QByteArray lk_Content = lk_File.readAll();
		mk_XhtmlOutStream << QString(lk_Content);
		lk_File.close();
	}
}


void k_Quantifier::handleScan(r_Scan& ar_Scan)
{
/*	if (QVariant(ar_Scan.ms_Id).toInt() != 3215)
		return;*/
	
	if (ar_Scan.mr_Spectrum.mi_PeaksCount == 0)
	{
		printf("Warning: Empty spectrum (scan #%s @ %1.2f minutes)!\n", ar_Scan.ms_Id.toStdString().c_str(), ar_Scan.md_RetentionTime);
		return;
	}
	
	updateFailureReason(ar_Scan.ms_Id, r_QuantitationFailureReason::NoMatchedTargetMass);
	
	QHash<QString, int> lk_TargetPeakStartIndex;
	QHash<QString, int> lk_TargetPeakEndIndex;
	
	// find all peaks in this spectrum
	QList<r_Peak> lk_AllPeaks = k_ScanIterator::findAllPeaks(ar_Scan.mr_Spectrum);
// 	printf("all peaks: %d\n", lk_AllPeaks.size());
	
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
	
/*	foreach (QString s, mk_TargetMzIndex.keys())
	{
		double ld_ScanMz = lk_AllPeaks[lk_PeakForTargetMz[mk_TargetMzIndex[s]]].md_PeakMz;
		double ld_TargetMz = mk_AllTargetMasses[mk_TargetMzIndex[s]];
		double ld_Error = fabs(ld_ScanMz - ld_TargetMz) / ld_TargetMz * 1000000.0;
		printf("%s: %1.6f - %1.6f (%1.2f ppm)\n", s.toStdString().c_str(), ld_ScanMz, ld_TargetMz, ld_Error);
	}*/
	
	// discard all target m/z matches that are bogus because
	// they are not within the specified mass accuracy
	QHash<int, int> lk_IncludePeakForTargetMz;
	QHash<int, int> lk_ExcludePeakForTargetMz;
	foreach (int li_TargetMzIndex, lk_PeakForTargetMz.keys())
	{
		double ld_TargetMz = mk_AllTargetMasses[li_TargetMzIndex];
		r_Peak lr_Peak = lk_AllPeaks[lk_PeakForTargetMz[li_TargetMzIndex]];
		double ld_PeakMz = lr_Peak.md_PeakMz;
		double ld_MzError = fabs(ld_TargetMz - ld_PeakMz);
		double ld_MaxMzError = ld_TargetMz * md_MassAccuracy / 1000000.0;
		if (ld_MzError <= ld_MaxMzError)
		{
			lk_IncludePeakForTargetMz[li_TargetMzIndex] = lk_PeakForTargetMz[li_TargetMzIndex];
			updateFailureReason(ar_Scan.ms_Id, r_QuantitationFailureReason::IsotopePeaksMissing);
		}
		ld_MaxMzError = ld_TargetMz * md_ExcludeMassAccuracy / 1000000.0;
		if (ld_MzError <= ld_MaxMzError)
			lk_ExcludePeakForTargetMz[li_TargetMzIndex] = lk_PeakForTargetMz[li_TargetMzIndex];
	}
	
	// now check for each peptide/charge whether all peaks are there
	foreach (QString ls_Peptide, mk_Peptides)
	{
		for (int li_Charge = mi_MinCharge; li_Charge <= mi_MaxCharge; ++li_Charge)
		{
			QHash<int, r_Peak> lk_LightPeaksInclude;
			QHash<int, r_Peak> lk_HeavyPeaksInclude;
			QHash<int, r_Peak> lk_LightPeaksExclude;
			QHash<int, r_Peak> lk_HeavyPeaksExclude;
			QList<double> lk_UnlabeledTargetMz;
			QList<double> lk_LabeledTargetMz;
			for (int li_Isotope = 0; li_Isotope < mi_WatchIsotopesCount; ++li_Isotope)
			{
				QString ls_Key = QString("%1-%2-unlabeled-%3").
					arg(ls_Peptide).arg(li_Charge).arg(li_Isotope);
				int li_TargetMzIndex = mk_TargetMzIndex[ls_Key];
				lk_UnlabeledTargetMz.push_back(mk_AllTargetMasses[li_TargetMzIndex]);
				if (lk_IncludePeakForTargetMz.contains(li_TargetMzIndex))
					lk_LightPeaksInclude[li_Isotope] = lk_AllPeaks[lk_IncludePeakForTargetMz[li_TargetMzIndex]];
				if (lk_ExcludePeakForTargetMz.contains(li_TargetMzIndex))
					lk_LightPeaksExclude[li_Isotope] = lk_AllPeaks[lk_ExcludePeakForTargetMz[li_TargetMzIndex]];
			} 
			for (int k = 0; k < mk_LabeledEnvelopeCountForPeptide[ls_Peptide]; ++k)
			{
				for (int li_Isotope = 0; li_Isotope < mi_WatchIsotopesCount; ++li_Isotope)
				{
					QString ls_Key = QString("%1-%2-labeled-%3-%4").
						arg(ls_Peptide).arg(li_Charge).arg(k).arg(li_Isotope);
					int li_TargetMzIndex = mk_TargetMzIndex[ls_Key];
					lk_LabeledTargetMz.push_back(mk_AllTargetMasses[li_TargetMzIndex]);
					if (lk_IncludePeakForTargetMz.contains(li_TargetMzIndex))
						lk_HeavyPeaksInclude[k * mi_WatchIsotopesCount + li_Isotope] = lk_AllPeaks[lk_IncludePeakForTargetMz[li_TargetMzIndex]];
					if (lk_ExcludePeakForTargetMz.contains(li_TargetMzIndex))
						lk_HeavyPeaksExclude[k * mi_WatchIsotopesCount + li_Isotope] = lk_AllPeaks[lk_ExcludePeakForTargetMz[li_TargetMzIndex]];
				}
				QString ls_Key = QString("%1-%2-forbidden-labeled-%3").arg(ls_Peptide).arg(li_Charge).arg(k);
				int li_TargetMzIndex = mk_TargetMzIndex[ls_Key];
				if (lk_ExcludePeakForTargetMz.contains(li_TargetMzIndex))
					lk_HeavyPeaksExclude[-1 - k] = lk_AllPeaks[lk_ExcludePeakForTargetMz[li_TargetMzIndex]];
			}
			QString ls_Key = QString("%1-%2-forbidden").arg(ls_Peptide).arg(li_Charge);
			int li_TargetMzIndex = mk_TargetMzIndex[ls_Key];
			if (lk_ExcludePeakForTargetMz.contains(li_TargetMzIndex))
				lk_LightPeaksExclude[-1] = lk_AllPeaks[lk_ExcludePeakForTargetMz[li_TargetMzIndex]];
			
			r_ScanQuantitationResult lr_ScanResult = 
				this->checkResult(lk_LightPeaksInclude, lk_HeavyPeaksInclude,
								   lk_LightPeaksExclude, lk_HeavyPeaksExclude,
								   ar_Scan, ls_Peptide, li_Charge, 
								   lk_UnlabeledTargetMz, lk_LabeledTargetMz);
			if (lr_ScanResult.mb_IsGood)
			{
				lr_ScanResult.ms_ScanHashKey = QString("%1.%2.%3.%4").arg(ms_CurrentSpot).arg(ar_Scan.ms_Id).arg(ls_Peptide).arg(li_Charge);
				lr_ScanResult.md_RetentionTime = ar_Scan.md_RetentionTime;
				lr_ScanResult.mi_Charge = li_Charge;
				lr_ScanResult.ms_ScanId = ar_Scan.ms_Id;
				/*
				mk_ScanHash.insert(lr_ScanResult.ms_ScanHashKey, r_Scan(ar_Scan));
				if (!mk_SpotResults.contains(ls_Peptide))
					mk_SpotResults[ls_Peptide] = QList<r_ScanQuantitationResult>();
				mk_SpotResults[ls_Peptide].push_back(lr_ScanResult);

				*/
				++mui_QuantitationResultCount;
			
				if (mk_CsvOutStream.device())
				{
// 						mk_CsvOutStream << "id,filename,scan id,peptide,amount light,amount heavy,retention time,charge,filter line,snr" << endl;
					mk_CsvOutStream << mui_QuantitationResultCount
						<< ",\"" << ms_CurrentSpot << "\""
						<< ",\"" << ar_Scan.ms_Id << "\""
						<< ",\"" << ls_Peptide << "\""
						<< "," << lr_ScanResult.md_AmountUnlabeled
						<< "," << lr_ScanResult.md_AmountLabeled
						<< "," << ar_Scan.md_RetentionTime
						<< "," << lr_ScanResult.mi_Charge
						<< ",\"" << ar_Scan.ms_FilterLine << "\""
						<< "," << lr_ScanResult.md_Snr
						<< endl;
				}

/*
				if (mk_XhtmlOutStream.device())
				{
					mk_XhtmlOutStream << QString("\n<!-- BEGIN PEPTIDE %1 -->\n").arg(ls_Peptide); 
					mk_XhtmlOutStream << "<tr style='cursor: pointer;' onmouseover='toggle(\"s" << mui_QuantitationResultCount << "\", \"block\")' onmouseout='toggle(\"s" << mui_QuantitationResultCount << "\", \"none\")'><td>" << mui_QuantitationResultCount << "</td>"
						<< "<td>" << ms_CurrentSpot << "</td>"
						<< "<td>" << ar_Scan.ms_Id << "</td>"
						<< "<td>" << ar_Scan.md_RetentionTime << "</td>"
						<< "<td>" << ls_Peptide << "</td>"
						<< "<td>" << lr_ScanResult.mi_Charge << "</td>"
						<< "<td>" << ar_Scan.ms_FilterLine << "</td>"
						<< "<td>" << lr_ScanResult.md_Snr << "</td>"
						<< "<td>" << lr_ScanResult.md_AmountUnlabeled << "</td>"
						<< "<td>" << lr_ScanResult.md_AmountLabeled << "</td>"
						<< "</tr>"
						<< endl;
					QString ls_Svg = this->renderScanAsSvg(ar_Scan, lr_ScanResult);
					ls_Svg.remove(QRegExp("<\\?xml.+\\?>"));
					ls_Svg.replace(QRegExp("width=\\\"[^\\\"]*\\\"\\s+height=\\\"[^\\\"]*\\\""), "width='950'");
					mk_XhtmlOutStream << "<div id='s" << mui_QuantitationResultCount << "' style='display: none; position: absolute; left: 32px; border: 1px solid #000; padding: 8px; background-color: #fff;'>";
					mk_XhtmlOutStream << ls_Svg;
					mk_XhtmlOutStream << "</div>" << endl;
					mk_XhtmlOutStream << QString("\n<!-- END PEPTIDE %1 -->\n").arg(ls_Peptide); 
				}
				*/
			}
		}
	}
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


QString k_Quantifier::renderScanAsSvg(r_Scan& ar_Scan, r_ScanQuantitationResult ar_QuantitationResult)
{
	/*
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
		foreach (double ld_Mz, ar_QuantitationResult.mk_UnlabeledTargetMz)
		{
			double x = (ld_Mz - xmin) * dx + x0;
			lk_Painter.drawLine(QPointF(x, y0), QPointF(x, y0 + dy));
		}
		foreach (double ld_Mz, ar_QuantitationResult.mk_LabeledTargetMz)
		{
			double x = (ld_Mz - xmin) * dx + x0;
			lk_Painter.drawLine(QPointF(x, y0), QPointF(x, y0 + dy));
		}
		
		// draw forbidden peak
		double ld_Mz = ar_QuantitationResult.mk_UnlabeledTargetMz.first() - NEUTRON / ar_QuantitationResult.mi_Charge;
		double x = (ld_Mz - xmin) * dx + x0;
		lk_Pen.setColor(QColor(128, 0, 0, 192));
		lk_Painter.setPen(lk_Pen);
		lk_Painter.drawLine(QPointF(x, y0), QPointF(x, y0 + dy));
		
		lk_Pen.setColor(QColor(192, 192, 192));
		lk_Painter.setPen(lk_Pen);
		
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
		
		foreach (r_Peak lr_Peak, (ar_QuantitationResult.mk_UnlabeledPeaks + ar_QuantitationResult.mk_LabeledPeaks))
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
	foreach (r_Peak lr_Peak, (ar_QuantitationResult.mk_UnlabeledPeaks + ar_QuantitationResult.mk_LabeledPeaks))
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
	*/
	return QString();
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
k_Quantifier::checkResult(QHash<int, r_Peak> ak_LightPeaksInclude, 
						   QHash<int, r_Peak> ak_HeavyPeaksInclude,
						   QHash<int, r_Peak> ak_LightPeaksExclude, 
						   QHash<int, r_Peak> ak_HeavyPeaksExclude,
						   r_Scan& ar_Scan, QString as_Peptide, int ai_Charge,
						   QList<double> ak_UnlabeledTargetMz,
						   QList<double> ak_LabeledTargetMz)
{
	r_ScanQuantitationResult lr_Result;
	lr_Result.mb_IsGood = false;
	
	updateFailureReason(ar_Scan.ms_Id, r_QuantitationFailureReason::LowSnr);

	// check whether all inclusion peaks have a good SNR
	QList<r_Peak> lk_Peaks = ak_LightPeaksInclude.values() + ak_HeavyPeaksInclude.values();
	
	// return if no inclusion peaks have been found
	if (lk_Peaks.empty())
		return lr_Result;
	
	double ld_MinSnr = lk_Peaks.first().md_Snr;
	for (int i = 1; i < lk_Peaks.size(); ++i)
		ld_MinSnr = std::min<double>(ld_MinSnr, lk_Peaks[i].md_Snr);
	
	if (ld_MinSnr < md_MinSnr)
		return lr_Result;
	
	updateFailureReason(ar_Scan.ms_Id, r_QuantitationFailureReason::Success);
		
	lr_Result.md_Snr = ld_MinSnr;
	
	double ld_LightSum = 0.0;
	double ld_HeavySum = 0.0;
	
	int li_LightPeakIncludeCount = 0;
	int li_HeavyPeakIncludeCount = 0;
	int li_LightPeakExcludeCount = 0;
	int li_HeavyPeakExcludeCount = 0;
	
	for (int li_Isotope = 0; li_Isotope < mi_WatchIsotopesCount; ++li_Isotope)
	{
		if (ak_LightPeaksInclude.contains(li_Isotope))
		{
			ld_LightSum += ak_LightPeaksInclude[li_Isotope].md_PeakArea;
			++li_LightPeakIncludeCount;
		}
		if (ak_LightPeaksExclude.contains(li_Isotope))
			++li_LightPeakExcludeCount;
		for (int k = 0; k < mk_LabeledEnvelopeCountForPeptide[as_Peptide]; ++k)
		{
			if (ak_HeavyPeaksInclude.contains(k * mi_WatchIsotopesCount + li_Isotope))
			{
				ld_HeavySum += ak_HeavyPeaksInclude[k * mi_WatchIsotopesCount + li_Isotope].md_PeakArea;
				++li_HeavyPeakIncludeCount;
			}
			if (ak_HeavyPeaksExclude.contains(k * mi_WatchIsotopesCount + li_Isotope))
				++li_HeavyPeakExcludeCount;
		}
	}
	
	if (li_LightPeakIncludeCount == mi_WatchIsotopesCount && 
		li_HeavyPeakIncludeCount == mi_WatchIsotopesCount * mk_LabeledEnvelopeCountForPeptide[as_Peptide])
	{
		// both isotope envelopes are complete, we get a ratio!
		// don't quantify if forbidden peak is present!
		if (ak_LightPeaksExclude.contains(-1))
			return lr_Result;
		lr_Result.md_AmountUnlabeled = ld_LightSum;
		lr_Result.md_AmountLabeled = ld_HeavySum;
		lr_Result.md_Ratio = ld_LightSum / ld_HeavySum;
	}
	else if (li_LightPeakIncludeCount == mi_WatchIsotopesCount && 
			 li_HeavyPeakExcludeCount == 0)
	{
		// light state only!
		if (ak_LightPeaksExclude.contains(-1))
			return lr_Result;
		lr_Result.md_AmountUnlabeled = ld_LightSum;
		lr_Result.md_AmountLabeled = 0.0;
		lr_Result.md_Ratio = std::numeric_limits<double>::infinity(); 
	}
	else if (li_LightPeakExcludeCount == 0 && 
			 li_HeavyPeakIncludeCount == mi_WatchIsotopesCount * mk_LabeledEnvelopeCountForPeptide[as_Peptide])
	{
		// heavy state only
		for (int k = 0; k < mk_LabeledEnvelopeCountForPeptide[as_Peptide]; ++k)
			if (ak_HeavyPeaksExclude.contains(-1 - k))
				return lr_Result;
		lr_Result.md_AmountUnlabeled = 0.0;
		lr_Result.md_AmountLabeled = ld_HeavySum;
		lr_Result.md_Ratio = 0.0;
	}
	else
	{
		// no complete isotope envelope
		return lr_Result;
	}
	
	lr_Result.mk_UnlabeledTargetMz = ak_UnlabeledTargetMz;
	lr_Result.mk_LabeledTargetMz = ak_LabeledTargetMz;
	
	QList<double> lk_TargetMz = ak_UnlabeledTargetMz + ak_LabeledTargetMz;
	qSort(lk_TargetMz);
		
	lr_Result.md_MinMz = lk_TargetMz.first() - 2.0 / ai_Charge - 0.5;
	lr_Result.md_MaxMz = lk_TargetMz.last() + 0.5;
	
	lr_Result.mk_UnlabeledPeaks = ak_LightPeaksInclude.values();
	lr_Result.mk_LabeledPeaks = ak_HeavyPeaksInclude.values();
	lr_Result.mb_IsGood = true;
	
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
